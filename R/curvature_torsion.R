.vc_curve_third <- function(curve,t,h3) {
  steps <- h3/c(1,2,4)
  estimates <- roundoff <- matrix(NA_real_,3,3)
  for (i in seq_along(steps)) {
    nodes <- t+c(-2,-1,1,2)*steps[i]
    widths <- diff(nodes)
    spans <- nodes[3:4]-nodes[1:2]
    if (any(!is.finite(c(nodes,widths,spans))) || any(widths<=0) ||
        nodes[2]>=t || nodes[3]<=t)
      stop("Third-derivative stencil collapsed or overflowed; change 'h3'.")
    values <- vapply(nodes,curve,numeric(3))
    slopes <- sweep(values[,2:4]-values[,1:3],2,widths,"/")
    left <- (slopes[,2]-slopes[,1])/spans[1]
    right <- (slopes[,3]-slopes[,2])/spans[2]
    estimates[i,] <- (right-left)/((nodes[4]/2-nodes[1]/2)/3)
    if (any(!is.finite(estimates[i,])))
      stop("Nonfinite third difference; change the scale or 'h3'.")
    width <- min(widths)
    roundoff[i,] <- (((16*.Machine$double.eps)*apply(abs(values),1,max))/width)/width/width
  }
  list(estimate=estimates[3,],steps=steps,estimates=estimates,
       roundoff_indicator=roundoff,
       indicator=abs(estimates[3,]-estimates[2,])+roundoff[3,])
}

.vc_curve_ratio <- function(numerator,denominator) {
  value <- exp(sum(log(numerator))-sum(log(denominator)))
  if (is.finite(value) && value>0) value else NA_real_
}

.vc_curve_bending <- function(frame) {
  curvature <- NA_real_
  status <- frame$status
  speed <- perpendicular <- NULL
  v <- frame$velocity
  a <- frame$acceleration
  a_scale <- max(abs(a))
  if (status=="unresolved_curvature" &&
      (a_scale==0 || all(.vc_curve_cross(v/max(abs(v)),a/a_scale)==0))) {
    curvature <- 0
    if (all(frame$derivative_sources=="supplied")) status <- "zero_curvature"
  }
  if (status=="ok") {
    v_scale <- max(abs(v))
    speed <- c(v_scale,.vc_norm(v/v_scale))
    cross_size <- .vc_norm(.vc_curve_cross(frame$tangent,a/a_scale))
    perpendicular <- c(a_scale,cross_size)
    curvature <- .vc_curve_ratio(perpendicular,rep(speed,2))
    if (is.na(curvature)) status <- "unrepresentable"
  }
  list(curvature=curvature,status=status,speed=speed,perpendicular=perpendicular)
}

#' Evaluate curvature and signed torsion of a spatial curve
#' @inheritParams vc_frenet_frame
#' @param jerk Optional function giving the third derivative as a finite real
#'   Cartesian three-vector. Evaluated only when the frame is resolved.
#' @param h3 Optional positive initial step for the direct third difference.
#'   The default is max(1, abs(t)) times machine epsilon to the power 1/5.
#' @return List with type curvature_torsion, frame, curvature, torsion, their
#'   separate statuses, the third derivative, its source and numerical
#'   diagnostics. Unavailable quantities are NA.
#' @details Curvature is the norm of the velocity crossed with acceleration
#'   divided by speed cubed. Torsion is the signed scalar triple product of the
#'   three derivatives divided by the squared norm of that cross product.
#'   Torsion requires resolved nonzero curvature. A vanishing computed cross
#'   product at resolved speed gives zero estimated curvature. Numerically
#'   unresolved curvature and speed are labeled. A torsion estimate with a
#'   third-derivative indicator at least as large as the binormal projection
#'   is labeled unresolved_torsion. A finite, nonzero magnitude outside
#'   the representable range is reported as NA with status unrepresentable.
#'   Numerical jerk uses direct four-point stencils at three scales. Stencils
#'   may leave the curve domain. The indicators are heuristic and do not
#'   establish smoothness or provide certified error bounds.
#' @examples
#' helix <- function(t)c(cos(t),sin(t),t/2)
#' measures <- vc_curvature_torsion(helix,0)
#' c(curvature=measures$curvature,torsion=measures$torsion)
#' @export
vc_curvature_torsion <- function(curve,t,velocity=NULL,acceleration=NULL,
                                jerk=NULL,h=NULL,h2=NULL,h3=NULL) {
  if (!is.null(jerk) && !is.function(jerk)) stop("'jerk' must be a function or NULL.")
  if (!is.null(h3)) .vc_positive(h3,"h3")
  frame <- vc_frenet_frame(curve,t,velocity,acceleration,h,h2)
  t <- frame$t
  bending <- .vc_curve_bending(frame)
  torsion <- NA_real_
  torsion_status <- "unavailable"
  third <- rep(NA_real_,3)
  diagnostics <- NULL
  source <- "not_evaluated"

  if (frame$status=="ok") {
    if (is.null(jerk)) {
      if (is.null(h3)) h3 <- max(1,abs(t))*.Machine$double.eps^(1/5)
      field <- .vc_cached_field(function(s).vc_curve_vector(curve,s),t,vector=TRUE)
      diagnostics <- .vc_curve_third(field$get,t,h3)
      third <- diagnostics$estimate
      source <- "numerical"
    } else {
      third <- .vc_curve_vector(jerk,t)
      source <- "supplied"
    }
    third_scale <- max(abs(third))
    projection <- if (third_scale==0) 0 else sum(frame$binormal*(third/third_scale))
    torsion <- if (projection==0) 0 else sign(projection)*
      .vc_curve_ratio(c(third_scale,abs(projection)),c(bending$speed,bending$perpendicular))
    torsion_status <- if (is.na(torsion)) "unrepresentable" else "ok"
    if (!is.null(diagnostics) && !is.na(torsion)) {
      active <- frame$binormal!=0
      indicator <- sum(abs(frame$binormal[active])*diagnostics$indicator[active])
      unresolved <- if (third_scale==0) indicator>0 else
        indicator>0 && indicator/third_scale>=abs(projection)
      if (unresolved) torsion_status <- "unresolved_torsion"
    }
  }
  list(type="curvature_torsion",frame=frame,curvature=bending$curvature,torsion=torsion,
       curvature_status=bending$status,torsion_status=torsion_status,
       jerk=third,jerk_source=source,diagnostics=diagnostics)
}

#' Display stored curvature and torsion on a sampled spatial curve
#' @inheritParams vc_plot_frenet_frame
#' @param measures Result of vc_curvature_torsion() for the same curve.
#' @return Plotly widget with the stored frame and a caption for the scalar
#'   measures and their statuses.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   helix <- function(t)c(cos(t),sin(t),t/2)
#'   vc_plot_curvature_torsion(vc_curve_sample(helix,c(0,2*pi)),
#'     vc_curvature_torsion(helix,pi/4),scale=.7)
#' }
#' @export
vc_plot_curvature_torsion <- function(result,measures,scale=1,title=NULL) {
  if (!is.list(measures) || !identical(measures$type,"curvature_torsion"))
    stop("Invalid curvature and torsion result.")
  p <- vc_plot_frenet_frame(result,measures$frame,scale=scale,
    title=if(is.null(title))"Curvature and torsion" else title)
  label <- function(value,status) {
    if (is.na(value)) return(paste("unavailable",status,sep=" / "))
    paste(format(value,digits=6),status,sep=" / ")
  }
  caption <- paste("Curvature:",label(measures$curvature,measures$curvature_status),
                   "| Torsion:",label(measures$torsion,measures$torsion_status))
  plotly::layout(p,annotations=list(list(text=caption,showarrow=FALSE,
    x=.5,y=0,xref="paper",yref="paper",yanchor="top")))
}
