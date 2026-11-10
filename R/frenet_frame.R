.vc_curve_cross <- function(a,b) {
  c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
}

#' Evaluate the Frenet frame of a spatial curve
#' @inheritParams vc_curve_sample
#' @param t One finite real parameter value.
#' @param velocity,acceleration Optional functions giving the first and second
#'   derivatives in Cartesian order. Each supplied function must return a
#'   finite real three-vector. The caller supplies its mathematical validity.
#' @param h,h2 Optional positive initial central steps for velocity and
#'   acceleration, respectively. Defaults use the Jacobian and Hessian scales.
#' @return List with type frenet_frame, t, point, velocity, acceleration,
#'   unit tangent, principal normal, binormal, status, derivative sources and
#'   numerical diagnostics. Each unavailable direction is a three-vector of NA.
#' @details The frame follows increasing parameter t. It uses B in the direction
#'   of the cross product of velocity and acceleration and N in the direction
#'   of B crossed with T. Zero supplied velocity gives zero_speed. A zero or
#'   numerically unresolved normal component leaves only the tangent available
#'   (unresolved_curvature). Numerical velocity can give unresolved_speed.
#'   Resolution uses the last refinement change and a roundoff indicator for
#'   each numerical derivative, plus machine precision for the cross product.
#'   These indicators and the stored sensitivity diagnostics are heuristic.
#'   Central stencils may leave the curve's domain or cross a nonsmooth point.
#'   Supply derivatives at such points when the required derivatives exist.
#'   Use lapply() to evaluate the frame at several parameters.
#' @examples
#' helix <- function(t)c(cos(t),sin(t),t/2)
#' vc_frenet_frame(helix,0,velocity=function(t)c(-sin(t),cos(t),.5),
#'   acceleration=function(t)c(-cos(t),-sin(t),0))
#' @export
vc_frenet_frame <- function(curve,t,velocity=NULL,acceleration=NULL,h=NULL,h2=NULL) {
  if (!is.function(curve)) stop("'curve' must be a function.")
  if (!.vc_scalar(t) || !is.null(dim(t))) stop("'t' must be a finite real scalar.")
  t <- as.double(t)
  if (!is.null(velocity) && !is.function(velocity)) stop("'velocity' must be a function or NULL.")
  if (!is.null(acceleration) && !is.function(acceleration)) stop("'acceleration' must be a function or NULL.")
  if (!is.null(h)) .vc_positive(h,"h")
  h2 <- .vc_matrix_step(t,h2,second=TRUE)
  field <- .vc_cached_field(function(s).vc_curve_vector(curve,s),t,vector=TRUE)
  point <- field$value
  v_diagnostics <- a_diagnostics <- NULL
  v_error <- a_error <- numeric(3)
  if (is.null(velocity)) {
    v_diagnostics <- vc_jacobian(field$get,t,h=h,warn=FALSE)
    v <- unname(v_diagnostics$estimate[,1])
    v_error <- vapply(v_diagnostics$results,function(entry) {
      d <- entry$derivative
      step <- d$diagnostics$step[3]
      left <- t-step
      right <- t+step
      magnitude <- max(abs(c(point[entry$row],field$get(left)[entry$row],
                            field$get(right)[entry$row])))
      d$changes[2]+(4*.Machine$double.eps*magnitude)/min(t-left,right-t)
    },numeric(1))
  } else v <- .vc_curve_vector(velocity,t)
  if (is.null(acceleration)) {
    a_diagnostics <- lapply(seq_len(3),function(i)
      .vc_second_difference(function(s)field$get(s)[i],t,h2,1e-8,1e-4))
    a <- vapply(a_diagnostics,function(d)d$estimate,numeric(1))
    a_error <- vapply(a_diagnostics,function(d)
      d$changes[2]+d$diagnostics$roundoff_indicator[3],numeric(1))
  } else a <- .vc_curve_vector(acceleration,t)

  tangent <- normal <- binormal <- rep(NA_real_,3)
  v_scale <- max(abs(v))
  a_scale <- max(abs(a))
  status <- if (v_scale==0 && all(v_error==0)) "zero_speed" else "unresolved_speed"
  if (v_scale>0) {
    v_unit_scale <- v/v_scale
    v_size <- .vc_norm(v_unit_scale)
    v_uncertainty <- v_error/v_scale/v_size
    if (all(is.finite(v_uncertainty)) && .vc_norm(v_uncertainty)<1) {
      tangent <- v_unit_scale/v_size
      status <- "unresolved_curvature"
      if (a_scale>0) {
        a_unit_scale <- a/a_scale
        cross <- .vc_curve_cross(tangent,a_unit_scale)
        cross_size <- .vc_norm(cross)
        a_uncertainty <- a_error/a_scale
        if (all(is.finite(a_uncertainty))) {
          threshold <- .vc_norm(a_uncertainty)+
            (2*.vc_norm(v_uncertainty)+8*.Machine$double.eps)*.vc_norm(a_unit_scale)
          if (cross_size>threshold) {
            binormal <- cross/cross_size
            normal <- .vc_curve_cross(binormal,tangent)
            status <- "ok"
          }
        }
      }
    }
  }
  list(type="frenet_frame",t=t,point=point,velocity=v,acceleration=a,
       tangent=tangent,normal=normal,binormal=binormal,status=status,
       derivative_sources=c(velocity=if(is.null(velocity))"numerical" else "supplied",
                            acceleration=if(is.null(acceleration))"numerical" else "supplied"),
       diagnostics=list(velocity=v_diagnostics,acceleration=a_diagnostics,
                        velocity_indicator=v_error,acceleration_indicator=a_error))
}

#' Display a stored Frenet frame on a sampled spatial curve
#' @inheritParams vc_plot_curve
#' @param frame Result of vc_frenet_frame() for the same curve at a parameter
#'   within the sampled interval.
#' @param scale Positive length of each displayed unit direction, in the
#'   physical coordinate units of the curve.
#' @param components Character vector selecting tangent, normal and/or binormal.
#' @return Plotly widget with the sampled curve, frame anchor and available
#'   selected directions. No curve or derivative functions are evaluated.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   helix <- function(t)c(cos(t),sin(t),t/2)
#'   vc_plot_frenet_frame(vc_curve_sample(helix,c(0,2*pi)),
#'     vc_frenet_frame(helix,pi/4),scale=.7)
#' }
#' @export
vc_plot_frenet_frame <- function(result,frame,scale=1,
    components=c("tangent","normal","binormal"),title=NULL) {
  if (!is.list(frame) || !identical(frame$type,"frenet_frame")) stop("Invalid Frenet frame.")
  .vc_positive(scale,"scale")
  colors <- c(tangent="#d1492f",normal="#18825a",binormal="#7b429f")
  if (!is.character(components) || !length(components) ||
      any(!components%in%names(colors)) || anyDuplicated(components))
    stop("Invalid frame components.")
  p <- vc_plot_curve(result,title=if(is.null(title))"Frenet frame" else title)
  if (frame$t<min(result$interval) || frame$t>max(result$interval))
    stop("Frame parameter is outside the sampled interval.")
  xyz <- frame$point
  p <- plotly::add_trace(p,x=xyz[1],y=xyz[2],z=xyz[3],type="scatter3d",
    mode="markers",name="Frame point",text=paste("t =",format(frame$t),frame$status),
    hoverinfo="x+y+z+text",marker=list(color="#252525",size=5),inherit=FALSE)
  for (component in components) {
    direction <- frame[[component]]
    if (anyNA(direction)) next
    end <- xyz+scale*direction
    if (any(!is.finite(end)) || all(end==xyz)) stop("Frame endpoint is not representable; change 'scale'.")
    p <- plotly::add_trace(p,x=c(xyz[1],end[1]),y=c(xyz[2],end[2]),z=c(xyz[3],end[3]),
      type="scatter3d",mode="lines+markers",name=component,
      line=list(color=colors[[component]],width=6),
      marker=list(color=colors[[component]],size=c(0,4)),inherit=FALSE)
  }
  p
}
