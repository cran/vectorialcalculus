#' Compare outward flux and volume divergence in Gauss' theorem
#' @param field Function of one Cartesian vector returning a finite
#'   three-component vector.
#' @param x_range,y_range,z_range Finite strictly increasing endpoint pairs
#'   defining an axis-aligned box.
#' @param divergence Optional function of position giving the scalar
#'   divergence. Defaults to vc_divergence().
#' @param h Optional numerical divergence step or three coordinate steps.
#' @inheritParams vc_triple_integrate
#' @return List with total outward flux, volume divergence, signed difference
#'   (flux minus divergence integral), comparison tolerance, per-face flux
#'   and quadrature diagnostics, outward normals and derivative sensitivity.
#' @details The six outward unit normals are exact signed coordinate axes.
#'   The same coordinate-specific cell counts are used on each face and in
#'   the volume, with x outermost, then y and z. The caller supplies a C1
#'   field on a neighborhood of the box. Central stencils can leave the
#'   box. Numerical agreement does not establish theorem hypotheses and
#'   quadrature indicators do not include derivative error.
#' @examples
#' vc_gauss(identity,c(0,1),c(0,1),c(0,1),n_outer=2,n_middle=2,n_inner=2)$outward_flux
#' @export
vc_gauss <- function(field,x_range,y_range,z_range,divergence=NULL,h=NULL,
    method=c("simpson","trapezoid"),n_outer=4L,n_middle=4L,n_inner=4L,
    abs_tol=1e-8,rel_tol=1e-8) {
  if (!is.function(field)) stop("'field' must be a function.")
  if (!is.null(divergence) && !is.function(divergence))
    stop("'divergence' must be a function or NULL.")
  ranges <- list(x=.vc_integral_interval(x_range),y=.vc_integral_interval(y_range),
    z=.vc_integral_interval(z_range))
  widths <- vapply(ranges,function(r)r[2]-r[1],numeric(1))
  if (any(!is.finite(widths)) || any(widths<=0)) stop("Box ranges must be finite and strictly increasing.")
  method <- match.arg(method)
  h <- .vc_matrix_step(c(0,0,0),h)
  counts <- list(n_outer,n_middle,n_inner)
  evaluate <- function(p).vc_spatial_field(field(p))
  lower <- function(i)function(t)ranges[[i]][1]
  upper <- function(i)function(t)ranges[[i]][2]
  faces <- list()
  normals <- matrix(0,6,3,dimnames=list(NULL,names(ranges)))
  for (axis in seq_len(3)) {
    free <- setdiff(seq_len(3),axis)
    for (side in 1:2) {
      label <- paste(names(ranges)[axis],c("min","max")[side],sep="_")
      direction <- c(-1,1)[side]
      integrand <- function(p) {
        xyz <- numeric(3)
        xyz[axis] <- ranges[[axis]][side]
        xyz[free] <- p
        direction*evaluate(xyz)[axis]
      }
      faces[[label]] <- vc_double_integrate(integrand,ranges[[free[1]]],
        lower(free[2]),upper(free[2]),method=method,n_outer=counts[[free[1]]],
        n_inner=counts[[free[2]]],abs_tol=abs_tol,rel_tol=rel_tol)
      normals[length(faces),axis] <- direction
    }
  }
  rownames(normals) <- names(faces)
  face_flux <- vapply(faces,function(r)r$estimate,numeric(1))
  flux <- .vc_finite(sum(face_flux),"Total outward flux")
  flux_error <- .vc_finite(sum(vapply(faces,function(r)r$error_estimate,numeric(1))),
    "Boundary error indicator")
  cache <- new.env(parent=emptyenv())
  evaluations <- 0L
  reliable <- TRUE
  density <- if (!is.null(divergence)) divergence else function(p) {
    key <- paste(sprintf("%a",p),collapse=";")
    if (exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache))
    r <- vc_divergence(evaluate,p,h=h,warn=FALSE)
    evaluations <<- evaluations+1L
    reliable <<- reliable && r$reliable
    assign(key,r$estimate,envir=cache)
    r$estimate
  }
  volume <- vc_triple_integrate(density,ranges$x,lower(2),upper(2),lower(3),upper(3),
    method=method,n_outer=n_outer,n_middle=n_middle,n_inner=n_inner,
    abs_tol=abs_tol,rel_tol=rel_tol)
  difference <- .vc_finite(flux-volume$estimate,"Gauss discrepancy")
  tolerance <- .vc_finite(max(abs_tol,rel_tol*max(abs(flux),abs(volume$estimate))),
    "Comparison tolerance")
  list(type="gauss",ranges=ranges,outward_flux=flux,volume_divergence=volume$estimate,
       difference=difference,tolerance=tolerance,within_tolerance=abs(difference)<=tolerance,
       error_indicator=.vc_finite(flux_error+volume$error_estimate,"Combined error indicator"),
       face_flux=face_flux,normals=normals,faces=faces,volume=volume,
       derivative_source=if(is.null(divergence))"numerical" else "supplied",
       derivative_sensitivity_ok=if(evaluations>0L)reliable else NA,
       derivative_evaluations=evaluations,
       note="Agreement is a numerical comparison under caller-supplied divergence theorem hypotheses.")
}

#' Display the box and outward normals for a Gauss comparison
#' @param result Result of vc_gauss().
#' @param title Optional character title.
#' @return Plotly widget with outward arrows, signed face flux in hover text
#'   and both stored integral estimates. The field is not reevaluated.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   result <- vc_gauss(identity,c(0,1),c(0,1),c(0,1),divergence=function(p)3,
#'     n_outer=2,n_middle=2,n_inner=2)
#'   vc_plot_gauss(result)
#' }
#' @export
vc_plot_gauss <- function(result,title=NULL) {
  if (!is.list(result) || !identical(result$type,"gauss")) stop("Invalid Gauss comparison.")
  r <- result$ranges
  solid <- vc_solid_cartesian(function(x)r$y[1],function(x)r$y[2],
    function(p)r$z[1],function(p)r$z[2],r$x,n=1,n_u=1,n_v=1)
  p <- vc_plot_solid_cartesian(solid,opacity=.2,
    title=if(is.null(title))"Gauss divergence theorem" else title)
  center <- vapply(r,mean,numeric(1))
  anchors <- matrix(rep(center,each=6),6,3)
  for (axis in seq_len(3)) anchors[2*axis+c(-1,0),axis] <- r[[axis]]
  scale <- .2*min(vapply(r,diff,numeric(1)))
  vectors <- scale*result$normals
  label <- paste(names(result$face_flux),"/ flux =",format(result$face_flux,digits=6))
  p <- plotly::add_trace(p,x=anchors[,1],y=anchors[,2],z=anchors[,3],
    u=vectors[,1],v=vectors[,2],w=vectors[,3],type="cone",anchor="tail",
    sizemode="absolute",sizeref=1,showscale=FALSE,name="Outward normals",
    text=label,hoverinfo="text",inherit=FALSE)
  caption <- paste("Outward flux =",format(result$outward_flux,digits=6),
    "| Volume divergence =",format(result$volume_divergence,digits=6),
    "| Difference =",format(result$difference,digits=3))
  plotly::layout(p,annotations=list(list(text=caption,showarrow=FALSE,
    x=.5,y=0,xref="paper",yref="paper",yanchor="top")),margin=list(b=55))
}
