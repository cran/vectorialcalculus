#' Integrate a Cartesian scalar field over a polar region
#' @inheritParams vc_double_integrate
#' @param interval Angular bounds in radians, reversed or equal allowed;
#'   absolute span at most 2*pi to prevent repeated angular coverage.
#' @param lower,upper Radius functions of angle, nonnegative and ordered
#'   at sampled sections.
#' @return Full nested result with type polar_integral, parameter sections
#'   (outer is theta, lower/upper are radii), and nested engine diagnostics.
#' @details Applies x=r*cos(theta), y=r*sin(theta) and the Jacobian r, then
#'   reuses vc_double_integrate unchanged. f always receives Cartesian x,y.
#'   Requires finite field samples, including the origin when sampled; no
#'   automatic removable-singularity extension. Negative radii and spans above
#'   one revolution are rejected. Caller justifies boundary behavior between
#'   samples and change-of-variables/integrability hypotheses. Angular reversal
#'   changes the signed result. Error indicators are heuristic, not bounds.
#' @examples
#' (function(v)1) |> vc_polar_integrate(c(0,2*pi),function(t)0,function(t)1)
#' @export
vc_polar_integrate <- function(f,interval,lower,upper,
                               method=c("simpson","trapezoid"),
                               n_outer=16L,n_inner=16L,abs_tol=1e-8,rel_tol=1e-8) {
  if(!is.function(f)||!is.function(lower)||!is.function(upper))
    stop("f, lower and upper must be functions.")
  interval <- .vc_integral_interval(interval)
  if(abs(diff(interval))>2*pi) stop("Angular span must be at most 2*pi.")
  radial <- function(g)function(t) {
    r <- .vc_value(g,t)
    if(r<0) stop("Radius must be nonnegative at sampled sections.")
    r
  }
  transformed <- function(v) {
    t <- v[1]
    r <- v[2]
    .vc_finite(.vc_field_value(f,c(r*cos(t),r*sin(t)))*r,"Weighted polar field")
  }
  out <- vc_double_integrate(transformed,interval,radial(lower),radial(upper),
    method=match.arg(method),n_outer=n_outer,n_inner=n_inner,
    abs_tol=abs_tol,rel_tol=rel_tol)
  out$type <- "polar_integral"
  out$coordinates <- "theta,r; field evaluated in Cartesian x,y; Jacobian r"
  out$note <- paste(out$note,"Change-of-variables hypotheses are caller assumptions.")
  out
}

#' Display a stored polar integration domain in Cartesian coordinates
#' @param result Result of vc_polar_integrate().
#' @param title Optional character title.
#' @return Plotly widget with sampled inner/outer boundary and radial sections.
#' @details Only stored sections are transformed. Curved boundaries are joined
#'   by straight chords. No quadrature or field/boundary evaluation is repeated.
#'   A full annulus is outlined by outer and reversed inner boundaries joined
#'   along the angular seam. This sampled domain is not a signed volume plot.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   r <- vc_polar_integrate(function(v)1,c(0,pi),function(t)1,function(t)2)
#'   vc_plot_polar_integral(r)
#' }
#' @export
vc_plot_polar_integral <- function(result,title=NULL) {
  if(!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"polar_integral")) stop("Invalid polar result.")
  if(!is.null(title)&&(!is.character(title)||length(title)!=1||is.na(title))) stop("Invalid title.")
  d <- result$sections
  if(!is.data.frame(d)||!nrow(d)) stop("No stored sections to display.")
  d <- d[order(d$outer),,drop=FALSE]
  t <- c(d$outer,rev(d$outer),d$outer[1])
  r <- c(d$upper,rev(d$lower),d$upper[1])
  p <- plotly::plot_ly(x=r*cos(t),y=r*sin(t),type="scatter",mode="lines",
    fill="toself",name="Sampled polar domain")
  t <- as.vector(rbind(d$outer,d$outer,NA_real_))
  r <- as.vector(rbind(d$lower,d$upper,NA_real_))
  p <- plotly::add_trace(p,x=r*cos(t),y=r*sin(t),type="scatter",mode="lines",
    name="Stored radial sections",inherit=FALSE)
  plotly::layout(p,title=if(is.null(title))"Polar domain in Cartesian coordinates" else title,
    xaxis=list(title="x"),yaxis=list(title="y",scaleanchor="x",scaleratio=1))
}
