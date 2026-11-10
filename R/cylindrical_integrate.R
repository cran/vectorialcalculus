#' Integrate a Cartesian field over a cylindrical-coordinate solid
#' @inheritParams vc_triple_integrate
#' @param interval Angular bounds in radians, reversed/equal allowed; span
#'   at most 2*pi to avoid repeated angular coverage.
#' @param radial_lower,radial_upper Nonnegative ordered radius functions of theta.
#' @param height_lower,height_upper Functions of one vector in theta,r order,
#'   returning ordered finite vertical bounds.
#' @return Full nested triple result with type cylindrical_integral; stored
#'   engine coordinates are theta,r,z, while f receives physical x,y,z.
#' @details Uses x=r*cos(theta), y=r*sin(theta), z=z and Jacobian r, then
#'   reuses vc_triple_integrate unchanged. All nested diagnostics survive.
#'   Error indicators remain heuristic; caller justifies change-of-variables,
#'   integrability and bounds between samples. Negative sampled radii and
#'   spans exceeding one turn are rejected. Sampled fields must be finite,
#'   including at the axis; no automatic removable-singularity extension.
#'   Angular reversal changes sign. No spherical or adaptive nested rule.
#' @examples
#' (function(v)1) |> vc_cylindrical_integrate(c(0,2*pi),function(t)0,function(t)1,
#'   function(v)0,function(v)2,n_outer=2,n_middle=2,n_inner=2)
#' @export
vc_cylindrical_integrate <- function(f,interval,radial_lower,radial_upper,
                                     height_lower,height_upper,
                                     method=c("simpson","trapezoid"),
                                     n_outer=8L,n_middle=8L,n_inner=8L,
                                     abs_tol=1e-8,rel_tol=1e-8) {
  if(!all(vapply(list(f,radial_lower,radial_upper,height_lower,height_upper),
    is.function,logical(1)))) stop("Field and boundaries must be functions.")
  interval <- .vc_integral_interval(interval)
  if(abs(diff(interval))>2*pi) stop("Angular span must be at most 2*pi.")
  radial <- function(g)function(t) {
    r <- .vc_value(g,t)
    if(r<0) stop("Radius must be nonnegative at sampled sections.")
    r
  }
  transformed <- function(v) {
    t <- v[1]; r <- v[2]
    .vc_finite(.vc_field_value(f,c(r*cos(t),r*sin(t),v[3]))*r,"Weighted cylindrical field")
  }
  out <- vc_triple_integrate(transformed,interval,radial(radial_lower),radial(radial_upper),
    height_lower,height_upper,method=match.arg(method),n_outer=n_outer,
    n_middle=n_middle,n_inner=n_inner,abs_tol=abs_tol,rel_tol=rel_tol)
  out$type <- "cylindrical_integral"
  out$coordinates <- "theta,r,z; Cartesian field; Jacobian r"
  out$note <- paste(out$note,"Change-of-variables hypotheses are caller assumptions.")
  out
}

#' Display stored sections of a cylindrical integration region
#' @param result Result of vc_cylindrical_integrate().
#' @param view Vertical fibers in Cartesian space or polar base projection.
#' @param title Optional character title.
#' @return Plotly widget based only on stored sections, without new quadrature
#'   or field/boundary evaluations.
#' @details Fibers are not a watertight surface. The base envelope includes
#'   degenerate height sections and is not a certified projection of positive-
#'   volume support. Curved base boundaries are joined by sampled chords.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   r <- vc_cylindrical_integrate(function(v)1,c(0,pi),function(t)0,function(t)1,
#'     function(v)0,function(v)2,n_outer=2,n_middle=2,n_inner=2)
#'   vc_plot_cylindrical_integral(r)
#' }
#' @export
vc_plot_cylindrical_integral <- function(result,view=c("sections","projection"),title=NULL) {
  if(!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"cylindrical_integral")) stop("Invalid cylindrical result.")
  view <- match.arg(view)
  if(!is.null(title)&&(!is.character(title)||length(title)!=1||is.na(title))) stop("Invalid title.")
  if(!length(result$outer_nodes)) stop("No stored sections to display.")
  if(view=="projection") {
    base <- result$projection
    base$type <- "polar_integral"
    return(vc_plot_polar_integral(base,title))
  }
  xyz <- lapply(seq_along(result$slices),function(i) {
    d <- result$slices[[i]]$sections
    if(!nrow(d)) return(matrix(numeric(),0,3))
    t <- result$outer_nodes[i]
    x <- as.vector(rbind(d$outer*cos(t),d$outer*cos(t),NA_real_))
    y <- as.vector(rbind(d$outer*sin(t),d$outer*sin(t),NA_real_))
    z <- as.vector(rbind(d$lower,d$upper,NA_real_))
    cbind(x,y,z)
  })
  m <- do.call(rbind,xyz)
  if(!nrow(m)) stop("No inner fibers to display.")
  p <- plotly::plot_ly(x=m[,1],y=m[,2],z=m[,3],type="scatter3d",mode="lines",name="Stored vertical fibers")
  plotly::layout(p,title=if(is.null(title))"Cylindrical region: Cartesian fibers" else title,
    scene=list(xaxis=list(title="x"),yaxis=list(title="y"),zaxis=list(title="z"),aspectmode="data"))
}
