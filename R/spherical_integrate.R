.vc_spherical_xyz <- function(theta,phi,rho) {
  s <- sin(phi)
  s[which(phi==0 | phi==pi)] <- 0
  cbind(x=rho*s*cos(theta),y=rho*s*sin(theta),z=rho*cos(phi))
}

#' Integrate a Cartesian field over a spherical-coordinate region
#' @inheritParams vc_triple_integrate
#' @param interval Azimuth theta bounds in radians, reversed/equal allowed;
#'   absolute span at most 2*pi.
#' @param polar_lower,polar_upper Functions of theta giving ordered polar
#'   angles phi in `[0,pi]`, measured from the positive z axis.
#' @param radial_lower,radial_upper Functions of one theta,phi vector giving
#'   nonnegative ordered radial bounds rho.
#' @return Full nested triple result of type spherical_integral. Stored engine
#'   coordinates are theta,phi,rho; f receives physical Cartesian x,y,z.
#' @details Maps `x=rho*sin(phi)*cos(theta)`, `y=rho*sin(phi)*sin(theta)`,
#'   `z=rho*cos(phi)`, using Jacobian `rho^2*sin(phi)`. Exact endpoint poles use
#'   zero sine. Reuses vc_triple_integrate unchanged, preserving all diagnostics
#'   and heuristic error indicators. They are not certified bounds or evidence
#'   of change-of-variables hypotheses, integrability or ordering between nodes.
#'   All sampled fields must be finite, including origin and poles even where
#'   the Jacobian vanishes; supply justified continuous extensions explicitly.
#'   Azimuth reversal changes sign. Negative radii, sampled polar angles outside
#'   `[0,pi]` and spans above one azimuth turn are rejected. Equal widths skip
#'   deeper functions. No automatic singularity removal or adaptive nested rule.
#' @examples
#' (function(v)1) |> vc_spherical_integrate(c(0,2*pi),function(t)0,function(t)pi,
#'   function(v)0,function(v)1,n_outer=2,n_middle=16,n_inner=2)
#' @export
vc_spherical_integrate <- function(f,interval,polar_lower,polar_upper,
                                   radial_lower,radial_upper,
                                   method=c("simpson","trapezoid"),
                                   n_outer=8L,n_middle=16L,n_inner=8L,
                                   abs_tol=1e-8,rel_tol=1e-8) {
  if(!all(vapply(list(f,polar_lower,polar_upper,radial_lower,radial_upper),
    is.function,logical(1)))) stop("Field and boundaries must be functions.")
  interval <- .vc_integral_interval(interval)
  if(abs(diff(interval))>2*pi) stop("Azimuth span must be at most 2*pi.")
  polar <- function(g)function(t) {
    phi <- .vc_value(g,t)
    if(phi<0 || phi>pi) stop("Polar angle must lie in [0,pi] at sampled sections.")
    phi
  }
  radial <- function(g)function(v) {
    rho <- .vc_field_value(g,v)
    if(rho<0) stop("Radius must be nonnegative at sampled sections.")
    rho
  }
  transformed <- function(v) {
    xyz <- as.numeric(.vc_spherical_xyz(v[1],v[2],v[3]))
    s <- if(v[2]==0 || v[2]==pi)0 else sin(v[2])
    jacobian <- .vc_finite(v[3]*(v[3]*s),"Spherical Jacobian")
    .vc_finite(.vc_field_value(f,xyz)*jacobian,"Weighted spherical field")
  }
  out <- vc_triple_integrate(transformed,interval,polar(polar_lower),polar(polar_upper),
    radial(radial_lower),radial(radial_upper),method=match.arg(method),
    n_outer=n_outer,n_middle=n_middle,n_inner=n_inner,abs_tol=abs_tol,rel_tol=rel_tol)
  out$type <- "spherical_integral"
  out$coordinates <- "theta,phi,rho; phi from +z; Cartesian field; Jacobian rho^2*sin(phi)"
  out$note <- paste(out$note,"Change-of-variables hypotheses are caller assumptions.")
  out
}

#' Display stored radial fibers of a spherical integration region
#' @param result Result of vc_spherical_integrate().
#' @param title Optional character title.
#' @return Plotly widget of stored radial segments in physical x,y,z coordinates.
#' @details No field/boundary reevaluation or quadrature. Radial fibers are not
#'   a watertight surface or a complete boundary mesh. Shell fibers start at
#'   the inner radius and preserve the empty core. Degenerate poles may overlap.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   r <- vc_spherical_integrate(function(v)1,c(0,pi/2),function(t)0,function(t)pi/2,
#'     function(v)1,function(v)2,n_outer=2,n_middle=2,n_inner=2)
#'   vc_plot_spherical_integral(r)
#' }
#' @export
vc_plot_spherical_integral <- function(result,title=NULL) {
  if(!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"spherical_integral")) stop("Invalid spherical result.")
  if(!is.null(title)&&(!is.character(title)||length(title)!=1||is.na(title))) stop("Invalid title.")
  if(!length(result$outer_nodes)) stop("No stored sections to display.")
  xyz <- lapply(seq_along(result$slices),function(i) {
    d <- result$slices[[i]]$sections
    if(!nrow(d)) return(matrix(numeric(),0,3))
    phi <- as.vector(rbind(d$outer,d$outer,NA_real_))
    rho <- as.vector(rbind(d$lower,d$upper,NA_real_))
    .vc_spherical_xyz(result$outer_nodes[i],phi,rho)
  })
  m <- do.call(rbind,xyz)
  if(!nrow(m)) stop("No radial fibers to display.")
  p <- plotly::plot_ly(x=m[,1],y=m[,2],z=m[,3],type="scatter3d",mode="lines",name="Stored radial fibers")
  plotly::layout(p,title=if(is.null(title))"Spherical region: Cartesian radial fibers" else title,
    scene=list(xaxis=list(title="x"),yaxis=list(title="y"),zaxis=list(title="z"),aspectmode="data"))
}
