#' Sample spherical solids as stored parameter patches
#' @param polar_lower,polar_upper Ordered polar-angle functions of theta;
#'   sampled phi must lie in `[0,pi]`, measured from positive z.
#' @param radial_lower,radial_upper Ordered nonnegative radius functions of
#'   one theta,phi vector.
#' @param interval Increasing azimuth theta bounds in radians; span at most
#'   `2*pi`. NULL with partition supplied.
#' @param n,n_u,n_v Positive subdivisions in theta, polar blend and radius blend.
#' @param partition Optional vc_partition result for custom theta strips.
#' @return Plain list of type solid_spherical with parameter (unchanged
#'   vc_solid_cartesian result in theta,phi,rho coordinates), Cartesian faces
#'   theta_min,theta_max,polar_lower,polar_upper,radial_inner,radial_outer,
#'   full_turn metadata (exact numeric span equals `2*pi`) and note.
#' @details Reuses the Cartesian parameter sampler and accepted spherical
#'   coordinate map unchanged. No volume, area, normals or quadrature computed.
#'   Exact pole sines are zero. Angular cuts and pole/origin degeneracies remain
#'   stored; cuts are not necessarily exterior boundary patches. Endpoint
#'   periodicity and pole consistency are not inferred from the angular ranges.
#'   Sampled ordering/shared endpoints cannot certify between-node behavior,
#'   exact closure, outward orientation or a watertight manifold. Curved patches
#'   are interpolated approximations. Original radial surface examples use
#'   constant polar bounds; theta-dependent polar bounds are also supported.
#' @examples
#' vc_solid_spherical(function(t)0,function(t)pi,function(v)1,function(v)2,
#'   c(0,2*pi),n=12,n_u=8,n_v=3)
#' @export
vc_solid_spherical <- function(polar_lower,polar_upper,radial_lower,radial_upper,
                               interval=NULL,n=24L,n_u=16L,n_v=8L,partition=NULL) {
  if(!all(vapply(list(polar_lower,polar_upper,radial_lower,radial_upper),is.function,logical(1))))
    stop("All angular and radial boundaries must be functions.")
  if(!is.null(partition)) {
    if(!is.list(partition)||!identical(partition$type,"partition"))stop("Invalid theta partition.")
    p <- vc_partition(interval,breaks=partition$breaks)
  } else p <- vc_partition(interval,n)
  span <- diff(p$interval)
  if(span>2*pi)stop("Azimuth span must be at most 2*pi.")
  polar <- function(g)function(t) {
    phi <- .vc_value(g,t)
    if(phi<0||phi>pi)stop("Sampled polar angle must lie between zero and pi.")
    phi
  }
  radial <- function(g)function(v) {
    rho <- .vc_field_value(g,v)
    if(rho<0)stop("Sampled radius must be nonnegative.")
    rho
  }
  param <- vc_solid_cartesian(polar(polar_lower),polar(polar_upper),
    radial(radial_lower),radial(radial_upper),n_u=n_u,n_v=n_v,partition=p)
  faces <- lapply(param$faces,function(f) {
    xyz <- .vc_spherical_xyz(as.numeric(f$x),as.numeric(f$y),as.numeric(f$z))
    if(any(!is.finite(xyz)))stop("Mapped face coordinates are nonfinite.")
    list(x=matrix(xyz[,1],nrow=nrow(f$x)),y=matrix(xyz[,2],nrow=nrow(f$x)),
         z=matrix(xyz[,3],nrow=nrow(f$x)))
  })
  names(faces) <- c("theta_min","theta_max","polar_lower","polar_upper","radial_inner","radial_outer")
  list(type="solid_spherical",parameter=param,faces=faces,full_turn=span==2*pi,
    note="Sampled spherical parameter patches; seam, pole and origin degeneracies do not certify exterior boundary or closure.")
}

#' Display stored spherical shells and optional angular cuts
#' @param result Result of vc_solid_spherical().
#' @inheritParams vc_plot_solid_cartesian
#' @param faces Named selection; default inner and outer radial shells.
#'   Use all stored names explicitly to include angular cut patches.
#' @return Plotly widget using unchanged Cartesian face/wire rendering, without
#'   boundary reevaluation or integration.
#' @details For a full shell, angular cuts are normally omitted and pole patches
#'   collapse. A ball's inner radius-zero shell collapses to the origin; select
#'   radial_outer alone if desired. For a restricted sector, select appropriate
#'   angular cuts explicitly. No automatic exterior-boundary union is computed.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   s <- vc_solid_spherical(function(t)0,function(t)pi,function(v)1,function(v)2,
#'     c(0,2*pi),n=12,n_u=8,n_v=2)
#'   vc_plot_solid_spherical(s,mode="both")
#' }
#' @export
vc_plot_solid_spherical <- function(result,mode=c("faces","wireframe","both"),
    faces=c("radial_inner","radial_outer"),opacity=.4,wire_step=2L,title=NULL) {
  if(!is.list(result)||!identical(result$type,"solid_spherical"))stop("Invalid spherical solid.")
  if(is.null(title))title <- "Spherical solid: sampled shells and cuts"
  vc_plot_solid_cartesian(list(type="solid_cartesian",faces=result$faces),
    mode=match.arg(mode),faces=faces,opacity=opacity,wire_step=wire_step,title=title)
}
