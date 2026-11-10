#' Sample a cylindrical solid using stored parameter-box faces
#' @param radial_lower,radial_upper Nonnegative ordered radius functions of theta.
#' @param bottom,top Height functions of one theta,r vector, ordered at nodes.
#' @param interval Increasing azimuth bounds in radians, span at most `2*pi`;
#'   NULL with partition supplied. Geometry does not use reversed orientation.
#' @param n,n_u,n_v Positive subdivisions in theta, radial blend and height blend.
#' @param partition Optional vc_partition result for custom theta strips.
#' @return Plain list of type solid_cylindrical with parameter (full unchanged
#'   vc_solid_cartesian result in theta,r,z coordinates), physical Cartesian
#'   faces named theta_min, theta_max, radial_inner, radial_outer, z_lower,
#'   z_upper, full_turn flag (exact numeric span equals `2*pi`), and note.
#' @details Reuses vc_solid_cartesian unchanged, then maps stored theta,r,z to
#'   `x=r*cos(theta)`, `y=r*sin(theta)`, `z=z`. No boundary reevaluation is needed
#'   for this mapping; no area, normals, volume or quadrature are computed.
#'   Sampled boundary order does not certify between-node order or topology.
#'   At the axis patches can collapse. At a full turn angular cut patches can
#'   overlap and are not necessarily exterior boundary: all parameter patches
#'   remain stored and caller selects which to draw. Periodic boundary agreement
#'   is not inferred; nonperiodic endpoints can describe different seam sections.
#'   Coincident coordinates can differ by trigonometric roundoff. No exact
#'   closure, outward orientation or watertight manifold is certified.
#' @examples
#' vc_solid_cylindrical(function(t)0,function(t)1,function(v)0,function(v)2,
#'   c(0,pi/2),n=4,n_u=3,n_v=3)
#' @export
vc_solid_cylindrical <- function(radial_lower,radial_upper,bottom,top,
                                 interval=NULL,n=24L,n_u=12L,n_v=12L,partition=NULL) {
  if(!is.function(radial_lower)||!is.function(radial_upper))stop("Radius bounds must be functions.")
  if(!is.null(partition)) {
    if(!is.list(partition)||!identical(partition$type,"partition"))stop("Invalid theta partition.")
    p <- vc_partition(interval,breaks=partition$breaks)
  } else p <- vc_partition(interval,n)
  span <- diff(p$interval)
  if(span>2*pi)stop("Azimuth span must be at most 2*pi.")
  radial <- function(g)function(t) {
    r <- .vc_value(g,t)
    if(r<0)stop("Radius must be nonnegative at sampled nodes.")
    r
  }
  param <- vc_solid_cartesian(radial(radial_lower),radial(radial_upper),bottom,top,
    n_u=n_u,n_v=n_v,partition=p)
  faces <- lapply(param$faces,function(f) {
    x <- f$y*cos(f$x);y <- f$y*sin(f$x)
    if(any(!is.finite(c(x,y))))stop("Mapped face coordinates are nonfinite.")
    list(x=x,y=y,z=f$z)
  })
  names(faces) <- c("theta_min","theta_max","radial_inner","radial_outer","z_lower","z_upper")
  list(type="solid_cylindrical",parameter=param,faces=faces,full_turn=span==2*pi,
    note="Sampled parameter patches; angular seams/axis degeneracies do not certify exterior boundary or closure.")
}

#' Display stored cylindrical-solid patches or wireframes
#' @param result Result of vc_solid_cylindrical().
#' @inheritParams vc_plot_solid_cartesian
#' @param faces Nonempty unique selection of stored physical face names.
#'   Angular cut patches remain selectable even for a full turn.
#' @return Plotly widget using unchanged face/wire renderer and stored Cartesian
#'   coordinates only. No geometry sampling, derivative or integral evaluation.
#' @details For a periodic full cylinder one can omit theta_min/theta_max and
#'   a degenerate radial_inner axis patch. A periodic annulus normally uses
#'   radial_inner, radial_outer, z_lower and z_upper. Caller chooses patches;
#'   this is not automatic boundary-union or normal-orientation analysis.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   s <- vc_solid_cylindrical(function(t)1,function(t)2,function(v)0,function(v)1,
#'     c(0,2*pi),n=12,n_u=2,n_v=2)
#'   vc_plot_solid_cylindrical(s,faces=c("radial_inner","radial_outer","z_lower","z_upper"))
#' }
#' @export
vc_plot_solid_cylindrical <- function(result,mode=c("faces","wireframe","both"),
    faces=names(result$faces),opacity=.4,wire_step=2L,title=NULL) {
  if(!is.list(result)||!identical(result$type,"solid_cylindrical"))stop("Invalid cylindrical solid.")
  if(is.null(title))title <- "Cylindrical solid: sampled parameter patches"
  vc_plot_solid_cartesian(list(type="solid_cartesian",faces=result$faces),
    mode=match.arg(mode),faces=faces,opacity=opacity,wire_step=wire_step,title=title)
}
