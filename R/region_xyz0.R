#' Planar region \eqn{\{(x, y): a \leq x \leq b, H1(x) \leq y \leq H2(x)\}} drawn at height z0
#'
#' Builds the planar region
#' \deqn{\Omega = \{(x,y):\ a \le x \le b,\ H_1(x) \le y \le H_2(x)\}}
#' and renders it as a thin patch on the plane \eqn{z = z_0}. Optionally it
#' draws the boundary curves, partition lines along \eqn{x}, and a surface
#' patch.
#'
#' @param H1,H2 Functions \code{H_i(x)}; lower/upper y-boundaries.
#' @param a,b Numeric scalars; x-interval endpoints with \eqn{a < b}.
#' @param z0 Numeric scalar; height of the display plane \eqn{z = z_0}.
#' @param D Integer > 0; number of x-partitions (controls vertical "slices"
#'   and grid density).
#' @param plot Logical; if \code{TRUE}, draw the region with \pkg{plotly}.
#' @param n_curve Integer; number of samples used to trace each boundary
#'   curve \eqn{y = H_i(x)}.
#' @param show_surface Logical; if \code{TRUE}, draw a thin surface patch of
#'   the region at \eqn{z = z_0}.
#' @param surface_colorscale Character; Plotly colorscale for the surface
#'   patch (for example, \code{"Blues"}).
#' @param surface_opacity Numeric in \eqn{[0,1]}; opacity of the surface
#'   patch.
#' @param show_surface_grid Logical; show grid/contours on the surface patch.
#' @param surface_grid_color Character; color for surface grid lines.
#' @param surface_grid_width Numeric; width for surface grid lines.
#' @param boundary_line List; style for the boundary polylines (for example,
#'   \code{list(color = "blue", width = 2)}).
#' @param partition_line List; style for the partition lines at
#'   \eqn{x = a + k(b-a)/D}.
#' @param show_end_edges Logical; draw edges at \eqn{x = a} and \eqn{x = b}.
#' @param scene List; Plotly 3D scene options.
#' @param bg List with \code{paper} and \code{plot} background colors.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{data}}{List (or tibble, depending on the implementation)
#'     with the sampled curves and/or the grid used.}
#'   \item{\code{plot}}{A \pkg{plotly} object when \code{plot = TRUE},
#'     otherwise \code{NULL}.}
#' }
#'
#' @examples
#' H1 <- function(x) 0
#' H2 <- function(x) 1 - x
#' # Region under H2 and above H1 in [0,1], drawn at z = 0
#' # region_xyz0(H1, H2, a = 0, b = 1, z0 = 0,
#' #             D = 20, plot = TRUE, show_surface = TRUE)
#'
#' @export
region_xyz0 <- function(
    H1, H2,
    a, b,
    z0,
    D,
    plot = TRUE,
    n_curve = 800,
    show_surface = FALSE,
    surface_colorscale = "Blues",
    surface_opacity    = 0.30,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    boundary_line  = list(color = "blue", width = 2),
    partition_line = list(color = "blue", width = 1),
    show_end_edges = TRUE,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  # basic checks (region_xy will do more)
  if (!is.function(H1) || !is.function(H2)) {
    stop("'H1' and 'H2' must be functions of the form Hi(x).", call. = FALSE)
  }
  if (!is.numeric(a) || !is.numeric(b) ||
      length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || b <= a) {
    stop("'a' and 'b' must be finite numeric scalars with b > a.", call. = FALSE)
  }
  if (!is.numeric(z0) || length(z0) != 1L || !is.finite(z0)) {
    stop("'z0' must be a finite numeric scalar.", call. = FALSE)
  }
  D <- as.integer(D)
  if (!is.finite(D) || D <= 0L) {
    stop("'D' must be a positive integer.", call. = FALSE)
  }

  # Delegate to region_xy() in 3D mode over plane z = z0
  region_xy(
    H1 = H1, H2 = H2,
    a = a, b = b, D = D,
    plot = plot,
    n_curve = n_curve,
    fill = isTRUE(show_surface),
    fillcolor = "rgba(49,130,189,0.25)",  # kept for compatibility with 2D version
    boundary_line  = boundary_line,
    partition_line = partition_line,
    show_end_edges = show_end_edges,
    axis_equal = FALSE,  # not applicable in 3D
    as_3d = TRUE,
    plane_z = z0,
    n_u = 30,
    surface_colorscale = surface_colorscale,
    surface_opacity    = surface_opacity,
    show_surface_grid  = show_surface_grid,
    surface_grid_color = surface_grid_color,
    surface_grid_width = surface_grid_width,
    scene = scene,
    bg = bg
  )
}#' Planar region with y-bounds H1(x) and H2(x) drawn at height z0
#'
#' Builds the planar region defined by
#' \eqn{a \le x \le b} and \eqn{H_1(x) \le y \le H_2(x)}, and renders it as
#' a thin patch on the plane \eqn{z = z_0}. Optionally it draws the boundary
#' curves, partition lines along \eqn{x}, and a surface patch.
#'
#' @param H1,H2 Functions \code{H_i(x)}; lower and upper y-boundaries.
#' @param a,b Numeric scalars; x-interval endpoints with \eqn{a < b}.
#' @param z0 Numeric scalar; height of the display plane \eqn{z = z_0}.
#' @param D Integer > 0; number of x-partitions (controls vertical slices
#'   and grid density).
#' @param plot Logical; if \code{TRUE}, draw the region with \pkg{plotly}.
#' @param n_curve Integer; number of samples used to trace each boundary
#'   curve \eqn{y = H_i(x)}.
#' @param show_surface Logical; if \code{TRUE}, draw a thin surface patch of
#'   the region at \eqn{z = z_0}.
#' @param surface_colorscale Character; Plotly colorscale for the surface
#'   patch (for example, \code{"Blues"}).
#' @param surface_opacity Numeric in \eqn{[0,1]}; opacity of the surface
#'   patch.
#' @param show_surface_grid Logical; show grid/contours on the surface patch.
#' @param surface_grid_color Character; color for surface grid lines.
#' @param surface_grid_width Numeric; width for surface grid lines.
#' @param boundary_line List; style for the boundary polylines (for example,
#'   \code{list(color = "blue", width = 2)}).
#' @param partition_line List; style for the partition lines at
#'   \eqn{x = a + k (b - a) / D}.
#' @param show_end_edges Logical; draw edges at \eqn{x = a} and \eqn{x = b}.
#' @param scene List; Plotly 3D scene options.
#' @param bg List with \code{paper} and \code{plot} background colors.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{data}}{List (or tibble, depending on the implementation)
#'     with the sampled curves and/or the grid used.}
#'   \item{\code{plot}}{A \pkg{plotly} object when \code{plot = TRUE},
#'     otherwise \code{NULL}.}
#' }
#'
#' @examples
#' H1 <- function(x) 0
#' H2 <- function(x) 1 - x
#' # Region under H2 and above H1 in [0,1], drawn at z = 0
#' # region_xyz0(H1, H2, a = 0, b = 1, z0 = 0,
#' #             D = 20, plot = TRUE, show_surface = TRUE)
#'
#' @export
region_xyz0 <- function(
    H1, H2,
    a, b,
    z0,
    D,
    plot = TRUE,
    n_curve = 800,
    show_surface = FALSE,
    surface_colorscale = "Blues",
    surface_opacity    = 0.30,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    boundary_line  = list(color = "blue", width = 2),
    partition_line = list(color = "blue", width = 1),
    show_end_edges = TRUE,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  # Basic checks (xy_region performs more detailed checks internally)
  if (!is.function(H1) || !is.function(H2)) {
    stop("'H1' and 'H2' must be functions of the form Hi(x).", call. = FALSE)
  }
  if (!is.numeric(a) || !is.numeric(b) ||
      length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || b <= a) {
    stop("'a' and 'b' must be finite numeric scalars with b > a.", call. = FALSE)
  }
  if (!is.numeric(z0) || length(z0) != 1L || !is.finite(z0)) {
    stop("'z0' must be a finite numeric scalar.", call. = FALSE)
  }
  D <- as.integer(D)
  if (!is.finite(D) || D <= 0L) {
    stop("'D' must be a positive integer.", call. = FALSE)
  }

  # Delegate to xy_region() in 3D mode over plane z = z0
  xy_region(
    H1 = H1, H2 = H2,
    a = a, b = b, D = D,
    plot = plot,
    n_curve = n_curve,
    fill = isTRUE(show_surface),
    fillcolor = "rgba(49,130,189,0.25)",  # kept for compatibility with 2D version
    boundary_line  = boundary_line,
    partition_line = partition_line,
    show_end_edges = show_end_edges,
    axis_equal = FALSE,  # not applicable in 3D
    as_3d = TRUE,
    plane_z = z0,
    n_u = 30,
    surface_colorscale = surface_colorscale,
    surface_opacity    = surface_opacity,
    show_surface_grid  = show_surface_grid,
    surface_grid_color = surface_grid_color,
    surface_grid_width = surface_grid_width,
    scene = scene,
    bg = bg
  )
}
