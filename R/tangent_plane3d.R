#' Tangent plane and normal vector to a surface z = f(x, y)
#'
#' @description
#' Computes the tangent plane to the graph of a scalar field \code{z = f(x, y)}
#' at a given point \code{(x0, y0)}, together with the associated normal
#' vector. Optionally, it displays the surface, the tangent plane, two
#' orthogonal cross-sections and the normal vector using \pkg{plotly}.
#'
#' @details
#' Given a differentiable function \code{f(x, y)} and a point
#' \code{(x0, y0)}, the function:
#' \itemize{
#'   \item approximates the partial derivatives \code{f_x(x0, y0)} and
#'         \code{f_y(x0, y0)} using centered finite differences with step
#'         \code{h},
#'   \item builds the tangent plane
#'         \code{g(x, y) = f(x0, y0) + f_x(x0, y0) * (x - x0) + f_y(x0, y0) * (y - y0)},
#'   \item constructs a normal vector \code{n = (-f_x, -f_y, 1)} and its
#'         unit version,
#'   \item encodes the plane in the form \code{a x + b y + c z + d = 0},
#'         with the coefficients returned in \code{plane_coeff}.
#' }
#'
#' When \code{plot = TRUE}, the function produces an interactive figure
#' containing:
#' \itemize{
#'   \item a patch of the original surface,
#'   \item a patch of the tangent plane centered at \code{(x0, y0)},
#'   \item two intersection curves of the surface along \code{x = x0} and
#'         \code{y = y0},
#'   \item the point of tangency and a segment in the direction of the
#'         normal vector.
#' }
#'
#' @param f Function of two variables \code{f(x, y)} returning a numeric
#'   scalar, representing the height \code{z}.
#' @param point Numeric vector of length 2 giving the point of tangency
#'   \code{c(x0, y0)}.
#' @param h Numeric step used in the centered finite-difference
#'   approximation of the partial derivatives. Must be strictly positive.
#' @param plot Logical; if \code{TRUE}, constructs a \pkg{plotly} figure
#'   with the surface, tangent plane, cross-sections and normal vector.
#' @param x_window Numeric half-width of the window in the \code{x}
#'   direction used to draw the surface patch.
#' @param y_window Numeric half-width of the window in the \code{y}
#'   direction used to draw the surface patch.
#' @param z_window Numeric half-height for the visible \code{z} range
#'   around \code{f(x0, y0)}.
#' @param grid Integer number of grid points used to sample the surface
#'   in each horizontal direction.
#' @param plane_window Numeric half-width of the square patch of the
#'   tangent plane drawn around \code{(x0, y0)}.
#' @param vec_N_factor Numeric scale factor applied to the unit normal
#'   vector when drawing the segment that represents the normal.
#' @param surface_opacity Numeric value between 0 and 1 controlling the
#'   opacity of the surface patch.
#' @param plane_opacity Numeric value between 0 and 1 controlling the
#'   opacity of the tangent-plane patch.
#' @param colors List with named entries that control colors in the plot:
#'   \itemize{
#'     \item \code{surface}: colorscale for the original surface,
#'     \item \code{plane}: colorscale for the tangent plane,
#'     \item \code{xcut}: color for the intersection curve along
#'           \code{y = y0},
#'     \item \code{ycut}: color for the intersection curve along
#'           \code{x = x0},
#'     \item \code{point}: color for the tangency point marker,
#'     \item \code{normal}: color for the normal vector segment.
#'   }
#' @param show_surface_grid Logical; if \code{TRUE}, draws a grid on the
#'   surface patch.
#' @param surface_grid_color,surface_grid_width Color and width of the
#'   surface grid lines.
#' @param show_axis_grid Logical; if \code{TRUE}, draws grid lines on the
#'   coordinate axes in the 3D scene.
#' @param axis_grid_color,axis_grid_width Color and width of the axis
#'   grid lines.
#' @param scene List with 3D scene options passed to \code{plotly::layout},
#'   typically including axis titles and \code{aspectmode}.
#' @param bg List with background colors for the figure, with fields
#'   \code{paper} and \code{plot}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{fx}, \code{fy}: approximations of the partial derivatives
#'         \code{f_x(x0, y0)} and \code{f_y(x0, y0)},
#'   \item \code{f0}: value \code{f(x0, y0)},
#'   \item \code{normal_unit}: unit normal vector at the point of tangency,
#'   \item \code{normal_raw}: non-normalized normal vector
#'         \code{(-fx, -fy, 1)},
#'   \item \code{plane_fun}: function \code{g(x, y)} for the tangent plane,
#'   \item \code{plane_coeff}: numeric vector \code{c(a, b, c, d)} such that
#'         \code{a x + b y + c z + d = 0} is the tangent-plane equation,
#'   \item \code{fig}: a \pkg{plotly} figure when \code{plot = TRUE},
#'         otherwise \code{NULL}.
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' f <- function(x, y) x^2 + y^2
#' tp <- tangent_plane3d(
#'   f,
#'   point = c(1, 1),
#'   plot = FALSE
#' )
#' tp$plane_coeff
#' \dontshow{\}}
#'
#' @export
tangent_plane3d <- function(
    f,
    point,
    h = 1e-4,
    plot = FALSE,
    x_window = 4, y_window = 4, z_window = 4,
    grid = 50,
    plane_window = 1,
    vec_N_factor = 1,
    surface_opacity = 0.85,
    plane_opacity   = 0.7,
    colors = list(
      surface = "Viridis",
      plane   = "Reds",
      xcut    = "black",
      ycut    = "black",
      point   = "blue",
      normal  = "red"
    ),
    show_surface_grid = FALSE,
    surface_grid_color = "rgba(0,0,0,0.35)",
    surface_grid_width = 1,
    show_axis_grid = TRUE,
    axis_grid_color = "#e5e5e5",
    axis_grid_width = 1,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  # Basic validation
  if (!is.function(f)) {
    stop("'f' must be a function of two arguments (x, y).", call. = FALSE)
  }
  if (!is.numeric(point) || length(point) != 2L || any(!is.finite(point))) {
    stop("'point' must be a finite numeric vector c(x0, y0).", call. = FALSE)
  }
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("'h' must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(grid) || length(grid) != 1L || !is.finite(grid) ||
      grid < 2 || abs(grid - round(grid)) > .Machine$double.eps^0.5) {
    stop("'grid' must be an integer >= 2.", call. = FALSE)
  }

  f <- match.fun(f)
  a <- point[1]; b <- point[2]

  # Local copies
  ff <- f
  hh <- h

  # Partial derivatives (central differences)
  d1x <- function(x, y) (ff(x + hh, y) - ff(x - hh, y)) / (2 * hh)
  d1y <- function(x, y) (ff(x, y + hh) - ff(x, y - hh)) / (2 * hh)

  f0 <- ff(a, b)
  fx <- d1x(a, b)
  fy <- d1y(a, b)

  # Tangent plane g(x,y) and normals
  g <- function(x, y) f0 + fx * (x - a) + fy * (y - b)
  n_raw  <- c(-fx, -fy, 1)
  n_unit <- n_raw / sqrt(sum(n_raw^2))

  # Coefficients a x + b y + c z + d = 0 (normal = n_raw, passes through (a,b,f0))
  dcoef <- fx * a + fy * b - f0
  plane_coeff <- c(a = -fx, b = -fy, c = 1, d = dcoef)

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need the 'plotly' package installed.", call. = FALSE)
    } else {
      # Grids for surface and plane
      xs <- seq(a - x_window, a + x_window, length.out = grid)
      ys <- seq(b - y_window, b + y_window, length.out = grid)
      Zf <- outer(ys, xs, function(yy, xx) ff(xx, yy))

      xg <- seq(a - plane_window, a + plane_window, length.out = max(20L, grid %/% 2L))
      yg <- seq(b - plane_window, b + plane_window, length.out = max(20L, grid %/% 2L))
      Zg <- outer(yg, xg, function(yy, xx) g(xx, yy))

      # Intersection curves x = a and y = b
      xcut <- seq(a - plane_window, a + plane_window, length.out = 200L)
      ycut <- seq(b - plane_window, b + plane_window, length.out = 200L)
      z_xcut <- vapply(xcut, function(x) ff(x, b), numeric(1))
      z_ycut <- vapply(ycut, function(y) ff(a, y), numeric(1))

      # Point and normal vector
      p0 <- c(a, b, f0)
      p1 <- p0 + vec_N_factor * n_unit

      contours_arg <- NULL
      if (isTRUE(show_surface_grid)) {
        contours_arg <- list(
          x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
          y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
          z = list(show = FALSE)
        )
      }

      plt <- plotly::plot_ly()
      plt <- plt |>
        plotly::add_surface(
          x = xs, y = ys, z = Zf,
          colorscale = colors$surface,
          showscale = FALSE, opacity = surface_opacity,
          contours = contours_arg
        ) |>
        plotly::add_surface(
          x = xg, y = yg, z = Zg,
          colorscale = colors$plane,
          showscale = FALSE, opacity = plane_opacity
        ) |>
        plotly::add_trace(
          x = xcut,
          y = rep(b, length(xcut)),
          z = z_xcut,
          type = "scatter3d",
          mode = "lines",
          line = list(color = colors$xcut, width = 2),
          showlegend = FALSE,
          hoverinfo = "none"
        ) |>
        plotly::add_trace(
          x = rep(a, length(ycut)),
          y = ycut,
          z = z_ycut,
          type = "scatter3d",
          mode = "lines",
          line = list(color = colors$ycut, width = 2),
          showlegend = FALSE,
          hoverinfo = "none"
        ) |>
        plotly::add_trace(
          x = p0[1],
          y = p0[2],
          z = p0[3],
          type = "scatter3d",
          mode = "markers",
          marker = list(color = colors$point, size = 4),
          showlegend = FALSE,
          hoverinfo = "none"
        ) |>
        plotly::add_trace(
          x = c(p0[1], p1[1]),
          y = c(p0[2], p1[2]),
          z = c(p0[3], p1[3]),
          type = "scatter3d",
          mode = "lines",
          line = list(color = colors$normal, width = 3),
          showlegend = FALSE,
          hoverinfo = "none"
        )

      # z range around f0
      zmin <- f0 - z_window
      zmax <- f0 + z_window

      # Scene adjustments without overwriting user settings
      scene_final <- scene
      if (is.null(scene_final$xaxis)) scene_final$xaxis <- list()
      if (is.null(scene_final$yaxis)) scene_final$yaxis <- list()
      if (is.null(scene_final$zaxis)) scene_final$zaxis <- list()

      scene_final$xaxis <- utils::modifyList(
        scene_final$xaxis,
        list(
          showgrid = show_axis_grid,
          gridcolor = axis_grid_color,
          gridwidth = axis_grid_width
        )
      )
      scene_final$yaxis <- utils::modifyList(
        scene_final$yaxis,
        list(
          showgrid = show_axis_grid,
          gridcolor = axis_grid_color,
          gridwidth = axis_grid_width
        )
      )
      scene_final$zaxis <- utils::modifyList(
        scene_final$zaxis,
        list(
          showgrid = show_axis_grid,
          gridcolor = axis_grid_color,
          gridwidth = axis_grid_width,
          range = c(zmin, zmax)
        )
      )

      plt <- plt |>
        plotly::layout(
          title = "Tangent plane and normal at (x0, y0)",
          scene = scene_final,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      fig <- plt
      print(plt)
    }
  }

  list(
    fx = fx,
    fy = fy,
    f0 = f0,
    normal_unit = n_unit,
    normal_raw  = n_raw,
    plane_fun   = g,
    plane_coeff = plane_coeff,
    fig = fig
  )
}
