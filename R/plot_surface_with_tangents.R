#' Surface with tangent lines at a point
#'
#' @description
#' Draws the surface \eqn{z = f(x, y)} on a rectangular domain and overlays
#' two tangent line segments at a given point \eqn{(x_0, y_0)}:
#' one tangent in the direction of the \eqn{x}-axis and one tangent in the
#' direction of the \eqn{y}-axis. The partial derivatives are approximated
#' numerically by central finite differences.
#'
#' @param f Scalar field, given as \code{function(x, y)} returning a numeric
#'   value.
#' @param x0,y0 Numeric scalars with the coordinates of the point where
#'   the tangent lines are drawn.
#' @param xlim Numeric vector \code{c(x_min, x_max)} giving the range of the
#'   \eqn{x}-axis used to draw the surface.
#' @param ylim Numeric vector \code{c(y_min, y_max)} giving the range of the
#'   \eqn{y}-axis used to draw the surface.
#' @param n Integer number of grid points per axis used to discretize the
#'   surface. Must be at least 20.
#' @param h Numeric step used for the central finite–difference
#'   approximation of the partial derivatives \eqn{f_x} and \eqn{f_y}.
#' @param t_len Numeric scalar giving half the length of the tangent
#'   segments along the \eqn{x} and \eqn{y} directions.
#' @param title_prefix Optional character string used as a prefix in the
#'   plot title (for example, the name of the function \eqn{f}).
#'
#' @return
#' A \pkg{plotly} object representing the surface \eqn{z = f(x, y)} together
#' with the point \eqn{(x_0, y_0, f(x_0, y_0))} and the two tangent line
#' segments. The object can be further modified with usual \pkg{plotly}
#' tools.
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' f <- function(x, y) sin(x) * cos(y)
#' p <- plot_surface_with_tangents(
#'   f,
#'   x0 = 1, y0 = 1,
#'   xlim = c(-3, 3),
#'   ylim = c(-3, 3),
#'   n = 80
#' )
#' # p
#' \dontshow{\}}
#'
#' @export
plot_surface_with_tangents <- function(
    f,
    x0, y0,
    xlim = c(-3, 3),
    ylim = c(-3, 3),
    n = 120,
    h = 1e-5,
    t_len = 0.75,
    title_prefix = "f"
) {
  # basic checks
  if (!is.function(f)) {
    stop("'f' must be a function of the form f(x, y).", call. = FALSE)
  }
  if (!is.numeric(x0) || !is.numeric(y0) ||
      length(x0) != 1L || length(y0) != 1L ||
      !is.finite(x0) || !is.finite(y0)) {
    stop("'x0' and 'y0' must be finite numeric scalars.", call. = FALSE)
  }
  if (!is.numeric(xlim) || length(xlim) != 2L ||
      !all(is.finite(xlim)) || xlim[2] <= xlim[1]) {
    stop("'xlim' must be numeric c(x_min, x_max) with x_max > x_min.", call. = FALSE)
  }
  if (!is.numeric(ylim) || length(ylim) != 2L ||
      !all(is.finite(ylim)) || ylim[2] <= ylim[1]) {
    stop("'ylim' must be numeric c(y_min, y_max) with y_max > y_min.", call. = FALSE)
  }
  n <- as.integer(n)
  if (!is.finite(n) || n < 20L) {
    stop("'n' must be an integer greater or equal to 20.", call. = FALSE)
  }
  h <- as.numeric(h)
  if (!is.finite(h) || h <= 0) {
    stop("'h' must be a positive finite numeric scalar.", call. = FALSE)
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function.", call. = FALSE)
  }

  # grid and surface
  x <- seq(xlim[1], xlim[2], length.out = n)
  y <- seq(ylim[1], ylim[2], length.out = n)
  XY <- expand.grid(x = x, y = y)
  f_vec <- Vectorize(function(xx, yy) f(xx, yy))
  Z <- matrix(
    f_vec(XY$x, XY$y),
    nrow = length(y), ncol = length(x),
    byrow = FALSE
  )

  # point and partials (central differences at (x0, y0))
  z0  <- f(x0, y0)
  fx0 <- (f(x0 + h, y0) - f(x0 - h, y0)) / (2 * h)
  fy0 <- (f(x0, y0 + h) - f(x0, y0 - h)) / (2 * h)

  # tangent parameters (along x and y)
  t <- seq(-t_len, t_len, length.out = 50L)

  x_tan_x <- x0 + t
  y_tan_x <- rep(y0, length(t))
  z_tan_x <- z0 + fx0 * t

  y_tan_y <- y0 + t
  x_tan_y <- rep(x0, length(t))
  z_tan_y <- z0 + fy0 * t

  # build plotly figure
  p <- plotly::plot_ly() |>
    # surface
    plotly::add_surface(
      x = x, y = y, z = Z,
      name = "Surface",
      showscale = FALSE
    ) |>
    # base point
    plotly::add_markers(
      x = x0, y = y0, z = z0,
      name = "Point (x0, y0)",
      hoverinfo = "text",
      text = paste0(
        "(x0, y0, z0) = (",
        signif(x0, 6), ", ",
        signif(y0, 6), ", ",
        signif(z0, 6), ")"
      )
    ) |>
    # tangent along x
    plotly::add_trace(
      x = x_tan_x, y = y_tan_x, z = z_tan_x,
      type = "scatter3d", mode = "lines",
      name = "Tangent along x",
      hoverinfo = "text",
      text = paste0(
        "x = ", round(x_tan_x, 3),
        "<br>y = ", round(y_tan_x, 3),
        "<br>z = ", signif(z_tan_x, 6)
      )
    ) |>
    # tangent along y
    plotly::add_trace(
      x = x_tan_y, y = y_tan_y, z = z_tan_y,
      type = "scatter3d", mode = "lines",
      name = "Tangent along y",
      hoverinfo = "text",
      text = paste0(
        "x = ", round(x_tan_y, 3),
        "<br>y = ", round(y_tan_y, 3),
        "<br>z = ", signif(z_tan_y, 6)
      )
    ) |>
    plotly::layout(
      title = sprintf(
        "%s: surface and tangent lines at (%.3f, %.3f)",
        title_prefix, x0, y0
      ),
      scene = list(
        xaxis = list(title = "x"),
        yaxis = list(title = "y"),
        zaxis = list(title = "z")
      ),
      legend = list(orientation = "h", x = 0.05, y = 1.02)
    )

  p
}
