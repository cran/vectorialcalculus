#' Partial derivatives of z = f(x, y) at a point with 3D visualization
#'
#' Numerically approximates the partial derivatives
#' \eqn{f_x(x_0, y_0)} and \eqn{f_y(x_0, y_0)} of a scalar field
#' \eqn{z = f(x, y)} at a given point \eqn{(x_0, y_0)} using central
#' finite differences.
#'
#' Optionally, it builds a 3D \pkg{plotly} surface for \eqn{z = f(x, y)}
#' on a rectangular window around \eqn{(x_0, y_0)} and overlays:
#' \itemize{
#'   \item the intersection curve of the surface with the plane \eqn{y = y_0}
#'         and its tangent line given by \eqn{f_x(x_0, y_0)};
#'   \item the intersection curve with the plane \eqn{x = x_0} and its
#'         tangent line given by \eqn{f_y(x_0, y_0)};
#'   \item the base point \eqn{(x_0, y_0, f(x_0, y_0))}.
#' }
#'
#' @param f Function \code{function(x, y)} returning a numeric scalar \code{f(x, y)}.
#' @param x0,y0 Numeric scalars; coordinates of the point where the partial
#'   derivatives are evaluated.
#' @param h Numeric step for the central finite differences. If \code{NULL},
#'   a default value is chosen as \code{1e-4 * (1 + max(abs(x0), abs(y0)))}.
#' @param xlim Numeric length-2 vector \code{c(x_min, x_max)}. Range for the
#'   \eqn{x}-axis used to draw the surface. If \code{NULL}, a symmetric window
#'   around \code{x0} is used.
#' @param ylim Numeric length-2 vector \code{c(y_min, y_max)}. Range for the
#'   \eqn{y}-axis used to draw the surface. If \code{NULL}, a symmetric window
#'   around \code{y0} is used.
#' @param nx,ny Integer grid sizes (number of points) along \code{x} and \code{y}
#'   for the surface plot. Recommended values are at least 20.
#' @param plot Logical; if \code{TRUE}, builds and returns a \pkg{plotly}
#'   surface plot. If \code{FALSE}, only the numeric derivatives are returned.
#' @param scene List with \pkg{plotly} 3D scene options (axis titles, aspect
#'   mode, and so on) passed to \code{plotly::layout()}.
#' @param bg List with background colors for \pkg{plotly}, with components
#'   \code{paper} and \code{plot}.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{fx}: numeric scalar, approximation of \eqn{f_x(x_0, y_0)}.
#'   \item \code{fy}: numeric scalar, approximation of \eqn{f_y(x_0, y_0)}.
#'   \item \code{f0}: numeric scalar, \eqn{f(x_0, y_0)}.
#'   \item \code{fig}: a \pkg{plotly} object if \code{plot = TRUE} and
#'         \pkg{plotly} is available; otherwise \code{NULL}.
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' f <- function(x, y) x^2 + 3 * x * y - y^2
#' res <- partial_derivatives_surface(
#'   f,
#'   x0 = 1, y0 = -1,
#'   xlim = c(-1, 3),
#'   ylim = c(-3, 1),
#'   nx = 60, ny = 60
#' )
#' res$fx
#' res$fy
#' \dontshow{\}}
#'
#' @export
partial_derivatives_surface <- function(
    f,
    x0, y0,
    h    = NULL,
    xlim = NULL,
    ylim = NULL,
    nx   = 60L,
    ny   = 60L,
    plot = TRUE,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  if (!is.function(f)) {
    stop("'f' must be a function of the form f(x, y).", call. = FALSE)
  }
  stopifnot(
    length(x0) == 1L, length(y0) == 1L,
    is.numeric(x0), is.numeric(y0),
    is.numeric(nx), is.numeric(ny),
    nx >= 20L, ny >= 20L
  )

  # Step for finite differences
  if (is.null(h)) {
    h <- 1e-4 * (1 + max(abs(x0), abs(y0)))
  }
  h <- as.numeric(h)
  if (!is.finite(h) || h <= 0) {
    stop("'h' must be a positive finite numeric scalar.", call. = FALSE)
  }

  # Central finite differences at (x0, y0)
  f0 <- f(x0, y0)
  fx <- (f(x0 + h, y0) - f(x0 - h, y0)) / (2 * h)
  fy <- (f(x0, y0 + h) - f(x0, y0 - h)) / (2 * h)

  # If no plot requested, return numeric values only
  if (!isTRUE(plot)) {
    return(list(fx = fx, fy = fy, f0 = f0, fig = NULL))
  }

  # Plot requested: check for plotly
  if (!requireNamespace("plotly", quietly = TRUE)) {
    warning(
      "Package 'plotly' is required for plotting. Returning numeric values only.",
      call. = FALSE
    )
    return(list(fx = fx, fy = fy, f0 = f0, fig = NULL))
  }

  # Default ranges around (x0, y0) if not provided
  if (is.null(xlim)) {
    xlim <- c(x0 - 1, x0 + 1)
  }
  if (is.null(ylim)) {
    ylim <- c(y0 - 1, y0 + 1)
  }
  if (length(xlim) != 2L || length(ylim) != 2L ||
      !is.finite(xlim[1]) || !is.finite(xlim[2]) ||
      !is.finite(ylim[1]) || !is.finite(ylim[2]) ||
      xlim[2] <= xlim[1] || ylim[2] <= ylim[1]) {
    stop("'xlim' and 'ylim' must be numeric c(min, max) with max > min.", call. = FALSE)
  }

  # Grid for the surface
  xs <- seq(xlim[1], xlim[2], length.out = nx)
  ys <- seq(ylim[1], ylim[2], length.out = ny)

  f_vec <- function(xx, yy) f(xx, yy)
  f_vec_v <- Vectorize(f_vec)

  X <- matrix(rep(xs, each = ny), nrow = ny, ncol = nx)
  Y <- matrix(rep(ys, times = nx), nrow = ny, ncol = nx)
  Z <- matrix(f_vec_v(X, Y), nrow = ny, ncol = nx)

  # Slice and tangent in x-direction (y = y0)
  xs_line <- seq(xlim[1], xlim[2], length.out = max(200L, nx))
  z_slice_x <- vapply(xs_line, function(xx) f(xx, y0), numeric(1))
  z_tan_x   <- f0 + fx * (xs_line - x0)

  # Slice and tangent in y-direction (x = x0)
  ys_line <- seq(ylim[1], ylim[2], length.out = max(200L, ny))
  z_slice_y <- vapply(ys_line, function(yy) f(x0, yy), numeric(1))
  z_tan_y   <- f0 + fy * (ys_line - y0)

  # Build plotly figure
  p <- plotly::plot_ly()

  # Surface z = f(x, y)
  p <- plotly::add_surface(
    p,
    x = xs, y = ys, z = Z,
    showscale = FALSE,
    opacity = 0.85,
    colorscale = list(
      list(0, "#deebf7"),
      list(1, "#3182bd")
    ),
    name = "z = f(x, y)"
  )

  # Slice at y = y0
  p <- plotly::add_trace(
    p,
    x = xs_line,
    y = rep(y0, length(xs_line)),
    z = z_slice_x,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 3),
    name = "slice y = y0"
  )

  # Tangent line in x-direction at (x0, y0)
  p <- plotly::add_trace(
    p,
    x = xs_line,
    y = rep(y0, length(xs_line)),
    z = z_tan_x,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "red", width = 4),
    name = "tangent in x"
  )

  # Slice at x = x0
  p <- plotly::add_trace(
    p,
    x = rep(x0, length(ys_line)),
    y = ys_line,
    z = z_slice_y,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 3, dash = "dot"),
    name = "slice x = x0"
  )

  # Tangent line in y-direction at (x0, y0)
  p <- plotly::add_trace(
    p,
    x = rep(x0, length(ys_line)),
    y = ys_line,
    z = z_tan_y,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "orange", width = 4),
    name = "tangent in y"
  )

  # Base point
  p <- plotly::add_markers(
    p,
    x = x0, y = y0, z = f0,
    marker = list(size = 6, color = "green"),
    name = "point (x0, y0, f(x0, y0))"
  )

  p <- plotly::layout(
    p,
    title = sprintf(
      "Partial derivatives at (x0, y0) = (%.3f, %.3f): fx = %.4g, fy = %.4g",
      x0, y0, fx, fy
    ),
    scene = scene,
    paper_bgcolor = bg$paper,
    plot_bgcolor  = bg$plot
  )

  print(p)

  list(
    fx  = fx,
    fy  = fy,
    f0  = f0,
    fig = p
  )
}
