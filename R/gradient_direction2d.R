#' Animate gradient and directional derivative on level curves (2D)
#'
#' Produces a Plotly animation showing level curves of a scalar field together
#' with the gradient direction at a point and a rotating unit direction vector.
#' The directional derivative value is displayed on screen for each frame.
#' A highlight is shown when the rotating direction aligns with the gradient,
#' which corresponds to the maximum directional derivative.
#'
#' The scalar field is
#' \deqn{
#'   z = f(x,y).
#' }
#' At the point
#' \deqn{
#'   (x_0,y_0),
#' }
#' the gradient vector is
#' \deqn{
#'   \nabla f(x_0,y_0) = \left(\frac{\partial f}{\partial x}(x_0,y_0),
#'   \frac{\partial f}{\partial y}(x_0,y_0)\right).
#' }
#' For a unit direction
#' \deqn{
#'   \mathbf{u}(\theta) = (\cos\theta,\sin\theta),
#' }
#' the directional derivative is
#' \deqn{
#'   D_{\mathbf{u}} f(x_0,y_0) = \nabla f(x_0,y_0)\cdot \mathbf{u}(\theta).
#' }
#' The maximum value over unit directions is
#' \deqn{
#'   \max_{\|\mathbf{u}\|=1} D_{\mathbf{u}} f(x_0,y_0) = \|\nabla f(x_0,y_0)\|,
#' }
#' and it occurs when
#' \deqn{
#'   \mathbf{u}(\theta)
#' }
#' points in the same direction as
#' \deqn{
#'   \nabla f(x_0,y_0).
#' }
#'
#' Partial derivatives are approximated numerically by central differences.
#'
#' @param f Function. A real-valued function f(x,y). It must accept two numeric
#'   arguments and return numeric values.
#' @param x0 Numeric scalar. x-coordinate of the base point.
#' @param y0 Numeric scalar. y-coordinate of the base point.
#' @param xlim Numeric vector of length 2. Range for x in the contour plot.
#' @param ylim Numeric vector of length 2. Range for y in the contour plot.
#' @param n_grid Integer. Grid size per axis for the contour computation.
#' @param theta_vals Numeric vector. Angles (radians) used as frames. If NULL,
#'   a default sequence from 0 to 2*pi is used.
#' @param h Numeric scalar. Step size for central differences.
#' @param arrow_scale Numeric scalar. Scale factor for drawing arrows. If NULL,
#'   an automatic scale based on the plot window is used.
#' @param frame_ms Integer. Frame duration in milliseconds.
#' @param transition_ms Integer. Transition duration in milliseconds.
#' @param title Character. Plot title. If NULL, a default title is used.
#' @param safe_mode Logical. If TRUE, use calmer animation defaults intended to
#'   reduce flicker and visual stress.
#' @param align_tol Numeric scalar. Angular tolerance (radians) used to decide
#'   when the rotating direction is considered aligned with the gradient.
#'
#' @return A plotly object (htmlwidget) with animation frames.
#'
#' @examples
#' \donttest{
#' library(plotly)
#'
#' f <- function(x, y) x^2 + 2*y^2
#' gradient_direction2d(
#'   f = f,
#'   x0 = 0.6,
#'   y0 = 0.4,
#'   xlim = c(-1.5, 1.5),
#'   ylim = c(-1.5, 1.5),
#'   safe_mode = TRUE,
#'   align_tol = 0.06
#' )
#' }
#'
#' @export
gradient_direction2d <- function(
    f,
    x0,
    y0,
    xlim,
    ylim,
    n_grid = 70L,
    theta_vals = NULL,
    h = 1e-4,
    arrow_scale = NULL,
    frame_ms = 220L,
    transition_ms = 220L,
    title = NULL,
    safe_mode = TRUE,
    align_tol = 0.08
) {

  if (!is.function(f)) stop("'f' must be a function.")
  if (!is.numeric(x0) || length(x0) != 1L) stop("'x0' must be a numeric scalar.")
  if (!is.numeric(y0) || length(y0) != 1L) stop("'y0' must be a numeric scalar.")
  if (!is.numeric(xlim) || length(xlim) != 2L) stop("'xlim' must be a numeric vector of length 2.")
  if (!is.numeric(ylim) || length(ylim) != 2L) stop("'ylim' must be a numeric vector of length 2.")
  if (!is.numeric(n_grid) || length(n_grid) != 1L || n_grid < 20L) stop("'n_grid' must be an integer >= 20.")
  if (!is.numeric(h) || length(h) != 1L || h <= 0) stop("'h' must be a positive numeric scalar.")
  if (!is.numeric(align_tol) || length(align_tol) != 1L || align_tol <= 0) stop("'align_tol' must be positive.")

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Please install it.")
  }

  if (is.null(theta_vals)) {
    theta_vals <- seq(0, 2*pi, length.out = 72)
  }
  if (!is.numeric(theta_vals) || length(theta_vals) < 2L) {
    stop("'theta_vals' must be a numeric vector with length >= 2.")
  }

  if (isTRUE(safe_mode)) {
    frame_ms <- max(as.integer(frame_ms), 200L)
    transition_ms <- max(as.integer(transition_ms), 200L)
  }

  # Numerical gradient at (x0, y0): central differences
  fx <- (f(x0 + h, y0) - f(x0 - h, y0)) / (2*h)
  fy <- (f(x0, y0 + h) - f(x0, y0 - h)) / (2*h)
  if (!is.numeric(fx) || !is.numeric(fy) || length(fx) != 1L || length(fy) != 1L) {
    stop("Unable to compute a numeric gradient at (x0, y0). Check 'f'.")
  }
  grad_norm <- sqrt(fx^2 + fy^2)

  # Unit gradient direction and its angle
  if (grad_norm == 0) {
    gx_unit <- 0
    gy_unit <- 0
    theta_g <- 0
  } else {
    gx_unit <- fx / grad_norm
    gy_unit <- fy / grad_norm
    theta_g <- atan2(gy_unit, gx_unit)
  }

  # Arrow scaling
  if (is.null(arrow_scale)) {
    arrow_scale <- 0.18 * min(abs(diff(xlim)), abs(diff(ylim)))
  }
  if (!is.numeric(arrow_scale) || length(arrow_scale) != 1L || arrow_scale <= 0) {
    stop("'arrow_scale' must be a positive numeric scalar.")
  }

  # Contour grid
  xs <- seq(xlim[1], xlim[2], length.out = as.integer(n_grid))
  ys <- seq(ylim[1], ylim[2], length.out = as.integer(n_grid))
  zmat <- outer(xs, ys, Vectorize(function(x, y) f(x, y)))
  if (!is.numeric(zmat)) stop("Contour computation failed: 'f' did not return numeric values on the grid.")

  # Stable frame ids
  frame_id <- sprintf("t%03d", seq_along(theta_vals))

  # Helper: shortest angular distance
  ang_dist <- function(a, b) {
    d <- a - b
    (d + pi) %% (2*pi) - pi
  }

  # Rotating direction vector data
  df_u <- lapply(seq_along(theta_vals), function(k) {
    th <- theta_vals[k]
    ux <- cos(th)
    uy <- sin(th)
    duf <- fx*ux + fy*uy
    data.frame(
      x = c(x0, x0 + arrow_scale * ux, NA_real_),
      y = c(y0, y0 + arrow_scale * uy, NA_real_),
      frame = frame_id[k],
      theta = th,
      duf = duf,
      hover = paste0(
        "theta = ", formatC(th, digits = 4, format = "f"),
        "<br>D_u f = ", formatC(duf, digits = 6, format = "f")
      )
    )
  })
  df_u <- do.call(rbind, df_u)

  # Highlight when aligned with gradient direction
  df_u_max <- lapply(seq_along(theta_vals), function(k) {
    th <- theta_vals[k]
    ux <- cos(th)
    uy <- sin(th)
    aligned <- (grad_norm > 0) && (abs(ang_dist(th, theta_g)) <= align_tol)
    if (isTRUE(aligned)) {
      data.frame(
        x = c(x0, x0 + arrow_scale * ux, NA_real_),
        y = c(y0, y0 + arrow_scale * uy, NA_real_),
        frame = frame_id[k]
      )
    } else {
      data.frame(
        x = c(NA_real_, NA_real_, NA_real_),
        y = c(NA_real_, NA_real_, NA_real_),
        frame = frame_id[k]
      )
    }
  })
  df_u_max <- do.call(rbind, df_u_max)

  # Fixed gradient arrow
  df_g <- data.frame(
    x = c(x0, x0 + arrow_scale * gx_unit, NA_real_),
    y = c(y0, y0 + arrow_scale * gy_unit, NA_real_)
  )

  # Base point
  df_p <- data.frame(x = x0, y = y0)

  # Always-visible label per frame
  x_label <- xlim[1] + 0.04 * (xlim[2] - xlim[1])
  y_label <- ylim[2] - 0.05 * (ylim[2] - ylim[1])

  df_label <- lapply(seq_along(theta_vals), function(k) {
    th <- theta_vals[k]
    ux <- cos(th)
    uy <- sin(th)
    duf <- fx*ux + fy*uy
    aligned <- (grad_norm > 0) && (abs(ang_dist(th, theta_g)) <= align_tol)
    status <- if (aligned) "MAX: u aligned with grad" else " "
    data.frame(
      x = x_label,
      y = y_label,
      frame = frame_id[k],
      txt = paste0(
        "D_u f(x0,y0) = ", formatC(duf, digits = 6, format = "f"),
        "    |grad f| = ", formatC(grad_norm, digits = 6, format = "f"),
        "<br>", status
      )
    )
  })
  df_label <- do.call(rbind, df_label)

  if (is.null(title)) title <- "Gradient and directional derivative (2D)"

  plotly::plot_ly() |>
    plotly::add_trace(
      x = xs, y = ys, z = zmat,
      type = "contour",
      showscale = FALSE,
      contours = list(coloring = "lines"),
      colorscale = "Greys",
      name = "Level curves",
      hoverinfo = "skip"
    ) |>
    plotly::add_trace(
      data = df_g,
      x = ~x, y = ~y,
      type = "scatter", mode = "lines",
      name = "grad f",
      hoverinfo = "skip"
    ) |>
    plotly::add_trace(
      data = df_p,
      x = ~x, y = ~y,
      type = "scatter", mode = "markers",
      name = "(x0,y0)",
      hoverinfo = "skip"
    ) |>
    plotly::add_trace(
      data = df_u,
      x = ~x, y = ~y,
      frame = ~frame,
      type = "scatter", mode = "lines",
      name = "u(theta)",
      text = ~hover,
      hoverinfo = "text"
    ) |>
    plotly::add_trace(
      data = df_u_max,
      x = ~x, y = ~y,
      frame = ~frame,
      type = "scatter", mode = "lines",
      name = "MAX alignment",
      line = list(width = 6),
      hoverinfo = "skip",
      showlegend = FALSE
    ) |>
    plotly::add_trace(
      data = df_label,
      x = ~x, y = ~y,
      frame = ~frame,
      type = "scatter", mode = "text",
      text = ~txt,
      textposition = "top left",
      textfont = list(size = 14, color = "black"),
      hoverinfo = "skip",
      showlegend = FALSE
    ) |>
    plotly::layout(
      title = title,
      xaxis = list(title = "x", range = xlim),
      yaxis = list(title = "y", range = ylim)
    ) |>
    plotly::animation_opts(
      frame = as.integer(frame_ms),
      transition = as.integer(transition_ms),
      easing = "linear",
      redraw = FALSE
    ) |>
    plotly::animation_button(
      x = 1, xanchor = "right",
      y = 1, yanchor = "top"
    ) |>
    plotly::animation_slider(
      currentvalue = list(prefix = "frame: ")
    )
}



