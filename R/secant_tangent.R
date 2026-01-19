#' Secant lines converge to the tangent line (Plotly)
#'
#' Approximates the derivative of a function at a point numerically and builds
#' an interactive Plotly animation showing how secant (incremental quotient)
#' lines converge to the tangent line as the step size decreases. The secant
#' point(s) used for the slope computation are also animated.
#'
#' The forward incremental quotient is
#' \deqn{
#'   m_h = \frac{f(x_0+h)-f(x_0)}{h}.
#' }
#' The central difference approximation is
#' \deqn{
#'   m_h = \frac{f(x_0+h)-f(x_0-h)}{2h}.
#' }
#' The tangent line model at x0 is
#' \deqn{
#'   y = f(x_0) + f'(x_0)\,(x-x_0).
#' }
#'
#' @param f Function. A real-valued function f(x). It must accept a numeric
#'   vector and return a numeric vector of the same length.
#' @param x0 Numeric scalar. Point where the derivative is approximated.
#' @param h_vals Numeric vector. Positive step sizes used as animation frames.
#'   If NULL, a default decreasing sequence is used.
#' @param method Character. Derivative approximation method:
#'   "forward" (default) or "central".
#' @param xlim Numeric vector of length 2. Plot range for x. If NULL, it is
#'   chosen automatically from x0 and h_vals.
#' @param n_curve Integer. Number of points used to draw the curve and lines.
#' @param frame_ms Integer. Frame duration in milliseconds.
#' @param transition_ms Integer. Transition duration in milliseconds.
#' @param title Character. Plot title. If NULL, a default title is used.
#' @param safe_mode Logical. If TRUE, use calmer animation defaults intended to
#'   reduce flicker and visual stress.
#'
#' @return A list with components:
#' \describe{
#' \item{plot}{A plotly object (htmlwidget) with animation frames.}
#' \item{derivative}{Numeric scalar. Derivative estimate using the smallest h.}
#' \item{data}{Data frame used for the animated secant lines (useful for debugging).}
#' }
#'
#' @examples
#' \donttest{
#' library(plotly)
#'
#' f <- function(x) x^2
#' out <- secant_tangent(f, x0 = 1)
#' out$plot
#' out$derivative
#'
#' g <- function(x) sin(x)
#' secant_tangent(g, x0 = 0.7, method = "central", h_vals = 2^(-(1:7)))
#' }
#'
#' @export
secant_tangent <- function(
    f,
    x0,
    h_vals = NULL,
    method = c("forward", "central"),
    xlim = NULL,
    n_curve = 400L,
    frame_ms = 220L,
    transition_ms = 220L,
    title = NULL,
    safe_mode = TRUE
) {

  if (!is.function(f)) stop("'f' must be a function.")
  if (!is.numeric(x0) || length(x0) != 1L) stop("'x0' must be a numeric scalar.")

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Please install it.")
  }

  method <- match.arg(method)

  if (is.null(h_vals)) {
    h_vals <- 2^(-(1:8))  # 1/2, 1/4, ..., 1/256
  }
  h_vals <- as.numeric(h_vals)
  if (any(!is.finite(h_vals)) || any(h_vals <= 0)) stop("'h_vals' must be positive finite numbers.")
  h_vals <- sort(unique(h_vals), decreasing = TRUE)

  n_curve <- as.integer(n_curve)
  if (is.na(n_curve) || n_curve < 80L) stop("'n_curve' must be an integer >= 80.")

  if (isTRUE(safe_mode)) {
    frame_ms <- max(as.integer(frame_ms), 200L)
    transition_ms <- max(as.integer(transition_ms), 200L)
  }

  f0 <- f(x0)
  if (!is.numeric(f0) || length(f0) != 1L || !is.finite(f0)) {
    stop("'f(x0)' must be a finite numeric scalar.")
  }

  # Choose x-limits automatically if not provided
  if (is.null(xlim)) {
    h_max <- max(h_vals)
    xlim <- c(x0 - 4*h_max, x0 + 4*h_max)
  }
  if (!is.numeric(xlim) || length(xlim) != 2L || !(xlim[2] > xlim[1])) {
    stop("'xlim' must be a numeric vector of length 2 with xlim[2] > xlim[1].")
  }

  # Curve data
  x_grid <- seq(xlim[1], xlim[2], length.out = n_curve)
  y_curve <- f(x_grid)
  if (!is.numeric(y_curve) || length(y_curve) != length(x_grid)) {
    stop("'f' must return a numeric vector of the same length as its input.")
  }
  df_curve <- data.frame(x = x_grid, y = y_curve)

  # Slope function
  slope_for_h <- function(h) {
    if (method == "forward") {
      (f(x0 + h) - f0) / h
    } else {
      (f(x0 + h) - f(x0 - h)) / (2*h)
    }
  }

  slopes <- vapply(h_vals, slope_for_h, numeric(1))
  if (any(!is.finite(slopes))) stop("Some slope approximations are not finite. Check 'f' near x0.")

  # Derivative estimate: use smallest h (last after decreasing sort)
  derivative_hat <- slopes[length(slopes)]

  # Tangent line (static) using derivative_hat
  df_tangent <- data.frame(
    x = x_grid,
    y = f0 + derivative_hat * (x_grid - x0)
  )

  # Stable frame ids
  frame_id <- sprintf("h%03d", seq_along(h_vals))

  # Animated secant lines: one trace per frame (x,y with NA separator)
  df_secant <- lapply(seq_along(h_vals), function(k) {
    h <- h_vals[k]
    m <- slopes[k]
    data.frame(
      x = c(x_grid, NA_real_),
      y = c(f0 + m * (x_grid - x0), NA_real_),
      frame = frame_id[k],
      h = h,
      slope = m,
      txt = paste0(
        "h = ", formatC(h, digits = 6, format = "g"),
        "<br>m_h = ", formatC(m, digits = 8, format = "f"),
        "<br>tangent slope = ", formatC(derivative_hat, digits = 8, format = "f")
      )
    )
  })
  df_secant <- do.call(rbind, df_secant)

  # Animated secant point(s)
  df_pts <- lapply(seq_along(h_vals), function(k) {
    h <- h_vals[k]
    if (method == "forward") {
      x1 <- x0 + h
      y1 <- f(x1)
      data.frame(
        x = c(x0, x1),
        y = c(f0, y1),
        which = c("x0", "x0+h"),
        frame = frame_id[k]
      )
    } else {
      x1 <- x0 - h
      x2 <- x0 + h
      y1 <- f(x1)
      y2 <- f(x2)
      data.frame(
        x = c(x0, x1, x2),
        y = c(f0, y1, y2),
        which = c("x0", "x0-h", "x0+h"),
        frame = frame_id[k]
      )
    }
  })
  df_pts <- do.call(rbind, df_pts)

  # Always-visible label per frame (safe location)
  y_all <- c(y_curve, df_tangent$y)
  y_top <- max(y_all, na.rm = TRUE)
  y_bot <- min(y_all, na.rm = TRUE)
  x_label <- xlim[1] + 0.03 * (xlim[2] - xlim[1])
  y_label <- y_top - 0.06 * (y_top - y_bot)

  df_label <- lapply(seq_along(h_vals), function(k) {
    h <- h_vals[k]
    m <- slopes[k]
    data.frame(
      x = x_label,
      y = y_label,
      frame = frame_id[k],
      txt = paste0(
        "h = ", formatC(h, digits = 6, format = "g"),
        "    m_h = ", formatC(m, digits = 6, format = "f"),
        "    f'(x0) approx = ", formatC(derivative_hat, digits = 6, format = "f")
      )
    )
  })
  df_label <- do.call(rbind, df_label)

  if (is.null(title)) {
    title <- "Secant lines converge to the tangent line"
  }

  plot <- plotly::plot_ly() |>
    # Base curve
    plotly::add_trace(
      data = df_curve,
      x = ~x, y = ~y,
      type = "scatter", mode = "lines",
      name = "y = f(x)",
      hoverinfo = "skip"
    ) |>
    # Static tangent line
    plotly::add_trace(
      data = df_tangent,
      x = ~x, y = ~y,
      type = "scatter", mode = "lines",
      name = "tangent (approx)",
      hoverinfo = "skip",
      line = list(width = 4)
    ) |>
    # Animated secant line
    plotly::add_trace(
      data = df_secant,
      x = ~x, y = ~y,
      frame = ~frame,
      type = "scatter", mode = "lines",
      name = "secant (m_h)",
      text = ~txt,
      hoverinfo = "text"
    ) |>
    # Animated points used to build the secant slope
    plotly::add_trace(
      data = df_pts,
      x = ~x, y = ~y,
      frame = ~frame,
      type = "scatter", mode = "markers",
      name = "secant points",
      text = ~which,
      hoverinfo = "text",
      marker = list(size = 9)
    ) |>
    # On-screen label
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
      yaxis = list(title = "y")
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
    plotly::animation_slider(currentvalue = list(prefix = "frame: "))

  list(plot = plot, derivative = derivative_hat, data = df_secant)
}

