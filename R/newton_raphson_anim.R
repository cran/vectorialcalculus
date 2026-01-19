#' Newton-Raphson root finding with tangent animation (Plotly)
#'
#' Builds a Plotly animation of the Newton-Raphson method for finding a root of
#' a real function. Each frame shows the tangent line at the current iterate
#' and how its x-intercept defines the next iterate.
#'
#' The Newton-Raphson update is
#' \deqn{
#'   x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}.
#' }
#' If a derivative function is not provided, the derivative is approximated
#' numerically by the central difference
#' \deqn{
#'   f'(x) \approx \frac{f(x+h)-f(x-h)}{2h}.
#' }
#'
#' @param f Function. A real-valued function f(x). Must accept a numeric vector
#'   and return a numeric vector of the same length.
#' @param x0 Numeric scalar. Initial guess.
#' @param df Optional function. Derivative f'(x). If NULL, a numerical derivative
#'   is used.
#' @param h Numeric scalar. Step size for numerical derivative (when df is NULL).
#' @param max_iter Integer. Maximum number of iterations.
#' @param tol Numeric scalar. Stopping tolerance based on |f(x_n)|.
#' @param xlim Numeric vector of length 2. Plot range for x. If NULL, it is chosen
#'   around the iterates.
#' @param n_curve Integer. Number of points used to draw the curve.
#' @param frame_ms Integer. Frame duration in milliseconds.
#' @param transition_ms Integer. Transition duration in milliseconds.
#' @param title Character. Plot title. If NULL, a default title is used.
#' @param safe_mode Logical. If TRUE, use calmer animation defaults intended to
#'   reduce flicker and visual stress.
#'
#' @return A list with components:
#' \describe{
#' \item{plot}{A plotly object (htmlwidget) with animation frames.}
#' \item{iterates}{Data frame with iterations (n, x, fx, dfx).}
#' \item{root}{Last iterate (approximate root).}
#' \item{converged}{Logical. TRUE if convergence was detected within max_iter.}
#' }
#'
#' @examples
#' \donttest{
#' library(plotly)
#'
#' f <- function(x) x^3 - 2*x - 5
#' out <- newton_raphson_anim(f, x0 = 2)
#' out$plot
#' out$root
#'
#' g <- function(x) cos(x) - x
#' newton_raphson_anim(g, x0 = 1)
#' }
#'
#' @export
newton_raphson_anim <- function(
    f,
    x0,
    df = NULL,
    h = 1e-4,
    max_iter = 10L,
    tol = 1e-8,
    xlim = NULL,
    n_curve = 600L,
    frame_ms = 600L,
    transition_ms = 400L,
    title = NULL,
    safe_mode = TRUE
) {

  if (!is.function(f)) stop("'f' must be a function.")
  if (!is.numeric(x0) || length(x0) != 1L) stop("'x0' must be a numeric scalar.")
  if (!is.null(df) && !is.function(df)) stop("'df' must be NULL or a function.")
  if (!is.numeric(h) || length(h) != 1L || h <= 0) stop("'h' must be a positive numeric scalar.")
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter < 1) stop("'max_iter' must be an integer >= 1.")
  if (!is.numeric(tol) || length(tol) != 1L || tol <= 0) stop("'tol' must be a positive numeric scalar.")

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Please install it.")
  }

  max_iter <- as.integer(max_iter)
  n_curve <- as.integer(n_curve)
  if (is.na(n_curve) || n_curve < 200L) stop("'n_curve' must be an integer >= 200.")

  if (isTRUE(safe_mode)) {
    frame_ms <- max(as.integer(frame_ms), 500L)
    transition_ms <- max(as.integer(transition_ms), 300L)
  }

  # Numerical derivative if df is NULL
  d_num <- function(x) (f(x + h) - f(x - h)) / (2*h)

  # Run Newton iterations
  xs <- numeric(max_iter + 1)
  fxs <- numeric(max_iter + 1)
  dfxs <- numeric(max_iter + 1)

  xs[1] <- x0
  fxs[1] <- f(x0)
  dfxs[1] <- if (is.null(df)) d_num(x0) else df(x0)

  converged <- FALSE
  n_stop <- 0L

  for (k in 1:max_iter) {
    dk <- dfxs[k]
    if (!is.finite(dk) || dk == 0) {
      # cannot continue
      n_stop <- k - 1L
      break
    }

    x_next <- xs[k] - fxs[k] / dk

    xs[k + 1] <- x_next
    fxs[k + 1] <- f(x_next)
    dfxs[k + 1] <- if (is.null(df)) d_num(x_next) else df(x_next)

    if (is.finite(fxs[k + 1]) && abs(fxs[k + 1]) <= tol) {
      converged <- TRUE
      n_stop <- k
      break
    }
    n_stop <- k
  }

  iter_df <- data.frame(
    n = 0:n_stop,
    x = xs[1:(n_stop + 1)],
    fx = fxs[1:(n_stop + 1)],
    dfx = dfxs[1:(n_stop + 1)]
  )

  # Choose xlim if NULL: cover iterates with margin
  if (is.null(xlim)) {
    xr <- range(iter_df$x, finite = TRUE)
    span <- xr[2] - xr[1]
    if (!is.finite(span) || span == 0) span <- 1
    xlim <- c(xr[1] - 0.8*span, xr[2] + 0.8*span)
  }
  if (!is.numeric(xlim) || length(xlim) != 2L || !(xlim[2] > xlim[1])) {
    stop("'xlim' must be a numeric vector of length 2 with xlim[2] > xlim[1].")
  }

  # Base curve
  x_grid <- seq(xlim[1], xlim[2], length.out = n_curve)
  y_curve <- f(x_grid)
  if (!is.numeric(y_curve) || length(y_curve) != length(x_grid)) {
    stop("'f' must return a numeric vector of the same length as its input.")
  }
  df_curve <- data.frame(x = x_grid, y = y_curve)
  df_axis <- data.frame(x = x_grid, y = rep(0, length(x_grid)))

  # Frame ids
  frame_id <- sprintf("it%02d", iter_df$n + 1L)

  # Tangent line per iteration + points
  df_tan <- lapply(seq_len(nrow(iter_df)), function(i) {
    xk <- iter_df$x[i]
    yk <- iter_df$fx[i]
    mk <- iter_df$dfx[i]

    # tangent line y = yk + mk (x - xk)
    yy <- yk + mk * (x_grid - xk)

    data.frame(
      x = c(x_grid, NA_real_),
      y = c(yy, NA_real_),
      frame = frame_id[i],
      txt = paste0(
        "n = ", iter_df$n[i],
        "<br>x_n = ", formatC(xk, digits = 10, format = "f"),
        "<br>f(x_n) = ", formatC(yk, digits = 6, format = "e"),
        "<br>f'(x_n) = ", formatC(mk, digits = 6, format = "e")
      )
    )
  })
  df_tan <- do.call(rbind, df_tan)

  # Points: (x_n, f(x_n)) and (x_{n+1}, 0)
  df_pts <- lapply(seq_len(nrow(iter_df)), function(i) {
    xk <- iter_df$x[i]
    yk <- iter_df$fx[i]

    if (i < nrow(iter_df)) {
      x_next <- iter_df$x[i + 1]
    } else {
      x_next <- xk
    }

    data.frame(
      x = c(xk, x_next),
      y = c(yk, 0),
      which = c("current", "next"),
      frame = frame_id[i]
    )
  })
  df_pts <- do.call(rbind, df_pts)

  # On-screen label per frame
  y_all <- c(y_curve, 0)
  y_top <- max(y_all, na.rm = TRUE)
  y_bot <- min(y_all, na.rm = TRUE)
  x_label <- xlim[1] + 0.03 * (xlim[2] - xlim[1])
  y_label <- y_top - 0.06 * (y_top - y_bot)

  df_label <- lapply(seq_len(nrow(iter_df)), function(i) {
    xk <- iter_df$x[i]
    yk <- iter_df$fx[i]
    mk <- iter_df$dfx[i]
    data.frame(
      x = x_label,
      y = y_label,
      frame = frame_id[i],
      txt = paste0(
        "n = ", iter_df$n[i],
        "    x_n = ", formatC(xk, digits = 8, format = "f"),
        "    f(x_n) = ", formatC(yk, digits = 3, format = "e"),
        "    f'(x_n) = ", formatC(mk, digits = 3, format = "e")
      )
    )
  })
  df_label <- do.call(rbind, df_label)

  if (is.null(title)) {
    title <- "Newton-Raphson method (tangent iteration)"
  }

  plot <- plotly::plot_ly() |>
    plotly::add_trace(
      data = df_curve,
      x = ~x, y = ~y,
      type = "scatter", mode = "lines",
      name = "y = f(x)",
      hoverinfo = "skip"
    ) |>
    plotly::add_trace(
      data = df_axis,
      x = ~x, y = ~y,
      type = "scatter", mode = "lines",
      name = "y = 0",
      hoverinfo = "skip"
    ) |>
    plotly::add_trace(
      data = df_tan,
      x = ~x, y = ~y,
      frame = ~frame,
      type = "scatter", mode = "lines",
      name = "tangent",
      text = ~txt,
      hoverinfo = "text"
    ) |>
    plotly::add_trace(
      data = df_pts,
      x = ~x, y = ~y,
      frame = ~frame,
      type = "scatter", mode = "markers",
      name = "points",
      text = ~which,
      hoverinfo = "text",
      marker = list(size = 9)
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
    plotly::animation_slider(currentvalue = list(prefix = "iter: "))

  list(plot = plot, iterates = iter_df, root = iter_df$x[nrow(iter_df)], converged = converged)
}
