#' Animate Riemann rectangles under a curve (2D)
#'
#' Builds an interactive Plotly animation of Riemann sums approximating the
#' area under a function on a closed interval.
#'
#' @param f Function. A real-valued function. It must accept a numeric vector
#'   and return a numeric vector of the same length.
#' @param a Numeric scalar. Left endpoint.
#' @param b Numeric scalar. Right endpoint. Must satisfy b > a.
#' @param n_vals Integer vector. Values of the number of subintervals used as
#'   animation frames. If NULL, a default increasing sequence is used.
#' @param method Character. Rule used for rectangle heights: "midpoint"
#'   (default), "left", or "right".
#' @param n_curve Integer. Number of points used to draw the base curve.
#' @param frame_ms Integer. Frame duration in milliseconds.
#' @param transition_ms Integer. Transition duration in milliseconds.
#' @param title Character. Plot title. If NULL, a default title is used.
#' @param show_sum Logical. If TRUE, show n and the value of the Riemann sum in
#'   hover text.
#' @param y0 Numeric scalar. Baseline for rectangles (default 0).
#'
#' @return A plotly object (htmlwidget) with animation frames.
#'
#' @examples
#' \donttest{
#' library(plotly)
#' f <- function(x) x^2
#' riemann_rectangles2d(f, 0, 1)
#' }
#'
#' @export
riemann_rectangles2d <- function(
    f,
    a,
    b,
    n_vals = NULL,
    method = c("midpoint", "left", "right"),
    n_curve = 400L,
    frame_ms = 900L,
    transition_ms = 0L,
    title = NULL,
    show_sum = TRUE,
    y0 = 0
) {

  if (!is.function(f)) stop("'f' must be a function.")
  if (!is.numeric(a) || length(a) != 1L) stop("'a' must be a numeric scalar.")
  if (!is.numeric(b) || length(b) != 1L) stop("'b' must be a numeric scalar.")
  if (!(b > a)) stop("'b' must be greater than 'a'.")
  if (!is.numeric(y0) || length(y0) != 1L) stop("'y0' must be a numeric scalar.")

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Please install it.")
  }

  method <- match.arg(method)

  if (is.null(n_vals)) n_vals <- c(4L, 8L, 16L, 32L, 64L)
  n_vals <- as.integer(n_vals)
  if (any(is.na(n_vals)) || any(n_vals <= 0L)) stop("'n_vals' must contain positive integers.")

  if (!is.numeric(n_curve) || length(n_curve) != 1L || n_curve < 10L) {
    stop("'n_curve' must be an integer >= 10.")
  }

  # Base curve
  x_curve <- seq(a, b, length.out = as.integer(n_curve))
  y_curve <- f(x_curve)
  if (!is.numeric(y_curve) || length(y_curve) != length(x_curve)) {
    stop("'f' must return a numeric vector of the same length as its input.")
  }
  df_curve <- data.frame(x = x_curve, y = y_curve)

  # One trace per frame: all rectangles concatenated, separated by NA
  make_rects_one_trace <- function(n) {
    n <- as.integer(n)

    breaks <- seq(a, b, length.out = n + 1L)
    dx <- diff(breaks)
    xL <- breaks[-length(breaks)]
    xR <- breaks[-1L]

    x_eval <- switch(
      method,
      left = xL,
      right = xR,
      midpoint = (xL + xR) / 2
    )

    h <- f(x_eval)
    if (!is.numeric(h) || length(h) != n) stop("Invalid output from 'f' while building rectangles.")

    Sn <- sum(h * dx)

    # Build rectangles (each rectangle is a closed polygon + NA separator)
    rects <- lapply(seq_len(n), function(i) {
      data.frame(
        x = c(xL[i], xR[i], xR[i], xL[i], xL[i], NA_real_),
        y = c(y0,    y0,    y0 + h[i], y0 + h[i], y0,    NA_real_),
        n = n,
        Sn = Sn
      )
    })

    df <- do.call(rbind, rects)

    if (isTRUE(show_sum)) {
      df$text <- paste0("n = ", n, "<br>S_n = ", formatC(Sn, digits = 8, format = "f"))
      df$text[is.na(df$x)] <- NA_character_
    } else {
      df$text <- "rectangle"
      df$text[is.na(df$x)] <- NA_character_
    }

    df
  }

  df_rects <- do.call(rbind, lapply(n_vals, make_rects_one_trace))

  if (is.null(title)) title <- "Riemann rectangles under y = f(x)"

  plotly::plot_ly() |>
    plotly::add_trace(
      data = df_curve,
      x = ~x, y = ~y,
      type = "scatter", mode = "lines",
      name = "y = f(x)"
    ) |>
    plotly::add_trace(
      data = df_rects,
      x = ~x, y = ~y,
      frame = ~n,
      type = "scatter", mode = "lines",
      fill = "toself",
      name = paste0("Rectangles (", method, ")"),
      showlegend = FALSE,
      text = ~text,
      hoverinfo = "text"
    ) |>
    plotly::layout(
      title = title,
      xaxis = list(title = "x"),
      yaxis = list(title = "y")
    ) |>
    plotly::animation_opts(
      frame = as.integer(frame_ms),
      transition = as.integer(transition_ms)
    )
}



