#' 1D Riemann sums with optional plot
#'
#' Computes lower, upper, and midpoint Riemann sums for a scalar function
#' \code{f(x)} on an interval \code{[xmin, xmax]}. Optionally draws a
#' 2D plot with rectangles and, if requested, the true curve.
#'
#' @param f Function \code{function(x)} returning numeric values.
#' @param xlim Numeric vector \code{c(xmin, xmax)} with \code{xmax > xmin}.
#' @param n Integer. Number of subintervals.
#' @param methods Character vector with any of \code{"lower"}, \code{"upper"},
#'   \code{"mid"} indicating which rectangle types to draw.
#' @param show_curve Logical. If \code{TRUE}, overlays the curve \code{f(x)}.
#' @param curve_res Integer. Number of points used to draw the curve.
#' @param colors Named list specifying fill colors for
#'   \code{list(lower=..., upper=..., mid=...)}.
#' @param alpha Numeric in \code{[0,1]}. Fill opacity for rectangles.
#' @param edge_color Color for rectangle borders.
#' @param edge_width Border width.
#' @param curve_color Color for the curve.
#' @param curve_width Line width for the curve.
#' @param show_baseline Logical. If \code{TRUE}, draws a horizontal baseline.
#' @param baseline Numeric. Y-value for the baseline.
#' @param baseline_color Baseline color.
#' @param baseline_width Baseline width.
#' @param warn_heavy Logical. If \code{TRUE}, warns when \code{n} is very large.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{lower_sum} Lower Riemann sum.
#'     \item \code{upper_sum} Upper Riemann sum.
#'     \item \code{mid_sum} Midpoint Riemann sum.
#'     \item \code{dx} Subinterval width.
#'     \item \code{x_breaks} Partition points.
#'     \item \code{figure} A plotly object, or \code{NULL} if not available.
#'   }
#'
#' @examples
#' f <- function(x) sin(2*x)
#' out <- riemann_sum_1d_plot(
#'   f, xlim = c(0, pi), n = 10,
#'   methods = c("lower","upper","mid"),
#'   show_curve = TRUE
#' )
#' out$mid_sum
#'
#' @export
riemann_sum_1d_plot <- function(
    f,
    xlim,
    n = 12L,
    methods = c("lower","upper","mid"),
    show_curve = TRUE,
    curve_res = 400L,
    colors = list(lower = "#a1d99b", upper = "#fc9272", mid = "#9ecae1"),
    alpha = 0.8,
    edge_color = "black",
    edge_width = 1.2,
    curve_color = "black",
    curve_width = 2,
    show_baseline = TRUE,
    baseline = 0,
    baseline_color = "gray50",
    baseline_width = 1,
    warn_heavy = TRUE
){
  stopifnot(is.function(f))

  if (length(xlim) != 2L || xlim[2] <= xlim[1]) {
    stop("'xlim' must be c(min, max) with max > min.", call. = FALSE)
  }

  n <- as.integer(n)
  curve_res <- as.integer(curve_res)

  if (!is.finite(n) || n < 1L) {
    stop("'n' must be a positive integer.", call. = FALSE)
  }

  methods <- intersect(tolower(methods), c("lower","upper","mid"))
  if (length(methods) == 0L) {
    stop("Choose at least one of methods = c('lower','upper','mid').", call. = FALSE)
  }

  have_plotly <- requireNamespace("plotly", quietly = TRUE)
  if (!have_plotly) {
    warning("Package 'plotly' not found; returning estimates only.", call. = FALSE)
  }
  if (isTRUE(warn_heavy) && n > 800L) {
    warning("High n: many polygons to draw; consider reducing 'n'.", call. = FALSE)
  }

  # partition
  x_breaks <- seq(xlim[1], xlim[2], length.out = n + 1L)
  dx <- diff(x_breaks)[1L]

  # endpoint & midpoint samples
  f_left  <- f(x_breaks[-(n + 1L)])
  f_right <- f(x_breaks[-1L])
  f_mid   <- f(0.5 * (x_breaks[-(n + 1L)] + x_breaks[-1L]))

  # per-interval heights
  h_lower <- pmin(f_left, f_right)
  h_upper <- pmax(f_left, f_right)
  h_mid   <- f_mid

  # sums
  lower_sum <- sum(h_lower) * dx
  upper_sum <- sum(h_upper) * dx
  mid_sum   <- sum(h_mid)   * dx

  # --- Plot
  fig <- NULL
  if (have_plotly) {
    add_alpha <- function(col, a) {
      if (grepl("^rgba\\(", col)) return(col)
      rgb <- grDevices::col2rgb(col)
      sprintf(
        "rgba(%d,%d,%d,%.3f)",
        rgb[1], rgb[2], rgb[3],
        max(0, min(1, a))
      )
    }

    p <- plotly::plot_ly()

    # baseline
    if (isTRUE(show_baseline)) {
      p <- plotly::add_trace(
        p, x = xlim, y = c(baseline, baseline),
        type = "scatter", mode = "lines",
        line = list(color = baseline_color, width = baseline_width),
        hoverinfo = "none", showlegend = FALSE
      )
    }

    # true curve
    if (isTRUE(show_curve)) {
      xs <- seq(xlim[1], xlim[2], length.out = max(50L, curve_res))
      ys <- f(xs)
      p <- plotly::add_trace(
        p, x = xs, y = ys,
        type = "scatter", mode = "lines",
        line = list(color = curve_color, width = curve_width),
        name = "f(x)"
      )
    }

    # helper to add all rectangles of one method as a single polygons trace
    add_rects <- function(p, heights, method_name, color_hex) {
      xpoly <- ypoly <- numeric(0L)
      for (i in seq_len(n)) {
        x0 <- x_breaks[i]
        x1 <- x_breaks[i + 1L]
        h  <- heights[i]
        # polygon: (x0,baseline) -> (x1,baseline) -> (x1,h) -> (x0,h) -> close
        xpoly <- c(xpoly, x0, x1, x1, x0, NA_real_)
        ypoly <- c(ypoly, baseline, baseline, h,  h,  NA_real_)
      }

      method_name <- tolower(method_name)

      plotly::add_trace(
        p,
        x = xpoly, y = ypoly,
        type = "scatter", mode = "lines",
        fill = "toself",
        fillcolor = add_alpha(color_hex, alpha),
        line = list(color = edge_color, width = edge_width),
        name = switch(
          method_name,
          lower = sprintf("Lower (sum = %.6g)", lower_sum),
          upper = sprintf("Upper (sum = %.6g)", upper_sum),
          mid   = sprintf("Midpoint (sum = %.6g)", mid_sum),
          method_name
        ),
        hoverinfo = "skip",
        inherit = FALSE
      )
    }

    # colores con respaldo por si falta alguno en 'colors'
    col_lower <- if (!is.null(colors$lower)) colors$lower else "#a1d99b"
    col_mid   <- if (!is.null(colors$mid))   colors$mid   else "#9ecae1"
    col_upper <- if (!is.null(colors$upper)) colors$upper else "#fc9272"

    # draw in a logical order (lower, mid, upper)
    if ("lower" %in% methods) {
      p <- add_rects(p, h_lower, "lower", col_lower)
    }
    if ("mid" %in% methods) {
      p <- add_rects(p, h_mid,   "mid",   col_mid)
    }
    if ("upper" %in% methods) {
      p <- add_rects(p, h_upper, "upper", col_upper)
    }

    p <- plotly::layout(
      p,
      title = "1D Riemann sums: upper, lower, midpoint",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      showlegend = TRUE
    )
    fig <- p
    print(fig)
  }

  list(
    lower_sum = lower_sum,
    upper_sum = upper_sum,
    mid_sum   = mid_sum,
    dx        = dx,
    x_breaks  = x_breaks,
    figure    = fig
  )
}
