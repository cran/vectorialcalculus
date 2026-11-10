#' Compare a function and its numerical Taylor polynomial interactively
#'
#' Show the sampled function and polynomial above a signed/absolute error
#' panel. Hover reports distance from the expansion point; zoom, pan and
#' the horizontal range slider explore behavior away from that point.
#' @param result Plain result from vc_taylor().
#' @inheritParams vc_plot_secant
#' @param xlim Optional increasing finite pair containing the expansion point.
#'   NULL uses the range of result$comparison$x together with x0, padded by
#'   ten percent; without a nonzero comparison range it uses x0 +/- 1.
#' @param colors Named character vector with curve, polynomial, error and
#'   point entries specifying Plotly/CSS colors.
#' @return A plotly htmlwidget with function, polynomial, expansion point,
#'   signed error and absolute error traces. Derivatives are not recalculated.
#' @details Error means the sampled difference f(x)-P(x), not a rigorous
#'   Taylor remainder bound. Undefined/erroring/nonfinite function samples
#'   create gaps and NA errors. Warnings during sampling are suppressed.
#'   Polynomial overflow also creates gaps; at least one polynomial sample
#'   must be finite. A grid can miss singularities or oscillations. The
#'   subtitle displays diagnostic status and the number of evaluable errors.
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   exp |> vc_taylor(x0 = 0, x = seq(-1,1,length.out=21)) |>
#'     vc_plot_taylor(f = exp)
#' }
#' @export
vc_plot_taylor <- function(result, f, xlim = NULL, n = 401L,
                           colors = c(curve = "#24476b", polynomial = "#c44736",
                                      error = "#188367", point = "#222222"),
                           title = NULL, xlab = "x", ylab = "y") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required only for visualization.")
  .vc_taylor_plot_validate(result, f, n, colors, title, xlab, ylab)
  if (is.null(xlim)) {
    points <- c(result$x0, result$comparison$x)
    bounds <- range(points)
    span <- bounds[2]-bounds[1]
    if (!is.finite(span)) stop("Comparison span overflowed; supply finite 'xlim'.")
    pad <- if (span > 0) span/10 else 1
    xlim <- bounds + c(-pad,pad)
  }
  if (!is.numeric(xlim) || is.complex(xlim) || length(xlim) != 2L ||
      any(!is.finite(xlim)) || xlim[1] >= xlim[2] ||
      result$x0 < xlim[1] || result$x0 > xlim[2])
    stop("'xlim' must be finite increasing limits containing x0.")
  grid <- seq(xlim[1], xlim[2], length.out = n)
  comparison <- .vc_taylor_compare(f, result$polynomial, grid, result$x0)
  if (!any(is.finite(comparison$approximation))) stop("No finite polynomial samples.")
  p <- .vc_taylor_curve_panel(comparison, result, colors)
  e <- .vc_taylor_error_panel(comparison, colors)
  plot <- plotly::subplot(p, e, nrows = 2L, shareX = TRUE, titleY = TRUE,
                          heights = c(.65,.35), margin = .08)
  if (is.null(title)) title <- paste("Taylor approximation, order", result$order)
  status <- if (result$reliable) "passed" else "FAILED"
  subtitle <- paste0("Sensitivity ", status, "; not a smoothness proof. Evaluable errors: ",
                      sum(comparison$status == "ok"), "/", n)
  plot <- plotly::layout(plot,
    title = list(text = paste0(title, "<br><sup>", subtitle, "</sup>")),
    xaxis = list(title = xlab, range = xlim,
                 rangeslider = list(visible = TRUE, thickness = .1)),
    yaxis = list(title = ylab), yaxis2 = list(title = "Sampled error", zeroline = TRUE),
    paper_bgcolor = "white", plot_bgcolor = "white", hovermode = "x unified",
    margin = list(t = 110, b = 130),
    legend = list(orientation = "h", y = 1.04, yanchor = "bottom"))
  plot$height <- 800
  plotly::config(plot, displaylogo = FALSE, scrollZoom = TRUE)
}

.vc_taylor_plot_validate <- function(result, f, n, colors, title, xlab, ylab) {
  if (!is.list(result) || !identical(result$type, "taylor") ||
      !.vc_scalar(result$x0) || !.vc_scalar(result$y0) ||
      !.vc_scalar(result$order) || !(result$order %in% c(1,2)) ||
      !is.function(result$polynomial) || !is.logical(result$reliable) ||
      length(result$reliable) != 1L || is.na(result$reliable))
    stop("'result' must be a valid vc_taylor() result.")
  if (!identical(.vc_value(f, result$x0), as.double(result$y0)))
    stop("'f' does not match the stored expansion value.")
  if (!.vc_scalar(n) || n < 2 || n != floor(n) || n > .Machine$integer.max)
    stop("'n' must be an integer >= 2 within R's integer range.")
  keys <- c("curve", "polynomial", "error", "point")
  if (!is.character(colors) || !all(keys %in% names(colors)) ||
      anyDuplicated(names(colors)) || anyNA(colors) || any(!nzchar(colors)))
    stop("'colors' must have unique curve, polynomial, error and point entries.")
  label <- function(z) is.character(z) && length(z) == 1L && !is.na(z)
  if ((!is.null(title) && !label(title)) || !label(xlab) || !label(ylab))
    stop("Plot labels must be scalar character values.")
}

.vc_taylor_curve_panel <- function(data, result, colors) {
  hover <- paste0("|x-x0|=", format(data$distance, digits = 6),
                  "<br>f-P=", ifelse(is.na(data$error), "not evaluable", format(data$error, digits = 6)),
                  "<br>|f-P|=", ifelse(is.na(data$abs_error), "not evaluable", format(data$abs_error, digits = 6)))
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p, x = data$x, y = data$value, type = "scatter", mode = "lines",
    name = "f(x) (sampled)", connectgaps = FALSE,
    line = list(color = colors[["curve"]], width = 2))
  p <- plotly::add_trace(p, x = data$x, y = data$approximation, text = hover,
    type = "scatter", mode = "lines", name = "Numerical Taylor polynomial",
    connectgaps = FALSE, line = list(color = colors[["polynomial"]], width = 2, dash = "dash"),
    hovertemplate = "x=%{x}<br>P=%{y}<br>%{text}<extra>Taylor</extra>")
  plotly::add_trace(p, x = result$x0, y = result$y0, type = "scatter", mode = "markers",
    name = "Expansion point", marker = list(color = colors[["point"]], size = 9))
}

.vc_taylor_error_panel <- function(data, colors) {
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p, x = data$x, y = data$error, type = "scatter", mode = "lines",
    name = "Signed sampled error f-P", connectgaps = FALSE,
    line = list(color = colors[["error"]], width = 2))
  plotly::add_trace(p, x = data$x, y = data$abs_error, type = "scatter", mode = "lines",
    name = "Absolute sampled error", connectgaps = FALSE,
    line = list(color = colors[["polynomial"]], width = 2, dash = "dot"))
}
