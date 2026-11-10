#' Plot an interactive secant and its function
#'
#' Display sampled function values and a local line with Plotly. Computation
#' is performed separately by vc_secant(), vc_tangent(), or vc_differential().
#' @param result Plain result list from the corresponding calculation function.
#' @param f The same deterministic scalar function used for the calculation.
#'   Function values at the stored anchor points must match the result.
#' @param xlim Optional increasing pair of finite real plot limits, containing
#'   the anchor points. NULL derives limits around the anchors.
#' @param n Integer number of curve samples, at least 2. The curve is a sampled
#'   visualization, not an exact representation or a continuity test.
#' @param colors Named character vector with entries curve, line, point,
#'   differential, and change; values are Plotly/CSS colors.
#' @param title Optional scalar plot title. NULL uses a descriptive default.
#' @param xlab,ylab Scalar axis labels.
#' @return An interactive plotly htmlwidget. Zoom, pan, hover and legend
#'   controls are available. No derivative calculation is performed here.
#' @details Function samples that raise errors or produce nonfinite/nonreal
#'   scalar values are gaps; warnings during curve sampling are suppressed.
#'   The plot fails if no finite sample exists. Finite samples across an
#'   unresolved discontinuity can still be joined: users must choose a suitable
#'   domain and sampling resolution. The anchor points are checked strictly.
#'   Graphical limits and line evaluations must be finite.
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   f <- function(x) x^2
#'   vc_secant(f, 1, 2) |> vc_plot_secant(f = f)
#' }
#' @export
vc_plot_secant <- function(result, f, xlim = NULL, n = 401L,
                           colors = c(curve = "#24476b", line = "#c44736",
                                      point = "#222222", differential = "#188367",
                                      change = "#8250a3"),
                           title = NULL, xlab = "x", ylab = "y") {
  .vc_plot_model(result, f, "secant", xlim, n, colors, title, xlab, ylab)
}

#' Plot an approximate tangent and its function
#'
#' Display the approximation and its sensitivity status. A plotted line does
#' not certify that the mathematical tangent exists.
#' @inheritParams vc_plot_secant
#' @return An interactive plotly htmlwidget with function, line and anchor.
#' @details All sampling limitations of vc_plot_secant() apply. The legend
#'   identifies an approximate tangent and the subtitle reports numerical
#'   sensitivity status. One-sided methods are identified explicitly.
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   vc_tangent(sin, 0) |> vc_plot_tangent(f = sin)
#' }
#' @export
vc_plot_tangent <- function(result, f, xlim = NULL, n = 401L,
                            colors = c(curve = "#24476b", line = "#c44736",
                                       point = "#222222", differential = "#188367",
                                       change = "#8250a3"),
                            title = NULL, xlab = "x", ylab = "y") {
  .vc_plot_model(result, f, "tangent", xlim, n, colors, title, xlab, ylab)
}

#' Plot a differential and the actual finite change
#'
#' Show dx, the approximate differential dy, and the finite change delta_y.
#' Dy and delta_y are shown at the same endpoint: use hover and legend
#' controls to compare overlapping segments.
#' @inheritParams vc_plot_secant
#' @return An interactive plotly htmlwidget with the function, approximate
#'   tangent, anchors, and separate labeled increment segments.
#' @details All sampling limitations of vc_plot_secant() apply. The subtitle
#'   identifies derivative sensitivity status and the selected method. The
#'   differential endpoint uses y0 + dy at the representable x1; dx is the
#'   requested increment, which can differ slightly from x1 - x0 in doubles.
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   f <- function(x) x^2
#'   vc_differential(f, 1, 0.2) |> vc_plot_differential(f = f)
#' }
#' @export
vc_plot_differential <- function(result, f, xlim = NULL, n = 401L,
                                 colors = c(curve = "#24476b", line = "#c44736",
                                            point = "#222222", differential = "#188367",
                                            change = "#8250a3"),
                                 title = NULL, xlab = "x", ylab = "y") {
  .vc_plot_model(result, f, "differential", xlim, n, colors, title, xlab, ylab)
}

.vc_plot_model <- function(result, f, type, xlim, n, colors, title, xlab, ylab) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required only for visualization.")
  }
  if (!is.list(result) || !identical(result$type, type)) {
    stop("'result' must come from vc_", type, "().")
  }
  fields <- switch(type, secant = c("x0", "y0", "x1", "y1", "slope"),
                   tangent = c("x0", "y0", "slope"),
                   differential = c("x0", "y0", "x1", "y1", "dx", "dy", "delta_y"))
  if (!all(vapply(result[fields], .vc_scalar, logical(1)))) {
    stop("Invalid or nonfinite result fields.")
  }
  tangent <- if (type == "differential") result$tangent else result
  if (type != "secant" &&
      (!is.list(tangent$derivative) ||
       !identical(tangent$derivative$method %in% c("central", "forward", "backward"), TRUE) ||
       !is.logical(tangent$derivative$reliable) ||
       length(tangent$derivative$reliable) != 1L || is.na(tangent$derivative$reliable))) {
    stop("Missing derivative diagnostics in 'result'.")
  }
  if (!.vc_scalar(tangent$slope)) stop("Invalid line slope.")
  if (!is.function(f)) stop("'f' must be a function.")
  anchors_x <- c(result$x0, if (type != "tangent") result$x1)
  anchors_y <- c(result$y0, if (type != "tangent") result$y1)
  actual_y <- vapply(anchors_x, function(x) .vc_value(f, x), numeric(1))
  if (!identical(as.double(actual_y), as.double(anchors_y))) {
    stop("'f' does not match the stored anchor values.")
  }
  if (!.vc_scalar(n) || n < 2 || n != floor(n) || n > .Machine$integer.max) {
    stop("'n' must be an integer >= 2 within R's integer range.")
  }
  required <- c("curve", "line", "point", "differential", "change")
  if (!is.character(colors) || !all(required %in% names(colors)) ||
      anyDuplicated(names(colors)) || anyNA(colors) || any(!nzchar(colors))) {
    stop("'colors' must provide unique named curve, line, point, differential and change colors.")
  }
  label <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
  if ((!is.null(title) && !label(title)) || !label(xlab) || !label(ylab)) {
    stop("'title', 'xlab' and 'ylab' must be scalar character labels.")
  }
  if (is.null(xlim)) {
    span <- max(anchors_x) - min(anchors_x)
    if (!is.finite(span)) stop("Anchor span overflowed; choose a smaller scale.")
    pad <- if (span > 0) span / 2 else max(1, abs(result$x0) / 10)
    xlim <- c(min(anchors_x) - pad, max(anchors_x) + pad)
  }
  if (!is.numeric(xlim) || is.complex(xlim) || length(xlim) != 2L ||
      any(!is.finite(xlim)) || xlim[1] >= xlim[2] ||
      min(anchors_x) < xlim[1] || max(anchors_x) > xlim[2]) {
    stop("'xlim' must be finite increasing limits containing the anchors.")
  }
  xs <- seq(xlim[1], xlim[2], length.out = n)
  ys <- vapply(xs, function(x) {
    value <- tryCatch(suppressWarnings(f(x)), error = function(e) NA_real_)
    if (.vc_scalar(value)) as.double(value) else NA_real_
  }, numeric(1))
  if (!any(is.finite(ys))) stop("No finite curve samples; change 'xlim' or 'n'.")
  line_y <- result$y0 + tangent$slope * (xlim - result$x0)
  if (any(!is.finite(line_y))) stop("Line values overflowed; narrow 'xlim'.")
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p, x = xs, y = ys, type = "scatter", mode = "lines",
                        name = "f(x) (sampled)", connectgaps = FALSE,
                        line = list(color = colors[["curve"]], width = 2),
                        hovertemplate = "x=%{x}<br>f(x)=%{y}<extra>Function</extra>")
  line_name <- if (type == "secant") "Secant" else "Approximate tangent"
  p <- plotly::add_trace(p, x = xlim, y = line_y, type = "scatter", mode = "lines",
                        name = line_name,
                        line = list(color = colors[["line"]], width = 2, dash = "dash"),
                        hovertemplate = "x=%{x}<br>line=%{y}<extra>Local line</extra>")
  p <- plotly::add_trace(p, x = anchors_x, y = anchors_y, type = "scatter", mode = "markers",
                        name = "Function anchors", marker = list(color = colors[["point"]], size = 9),
                        hovertemplate = "x=%{x}<br>y=%{y}<extra>Anchor</extra>")
  subtitle <- "Secant slope; no derivative or continuity claim"
  if (type != "secant") {
    status <- if (tangent$derivative$reliable) "passed" else "FAILED"
    subtitle <- paste(tangent$derivative$method, "differences; sensitivity", status,
                      "- not proof of differentiability")
  }
  if (type == "differential") {
    linear_y <- .vc_finite(result$y0 + result$dy, "Differential endpoint")
    p <- plotly::add_trace(p, x = c(result$x0, result$x1),
                          y = rep(result$y0, 2), type = "scatter", mode = "lines",
                          name = paste0("dx = ", format(result$dx, digits = 6)),
                          line = list(color = colors[["point"]], width = 3))
    p <- plotly::add_trace(p, x = rep(result$x1, 2), y = c(result$y0, linear_y),
                          type = "scatter", mode = "lines+markers",
                          name = paste0("dy (approx.) = ", format(result$dy, digits = 6)),
                          line = list(color = colors[["differential"]], width = 5),
                          marker = list(color = colors[["differential"]], size = 6))
    p <- plotly::add_trace(p, x = rep(result$x1, 2), y = c(result$y0, result$y1),
                          type = "scatter", mode = "lines",
                          name = paste0("delta_y = ", format(result$delta_y, digits = 6)),
                          line = list(color = colors[["change"]], width = 2, dash = "dot"))
  }
  if (is.null(title)) title <- switch(type, secant = "Secant line",
                                    tangent = "Approximate tangent line",
                                    differential = "Differential and finite change")
  p <- plotly::layout(p, title = list(text = paste0(title, "<br><sup>", subtitle, "</sup>")),
                     xaxis = list(title = xlab, range = xlim, zeroline = FALSE),
                     yaxis = list(title = ylab, zeroline = FALSE),
                     hovermode = "closest", paper_bgcolor = "white", plot_bgcolor = "white",
                     margin = list(t = 90, b = 60),
                     legend = list(orientation = "h", y = -0.2))
  plotly::config(p, displaylogo = FALSE, scrollZoom = TRUE)
}
