#' Unit tangent vectors along a 3D parametric curve
#'
#' @description
#' Computes numerical unit tangent vectors of a three-dimensional
#' parametric curve at selected values of the parameter. The curve is
#' defined by three functions that give its coordinate components. For
#' each evaluation point, the first derivative of the curve is
#' approximated numerically and then normalized to obtain a unit tangent
#' direction.
#'
#' @details
#' For every element of \code{t_points}, the function:
#' \itemize{
#'   \item computes a centered finite-difference approximation of the
#'         first derivative of the curve,
#'   \item evaluates the magnitude of that derivative,
#'   \item divides the derivative by its magnitude to obtain a unit vector
#'         pointing in the direction of motion of the curve at that point.
#' }
#'
#' If the magnitude of the first derivative is extremely small at a given
#' parameter value, the tangent direction becomes numerically unstable; in
#' such cases, the function returns \code{NA} for the corresponding
#' components and may emit a diagnostic message.
#'
#' Optionally, the curve and the associated tangent directions can be shown
#' in an interactive 3D plot using \pkg{plotly}. Short line segments
#' representing the tangent direction can be anchored at each evaluation
#' point. The sampled curve, the reference points and the tangent segments
#' may be displayed or hidden independently.
#'
#' @param X Function returning the \code{x} coordinate of the curve as a
#' function of the parameter \code{t}.
#' @param Y Function returning the \code{y} coordinate of the curve as a
#' function of the parameter \code{t}.
#' @param Z Function returning the \code{z} coordinate of the curve as a
#' function of the parameter \code{t}.
#' @param a Lower endpoint of the parameter interval.
#' @param b Upper endpoint of the parameter interval.
#' @param t_points Numeric vector of parameter values at which the tangent
#' direction is evaluated and optionally plotted.
#' @param h Step size for centered finite-difference approximations.
#' @param plot Logical; if \code{TRUE}, shows a 3D \pkg{plotly} visualization
#' of the curve and tangent segments.
#' @param n_samples Number of points used to sample and display the curve in
#' the 3D plot.
#' @param vec_scale Base length used for the tangent segments. If
#' \code{NULL}, it is estimated as a small fraction of the overall size of
#' the sampled curve.
#' @param vec_factor Multiplicative factor applied to \code{vec_scale} to
#' control the visual size of the tangent segments.
#' @param curve_line List with \pkg{plotly} style options for drawing the
#' base curve.
#' @param T_line List with \pkg{plotly} style options for the tangent
#' segments.
#' @param show_curve Logical; if \code{TRUE}, the base curve is included in
#' the plot.
#' @param show_points Logical; if \code{TRUE}, the evaluation points are
#' marked in the plot.
#' @param point_marker List with \pkg{plotly} marker options for the
#' evaluation points.
#' @param scene List with 3D scene settings for \pkg{plotly}.
#' @param bg Background colors for the figure, usually a list with entries
#' such as \code{paper} and \code{plot}.
#' @param tol Numeric tolerance for detecting situations in which the first
#' derivative is too small to define a stable tangent direction.
#'
#' @return
#' A tibble with columns \code{t}, \code{x}, \code{y}, \code{z},
#' \code{Tx}, \code{Ty} and \code{Tz}, where the last three columns contain
#' the components of the unit tangent vector at each evaluation point.
#'
#' @examples
#' X <- function(t) t*cos(t)
#' Y <- function(t) t*sin(3*t)
#' Z <- function(t) t
#' tangent3d(X, Y, Z, a = 0, b = 2*pi, t_points = c(pi/3, pi, 5*pi/3))
#'
#' @export
tangent3d <- function(
    X, Y, Z,
    a, b,
    t_points,
    h = 1e-4,
    plot = FALSE,
    n_samples = 400,
    vec_scale = NULL,
    vec_factor = 1,
    curve_line = list(color = "blue", width = 2, dash = "solid"),
    T_line    = list(color = "red",  width = 5, dash = "solid"),
    show_curve  = TRUE,
    show_points = TRUE,
    point_marker = list(color = "black", size = 3, symbol = "circle"),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    bg = list(paper = "white", plot = "white"),
    tol = 1e-10
) {
  # basic checks
  if (!is.function(X) || !is.function(Y) || !is.function(Z)) {
    stop("'X', 'Y' and 'Z' must be functions of t.", call. = FALSE)
  }
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || b < a) {
    stop("'a' and 'b' must be finite numeric scalars with b >= a.", call. = FALSE)
  }
  if (!is.numeric(t_points) || any(!is.finite(t_points))) {
    stop("'t_points' must be finite numeric.", call. = FALSE)
  }
  if (any(t_points < a | t_points > b)) {
    stop("All 't_points' must lie within [a, b].", call. = FALSE)
  }
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("'h' must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(vec_factor) || length(vec_factor) != 1L ||
      !is.finite(vec_factor) || vec_factor <= 0) {
    stop("'vec_factor' must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) ||
      n_samples < 2 || abs(n_samples - round(n_samples)) > .Machine$double.eps^0.5) {
    stop("'n_samples' must be an integer >= 2.", call. = FALSE)
  }

  # internal helpers ---------------------------------------------------------
  d1 <- function(f, t, h) (f(t + h) - f(t - h)) / (2 * h)
  r  <- function(t) c(X(t), Y(t), Z(t))
  dot <- function(a, b) sum(a * b)
  nrm <- function(a) sqrt(dot(a, a))

  compute_one <- function(t0) {
    r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
    n1 <- nrm(r1)
    if (n1 < tol) {
      c(t = t0, r(t0), Tx = NA_real_, Ty = NA_real_, Tz = NA_real_)
    } else {
      c(t = t0, r(t0), r1 / n1)
    }
  }

  # compute tangent data -----------------------------------------------------
  M <- t(vapply(t_points, compute_one, numeric(1 + 3 + 3)))
  colnames(M) <- c("t", "x", "y", "z", "Tx", "Ty", "Tz")
  out <- tibble::as_tibble(M)

  # optional plot ------------------------------------------------------------
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      # if curve not sampled yet (or vec_scale needed), sample it
      data_curve <- if (isTRUE(show_curve) || is.null(vec_scale)) {
        curve_sample3d(X, Y, Z, a, b, n_samples = n_samples)
      } else {
        NULL
      }

      if (is.null(vec_scale)) {
        if (!is.null(data_curve) && nrow(data_curve) > 1L) {
          rx <- diff(range(data_curve$x))
          ry <- diff(range(data_curve$y))
          rz <- diff(range(data_curve$z))
          vec_scale <- 0.05 * max(rx, ry, rz)
        } else {
          vec_scale <- 1
        }
      }
      vec_scale <- vec_scale * vec_factor

      make_segment_df <- function(p0, v, s) {
        tibble::tibble(
          x = c(p0[1], p0[1] + s * v[1]),
          y = c(p0[2], p0[2] + s * v[2]),
          z = c(p0[3], p0[3] + s * v[3])
        )
      }

      plt <- if (isTRUE(show_curve) && !is.null(data_curve)) {
        plotly::plot_ly(
          data_curve, x = ~x, y = ~y, z = ~z,
          type = "scatter3d", mode = "lines",
          line = curve_line, hoverinfo = "none", showlegend = FALSE
        )
      } else {
        plotly::plot_ly()
      }

      if (isTRUE(show_points)) {
        plt <- plt |>
          plotly::add_trace(
            x = out$x, y = out$y, z = out$z,
            type = "scatter3d", mode = "markers",
            marker = point_marker, showlegend = FALSE, hoverinfo = "none"
          )
      }

      for (i in seq_len(nrow(out))) {
        v <- as.numeric(out[i, c("Tx", "Ty", "Tz")])
        if (any(is.na(v))) next
        p0 <- as.numeric(out[i, c("x", "y", "z")])
        seg <- make_segment_df(p0, v, vec_scale)
        plt <- plt |>
          plotly::add_trace(
            data = seg, x = ~x, y = ~y, z = ~z,
            type = "scatter3d", mode = "lines",
            line = T_line, hoverinfo = "none", showlegend = FALSE
          )
      }

      plt <- plt |>
        plotly::layout(
          title = "Unit tangent vectors T",
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )
      print(plt)
    }
  }

  out
}
