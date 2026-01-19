#' Binormal vectors along a 3D parametric curve
#'
#' @description
#' Computes numerical binormal vectors of a three-dimensional parametric
#' curve at selected parameter values. The curve is given by three coordinate
#' functions. At each evaluation point, the first and second derivatives of
#' the curve are approximated numerically, and their cross-product direction
#' is normalized to obtain the binormal vector.
#'
#' @details
#' For every value in \code{t_points}, the function:
#' \itemize{
#'   \item computes centered finite-difference approximations of the first
#'         and second derivatives,
#'   \item forms a direction perpendicular to both derivatives,
#'   \item normalizes that direction to obtain a unit binormal vector.
#' }
#'
#' If the first derivative is extremely small or if the first and second
#' derivative vectors are nearly parallel, the binormal direction cannot be
#' determined reliably. In these cases, the function returns \code{NA} for
#' the affected components.
#'
#' Optionally, the function can display the curve and the associated binormal
#' segments in an interactive 3D plot using \pkg{plotly}. The sampled curve,
#' evaluation points and binormal segments can be shown or hidden
#' independently.
#'
#' @param X Function returning the \code{x} coordinate of the curve as a
#' function of the parameter \code{t}.
#' @param Y Function returning the \code{y} coordinate of the curve as a
#' function of the parameter \code{t}.
#' @param Z Function returning the \code{z} coordinate of the curve as a
#' function of the parameter \code{t}.
#' @param a Lower endpoint of the parameter interval.
#' @param b Upper endpoint of the parameter interval.
#' @param t_points Numeric vector of parameter values for evaluation and
#' optional plotting.
#' @param h Step size for centered finite-difference approximations.
#' @param plot Logical; if \code{TRUE}, produces a \pkg{plotly} 3D
#' visualization showing the curve and the binormal vectors.
#' @param n_samples Number of points used to sample and display the curve in
#' the plot.
#' @param vec_scale Base length used for binormal segments. If \code{NULL},
#' it is estimated as a small proportion of the size of the sampled curve.
#' @param vec_factor Multiplicative factor applied to \code{vec_scale} to
#' adjust segment length.
#' @param curve_line List with \pkg{plotly} style options for drawing the
#' base curve.
#' @param B_line List with \pkg{plotly} style options for the binormal
#' segments.
#' @param show_curve Logical; if \code{TRUE}, the base curve is included in
#' the plot.
#' @param show_points Logical; if \code{TRUE}, the evaluation points are
#' marked in the plot.
#' @param point_marker List with \pkg{plotly} marker options for the
#' evaluation points.
#' @param scene List with 3D scene settings for the \pkg{plotly} figure.
#' @param bg Background color configuration for the \pkg{plotly} figure.
#' @param tol Numeric tolerance for detecting situations in which the
#' derivative information is too weak to determine a stable binormal
#' direction.
#'
#' @return
#' A tibble with columns \code{t}, \code{x}, \code{y}, \code{z},
#' \code{Bx}, \code{By} and \code{Bz}, containing the components of the
#' binormal vector at each evaluation point.
#'
#' @examples
#' X <- function(t) t * cos(t)
#' Y <- function(t) t * sin(3 * t)
#' Z <- function(t) t
#' binormal3d(X, Y, Z, a = 0, b = 2 * pi, t_points = c(pi / 3, pi, 5 * pi / 3))
#'
#' @importFrom tibble as_tibble
#' @export
binormal3d <- function(
    X, Y, Z,
    a, b,
    t_points,
    h = 1e-4,
    plot = FALSE,
    n_samples = 400,
    vec_scale = NULL,
    vec_factor = 1,
    curve_line = list(color = "blue",  width = 2, dash = "solid"),
    B_line    = list(color = "black", width = 5, dash = "solid"),
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

  # --- Argument validation ---------------------------------------------------
  if (b < a) stop("'b' must be >= 'a'.", call. = FALSE)
  if (!is.numeric(t_points) || any(!is.finite(t_points)))
    stop("'t_points' must be finite numeric.", call. = FALSE)
  if (any(t_points < a | t_points > b))
    stop("All 't_points' must lie within [a, b].", call. = FALSE)

  if (vec_factor <= 0) stop("'vec_factor' must be > 0.", call. = FALSE)
  if (!is.numeric(h) || h <= 0)
    stop("'h' must be positive.", call. = FALSE)

  # --- Internal helpers ------------------------------------------------------
  d1 <- function(f, t, h) (f(t + h) - f(t - h)) / (2 * h)
  d2 <- function(f, t, h) (f(t + h) - 2 * f(t) + f(t - h)) / (h * h)

  r  <- function(t) c(X(t), Y(t), Z(t))
  dot <- function(a, b) sum(a * b)
  norm <- function(a) sqrt(dot(a, a))
  cross <- function(a, b) c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )

  # --- Pointwise binormal computation ---------------------------------------
  compute_one <- function(t0) {
    r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
    r2 <- c(d2(X, t0, h), d2(Y, t0, h), d2(Z, t0, h))

    n1 <- norm(r1)
    if (n1 < tol)
      return(c(t = t0, r(t0), Bx = NA_real_, By = NA_real_, Bz = NA_real_))

    c12 <- cross(r1, r2)
    n12 <- norm(c12)
    if (n12 < tol)
      return(c(t = t0, r(t0), Bx = NA_real_, By = NA_real_, Bz = NA_real_))

    B <- c12 / n12
    c(t = t0, r(t0), B)
  }

  M <- t(vapply(t_points, compute_one, numeric(1 + 3 + 3)))
  colnames(M) <- c("t", "x", "y", "z", "Bx", "By", "Bz")
  out <- tibble::as_tibble(M)

  # --- Optional Plot ---------------------------------------------------------
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      data_curve <- if (isTRUE(show_curve) || is.null(vec_scale)) {
        curve_sample3d(X, Y, Z, a, b, n_samples = n_samples)
      } else {
        NULL
      }

      if (is.null(vec_scale)) {
        rx <- if (!is.null(data_curve)) diff(range(data_curve$x)) else 1
        ry <- if (!is.null(data_curve)) diff(range(data_curve$y)) else 1
        rz <- if (!is.null(data_curve)) diff(range(data_curve$z)) else 1
        vec_scale <- 0.05 * max(rx, ry, rz)
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
            marker = point_marker, showlegend = FALSE,
            hoverinfo = "none"
          )
      }

      for (i in seq_len(nrow(out))) {
        v  <- as.numeric(out[i, c("Bx", "By", "Bz")])
        if (any(is.na(v))) next
        p0 <- as.numeric(out[i, c("x", "y", "z")])
        seg <- make_segment_df(p0, v, vec_scale)

        plt <- plt |>
          plotly::add_trace(
            data = seg, x = ~x, y = ~y, z = ~z,
            type = "scatter3d", mode = "lines",
            line = B_line, hoverinfo = "none", showlegend = FALSE
          )
      }

      plt <- plt |>
        plotly::layout(
          title = "Binormal vectors B",
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      print(plt)
    }
  }

  out
}
