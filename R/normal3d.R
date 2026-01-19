#' Principal normal vectors along a 3D curve
#'
#' @description
#' Computes numerical principal normal vectors of a three-dimensional
#' parametric curve at several parameter values. The curve is described by
#' three coordinate functions \code{X}, \code{Y} and \code{Z}. At each
#' evaluation point, the function approximates the first and second
#' derivatives of the curve, builds the unit tangent and binormal vectors,
#' and then obtains the principal normal as the unit vector orthogonal to
#' both of them.
#'
#' @details
#' For every parameter value in \code{t_points}, the function:
#' \itemize{
#'   \item approximates the first derivative of the curve with respect to the
#'         parameter,
#'   \item normalizes this derivative to obtain a unit tangent direction,
#'   \item uses the first and second derivative vectors to construct a
#'         direction orthogonal to the tangent and interprets it as a
#'         binormal direction,
#'   \item builds the principal normal direction as a unit vector orthogonal
#'         to both the tangent and the binormal.
#' }
#'
#' When the curvature of the curve at a given parameter value is extremely
#' small, the normal direction becomes poorly defined from a numerical point
#' of view. In such situations, the function marks the corresponding
#' components of the normal vector as \code{NA}.
#'
#' Optionally, the function can display the curve and the associated normal
#' segments in a 3D interactive plot using \pkg{plotly}. The base curve,
#' the evaluation points and the normal segments can be shown or hidden
#' independently.
#'
#' @param X Function giving the \code{x} coordinate of the curve as a function
#' of the parameter \code{t}.
#' @param Y Function giving the \code{y} coordinate of the curve as a function
#' of the parameter \code{t}.
#' @param Z Function giving the \code{z} coordinate of the curve as a function
#' of the parameter \code{t}.
#' @param a Lower endpoint of the parameter interval.
#' @param b Upper endpoint of the parameter interval.
#' @param t_points Numeric vector of parameter values at which the principal
#' normal is evaluated and, optionally, plotted.
#' @param h Step size for the centered finite-difference approximations used
#' to compute derivatives.
#' @param plot Logical; if \code{TRUE}, displays a 3D plot of the curve and
#' the corresponding normal segments using \pkg{plotly}.
#' @param n_samples Number of points used to sample and draw the curve for
#' plotting purposes.
#' @param vec_scale Base length used for the normal segments. If
#' \code{NULL}, it is estimated as a small fraction of the overall size
#' of the sampled curve.
#' @param vec_factor Multiplicative factor applied to \code{vec_scale} to
#' control the visual length of the normal segments.
#' @param curve_line List with \pkg{plotly} style options for the base curve.
#' @param N_line List with \pkg{plotly} style options for the normal segments.
#' @param show_curve Logical; if \code{TRUE}, the base curve is drawn.
#' @param show_points Logical; if \code{TRUE}, the evaluation points are
#' marked on the curve.
#' @param point_marker List with \pkg{plotly} marker options for the
#' evaluation points.
#' @param scene List with 3D scene settings for \pkg{plotly}.
#' @param bg Background colors for the figure, given as a list with entries
#' such as \code{paper} and \code{plot}.
#' @param tol Numeric tolerance used to detect singular or nearly singular
#' situations in which the normal direction cannot be computed reliably.
#'
#' @return
#' A tibble with columns \code{t}, \code{x}, \code{y}, \code{z}, \code{Nx},
#' \code{Ny} and \code{Nz}, where the last three columns contain the
#' components of the principal normal vector at each parameter value.
#'
#' @examples
#' X <- function(t) t*cos(t)
#' Y <- function(t) t*sin(3*t)
#' Z <- function(t) t
#' normal3d(X, Y, Z, a = 0, b = 2*pi, t_points = c(pi/3, pi, 5*pi/3))
#'
#' @export
normal3d <- function(
    X, Y, Z,
    a, b,
    t_points,
    h = 1e-4,
    plot = FALSE,
    n_samples = 400,
    vec_scale = NULL,
    vec_factor = 1,
    curve_line = list(color = "blue",  width = 2, dash = "solid"),
    N_line    = list(color = "green", width = 5, dash = "solid"),
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
  if (!is.function(X)) stop("'X' must be a function of t.", call. = FALSE)
  if (!is.function(Y)) stop("'Y' must be a function of t.", call. = FALSE)
  if (!is.function(Z)) stop("'Z' must be a function of t.", call. = FALSE)

  if (b < a) stop("'b' must be >= 'a'.", call. = FALSE)
  if (!is.numeric(t_points) || any(!is.finite(t_points)))
    stop("'t_points' must be finite numeric.", call. = FALSE)
  if (any(t_points < a | t_points > b))
    stop("All 't_points' must lie within [a, b].", call. = FALSE)
  if (!is.numeric(vec_factor) || length(vec_factor) != 1L || vec_factor <= 0)
    stop("'vec_factor' must be a positive numeric scalar.", call. = FALSE)
  if (!is.numeric(h) || length(h) != 1L || h <= 0)
    stop("'h' must be a positive numeric scalar.", call. = FALSE)

  # internal helpers ---------------------------------------------------------
  d1 <- function(f, t, h) (f(t + h) - f(t - h)) / (2 * h)
  d2 <- function(f, t, h) (f(t + h) - 2 * f(t) + f(t - h)) / (h * h)

  r  <- function(t) c(X(t), Y(t), Z(t))

  dot  <- function(a, b) sum(a * b)
  norm <- function(a) sqrt(dot(a, a))

  cross <- function(a, b) c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )

  compute_one <- function(t0) {
    r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
    r2 <- c(d2(X, t0, h), d2(Y, t0, h), d2(Z, t0, h))

    n1 <- norm(r1)
    if (n1 < tol) {
      return(c(t = t0, r(t0), Nx = NA_real_, Ny = NA_real_, Nz = NA_real_))
    }

    T  <- r1 / n1
    c12 <- cross(r1, r2)
    n12 <- norm(c12)

    if (n12 < tol) {
      return(c(t = t0, r(t0), Nx = NA_real_, Ny = NA_real_, Nz = NA_real_))
    }

    B <- c12 / n12
    N <- cross(B, T)
    N <- N / norm(N)

    c(t = t0, r(t0), N)
  }

  M <- t(vapply(t_points, compute_one, numeric(1 + 3 + 3)))
  colnames(M) <- c("t", "x", "y", "z", "Nx", "Ny", "Nz")
  out <- tibble::as_tibble(M)

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
            marker = point_marker, showlegend = FALSE, hoverinfo = "none"
          )
      }

      for (i in seq_len(nrow(out))) {
        v  <- as.numeric(out[i, c("Nx", "Ny", "Nz")])
        if (any(is.na(v))) next
        p0 <- as.numeric(out[i, c("x", "y", "z")])
        seg <- make_segment_df(p0, v, vec_scale)
        plt <- plt |>
          plotly::add_trace(
            data = seg, x = ~x, y = ~y, z = ~z,
            type = "scatter3d", mode = "lines",
            line = N_line, hoverinfo = "none", showlegend = FALSE
          )
      }

      plt <- plt |>
        plotly::layout(
          title = "Principal normal vectors N",
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )
      print(plt)
    }
  }

  out
}
