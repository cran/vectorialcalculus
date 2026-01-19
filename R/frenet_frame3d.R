#' Frenet-Serret frame for a 3D parametric curve
#'
#' @description
#' Computes the Frenet-Serret frame, that is, the tangent, normal and
#' binormal vectors of a three dimensional parametric curve at selected
#' values of the parameter. The frame is obtained from numerical
#' approximations of the first and second derivatives of the curve.
#' Optionally, the curve and the three vector fields can be displayed in a
#' 3D interactive visualization using \pkg{plotly}.
#'
#' @details
#' At each parameter value in \code{t_points}, the function:
#' \itemize{
#'   \item computes finite difference approximations of the first and second
#'         derivatives of the curve,
#'   \item normalizes the first derivative to obtain the unit tangent
#'         direction,
#'   \item uses the first and second derivatives to construct a principal
#'         normal direction,
#'   \item constructs the binormal direction as a unit vector orthogonal to
#'         both the tangent and the normal,
#'   \item evaluates a numerical estimate of the curvature using the same
#'         derivative information.
#' }
#'
#' If the derivative information is too small or nearly degenerate (for
#' example, when the tangent direction cannot be reliably obtained), some
#' components of the frame may be set to \code{NA}. The tolerance parameter
#' \code{tol} controls how these situations are detected.
#'
#' When \code{plot = TRUE}, the function displays:
#' \itemize{
#'   \item a sampled representation of the curve,
#'   \item the evaluation points,
#'   \item short line segments indicating the tangent, normal and binormal
#'         directions at each evaluation point.
#' }
#'
#' All visual elements can be styled or shown selectively through the
#' corresponding arguments.
#'
#' @param X Function returning the \code{x} coordinate of the curve as a
#' function of the parameter \code{t}.
#' @param Y Function returning the \code{y} coordinate of the curve.
#' @param Z Function returning the \code{z} coordinate of the curve.
#' @param a Lower endpoint of the parameter interval.
#' @param b Upper endpoint of the parameter interval.
#' @param t_points Numeric vector with the parameter values where the frame
#' is computed and optionally plotted.
#' @param h Step size for centered finite difference approximations.
#' @param plot Logical; if \code{TRUE}, shows a 3D \pkg{plotly} visualization
#' of the curve together with the three vector fields.
#' @param n_samples Number of points used to sample the curve for plotting.
#' @param vec_scale Base scaling factor for the vector segments. If
#' \code{NULL}, it is estimated from the overall size of the sampled curve.
#' @param curve_line Style options for drawing the base curve.
#' @param T_line Style options for tangent vector segments.
#' @param N_line Style options for normal vector segments.
#' @param B_line Style options for binormal vector segments.
#' @param show_curve Logical; if \code{TRUE}, the base curve appears in the
#' plot.
#' @param show_points Logical; if \code{TRUE}, the evaluation points are
#' marked on the curve.
#' @param point_marker Plotly marker style for the evaluation points.
#' @param scene Plotly 3D scene configuration.
#' @param bg Background settings for the \pkg{plotly} figure.
#' @param tol Numeric tolerance used to detect degenerate derivative
#' situations.
#'
#' @return
#' A tibble containing the parameter values and the coordinates of:
#' \itemize{
#'   \item the point on the curve,
#'   \item the tangent vector,
#'   \item the normal vector,
#'   \item the binormal vector,
#'   \item a numerical estimate of the curvature.
#' }
#' Columns are named
#' \code{t, x, y, z, Tx, Ty, Tz, Nx, Ny, Nz, Bx, By, Bz, kappa}.
#'
#' @examples
#' X <- function(t) t*cos(t)
#' Y <- function(t) t*sin(3*t)
#' Z <- function(t) t
#' frenet_frame3d(
#'   X, Y, Z, a = 0, b = 2*pi,
#'   t_points = c(pi/3, pi, 5*pi/3)
#' )
#'
#' @export
frenet_frame3d <- function(
    X, Y, Z,
    a, b,
    t_points,
    h = 1e-4,
    plot = FALSE,
    n_samples = 400,
    vec_scale = NULL,
    curve_line = list(color = "blue",  width = 2, dash = "solid"),
    T_line    = list(color = "red",    width = 4, dash = "solid"),
    N_line    = list(color = "green",  width = 4, dash = "solid"),
    B_line    = list(color = "black",  width = 4, dash = "solid"),
    show_curve = TRUE,
    show_points = TRUE,
    point_marker = list(color = "blue", size = 3, symbol = "circle"),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    bg = list(paper = "white", plot = "white"),
    tol = 1e-10
) {

  # sanity checks
  if (b < a) stop("'b' must be >= 'a'.", call. = FALSE)
  if (!is.numeric(t_points) || any(!is.finite(t_points))) {
    stop("'t_points' must be finite numeric.", call. = FALSE)
  }
  if (any(t_points < a | t_points > b)) {
    stop("All 't_points' must lie in [a, b].", call. = FALSE)
  }
  if (!is.numeric(h) || h <= 0) {
    stop("'h' must be positive.", call. = FALSE)
  }

  # derivatives & utilities
  d1 <- function(f, t, h) (f(t + h) - f(t - h)) / (2 * h)
  d2 <- function(f, t, h) (f(t + h) - 2 * f(t) + f(t - h)) / (h * h)

  r <- function(t) c(X(t), Y(t), Z(t))

  dot  <- function(a, b) sum(a * b)
  norm <- function(a) sqrt(dot(a, a))

  cross <- function(a, b) c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )

  # compute Frenet frame at one t
  compute_one <- function(t0) {
    r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
    r2 <- c(d2(X, t0, h), d2(Y, t0, h), d2(Z, t0, h))

    n1 <- norm(r1)
    if (n1 < tol) {
      return(c(t = t0, r(t0), rep(NA_real_, 9L), kappa = NA_real_))
    }

    T <- r1 / n1

    c12 <- cross(r1, r2)
    n12 <- norm(c12)

    # zero curvature
    if (n12 < tol) {
      return(c(t = t0, r(t0), T, rep(NA_real_, 6L), kappa = 0))
    }

    B <- c12 / n12
    N <- cross(B, T)
    N <- N / norm(N)

    kappa <- n12 / (n1^3)

    c(t = t0, r(t0), T, N, B, kappa = kappa)
  }

  # assemble output
  M <- t(vapply(t_points, compute_one,
                numeric(1 + 3 + 3 + 3 + 3 + 1)))
  colnames(M) <- c(
    "t", "x", "y", "z",
    "Tx", "Ty", "Tz",
    "Nx", "Ny", "Nz",
    "Bx", "By", "Bz",
    "kappa"
  )

  out <- tibble::as_tibble(M)

  # optional plot
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

      make_seg <- function(p0, v, s) {
        tibble::tibble(
          x = c(p0[1], p0[1] + s * v[1]),
          y = c(p0[2], p0[2] + s * v[2]),
          z = c(p0[3], p0[3] + s * v[3])
        )
      }

      plt <- if (!is.null(data_curve) && isTRUE(show_curve)) {
        plotly::plot_ly(
          data_curve, x = ~x, y = ~y, z = ~z,
          type = "scatter3d", mode = "lines",
          line = curve_line,
          hoverinfo = "none",
          showlegend = FALSE
        )
      } else {
        plotly::plot_ly()
      }

      if (isTRUE(show_points)) {
        plt <- plt |>
          plotly::add_trace(
            x = out$x, y = out$y, z = out$z,
            type = "scatter3d", mode = "markers",
            marker = point_marker,
            hoverinfo = "none",
            showlegend = FALSE
          )
      }

      add_vec <- function(plt, cols, style) {
        for (i in seq_len(nrow(out))) {
          v <- as.numeric(out[i, cols])
          if (any(is.na(v))) next
          p0 <- as.numeric(out[i, c("x", "y", "z")])
          seg <- make_seg(p0, v, vec_scale)
          plt <- plt |>
            plotly::add_trace(
              data = seg,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = style,
              hoverinfo = "none",
              showlegend = FALSE
            )
        }
        plt
      }

      plt <- add_vec(plt, c("Tx", "Ty", "Tz"), T_line)
      plt <- add_vec(plt, c("Nx", "Ny", "Nz"), N_line)
      plt <- add_vec(plt, c("Bx", "By", "Bz"), B_line)

      plt <- plt |>
        plotly::layout(
          title = "Frenet-Serret frame (T, N, B)",
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      print(plt)
    }
  }

  out
}
