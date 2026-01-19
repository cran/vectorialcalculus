#' 2D line integral of a vector field with visualization
#'
#' This function computes a numerical approximation of the line integral
#' of a planar vector field along a parametric curve r(t) on the interval
#' from a to b. The derivative of the curve is approximated by finite
#' differences and the integral is evaluated either by an adaptive numerical
#' integrator or by a composite Simpson rule.
#'
#' Optionally, the function can build an interactive \pkg{plotly} figure that
#' shows a grid of arrows for the vector field and the curve colored by the
#' local power field(r(t)) * r'(t).
#'
#' @param field Vector field in the plane. A function \code{function(x, y)} that
#'   returns a numeric vector of length 2 \code{c(Fx, Fy)}.
#' @param r Parametric curve in the plane. A function \code{function(t)} that
#'   returns a numeric vector of length 2 \code{c(x, y)}.
#' @param a Numeric scalar. Left endpoint of the parameter interval.
#' @param b Numeric scalar. Right endpoint of the parameter interval.
#'   Must satisfy \code{b > a}.
#' @param plot Logical. If \code{TRUE}, an interactive \pkg{plotly} figure is
#'   created with the field and the curve.
#' @param n_curve Integer. Number of parameter values used to sample the curve.
#' @param grid_n Integer. Number of grid points per axis used to draw the
#'   vector field arrows.
#' @param padding Numeric scalar. Relative margin added around the bounding
#'   box of the curve when building the field grid.
#' @param h Numeric scalar or \code{NULL}. Step size used in the finite
#'   difference approximation of the derivative of the curve. If \code{NULL},
#'   a small step is chosen automatically, based on the length of the interval
#'   \code{b - a}.
#' @param method Character string. Integration method for the line integral.
#'   One of \code{"adaptive"} (uses \code{stats::integrate}) or
#'   \code{"simpson"} (composite Simpson rule).
#' @param n_simpson Integer. Number of subintervals used when
#'   \code{method = "simpson"}. If it is odd, it is increased by one
#'   internally.
#' @param arrow_scale Numeric scalar. Controls the overall length of the field
#'   arrows as a fraction of the plot span.
#' @param normalize_bias Numeric scalar. Saturation parameter used to avoid
#'   extremely long arrows for large field magnitudes.
#' @param field_color Character string. Color used for the field arrows.
#' @param field_width Numeric scalar. Line width for the field arrows.
#' @param traj_palette Color scale used to represent the power along the
#'   trajectory. Passed to \pkg{plotly} as a colorscale name.
#' @param traj_width Numeric scalar. Line width for the trajectory.
#' @param show_markers Logical. If \code{TRUE}, markers are drawn along the
#'   trajectory in addition to the line.
#' @param scene List with \pkg{plotly} scene options (axis titles, aspect
#'   mode, etc.) passed to \code{plotly::layout()}.
#' @param bg List with background colors for \pkg{plotly}, with components
#'   \code{paper} and \code{plot}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{value}: numeric value of the line integral.
#'   \item \code{samples}: data frame with sampled points, velocities
#'         and power along the trajectory.
#'   \item \code{fig}: \pkg{plotly} object when \code{plot = TRUE},
#'         otherwise \code{NULL}.
#' }
#'
#' @examples
#' # Simple example:
#' # field(x, y) = (y, -x), r(t) = (cos t, sin t), t in [0, 2*pi]
#' line_integral_vector2d(
#'   field = function(x, y) c(y, -x),
#'   r = function(t) c(cos(t), sin(t)),
#'   a = 0, b = 2*pi, plot = FALSE
#' )
#'
#' @export
line_integral_vector2d <- function(
    field, r, a, b,
    plot = TRUE,
    n_curve = 600,
    grid_n  = 15,
    padding = 0.15,
    h = NULL,
    method = c("adaptive","simpson"),
    n_simpson = 1000,
    arrow_scale = 0.08,
    normalize_bias = 1,
    field_color = "rgba(0,0,0,0.55)",
    field_width = 1.8,
    traj_palette = "RdBu",
    traj_width = 5,
    show_markers = FALSE,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
){
  method <- match.arg(method)

  if (!is.function(field)) stop("'field' must be function(x, y) -> c(Fx, Fy).", call. = FALSE)
  if (!is.function(r)) stop("'r' must be function(t) -> c(x, y).", call. = FALSE)
  if (!is.numeric(a) || !is.numeric(b) ||
      length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || b <= a) {
    stop("'a' and 'b' must be finite scalars with b > a.", call. = FALSE)
  }

  # small helper to evaluate r(t) safely
  eval_r <- function(t) {
    out <- r(t)
    if (!is.numeric(out) || length(out) != 2L || any(!is.finite(out))) {
      stop("r(t) must return a finite numeric vector of length 2.", call. = FALSE)
    }
    as.numeric(out)
  }

  eval_field <- function(x, y) {
    out <- field(x, y)
    if (!is.numeric(out) || length(out) != 2L || any(!is.finite(out))) {
      stop("field(x, y) must return a finite numeric vector of length 2.", call. = FALSE)
    }
    as.numeric(out)
  }

  # Path samples (for plotting and domain box)
  tt  <- seq(a, b, length.out = n_curve)
  Rxy <- t(vapply(tt, eval_r, numeric(2)))
  xs  <- Rxy[, 1]
  ys  <- Rxy[, 2]

  # r'(t) via centered finite differences (clamped near ends)
  if (is.null(h)) {
    h <- (b - a) * 1e-4 + 1e-6
  } else {
    if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
      stop("'h' must be a positive finite numeric scalar or NULL.", call. = FALSE)
    }
  }

  rprime <- function(t) {
    t1 <- max(a, min(b, t - h))
    t2 <- max(a, min(b, t + h))
    if (t2 == t1) {
      stop("Step size 'h' is too small relative to [a, b].", call. = FALSE)
    }
    (eval_r(t2) - eval_r(t1)) / (t2 - t1)
  }

  # Power density p(t) = field(r(t)) · r'(t)
  Fr <- t(vapply(seq_along(tt), function(i) eval_field(xs[i], ys[i]), numeric(2)))
  rp <- t(vapply(tt, rprime, numeric(2)))
  speed <- sqrt(rowSums(rp^2))
  power <- rowSums(Fr * rp)

  # Integrand and integral value
  integrand <- function(t) {
    p  <- eval_r(t)
    v  <- rprime(t)
    fv <- eval_field(p[1], p[2])
    sum(fv * v)
  }
  integrand_vec <- Vectorize(integrand, "t")

  value <- if (method == "adaptive") {
    stats::integrate(integrand_vec, lower = a, upper = b, rel.tol = 1e-6)$value
  } else {
    n <- as.integer(n_simpson)
    if (n <= 0L) stop("'n_simpson' must be positive.", call. = FALSE)
    if (n %% 2L == 1L) n <- n + 1L
    tS <- seq(a, b, length.out = n + 1L)
    yS <- integrand_vec(tS)
    hS <- (b - a) / n
    hS * (yS[1] + yS[length(yS)] +
            4 * sum(yS[seq(2, n, by = 2)]) +
            2 * sum(yS[seq(3, n - 1, by = 2)])) / 3
  }

  # Field grid around path
  rx <- range(xs)
  ry <- range(ys)
  pad_x <- padding * diff(rx)
  pad_y <- padding * diff(ry)
  if (!is.finite(pad_x) || pad_x == 0) pad_x <- 1
  if (!is.finite(pad_y) || pad_y == 0) pad_y <- 1
  xr <- c(rx[1] - pad_x, rx[2] + pad_x)
  yr <- c(ry[1] - pad_y, ry[2] + pad_y)

  gx <- seq(xr[1], xr[2], length.out = grid_n)
  gy <- seq(yr[1], yr[2], length.out = grid_n)

  grid_pts <- expand.grid(x = gx, y = gy)
  Fgrid    <- t(vapply(seq_len(nrow(grid_pts)),
                       function(i) eval_field(grid_pts$x[i], grid_pts$y[i]),
                       numeric(2)))
  mag <- sqrt(rowSums(Fgrid^2))

  # Saturate and scale arrows
  Fsat <- Fgrid / sqrt(mag^2 + normalize_bias)
  mlen <- sqrt(rowSums(Fsat^2))
  span <- max(diff(xr), diff(yr))
  if (!is.finite(span) || span <= 0) span <- 1
  L <- arrow_scale * span * mlen
  diru <- Fsat
  nz <- mlen > 0
  diru[nz, ]  <- Fsat[nz, , drop = FALSE] / mlen[nz]
  diru[!nz, ] <- 0

  X0 <- grid_pts$x
  Y0 <- grid_pts$y
  X1 <- X0 + diru[, 1] * L
  Y1 <- Y0 + diru[, 2] * L

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need the 'plotly' package installed.", call. = FALSE)
    } else {
      plt <- plotly::plot_ly()

      # Field as line segments (with NA separators)
      xs_ar <- as.numeric(rbind(X0, X1, NA))
      ys_ar <- as.numeric(rbind(Y0, Y1, NA))

      plt <- plt |>
        plotly::add_trace(
          x = xs_ar, y = ys_ar, z = 0,
          type = "scatter3d", mode = "lines",
          line = list(color = field_color, width = field_width),
          hoverinfo = "none", showlegend = FALSE
        ) |>
        plotly::add_trace(
          x = xs, y = ys, z = rep(0, length(xs)),
          type = "scatter3d",
          mode = if (isTRUE(show_markers)) "lines+markers" else "lines",
          line = list(width = traj_width,
                      color = power,
                      colorscale = traj_palette),
          marker = list(
            size = 3,
            color = power,
            colorscale = traj_palette,
            showscale = TRUE
          ),
          showlegend = FALSE,
          name = "trajectory"
        ) |>
        plotly::layout(
          title = sprintf("Vector line integral: Work ~ %.6g", value),
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      fig <- plt
      print(fig)
    }
  }

  samples <- data.frame(
    t = tt,
    x = xs,
    y = ys,
    px = rp[, 1],
    py = rp[, 2],
    speed = speed,
    power = power
  )

  list(value = value, samples = samples, fig = fig)
}
