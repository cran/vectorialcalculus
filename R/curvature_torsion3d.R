# -------------------------------------------------------------------
#' Curvature and torsion of a 3D parametric curve
#'
#' @description
#' Computes numerical curvature and torsion of a three-dimensional
#' parametric curve at a specific value of the parameter. The curve is
#' described by three functions that give the coordinate components.
#' All derivatives are approximated using centered finite differences
#' of first, second and third order.
#'
#' @details
#' The curvature at the evaluation point measures how sharply the curve
#' bends at that location. It is computed from the first and second
#' derivative vectors. The torsion measures how the curve deviates from
#' being planar and is computed from the first, second and third
#' derivative vectors. If the first derivative vector is nearly zero, or
#' if the first and second derivative vectors are nearly parallel, the
#' torsion becomes undefined; in such cases the function returns
#' \code{NA} and provides a diagnostic message.
#'
#' Optionally, the function can display a small segment of the curve
#' around the evaluation point using \pkg{plotly}. The point where the
#' curvature and torsion are computed is highlighted in the 3D plot.
#'
#' @param X,Y,Z Functions of \code{t} returning the three coordinate
#' components of the parametric curve.
#' @param t0 Value of the parameter at which curvature and torsion are
#' evaluated.
#' @param h Step size for the finite-difference approximations.
#' Smaller values give more accuracy but may amplify numerical noise.
#' @param plot Logical; if \code{TRUE}, displays a 3D plot of a short
#' segment of the curve around the evaluation point.
#' @param window Length of the parameter interval shown when
#' \code{plot = TRUE}. The interval is centered at \code{t0}.
#' @param n_samples Number of points used to draw the curve segment in
#' the 3D plot.
#' @param line A list defining the visual style of the curve in the
#' 3D plot.
#' @param point A list defining the visual style of the marker placed
#' at the evaluation point.
#' @param scene A list with 3D axis settings for \pkg{plotly}.
#' @param bg Background colors for the \pkg{plotly} figure.
#' @param tol Numeric tolerance used to detect degenerate situations in
#' which curvature or torsion cannot be reliably computed.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{kappa}: numerical curvature at the evaluation point.
#'   \item \code{tau}: numerical torsion at the evaluation point, or
#'         \code{NA} if the computation is unstable.
#'   \item \code{t0}: the parameter value where the evaluation was made.
#'   \item \code{point}: a numeric vector containing the coordinates of
#'         the curve at \code{t0}.
#'   \item \code{r1}, \code{r2}, \code{r3}: numeric approximations to the
#'         first, second and third derivative vectors at \code{t0}.
#' }
#'
#' @examples
#' # Example curve
#' X <- function(t) t^2 * cos(t)
#' Y <- function(t) t^3 * sin(3 * t)
#' Z <- function(t) t
#' res <- curvature_torsion3d(X, Y, Z, t0 = pi)
#' res$kappa
#' res$tau
#'
#' # \donttest{ if (requireNamespace("plotly", quietly = TRUE)) {
#' #   curvature_torsion3d(
#' #     X, Y, Z, t0 = pi, plot = TRUE,
#' #     window = 1.0, n_samples = 200,
#' #     line  = list(color = "red",  width = 2),
#' #     point = list(color = "black", size = 5, symbol = "circle")
#' #   )
#' # } }
#'
#' @importFrom tibble tibble
#' @export
curvature_torsion3d <- function(
    X, Y, Z, t0,
    h = 1e-4,
    plot = FALSE,
    window = 1.0,
    n_samples = 200,
    line  = list(color = "red",   width = 2, dash = "solid"),
    point = list(color = "black", size  = 5, symbol = "circle"),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    bg = list(paper = "white", plot = "white"),
    tol = 1e-10
) {
  # --------- minimal validation ----------
  if (!is.function(X) || !is.function(Y) || !is.function(Z)) {
    stop("X, Y, Z must be functions of one variable t.", call. = FALSE)
  }
  if (!is.numeric(t0) || length(t0) != 1L || !is.finite(t0)) {
    stop("t0 must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(h) || h <= 0) {
    stop("h must be positive.", call. = FALSE)
  }

  # --------- numeric derivatives (centered) ----------
  d1 <- function(f, t, h) (f(t + h) - f(t - h)) / (2 * h)
  d2 <- function(f, t, h) (f(t + h) - 2 * f(t) + f(t - h)) / (h * h)
  # Third derivative, centered finite difference (step h):
  # f'''(t) \u2248 (f(t-2h) - 2f(t-h) + 2f(t+h) - f(t+2h)) / (2 h^3)
  d3 <- function(f, t, h) (f(t - 2 * h) - 2 * f(t - h) +
                             2 * f(t + h) - f(t + 2 * h)) / (2 * h^3)

  r  <- function(t) c(X(t), Y(t), Z(t))
  r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
  r2 <- c(d2(X, t0, h), d2(Y, t0, h), d2(Z, t0, h))
  r3 <- c(d3(X, t0, h), d3(Y, t0, h), d3(Z, t0, h))

  # --------- vector algebra ----------
  dot  <- function(a, b) sum(a * b)
  norm <- function(a) sqrt(dot(a, a))
  cross <- function(a, b) c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )

  r1_norm <- norm(r1)
  if (r1_norm < tol) {
    stop("||r'(t0)|| \u2248 0. Curvature/torsion undefined at t0 (possible singular point).", call. = FALSE)
  }

  c12   <- cross(r1, r2)
  n_c12 <- norm(c12)

  kappa <- n_c12 / (r1_norm^3)

  # torsion: triple product / ||r' x r''||^2
  # det([r1 r2 r3]) = r1 \u00B7 (r2 \u00D7 r3)
  triple <- dot(r1, cross(r2, r3))
  tau <- if (n_c12 < tol) {
    warning("||r'(t0) \u00D7 r''(t0)|| \u2248 0. Torsion undefined at t0; returning NA.", call. = FALSE)
    NA_real_
  } else {
    triple / (n_c12^2)
  }

  # --------- optional plot ----------
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      a <- t0 - window / 2
      b <- t0 + window / 2

      ts <- seq(a, b, length.out = n_samples)
      data <- tibble::tibble(
        t = ts,
        x = vapply(ts, X, numeric(1)),
        y = vapply(ts, Y, numeric(1)),
        z = vapply(ts, Z, numeric(1))
      )
      p0 <- r(t0)

      plt <- plotly::plot_ly(
        data = data,
        x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "lines",
        line = line
      ) |>
        plotly::add_trace(
          x = p0[1], y = p0[2], z = p0[3],
          type = "scatter3d", mode = "markers",
          marker = point,
          showlegend = FALSE,
          hoverinfo = "none"
        ) |>
        plotly::layout(
          title = paste0(
            "\u03BA \u2248 ", signif(kappa, 6),
            "   \u00B7   \u03C4 \u2248 ",
            ifelse(is.na(tau), "NA", signif(tau, 6))
          ),
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      print(plt)
    }
  }

  list(
    kappa = kappa,
    tau   = tau,
    t0    = t0,
    point = r(t0),
    r1 = r1, r2 = r2, r3 = r3
  )
}
