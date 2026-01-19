#' Sample a 3D parametric curve
#'
#' Generates a tibble with columns \code{t, x, y, z} by evaluating
#' the parametric curve \eqn{(X(t), Y(t), Z(t))} on the interval
#' \eqn{[a, b]} at a given number of sample points.
#'
#' @param X,Y,Z Functions of one variable \code{t}, e.g. \code{function(t) 2 * cos(t)}.
#' @param a,b Numeric parameter limits for \code{t}.
#' @param n_samples Integer. Number of sample points along the curve.
#'
#' @return
#' A tibble with columns \code{t}, \code{x}, \code{y}, \code{z},
#' where \code{x = X(t)}, \code{y = Y(t)}, \code{z = Z(t)}.
#'
#' @seealso [plot_curve3d()], [arc_length3d()]
#'
#' @examples
#' X <- function(t) 2 * cos(t)
#' Y <- function(t) 3 * sin(t)
#' Z <- function(t) t / 5
#' curve_sample3d(X, Y, Z, 0, 2 * pi, n_samples = 100)
#'
#' @importFrom tibble tibble
#' @export
curve_sample3d <- function(X, Y, Z, a, b, n_samples = 400) {
  ts <- seq(a, b, length.out = n_samples)
  tibble::tibble(
    t = ts,
    x = vapply(ts, X, numeric(1)),
    y = vapply(ts, Y, numeric(1)),
    z = vapply(ts, Z, numeric(1))
  )
}
