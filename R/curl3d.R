#' Numerical curl of a three-dimensional vector field
#'
#' @description
#' Computes the curl of a vector field in three dimensions at a given
#' point, using second-order central finite differences.
#' The field may optionally depend on a time parameter; if so, the curl
#' is evaluated at a fixed time value.
#'
#' @details
#' The vector field must be a function that returns a numeric vector of
#' length three representing the components of the field at a point.
#' The curl is obtained by approximating the partial derivatives of the
#' field components with respect to each coordinate direction using
#' symmetric finite differences.
#'
#' The step size for each coordinate can be:
#' \itemize{
#'   \item a single scalar used for all three axes,
#'   \item a numeric vector of length three providing separate steps for
#'         the x, y and z directions,
#'   \item or \code{NULL}, in which case automatic step sizes are chosen
#'         based on the evaluation point.
#' }
#'
#' The method currently implemented is the second-order central
#' differencing scheme. Smaller step sizes may provide more accurate
#' results for rapidly varying fields, at the cost of increased
#' sensitivity to floating-point error.
#'
#' @param field A function representing the vector field. It can be defined
#'        as \code{function(x, y, z)} or as \code{function(x, y, z, t)}.
#'        It must return a numeric vector \code{c(Fx, Fy, Fz)}.
#' @param x,y,z Numeric scalars giving the coordinates of the evaluation
#'        point.
#' @param h Step size for finite differences. It may be:
#'        \itemize{
#'          \item a single numeric value,
#'          \item a numeric vector of length three,
#'          \item or \code{NULL} (automatic selection).
#'        }
#' @param tval Time value used when the vector field depends on time.
#'        Default is \code{0}.
#' @param method Differencing scheme. Currently only \code{"central"} is
#'        supported.
#'
#' @return A named numeric vector of length three containing the curl
#'         components at the evaluation point. The components are named
#'         \code{omega_x}, \code{omega_y} and \code{omega_z}.
#'
#' @examples
#' # Simple rotating field: curl is constant in the third component
#' field1 <- function(x, y, z) c(-y, x, 0.6)
#' curl3d(field1, x = 0.1, y = -0.3, z = 2)
#'
#' # Time-dependent example (time does not affect the curl):
#' field2 <- function(x, y, z, t) c(-y, x + t, z)
#' curl3d(field2, x = 1, y = 2, z = 3, tval = 5)
#'
#' # Using a smaller step size for more precision:
#' curl3d(field1, x = 1, y = 1, z = 1, h = 1e-5)
#'
#' @export
curl3d <- function(field, x, y, z, h = NULL, tval = 0, method = c("central")) {
  method <- match.arg(method)
  if (!is.function(field)) {
    stop("'field' must be function(x,y,z) or function(x,y,z,t).", call. = FALSE)
  }
  if (!is.numeric(x) || !is.numeric(y) || !is.numeric(z) ||
      length(x) != 1L || length(y) != 1L || length(z) != 1L ||
      !is.finite(x) || !is.finite(y) || !is.finite(z)) {
    stop("'x','y','z' must be finite numeric scalars.", call. = FALSE)
  }

  # step sizes
  if (is.null(h)) {
    hx <- 1e-4 * (1 + abs(x))
    hy <- 1e-4 * (1 + abs(y))
    hz <- 1e-4 * (1 + abs(z))
  } else if (length(h) == 1L) {
    hx <- hy <- hz <- as.numeric(h)
  } else if (length(h) == 3L) {
    hx <- as.numeric(h[1])
    hy <- as.numeric(h[2])
    hz <- as.numeric(h[3])
  } else {
    stop("'h' must be a scalar, a length-3 vector, or NULL.", call. = FALSE)
  }
  if (any(!is.finite(c(hx, hy, hz))) || any(c(hx, hy, hz) <= 0)) {
    stop("Step sizes 'h' must be positive and finite.", call. = FALSE)
  }

  # wrapper: evaluate field with or without time
  field_eval <- function(xx, yy, zz, tt) {
    out <- if (length(formals(field)) >= 4L) field(xx, yy, zz, tt) else field(xx, yy, zz)
    if (!is.numeric(out) || length(out) != 3L || any(!is.finite(out))) {
      stop("field(x,y,z[,t]) must return a finite numeric vector of length 3.", call. = FALSE)
    }
    out
  }

  # central differences (2nd order)
  # dFz/dy
  Fp <- field_eval(x, y + hy, z, tval)[3]; Fm <- field_eval(x, y - hy, z, tval)[3]
  dFz_dy <- (Fp - Fm) / (2 * hy)

  # dFy/dz
  Fp <- field_eval(x, y, z + hz, tval)[2]; Fm <- field_eval(x, y, z - hz, tval)[2]
  dFy_dz <- (Fp - Fm) / (2 * hz)

  # dFx/dz
  Fp <- field_eval(x, y, z + hz, tval)[1]; Fm <- field_eval(x, y, z - hz, tval)[1]
  dFx_dz <- (Fp - Fm) / (2 * hz)

  # dFz/dx
  Fp <- field_eval(x + hx, y, z, tval)[3]; Fm <- field_eval(x - hx, y, z, tval)[3]
  dFz_dx <- (Fp - Fm) / (2 * hx)

  # dFy/dx
  Fp <- field_eval(x + hx, y, z, tval)[2]; Fm <- field_eval(x - hx, y, z, tval)[2]
  dFy_dx <- (Fp - Fm) / (2 * hx)

  # dFx/dy
  Fp <- field_eval(x, y + hy, z, tval)[1]; Fm <- field_eval(x, y - hy, z, tval)[1]
  dFx_dy <- (Fp - Fm) / (2 * hy)

  # curl
  omega_x <- dFz_dy - dFy_dz
  omega_y <- dFx_dz - dFz_dx
  omega_z <- dFy_dx - dFx_dy
  c(omega_x = omega_x, omega_y = omega_y, omega_z = omega_z)
}
