#' Total differential of a scalar field in R^n
#'
#' Computes the gradient of a scalar field \code{f(x)} at a point \code{x0}
#' using central finite differences. It also returns the total differential,
#' understood as the linear map \code{v -> grad f(x0) %*% v}.
#'
#' The function \code{f} must take a single numeric vector \code{x} and
#' return a single numeric value.
#'
#' @param f Function of one argument \code{x} (numeric vector) returning a
#'   numeric scalar.
#' @param x0 Numeric vector giving the evaluation point.
#' @param h Step size for finite differences. Can be \code{NULL} (automatic),
#'   a scalar, or a vector of the same length as \code{x0}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{point}: the numeric vector \code{x0},
#'   \item \code{value}: the numeric value \code{f(x0)},
#'   \item \code{gradient}: numeric vector with the gradient at \code{x0},
#'   \item \code{differential}: a function \code{d(v)} that returns
#'         \code{sum(gradient * v)}.
#' }
#'
#' @examples
#' f <- function(x) x[1]^2 + 3 * x[2]^2
#' out <- total_differential_nd(f, c(1, 2))
#' out$gradient
#' out$differential(c(1, 0))  # directional derivative in direction (1,0)
#'
#' @export
total_differential_nd <- function(f, x0, h = NULL) {

  if (!is.function(f)) {
    stop("'f' must be a function of one numeric vector.", call. = FALSE)
  }
  if (!is.numeric(x0) || length(x0) < 1L) {
    stop("'x0' must be a non-empty numeric vector.", call. = FALSE)
  }
  if (any(!is.finite(x0))) {
    stop("'x0' must contain only finite values.", call. = FALSE)
  }

  x0 <- as.numeric(x0)
  n  <- length(x0)

  # step sizes
  if (is.null(h)) {
    h_vec <- 1e-4 * (1 + abs(x0))
  } else {
    h_vec <- as.numeric(h)
    if (length(h_vec) == 1L) {
      h_vec <- rep(h_vec, n)
    } else if (length(h_vec) != n) {
      stop("'h' must be scalar or have the same length as 'x0'.", call. = FALSE)
    }
  }
  if (any(!is.finite(h_vec)) || any(h_vec <= 0)) {
    stop("'h' must be strictly positive and finite.", call. = FALSE)
  }

  f0 <- f(x0)
  if (!is.numeric(f0) || length(f0) != 1L || !is.finite(f0)) {
    stop("'f(x0)' must be a finite numeric scalar.", call. = FALSE)
  }

  # gradient (central finite differences)
  grad <- numeric(n)
  for (i in seq_len(n)) {
    e <- numeric(n)
    e[i] <- h_vec[i]
    fp <- f(x0 + e)
    fm <- f(x0 - e)

    if (!is.numeric(fp) || length(fp) != 1L || !is.finite(fp) ||
        !is.numeric(fm) || length(fm) != 1L || !is.finite(fm)) {
      stop("Values 'f(x0 +/- h e_i)' must be finite numeric scalars.", call. = FALSE)
    }

    grad[i] <- (fp - fm) / (2 * h_vec[i])
  }

  differential_fun <- function(v) {
    v <- as.numeric(v)
    if (length(v) != n) {
      stop("Direction 'v' must have the same length as 'x0'.", call. = FALSE)
    }
    if (any(!is.finite(v))) {
      stop("Direction 'v' must contain only finite values.", call. = FALSE)
    }
    sum(grad * v)
  }

  list(
    point       = x0,
    value       = f0,
    gradient    = grad,
    differential = differential_fun
  )
}
