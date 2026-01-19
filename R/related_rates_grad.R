#' Related rates via the gradient (implicit constraint)
#'
#' Computes related rates for an implicit constraint using the gradient.
#' Let
#' \deqn{
#'   g(\mathbf{x}) = 0,\quad \mathbf{x}=\mathbf{x}(t)\in\mathbb{R}^k.
#' }
#' Differentiating with respect to time yields
#' \deqn{
#'   \nabla g(\mathbf{x}(t))\cdot \mathbf{x}'(t) = 0.
#' }
#' If all components of \eqn{\mathbf{x}'(t)} are known except one, the missing
#' rate is determined by this orthogonality condition (velocity tangent to the
#' level set).
#'
#' @param g Function. A scalar function g(x1, x2, ..., xk) defining the implicit
#'   constraint g = 0. It must accept k numeric arguments and return a numeric
#'   scalar.
#' @param x Numeric vector of length k. Point where rates are evaluated.
#'   This point should satisfy g(x) = 0 (approximately).
#' @param known_rates Numeric vector of known rates. It can be:
#'   (1) a named numeric vector with names matching var_names (e.g. "x","y"),
#'   (2) a named numeric vector with names matching paste0("d", var_names)
#'       (e.g. "dx","dy"), or
#'   (3) an unnamed numeric vector of length k with NA for the unknown component.
#' @param solve_for Integer or character. Which rate to solve for. If integer,
#'   it is the position in var_names (1..k). If character, it may be either a
#'   variable name (e.g. "y") or a rate name (e.g. "dy").
#' @param var_names Character vector of length k. Variable names. If NULL,
#'   defaults to c("x1","x2",...,"xk").
#' @param h Numeric scalar. Step size for numerical partial derivatives.
#'
#' @return A list with components:
#' \describe{
#' \item{rate}{Numeric scalar. The solved rate (the requested component of x'(t)).}
#' \item{rates}{Named numeric vector length k. Full rate vector x'(t). Names are
#'   paste0("d", var_names).}
#' \item{grad}{Numeric vector length k. Gradient of g at x.}
#' \item{dot}{Numeric scalar. Dot product grad . rates (should be near 0).}
#' \item{gx}{Numeric scalar. g(x) value (should be near 0 if x is on the constraint).}
#' }
#'
#' @examples
#' # Ladder (implicit circle): x^2 + y^2 = L^2
#' # Suppose L = 5, at the instant (x,y) = (4,3) and dx/dt = -1.
#' # Then dy/dt = -(x/y) dx/dt = (4/3).
#' g <- function(x, y) x^2 + y^2 - 25
#'
#' out <- related_rates_grad(
#'   g = g,
#'   x = c(4, 3),
#'   known_rates = c(dx = -1, dy = NA),
#'   solve_for = "dy",
#'   var_names = c("x", "y")
#' )
#' out$rate
#' out$rates
#'
#' @export
related_rates_grad <- function(
    g,
    x,
    known_rates,
    solve_for,
    var_names = NULL,
    h = 1e-6
) {

  if (!is.function(g)) stop("'g' must be a function.")
  if (!is.numeric(x) || length(x) < 1L) stop("'x' must be a numeric vector.")
  if (any(!is.finite(x))) stop("'x' must be finite numeric values.")
  k <- length(x)

  if (is.null(var_names)) {
    var_names <- paste0("x", seq_len(k))
  }
  if (!is.character(var_names) || length(var_names) != k) {
    stop("'var_names' must be a character vector of length length(x).")
  }

  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("'h' must be a positive finite numeric scalar.")
  }

  rate_names <- paste0("d", var_names)

  # Normalize known_rates into a numeric vector length k with NA for unknown(s)
  rates <- rep(NA_real_, k)
  names(rates) <- rate_names

  if (is.numeric(known_rates) && length(known_rates) == k && is.null(names(known_rates))) {
    rates[] <- known_rates
  } else if (is.numeric(known_rates) && !is.null(names(known_rates))) {

    nm <- names(known_rates)

    # Accept either var names ("x","y") or rate names ("dx","dy")
    idx_var <- match(nm, var_names)
    idx_rate <- match(nm, rate_names)

    if (all(!is.na(idx_var))) {
      rates[idx_var] <- known_rates
    } else if (all(!is.na(idx_rate))) {
      rates[idx_rate] <- known_rates
    } else {
      stop("Names in 'known_rates' must match 'var_names' (e.g. 'x') or paste0('d', var_names) (e.g. 'dx').")
    }

  } else {
    stop("'known_rates' must be a numeric vector (named recommended).")
  }

  # Determine index to solve for
  j <- NA_integer_
  if (is.character(solve_for) && length(solve_for) == 1L) {
    # accept "y" or "dy"
    if (solve_for %in% var_names) {
      j <- match(solve_for, var_names)
    } else if (solve_for %in% rate_names) {
      j <- match(solve_for, rate_names)
    } else {
      j <- NA_integer_
    }
  } else if (is.numeric(solve_for) && length(solve_for) == 1L) {
    j <- as.integer(solve_for)
  } else {
    stop("'solve_for' must be a single integer or a single character name.")
  }
  if (is.na(j) || j < 1L || j > k) stop("'solve_for' must match one variable (e.g. 'y') or one rate (e.g. 'dy').")

  # Evaluate g(x)
  gx <- do.call(g, as.list(x))
  if (!is.numeric(gx) || length(gx) != 1L || !is.finite(gx)) {
    stop("'g' must return a finite numeric scalar at 'x'.")
  }

  # Numerical gradient by central differences
  grad <- numeric(k)
  for (i in seq_len(k)) {
    xp <- x; xm <- x
    xp[i] <- xp[i] + h
    xm[i] <- xm[i] - h
    gp <- do.call(g, as.list(xp))
    gm <- do.call(g, as.list(xm))
    if (!is.numeric(gp) || !is.numeric(gm) || length(gp) != 1L || length(gm) != 1L ||
        !is.finite(gp) || !is.finite(gm)) {
      stop("Non-finite values encountered while computing numerical gradient.")
    }
    grad[i] <- (gp - gm) / (2*h)
  }
  names(grad) <- var_names

  # Need exactly one unknown: the one we solve for
  na_pos <- which(is.na(rates))
  if (length(na_pos) == 0L) {
    dot_val <- sum(grad * unname(rates))
    return(list(rate = rates[j], rates = rates, grad = grad, dot = dot_val, gx = gx))
  }
  if (length(na_pos) > 1L || na_pos[1] != j) {
    stop("Provide all rates except the one indicated by 'solve_for'.")
  }

  gj <- grad[j]
  if (!is.finite(gj) || gj == 0) {
    stop("Cannot solve: partial derivative for the requested variable is zero (or non-finite) at 'x'.")
  }

  # Solve grad . rates = 0 for rates[j]
  sum_known <- sum(grad[-j] * rates[-j], na.rm = TRUE)
  rates[j] <- -sum_known / gj

  dot_val <- sum(grad * rates)

  list(rate = rates[j], rates = rates, grad = grad, dot = dot_val, gx = gx)
}

