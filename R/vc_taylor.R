#' Numerical Taylor approximations of order one or two
#'
#' Construct P(x) = f(x0) + d1*(x-x0) + d2*(x-x0)^2/2 (order two),
#' or omit the quadratic term for order one. Derivatives are numerical.
#' @inheritParams vc_tangent
#' @param order Integer 1 or 2. Higher orders are not implemented.
#' @param x Optional nonempty finite real vector at which to compare f and P.
#'   NULL leaves comparison empty; the returned polynomial is still evaluable.
#' @param h2 Positive finite initial step for the second derivative. NULL uses
#'   max(1, abs(x0)) * .Machine$double.eps^(1/4). Ignored for order one,
#'   but a supplied value must still be valid. h controls only the first derivative.
#' @details The first derivative reuses vc_derivative() with its selected
#'   method and diagnostics. Order two uses a direct three-point stencil,
#'   never repeated numerical differentiation. With actual representable
#'   displacements a = x0-left and b = right-x0, the second derivative is
#'   2 * ((f(right)-f(x0))/b - (f(x0)-f(left))/a) / (a+b).
#'   Equal spacing gives the usual centered second difference, with nominal
#'   second-order truncation for sufficiently smooth functions. Asymmetric
#'   rounded stencils can reduce that order.
#'
#'   Second-derivative steps h2, h2/2, h2/4 are compared using the same
#'   absolute/relative tolerance rule as vc_derivative(). An additional
#'   roundoff indicator is 4*epsilon*max(abs(stencil values))/min(a,b)^2.
#'   The finest indicator must not exceed atol + rtol*abs(d2). This is a
#'   conservative floating-point scale heuristic, not a rigorous error bound.
#'   The finest estimate is retained even if checks fail. Order two always
#'   requires finite values on both sides, even with a one-sided first derivative.
#'
#'   A genuine Taylor interpretation assumes the required derivatives exist
#'   near x0. Numerical checks do not establish this assumption. Comparison
#'   errors are sampled differences f(x)-P(x), not certified remainder bounds.
#'   Optional comparison points with failed/nonfinite function evaluations have
#'   value and error NA, with explicit status; evaluation warnings are suppressed.
#'   Nonfinite polynomial values are likewise NA in the comparison table.
#'   Large offsets, cancellation and unresolved oscillations can evade checks.
#'   Checks are conservative: the existing first-derivative side comparison
#'   can flag smooth stationary points, and the second-derivative roundoff
#'   indicator can flag a nonzero constant with exactly zero curvature.
#'   Such results are retained with failed diagnostics, not declared invalid.
#' @return Plain list with type, x0, y0, order, coefficients (in powers of
#'   x-x0, including the constant), polynomial (a vectorized function),
#'   derivative (full vc_derivative result), second_derivative (diagnostics
#'   list or NULL), comparison (data frame or NULL), reliable and note.
#'   The comparison columns are x, distance, value, approximation, error,
#'   abs_error and status. reliable means only that applicable numerical checks
#'   passed. The polynomial rejects invalid inputs and overflowing evaluations.
#'   second_derivative contains estimate, h, diagnostics (step, left_step,
#'   right_step, estimate, roundoff_indicator), changes, stable, resolved,
#'   and reliable. stable compares successive estimates; resolved checks the
#'   finest roundoff indicator; reliable is their conjunction.
#' @examples
#' a <- vc_taylor(exp, 0, order = 2, x = c(-0.2, 0, 0.2))
#' a$comparison
#' a$polynomial(c(0, 0.1))
#' sin |> vc_taylor(x0 = 0, order = 1, x = c(0, 0.2))
#' @export
vc_taylor <- function(f, x0, order = 2L, x = NULL,
                      method = c("central", "forward", "backward"),
                      h = NULL, h2 = NULL, atol = 1e-8, rtol = 1e-4, warn = TRUE) {
  if (!.vc_scalar(order) || !(order %in% c(1,2))) stop("'order' must be 1 or 2.")
  if (!is.null(h2) && (!.vc_scalar(h2) || h2 <= 0)) stop("'h2' must be positive and finite.")
  if (!is.null(x)) .vc_taylor_points(x)
  # Validate warn/tolerances via the approved first-derivative engine.
  d1 <- vc_derivative(f, x0, method = match.arg(method), h = h,
                       atol = atol, rtol = rtol, warn = warn)
  y0 <- .vc_value(f, x0)
  d2 <- NULL
  coefficients <- c(constant = y0, linear = d1$estimate)
  if (order == 2) {
    if (is.null(h2)) h2 <- max(1, abs(x0)) * .Machine$double.eps^(1/4)
    d2 <- .vc_second_difference(f, x0, h2, atol, rtol)
    if (warn && !d2$reliable) warning(
      "Second-derivative sensitivity or roundoff checks failed; no accuracy guarantee.", call. = FALSE)
    coefficients <- c(coefficients, quadratic = d2$estimate / 2)
  }
  polynomial <- .vc_taylor_polynomial(coefficients, x0)
  reliable <- d1$reliable && (order == 1 || d2$reliable)
  list(type = "taylor", x0 = x0, y0 = y0, order = as.integer(order),
       coefficients = coefficients, polynomial = polynomial, derivative = d1,
       second_derivative = d2,
       comparison = if (!is.null(x)) .vc_taylor_compare(f, polynomial, x, x0),
       reliable = reliable,
       note = "Numerical diagnostics do not prove smoothness or bound the Taylor remainder.")
}

.vc_taylor_points <- function(x) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || !length(x) || any(!is.finite(x))) {
    stop("Evaluation points must be a nonempty finite real vector.")
  }
  invisible(NULL)
}

.vc_taylor_polynomial <- function(coefficients, center) {
  force(coefficients)
  force(center)
  function(x) {
    .vc_taylor_points(x)
    z <- as.double(x) - center
    value <- rep(unname(coefficients[length(coefficients)]), length(z))
    for (i in rev(seq_len(length(coefficients)-1L))) value <- value*z + coefficients[i]
    if (any(!is.finite(z)) || any(!is.finite(value))) stop("Polynomial evaluation overflowed.")
    as.double(value)
  }
}

.vc_taylor_compare <- function(f, polynomial, x, x0) {
  .vc_taylor_points(x)
  approximation <- vapply(x, function(z) tryCatch(polynomial(z), error = function(e) NA_real_), numeric(1))
  value <- vapply(x, function(z) {
    y <- tryCatch(suppressWarnings(f(z)), error = function(e) NA_real_)
    if (.vc_scalar(y)) as.double(y) else NA_real_
  }, numeric(1))
  error <- value - approximation
  status <- rep("ok", length(x))
  status[!is.finite(value)] <- "function_not_finite"
  status[!is.finite(approximation)] <- "polynomial_not_finite"
  status[is.finite(value) & is.finite(approximation) & !is.finite(error)] <- "error_overflow"
  error[!is.finite(error)] <- NA_real_
  data.frame(x = as.double(x), distance = abs(as.double(x)-x0), value = value,
             approximation = approximation, error = error, abs_error = abs(error), status = status)
}

.vc_second_difference <- function(f, x0, h2, atol, rtol) {
  steps <- h2 / c(1,2,4)
  if (any(steps == 0)) stop("Second-derivative step underflow; increase 'h2'.")
  y0 <- .vc_value(f, x0)
  estimates <- indicators <- left_step <- right_step <- numeric(3L)
  for (i in seq_along(steps)) {
    left <- x0-steps[i]
    right <- x0+steps[i]
    if (!is.finite(left) || !is.finite(right) || left >= x0 || right <= x0)
      stop("Second-derivative stencil collapsed or overflowed; change 'h2'.")
    a <- x0-left
    b <- right-x0
    yl <- .vc_value(f, left)
    yr <- .vc_value(f, right)
    # a/2+b/2 avoids overflow in a+b. These are raw function differences,
    # not calls to the first-derivative engine.
    estimates[i] <- ((yr-y0)/b - (y0-yl)/a) / (a/2+b/2)
    if (!is.finite(estimates[i])) stop("Nonfinite second difference; change the scale or 'h2'.")
    indicators[i] <- (4*.Machine$double.eps*max(abs(c(yl,y0,yr))) / min(a,b)) / min(a,b)
    if (is.nan(indicators[i])) indicators[i] <- Inf
    left_step[i] <- a
    right_step[i] <- b
  }
  close <- function(a,b) {
    scale <- max(abs(a),abs(b),atol,.Machine$double.xmin)
    abs(a/scale-b/scale) <= atol/scale + rtol*max(abs(a/scale),abs(b/scale))
  }
  stable <- close(estimates[1], estimates[2]) && close(estimates[2], estimates[3])
  scale <- max(abs(estimates[3]), atol, .Machine$double.xmin)
  resolved <- is.finite(indicators[3]) &&
    indicators[3]/scale <= atol/scale + rtol*abs(estimates[3]/scale)
  list(estimate = estimates[3], h = h2,
       diagnostics = data.frame(step = steps, left_step = left_step,
                                right_step = right_step, estimate = estimates,
                                roundoff_indicator = indicators),
       changes = abs(diff(estimates)), stable = stable, resolved = resolved,
       reliable = stable && resolved)
}
