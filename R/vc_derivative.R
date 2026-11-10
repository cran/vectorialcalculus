#' Approximate a derivative with finite differences
#'
#' Evaluate finite difference quotients at three decreasing step sizes.
#' These are numerical approximations, never a proof of differentiability.
#'
#' @param f Function accepting a real scalar and returning one finite real scalar.
#'   It should be deterministic. Vectorization is not required.
#' @param x Finite real scalar at which to approximate the derivative.
#' @param method One of "central", "forward", or "backward".
#' @param h Positive finite real initial step. NULL uses
#'   max(1, abs(x)) times machine epsilon to the power 1/3 for central
#'   differences or 1/2 for one-sided differences. This heuristic is not
#'   optimal for every function or physical scale.
#' @param atol Nonnegative finite absolute tolerance for sensitivity comparisons.
#' @param rtol Nonnegative finite relative tolerance for sensitivity comparisons.
#' @param warn Logical scalar; warn when sensitivity diagnostics fail.
#'
#' @details For positive s the forward quotient is (f(x+s)-f(x))/s,
#'   the backward quotient is (f(x)-f(x-s))/s, and the central quotient is
#'   (f(x+s)-f(x-s))/(2*s). The actual representable displacements are used
#'   as denominators to account for floating-point rounding of x+s and x-s.
#'   For sufficiently smooth functions these have nominal truncation orders
#'   one, one, and two, respectively, on an equally spaced stencil.
#'
#'   Steps h, h/2 and h/4 are evaluated. The last approximation is returned.
#'   Stability requires both successive changes to be at most
#'   atol + rtol * max(abs(the two estimates)). For central differences,
#'   the finest left and right quotients must also agree by this criterion.
#'   Failure produces a diagnostic and, by default, a warning. Reported
#'   changes are sensitivity measures, not rigorous error bounds.
#'
#'   Central differences require finite values on both sides, and one-sided
#'   methods require only their selected side and x. Out-of-domain values,
#'   collapsed or overflowing stencils, and nonfinite quotients are errors.
#'   One-sided results need not equal a two-sided derivative, especially at
#'   domain boundaries. Stable results can occur at discontinuities, cusps,
#'   or unresolved oscillations; even agreement of all checks is not proof.
#'   Cancellation, underflow, and large additive offsets can also conceal error.
#'
#' @return A plain list with estimate, x, method, h, a diagnostics data frame
#'   (step, forward, backward, estimate; unused sides are NA), changes
#'   (absolute successive changes), stable, side_agreement (NA for one-sided
#'   methods), reliable (the conjunction of applicable checks, not a
#'   mathematical guarantee), and a character note explaining that limitation.
#' @examples
#' vc_derivative(function(x) x^3, 2)$estimate
#' sin |> vc_derivative(x = 0, method = "central")
#' vc_derivative(abs, 0, warn = FALSE)$reliable
#' vc_derivative(sqrt, 1, method = "forward")$estimate
#' @export
vc_derivative <- function(f, x, method = c("central", "forward", "backward"),
                          h = NULL, atol = 1e-8, rtol = 1e-4, warn = TRUE) {
  scalar <- function(z) {
    is.numeric(z) && !is.complex(z) && length(z) == 1L && is.finite(z)
  }
  if (!is.function(f)) stop("'f' must be a function.")
  if (!scalar(x)) stop("'x' must be a finite real scalar.")
  method <- match.arg(method)
  if (!scalar(atol) || atol < 0) stop("'atol' must be nonnegative and finite.")
  if (!scalar(rtol) || rtol < 0) stop("'rtol' must be nonnegative and finite.")
  if (!is.logical(warn) || length(warn) != 1L || is.na(warn)) {
    stop("'warn' must be TRUE or FALSE.")
  }
  if (is.null(h)) {
    h <- max(1, abs(x)) * .Machine$double.eps^
      (if (method == "central") 1/3 else 1/2)
  }
  if (!scalar(h) || h <= 0) stop("'h' must be a positive finite real scalar.")
  evaluate <- function(z) {
    value <- f(z)
    if (!scalar(value)) stop("'f' must return one finite real scalar at every stencil point.")
    as.double(value)
  }
  f0 <- evaluate(x)
  steps <- h / c(1, 2, 4)
  if (any(steps == 0)) stop("Step underflow: choose a larger 'h'.")
  forward <- backward <- rep(NA_real_, 3L)
  for (i in seq_along(steps)) {
    if (method != "backward") {
      right <- x + steps[i]
      if (!is.finite(right) || right <= x) stop("Right stencil collapsed or overflowed; change 'h'.")
      forward[i] <- (evaluate(right) - f0) / (right - x)
    }
    if (method != "forward") {
      left <- x - steps[i]
      if (!is.finite(left) || left >= x) stop("Left stencil collapsed or overflowed; change 'h'.")
      backward[i] <- (f0 - evaluate(left)) / (x - left)
    }
  }
  estimates <- switch(method, forward = forward, backward = backward,
                      central = forward / 2 + backward / 2)
  if (any(!is.finite(estimates)) ||
      any(!is.finite(forward[!is.na(forward)])) ||
      any(!is.finite(backward[!is.na(backward)]))) {
    stop("Nonfinite difference quotient; change function scale or 'h'.")
  }
  close <- function(a, b) {
    # Scale before subtraction to avoid overflow in comparisons.
    scale <- max(abs(a), abs(b), atol, .Machine$double.xmin)
    abs(a / scale - b / scale) <= atol / scale + rtol * max(abs(a / scale), abs(b / scale))
  }
  stable <- close(estimates[1], estimates[2]) && close(estimates[2], estimates[3])
  side_agreement <- if (method == "central") close(forward[3], backward[3]) else NA
  reliable <- stable && (method != "central" || side_agreement)
  note <- paste("Sensitivity checks are numerical evidence only;",
                "they do not establish differentiability or bound the error.")
  if (warn && !reliable) warning("Numerical sensitivity checks failed. ", note, call. = FALSE)
  list(estimate = estimates[3], x = x, method = method, h = h,
       diagnostics = data.frame(step = steps, forward = forward,
                                backward = backward, estimate = estimates),
       changes = abs(diff(estimates)), stable = stable,
       side_agreement = side_agreement, reliable = reliable, note = note)
}
