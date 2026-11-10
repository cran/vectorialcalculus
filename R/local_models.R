# Internal scalar evaluation shared by the local models and plots.
.vc_scalar <- function(x) {
  is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x)
}

.vc_value <- function(f, x) {
  if (!is.function(f)) stop("'f' must be a function.")
  if (!.vc_scalar(x)) stop("Evaluation points must be finite real scalars.")
  value <- f(x)
  if (!.vc_scalar(value)) stop("'f' must return one finite real scalar.")
  as.double(value)
}

.vc_finite <- function(value, label) {
  if (!.vc_scalar(value)) stop(label, " is not finite; change the scale or inputs.")
  as.double(value)
}

#' Compute the secant through two points of a real function
#'
#' Compute the finite slope (f(x1)-f(x0))/(x1-x0). This is a secant slope,
#' not a derivative or evidence of differentiability.
#' @param f Deterministic function accepting a real scalar and returning one
#'   finite real scalar. Vectorization is not required.
#' @param x0 Finite real scalar giving the first abscissa.
#' @param x1 Finite real scalar giving a distinct second abscissa.
#' @return Plain list with type ("secant"), x0, y0, x1, y1, dx, dy, slope,
#'   and note. The line is y0 + slope * (x - x0). No intercept is computed,
#'   avoiding unnecessary cancellation in an intercept representation.
#' @details Both function values, the displacement, the ordinate difference,
#'   and the quotient must be finite. Reversing the endpoints preserves the
#'   slope. No smoothness assumption is needed. Floating-point cancellation
#'   can affect the result; this operation does not estimate numerical error.
#' @examples
#' vc_secant(function(x) x^2, 1, 2)
#' sin |> vc_secant(x0 = 0, x1 = 0.1)
#' @export
vc_secant <- function(f, x0, x1) {
  y0 <- .vc_value(f, x0)
  y1 <- .vc_value(f, x1)
  if (x0 == x1) stop("'x0' and 'x1' must be distinct.")
  dx <- .vc_finite(as.double(x1) - as.double(x0), "Abscissa displacement")
  dy <- .vc_finite(y1 - y0, "Ordinate displacement")
  slope <- .vc_finite(dy / dx, "Secant slope")
  list(type = "secant", x0 = x0, y0 = y0, x1 = x1, y1 = y1,
       dx = dx, dy = dy, slope = slope,
       note = "A secant slope is not a derivative or a differentiability test.")
}

#' Construct an approximate tangent line
#'
#' Use vc_derivative() to construct the local line y0 + slope * (x - x0).
#' The numerical line is not proof that a mathematical tangent exists.
#' @inheritParams vc_secant
#' @inheritParams vc_derivative
#' @param x0 Finite real scalar evaluation point of the local model.
#' @param h Positive finite real initial derivative step, or NULL for the
#'   heuristic used by vc_derivative(), evaluated with x = x0.
#' @return Plain list with type ("tangent"), x0, y0, slope, derivative
#'   (the full vc_derivative result), and note. Inspect derivative$reliable
#'   and its diagnostics before interpreting the slope. Failed diagnostics
#'   retain the approximation, with a warning unless warn = FALSE.
#' @details A two-sided tangent requires differentiability at x0. One-sided
#'   methods instead construct a one-sided approximation, which may differ
#'   from a two-sided tangent. All limitations of vc_derivative apply.
#' @examples
#' vc_tangent(function(x) x^2, 1)
#' sin |> vc_tangent(x0 = 0)
#' vc_tangent(abs, 0, warn = FALSE)$derivative$reliable
#' @export
vc_tangent <- function(f, x0, method = c("central", "forward", "backward"),
                       h = NULL, atol = 1e-8, rtol = 1e-4, warn = TRUE) {
  derivative <- vc_derivative(f, x = x0, method = match.arg(method),
                              h = h, atol = atol, rtol = rtol, warn = warn)
  list(type = "tangent", x0 = x0, y0 = .vc_value(f, x0),
       slope = derivative$estimate, derivative = derivative,
       note = derivative$note)
}

#' Compare a differential with a finite change
#'
#' Approximate dy = f'(x0) * dx and compare it with the actual sampled
#' change delta_y = f(x0 + dx) - f(x0). The differential is a local linear
#' approximation; it is generally different from the finite change.
#' @inheritParams vc_tangent
#' @param dx Finite real scalar displacement. Zero and negative values are
#'   allowed. For nonzero dx, x0 + dx must be finite and distinct from x0.
#' @return Plain list with type ("differential"), x0, y0, x1, y1, dx, dy,
#'   delta_y, remainder (delta_y - dy), tangent (full vc_tangent result),
#'   and note. dx is the requested displacement; actual_dx is the representable
#'   x1 - x0, also returned. The difference can matter at floating-point scale.
#' @details In exact mathematics, differentiability implies
#'   delta_y = f'(x0)*dx + o(abs(dx)) as dx tends to zero. Here the slope
#'   is numerical: remainder includes derivative approximation and rounding
#'   errors, and is not a rigorous Taylor error bound. The function must be
#'   finite at x0 + dx, as well as at the derivative stencil. One-sided
#'   derivative choices do not certify a two-sided differential.
#' @examples
#' vc_differential(function(x) x^2, 1, dx = 0.1)
#' sin |> vc_differential(x0 = 0, dx = -0.1)
#' @export
vc_differential <- function(f, x0, dx,
                            method = c("central", "forward", "backward"),
                            h = NULL, atol = 1e-8, rtol = 1e-4, warn = TRUE) {
  if (!.vc_scalar(dx)) stop("'dx' must be a finite real scalar.")
  tangent <- vc_tangent(f, x0, method = match.arg(method), h = h,
                        atol = atol, rtol = rtol, warn = warn)
  x1 <- .vc_finite(as.double(x0) + as.double(dx), "Endpoint")
  if (dx != 0 && x1 == x0) stop("Displacement collapsed; change 'dx'.")
  actual_dx <- .vc_finite(x1 - as.double(x0), "Actual displacement")
  y1 <- .vc_value(f, x1)
  dy <- .vc_finite(tangent$slope * dx, "Differential")
  delta_y <- .vc_finite(y1 - tangent$y0, "Finite change")
  remainder <- .vc_finite(delta_y - dy, "Remainder")
  list(type = "differential", x0 = x0, y0 = tangent$y0, x1 = x1, y1 = y1,
       dx = dx, actual_dx = actual_dx, dy = dy, delta_y = delta_y,
       remainder = remainder, tangent = tangent, note = tangent$note)
}
