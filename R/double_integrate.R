#' Integrate over a Cartesian region with nested fixed quadrature
#' @param f Scalar field taking one numeric vector in x,y order.
#' @param interval Finite outer bounds, possibly reversed or equal.
#' @param lower,upper Functions giving finite ordered inner bounds.
#' @param order Inner then outer coordinate: yx (type I) or xy (type II).
#' @param method Simpson or trapezoid rule.
#' @param n_outer,n_inner Initial cell counts, refined at both levels.
#' @param abs_tol,rel_tol Nonnegative diagnostic tolerances, one positive.
#' @return List with estimate, error_estimate, tolerance_met, full outer and
#'   cached inner vc_integrate results, sections, order, interval and status.
#' @details Reuses vc_integrate unchanged. Adds its outer error indicator to
#'   the sum of absolute finest outer weights times inner error indicators.
#'   This is heuristic, never a certified bound. All inner flags are retained.
#'   Finite sampling does not prove integrability, Fubini hypotheses or boundary
#'   ordering between samples. Split known singularities and complex regions.
#'   Outer reversal changes sign; sampled inner reversal is rejected. Equal
#'   bounds contribute zero without field evaluations. Signed integrals are not
#'   geometric volumes. Polar and adaptive nested quadrature are not included.
#' @examples
#' (function(v)v[1]+2*v[2]) |>
#'   vc_double_integrate(c(0,1),function(x)0,function(x)x)
#' @export
vc_double_integrate <- function(f,interval,lower,upper,order=c("yx","xy"),
                                method=c("simpson","trapezoid"),
                                n_outer=16L,n_inner=16L,abs_tol=1e-8,rel_tol=1e-8) {
  if(!is.function(f)||!is.function(lower)||!is.function(upper))
    stop("f, lower and upper must be functions.")
  order <- match.arg(order)
  method <- match.arg(method)
  interval <- .vc_integral_interval(interval)
  for(n in list(n_outer,n_inner)) {
    .vc_integer(n,"cell count",if(method=="simpson")2 else 1)
    if(n>.Machine$integer.max/2 || (method=="simpson" && n%%2!=0))
      stop("Invalid cell count for refinement/rule.")
  }
  nodes <- numeric()
  inner <- bounds <- list()
  slice <- function(t) {
    i <- match(t,nodes)
    if(!is.na(i)) return(inner[[i]]$estimate)
    a <- .vc_value(lower,t)
    b <- .vc_value(upper,t)
    if(a>b) stop("Inner lower bound exceeds upper bound at sampled section.")
    z <- vc_integrate(function(s).vc_field_value(f,if(order=="yx")c(t,s) else c(s,t)),
      c(a,b),method,n_inner,abs_tol,rel_tol)
    nodes <<- c(nodes,t)
    inner[[length(inner)+1L]] <<- z
    bounds[[length(bounds)+1L]] <<- c(a,b)
    z$estimate
  }
  outer <- vc_integrate(slice,interval,method,n_outer,abs_tol,rel_tol)
  propagated <- 0
  if(!is.null(outer$table)) {
    ids <- match(outer$table$x,nodes)
    propagated <- .vc_finite(sum(abs(outer$table$weight)*
      vapply(inner[ids],function(z)z$error_estimate,numeric(1))),"Inner error indicator")
  }
  error <- .vc_finite(outer$error_estimate+propagated,"Combined error indicator")
  sections <- data.frame(outer=nodes,lower=vapply(bounds,function(z)z[1],numeric(1)),
    upper=vapply(bounds,function(z)z[2],numeric(1)),
    estimate=vapply(inner,function(z)z$estimate,numeric(1)),
    error_estimate=vapply(inner,function(z)z$error_estimate,numeric(1)),
    tolerance_met=vapply(inner,function(z)z$tolerance_met,logical(1)))
  list(type="double_integral",estimate=outer$estimate,error_estimate=error,
    outer_error=outer$error_estimate,inner_error=propagated,tolerance=outer$tolerance,
    tolerance_met=outer$tolerance_met && all(sections$tolerance_met) && error<=outer$tolerance,
    outer=outer,inner=inner,sections=sections,order=order,interval=interval,method=method,
    status=if(length(nodes))"nested_fixed_diagnostic" else "zero_interval",
    note="Signed integral; nested error indicators are not certified bounds.")
}
