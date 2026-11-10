.vc_accumulation_x <- function(x) {
  if(!is.numeric(x)||is.complex(x)||!is.null(dim(x))||!length(x)||any(!is.finite(x)))
    stop("'x' must be a nonempty finite real vector.")
  as.double(x)
}

#' Evaluate a signed accumulation function
#' @inheritParams vc_integrate
#' @param x Nonempty finite real vector of upper limits. Input order and
#'   duplicates are retained; limits below lower reverse orientation.
#' @param lower Finite real scalar fixed lower integration limit.
#' @return Plain list with type, lower, method, table (x, estimate,
#'   error_estimate, tolerance_met, status), results (full vc_integrate
#'   objects, one per x) and note.
#' @details Evaluates A(x) = integral from lower to x of f(t) dt independently
#'   at each upper limit; it does not accumulate rounding error along a grid.
#'   Equal limits use the integral engine's zero-interval convention. Each
#'   integral retains its diagnostics; a returned finite estimate need not
#'   satisfy the numerical tolerance. These are signed integrals, not geometric
#'   areas. Finite quadrature cannot prove integrability. Split known
#'   singularities/discontinuities as appropriate; infinite bounds are absent.
#' @examples
#' (function(x)x^2) |> vc_accumulation(x=c(-1,0,1,2),lower=0)
#' @export
vc_accumulation <- function(f,x,lower=0,method=c("adaptive","trapezoid","simpson"),
                            n=32L,abs_tol=1e-8,rel_tol=1e-8,subdivisions=100L) {
  x <- .vc_accumulation_x(x)
  .vc_finite(lower,"lower")
  method <- match.arg(method)
  results <- lapply(x,function(z)vc_integrate(f,c(lower,z),method,n,abs_tol,rel_tol,subdivisions))
  list(type="accumulation",lower=as.double(lower),method=method,
       table=data.frame(x=x,estimate=vapply(results,function(z)z$estimate,numeric(1)),
         error_estimate=vapply(results,function(z)z$error_estimate,numeric(1)),
         tolerance_met=vapply(results,function(z)z$tolerance_met,logical(1)),
         status=vapply(results,function(z)z$status,character(1))),
       results=results,note="Signed accumulation; quadrature diagnostics are not certified bounds.")
}

#' Explore the fundamental theorem of calculus numerically
#' @inheritParams vc_accumulation
#' @param h Positive finite initial difference step, or NULL for the
#'   vc_derivative() scale heuristic. Smaller is not always more accurate.
#' @param derivative_method Central, forward or backward differences applied
#'   to the numerical accumulation function. One-sided comparisons are
#'   one-sided evidence only and can be used near domain boundaries.
#' @param atol,rtol Nonnegative finite sensitivity/agreement tolerances.
#'   Agreement threshold is atol + rtol*abs(f(x)).
#' @return Plain list with type, lower, method, derivative_method, table
#'   (x, value, estimate, difference, agreement, derivative_reliable,
#'   integration_ok, quadrature_indicator, quadrature_resolved, status, message),
#'   results (derivative diagnostics, full sampled integrals and propagated
#'   quadrature diagnostics per point), and note. Unavailable estimates and
#'   differences are NA. Status is consistent, mismatch or inconclusive.
#' @details Reuses vc_integrate() and vc_derivative(), not an analytical
#'   substitution of f(x) for A'(x). At each of the three difference scales,
#'   quadrature error estimates are multiplied by the absolute finite-difference
#'   coefficients and summed. This indicator accounts for amplification by
#'   small steps, but remains heuristic because its input errors are estimates.
#'   A consistent/mismatch status requires reliable derivative diagnostics,
#'   all sampled integrals meeting their diagnostic tolerances, and propagated
#'   indicators no larger than the agreement threshold at all three scales.
#'   Otherwise status is inconclusive, even if the computed values agree.
#'   Integration/derivative evaluation errors around a finite f(x) are retained
#'   as inconclusive messages. Nonfinite f(x) and invalid controls are errors.
#'   Integrability on the intervening interval and continuity at x are standard
#'   hypotheses for A'(x)=f(x) at an interior point. They are not established by
#'   this procedure. A point discontinuity can produce a numerical mismatch;
#'   a symmetric quotient alone does not prove a derivative at a jump.
#' @examples
#' (function(x)x^2) |> vc_fundamental(x=c(0,1,2),lower=0)
#' @export
vc_fundamental <- function(f,x,lower=0,method=c("adaptive","trapezoid","simpson"),
                           n=32L,abs_tol=1e-8,rel_tol=1e-8,subdivisions=100L,
                           h=NULL,derivative_method=c("central","forward","backward"),
                           atol=1e-8,rtol=1e-4) {
  x <- .vc_accumulation_x(x)
  .vc_finite(lower,"lower")
  method <- match.arg(method)
  derivative_method <- match.arg(derivative_method)
  .vc_optional_step(h,"h")
  .vc_positive(atol,"atol",TRUE)
  .vc_positive(rtol,"rtol",TRUE)
  # Validate quadrature controls before retaining point-specific failures.
  vc_integrate(f,c(lower,lower),method,n,abs_tol,rel_tol,subdivisions)
  points <- lapply(x,function(z) {
    value <- .vc_value(f,z)
    threshold <- .vc_finite(atol+rtol*abs(value),"Agreement threshold")
    locations <- numeric()
    integrals <- list()
    sample <- function(t) {
      i <- match(t,locations)
      if(is.na(i)) {
        q <- vc_integrate(f,c(lower,t),method,n,abs_tol,rel_tol,subdivisions)
        locations <<- c(locations,t)
        integrals[[length(integrals)+1L]] <<- q
        q
      } else integrals[[i]]
    }
    row <- data.frame(x=z,value=value,estimate=NA_real_,difference=NA_real_,agreement=FALSE,
      derivative_reliable=FALSE,integration_ok=FALSE,quadrature_indicator=NA_real_,
      quadrature_resolved=FALSE,status="inconclusive",message="")
    d <- tryCatch(vc_derivative(function(t)sample(t)$estimate,z,
      method=derivative_method,h=h,atol=atol,rtol=rtol,warn=FALSE),error=function(e)e)
    diagnostic <- NULL
    if(inherits(d,"error")) {
      row$message <- conditionMessage(d)
      d <- NULL
    } else {
      steps <- d$diagnostics$step
      indicators <- vapply(steps,function(s) {
        q0 <- sample(z)$error_estimate
        if(derivative_method=="forward") {
          (sample(z+s)$error_estimate+q0)/((z+s)-z)
        } else if(derivative_method=="backward") {
          (sample(z-s)$error_estimate+q0)/(z-(z-s))
        } else {
          dr <- (z+s)-z
          dl <- z-(z-s)
          sample(z+s)$error_estimate/(2*dr)+sample(z-s)$error_estimate/(2*dl)+
            abs(.5/dl-.5/dr)*q0
        }
      },numeric(1))
      # Overflow in an uncertainty indicator is unresolved, never success.
      indicators[!is.finite(indicators)] <- Inf
      diagnostic <- data.frame(step=steps,quadrature_indicator=indicators)
      row$estimate <- d$estimate
      row$difference <- .vc_finite(d$estimate-value,"Derivative comparison")
      row$agreement <- abs(row$difference)<=threshold
      row$derivative_reliable <- d$reliable
      row$integration_ok <- all(vapply(integrals,function(q)q$tolerance_met,logical(1)))
      row$quadrature_indicator <- max(indicators)
      row$quadrature_resolved <- all(is.finite(indicators)) && all(indicators<=threshold)
      if(d$reliable && row$integration_ok && row$quadrature_resolved)
        row$status <- if(row$agreement)"consistent" else "mismatch"
      row$message <- if(row$status=="inconclusive")
        "Derivative sensitivity or quadrature diagnostics unresolved." else
        "Numerical comparison only; theorem hypotheses were not verified."
    }
    list(row=row,derivative=d,locations=locations,integrals=integrals,diagnostics=diagnostic)
  })
  list(type="fundamental",lower=as.double(lower),method=method,derivative_method=derivative_method,
    table=do.call(rbind,lapply(points,function(p)p$row)),results=points,
    note="Numerical evidence only; continuity, differentiability and theorem hypotheses are not proved.")
}
