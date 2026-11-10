.vc_integral_interval <- function(interval) {
  if(!is.numeric(interval)||is.complex(interval)||length(interval)!=2L||
     any(!is.finite(interval))) stop("'interval' must contain two finite real bounds.")
  interval <- as.double(interval)
  if(!is.finite(interval[2]-interval[1])) stop("Interval width must be finite.")
  interval
}

.vc_quadrature <- function(f,bounds,n,method) {
  p <- vc_partition(bounds,n)
  x <- p$breaks
  if(any(p$cells$width<=0)) stop("Unrepresentable quadrature grid.")
  y <- vapply(x,function(z).vc_value(f,z),numeric(1))
  if(method=="trapezoid") {
    w <- c(p$cells$width[1]/2,
           if(n>1) (p$cells$width[-n]+p$cells$width[-1])/2 else numeric(),
           p$cells$width[n]/2)
  } else {
    # Simpson's rule assumes an equally spaced grid, up to rounding.
    h <- (bounds[2]-bounds[1])/n
    if(any(abs(p$cells$width-h)>1e-6*h)) stop("Simpson grid is not numerically uniform.")
    w <- c(1,rep(c(4,2),length.out=n-1),1)*h/3
  }
  terms <- y*w
  if(any(!is.finite(terms))) stop("Nonfinite quadrature contribution.")
  list(estimate=.vc_finite(sum(terms),"Integral estimate"),
       table=data.frame(x=x,value=y,weight=w,contribution=terms))
}

#' Approximate a signed definite integral
#' @inheritParams vc_riemann
#' @param interval Two finite real bounds, allowing reversed or equal bounds.
#' @param method Quadrature rule: trapezoid, Simpson, or adaptive QUADPACK
#'   through stats::integrate(). Improper infinite bounds are not supported.
#' @param n Initial number of equal cells for fixed rules; positive integer,
#'   even and at least two for Simpson. Ignored by the adaptive rule.
#' @param abs_tol,rel_tol Nonnegative finite absolute and relative diagnostic
#'   tolerances; at least one must be positive. Adaptive integration also
#'   applies the tolerance restrictions of stats::integrate().
#' @param subdivisions Positive integer maximum adaptive subinterval count;
#'   ignored by fixed rules.
#' @return Plain list with type, interval, method, estimate, error_estimate,
#'   tolerance, tolerance_met, status, message, n (finest fixed cell count or
#'   NA), table (finest fixed nodes, values, signed weights/contributions),
#'   refinement (coarse/fine estimates), subdivisions and note. Adaptive
#'   internal nodes are unavailable and table/refinement are NULL.
#' @details Fixed rules evaluate n and 2*n cells and return the finer result.
#'   Richardson diagnostics use abs(fine-coarse)/3 for trapezoids and /15 for
#'   Simpson. They assume sufficiently smooth functions and an asymptotic
#'   error regime; aliasing or missed narrow features can make them misleading.
#'   tolerance_met is only a numerical diagnostic, never a certified bound.
#'   Adaptive status preserves QUADPACK failure messages without claiming
#'   convergence. Evaluation errors and nonfinite sampled values/results stop.
#'   Scalar-only functions are supported. Reversed bounds negate the integral;
#'   equal bounds return zero without evaluating f. Integrals are signed, not
#'   geometric areas. No method proves integrability or discovers every
#'   singularity; split known discontinuities/singularities explicitly.
#' @examples
#' (function(x)x^2) |> vc_integrate(c(0,2),method="simpson")
#' vc_integrate(sin,c(0,pi),method="adaptive")
#' @export
vc_integrate <- function(f,interval,method=c("trapezoid","simpson","adaptive"),
                         n=32L,abs_tol=1e-8,rel_tol=1e-8,subdivisions=100L) {
  if(!is.function(f)) stop("'f' must be a function.")
  interval <- .vc_integral_interval(interval)
  method <- match.arg(method)
  .vc_positive(abs_tol,"abs_tol",TRUE)
  .vc_positive(rel_tol,"rel_tol",TRUE)
  if(abs_tol==0 && rel_tol==0) stop("At least one tolerance must be positive.")
  if(method=="adaptive") {
    .vc_integer(subdivisions,"subdivisions",1)
    if(abs_tol==0 && rel_tol<max(50*.Machine$double.eps,5e-29))
      stop("Adaptive relative tolerance is too small when abs_tol is zero.")
  } else {
    .vc_integer(n,"n",if(method=="simpson")2 else 1)
    if(n>.Machine$integer.max/2) stop("'n' is too large for refinement.")
    if(method=="simpson" && n%%2!=0) stop("Simpson requires an even 'n'.")
  }
  note <- paste("Signed integral; error estimates are not certified bounds.",
                "Finite sampling can miss oscillations, singularities or narrow features.")
  out <- list(type="integral",interval=interval,method=method,estimate=0,
    error_estimate=0,tolerance=abs_tol,tolerance_met=TRUE,status="zero_interval",
    message="Equal bounds: f was not evaluated.",n=NA_integer_,table=NULL,
    refinement=NULL,subdivisions=NA_integer_,note=note)
  if(interval[1]==interval[2]) return(out)
  bounds <- sort(interval)
  direction <- sign(interval[2]-interval[1])
  if(method=="adaptive") {
    z <- stats::integrate(function(x)vapply(x,function(t).vc_value(f,t),numeric(1)),
       bounds[1],bounds[2],abs.tol=abs_tol,rel.tol=rel_tol,
       subdivisions=subdivisions,stop.on.error=FALSE)
    out$estimate <- .vc_finite(direction*z$value,"Integral estimate")
    out$error_estimate <- .vc_finite(z$abs.error,"Error estimate")
    out$message <- z$message
    out$subdivisions <- z$subdivisions
    out$status <- if(identical(z$message,"OK")) "adaptive_ok" else "adaptive_failure"
  } else {
    coarse <- .vc_quadrature(f,bounds,n,method)
    fine <- .vc_quadrature(f,bounds,2*n,method)
    delta <- .vc_finite(fine$estimate-coarse$estimate,"Refinement difference")
    out$estimate <- direction*fine$estimate
    out$error_estimate <- abs(delta)/if(method=="trapezoid")3 else 15
    out$n <- 2*n
    out$table <- fine$table
    out$table$weight <- direction*out$table$weight
    out$table$contribution <- direction*out$table$contribution
    out$refinement <- data.frame(n=c(n,2*n),estimate=direction*c(coarse$estimate,fine$estimate))
    out$status <- "fixed_rule_diagnostic"
    out$message <- "Returned refined grid; error estimate assumes an asymptotic smooth regime."
  }
  out$tolerance <- .vc_finite(max(abs_tol,rel_tol*abs(out$estimate)),"Tolerance")
  out$tolerance_met <- out$error_estimate<=out$tolerance && out$status!="adaptive_failure"
  out
}

#' Compare numerical quadrature rules
#' @inheritParams vc_integrate
#' @param methods Nonempty unique vector of quadrature method names.
#' @param reference Optional finite analytical/reference integral.
#' @return Plain list with type, interval, results (complete integration
#'   objects), table (method, estimate, error_estimate, tolerance_met, status,
#'   error, abs_error), reference and note. Actual reference errors are NA
#'   when reference is absent. No best method is asserted.
#' @examples
#' vc_integral_compare(sin,c(0,pi),reference=2)$table
#' @export
vc_integral_compare <- function(f,interval,methods=c("trapezoid","simpson","adaptive"),
                                n=32L,abs_tol=1e-8,rel_tol=1e-8,
                                subdivisions=100L,reference=NULL) {
  if(!is.character(methods)||!length(methods)||anyNA(methods)||anyDuplicated(methods)||
     any(!methods %in% c("trapezoid","simpson","adaptive"))) stop("Invalid 'methods'.")
  if(!is.null(reference)&&!.vc_scalar(reference)) stop("Invalid 'reference'.")
  results <- lapply(methods,function(m)vc_integrate(f,interval,m,n,abs_tol,rel_tol,subdivisions))
  names(results) <- methods
  estimates <- vapply(results,function(z)z$estimate,numeric(1))
  error <- if(is.null(reference)) rep(NA_real_,length(methods)) else estimates-reference
  if(any(!is.finite(error[!is.na(error)]))) stop("Reference error overflowed.")
  list(type="integral_comparison",interval=results[[1]]$interval,results=results,
       table=data.frame(method=methods,estimate=unname(estimates),
         error_estimate=vapply(results,function(z)z$error_estimate,numeric(1)),
         tolerance_met=vapply(results,function(z)z$tolerance_met,logical(1)),
         status=vapply(results,function(z)z$status,character(1)),error=unname(error),abs_error=unname(abs(error))),
       reference=reference,note="Estimated errors and reference errors are distinct; no certified ranking.")
}
