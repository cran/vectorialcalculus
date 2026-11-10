#' Newton-Raphson iteration with explicit termination diagnostics
#' @inheritParams vc_tangent
#' @param df Optional function returning a finite real derivative scalar.
#'   NULL uses vc_derivative(), retaining its diagnostics.
#' @param max_iter Nonnegative integer maximum number of accepted updates.
#' @param f_tol Nonnegative absolute numerical residual tolerance.
#' @param x_tol Nonnegative relative step tolerance, scaled by max(1,abs(x_next)).
#' @param derivative_tol Nonnegative absolute safety threshold for slope size.
#' @return Plain list with type, x0, estimate (last accepted iterate), converged,
#'   status, iterations (n, x, fx, dfx, step, residual, derivative_status),
#'   derivative_diagnostics, rejected_trial (or NULL), control (input numerical
#'   settings), and note. status is
#'   converged, max_iter, derivative_too_small, derivative_unreliable,
#'   nonfinite_derivative, nonfinite_evaluation, nonfinite_step, stagnation
#'   or cycle_detected. Unaccepted trials never replace the last valid iterate.
#' @details Uses x_next=x-f(x)/df(x), without damping or fallback. The initial
#'   residual is checked before evaluating the slope. Later acceptance as
#'   converged requires residual <=f_tol and a relative step <=x_tol, or an
#'   exactly zero computed residual. A small step with a large residual is
#'   stagnation. Repeated iterates within x_tol flag a cycle, not convergence.
#'   These finite numerical tests are not a proof of a mathematical root.
#'   Near a simple root, sufficient smoothness and a sufficiently close start
#'   give the usual local Newton theory; arbitrary starts and multiple roots
#'   have no such general guarantee. Finite-difference derivatives also do not
#'   promise the quadratic rate of exact-derivative Newton theory. Small slopes, scale, rounding, restricted
#'   domains and poor derivative estimates can prevent convergence. Numerical
#'   derivative diagnostics that fail terminate the iteration conservatively.
#'   A supplied df is checked for finite scalar output, not independently
#'   certified to be the derivative of f. Warnings during numerical derivative
#'   probes and trial evaluations are suppressed; invalid results are retained
#'   as failure diagnostics rather than treated as usable values.
#' @examples
#' f <- function(x) x^2-2
#' f |> vc_newton(x0=1)
#' vc_newton(function(x) x^2, 0)$status
#' vc_newton(function(x) x^3-2*x+2, 0)$status
#' @export
vc_newton <- function(f,x0,df=NULL,h=NULL,max_iter=50L,f_tol=1e-8,
                      x_tol=1e-8,derivative_tol=1e-12,
                      method=c("central","forward","backward")) {
  if (!is.null(df) && !is.function(df)) stop("'df' must be a function or NULL.")
  .vc_integer(max_iter,"max_iter",0)
  .vc_positive(f_tol,"f_tol",TRUE)
  .vc_positive(x_tol,"x_tol",TRUE)
  .vc_positive(derivative_tol,"derivative_tol",TRUE)
  .vc_optional_step(h,"h")
  method <- match.arg(method)
  fx <- .vc_value(f,x0)
  iterations <- data.frame(n=0L,x=as.double(x0),fx=fx,dfx=NA_real_,step=NA_real_,
                           residual=abs(fx),derivative_status="not_evaluated")
  diagnostics <- list()
  rejected <- NULL
  converged <- abs(fx)<=f_tol
  status <- if (converged) "converged" else "max_iter"
  if (!converged) for(k in seq_len(max_iter)) {
    x <- iterations$x[k]
    fx <- iterations$fx[k]
    d <- if (is.null(df)) tryCatch(suppressWarnings(vc_derivative(f,x,method=method,h=h,warn=FALSE)),error=function(e) NULL) else NULL
    slope <- if (!is.null(df)) .vc_safe_value(df,x) else if (is.null(d)) NA_real_ else d$estimate
    diagnostics[[k]] <- if (!is.null(df)) list(source="supplied",estimate=slope) else
      if (is.null(d)) list(source="numerical",error="invalid_stencil") else d
    iterations$dfx[k] <- slope
    iterations$derivative_status[k] <- if (!is.null(df)) "supplied" else
      if (!is.null(d) && d$reliable) "passed" else "failed"
    if (!is.finite(slope)) { status <- "nonfinite_derivative"; break }
    if (abs(slope)<=derivative_tol) { status <- "derivative_too_small"; break }
    if (is.null(df) && !d$reliable) { status <- "derivative_unreliable"; break }
    step <- -fx/slope
    next_x <- x+step
    if (!is.finite(step) || !is.finite(next_x)) {
      status <- "nonfinite_step"; rejected <- list(x=next_x,step=step); break
    }
    next_fx <- .vc_safe_value(f,next_x)
    if (!is.finite(next_fx)) {
      status <- "nonfinite_evaluation"; rejected <- list(x=next_x,step=step,value=next_fx); break
    }
    actual_step <- next_x-x
    iterations <- rbind(iterations,data.frame(n=k,x=next_x,fx=next_fx,dfx=NA_real_,
      step=actual_step,residual=abs(next_fx),derivative_status="not_evaluated"))
    small_step <- abs(actual_step)/max(1,abs(next_x))<=x_tol
    if (next_fx==0 || (abs(next_fx)<=f_tol && small_step)) {
      status <- "converged"; converged <- TRUE; break
    }
    if (small_step) { status <- "stagnation"; break }
    previous <- iterations$x[seq_len(k)]
    if (any(abs(previous-next_x)/pmax(1,abs(next_x))<=x_tol)) {
      status <- "cycle_detected"; break
    }
  }
  list(type="newton",x0=x0,estimate=iterations$x[nrow(iterations)],converged=converged,
       status=status,iterations=iterations,derivative_diagnostics=diagnostics,
       rejected_trial=rejected,
       control=list(max_iter=max_iter,f_tol=f_tol,x_tol=x_tol,derivative_tol=derivative_tol,method=method,h=h),
       note="Numerical stopping criteria do not prove a root or guarantee convergence from arbitrary starts.")
}
