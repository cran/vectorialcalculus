#' Newton iteration for a square nonlinear system
#' @param f Deterministic function accepting one coordinate vector and returning
#'   a finite real vector of the same length (at least two). Output dimension
#'   and optional unique/nonempty names must remain fixed. Wrap f(x,y)
#'   explicitly as `function(v)f(v[1],v[2])`; no signature guessing is performed.
#' @param x0 Finite initial coordinate vector of length at least two. Optional
#'   unique/nonempty names are passed to f and jacobian.
#' @param jacobian Optional function accepting the same vector and returning
#'   a finite real square matrix in output-row/input-column order. Entries are
#'   checked, but not certified to be derivatives; matrix names never reorder
#'   entries. NULL reuses vc_jacobian with complete sensitivity diagnostics.
#' @param h NULL or positive scalar/per-coordinate first-difference steps.
#' @param method One difference method or one per coordinate: central, forward
#'   or backward. Used only for the numerical Jacobian.
#' @param max_iter Nonnegative integer maximum accepted updates.
#' @param f_tol Nonnegative absolute infinity-norm residual tolerance.
#' @param x_tol Nonnegative relative step tolerance, using maximum coordinate
#'   change divided by max(1,max(abs(next_point))). A small step alone is not
#'   convergence.
#' @param jacobian_tol Nonnegative minimum reciprocal one-norm condition
#'   estimate. Matrices at or below this threshold are rejected before solving.
#' @param backtrack Logical; TRUE halves a Newton step until its finite
#'   residual satisfies residual_next <= (1-1e-4*damping)*residual_current,
#'   or residual_next <= f_tol. This safeguard does not guarantee convergence.
#'   FALSE accepts any finite trial, recovering undamped Newton iteration.
#' @param max_backtrack Nonnegative integer maximum step halvings per update.
#'   The full step is always attempted first.
#' @param atol,rtol Nonnegative absolute/relative sensitivity tolerances passed
#'   to vc_jacobian; failed numerical entry flags terminate conservatively.
#' @return Plain list with type, x0, estimate (last accepted point), value,
#'   converged, status, points and values (rows aligned with iterations),
#'   iterations (n, residual, step_norm, relative_step, damping, jacobian_rcond),
#'   jacobian_diagnostics, trials (all accepted/rejected attempted steps with
#'   point, value, residual, error and damping), control and note. A zero-update
#'   or failed search retains the initial row. Status is converged, max_iter,
#'   jacobian_unavailable, jacobian_unreliable, jacobian_singular,
#'   linear_solve_failed, nonfinite_step, invalid_trial, line_search_failed,
#'   stagnation or cycle_detected. Input/initial-field errors raise errors;
#'   failures after initialization return diagnostics without erasing history.
#' @details Solves J*s=-f by a linear solve, never by an explicit inverse.
#'   The initial residual is checked before asking for a Jacobian. Later
#'   numerical convergence requires residual <= f_tol and a small relative
#'   step, or an exactly zero computed residual. A small step with excessive
#'   residual is stagnation. Repeated accepted points within x_tol trigger a
#'   cycle diagnostic. Rejected trials never replace the last valid iterate.
#'
#'   This recovers original newton_raphson2d's square-system recurrence and
#'   iterate records, separates graphics, generalizes dimension and retains
#'   conditioning, stencil and domain failures. Exact-Newton local theory
#'   requires appropriate differentiability, a nonsingular true Jacobian and
#'   a sufficiently close initial point. Arbitrary starts, damping, numerical
#'   derivatives, multiple roots and singular systems have no promised global
#'   convergence or quadratic rate. Finite residual tests do not prove a root,
#'   regularity, uniqueness or that a supplied Jacobian is correct. Fixed
#'   tolerances are scale-dependent; rescale poorly conditioned systems.
#' @examples
#' f <- function(v)c(v[1]^2+v[2]^2-1,v[1]-v[2])
#' f |> vc_newton_system(c(.8,.2))
#' vc_newton_system(function(v)v,c(0,0))$status
#' @export
vc_newton_system <- function(f,x0,jacobian=NULL,h=NULL,method="central",
    max_iter=50L,f_tol=1e-8,x_tol=1e-8,jacobian_tol=1e-12,
    backtrack=TRUE,max_backtrack=20L,atol=1e-8,rtol=1e-4) {
  x0 <- .vc_point(x0,"x0");n <- length(x0)
  if(n<2) stop("Use vc_newton for one-variable equations.")
  if(!is.function(f)) stop("'f' must be a function.")
  if(!is.null(jacobian)&&!is.function(jacobian)) stop("Invalid 'jacobian'.")
  .vc_matrix_step(x0,h)
  if(!is.character(method)||!length(method)%in%c(1L,n)||anyNA(method)||
     any(!method%in%c("central","forward","backward"))) stop("Invalid 'method'.")
  .vc_integer(max_iter,"max_iter",0);.vc_integer(max_backtrack,"max_backtrack",0)
  .vc_positive(f_tol,"f_tol",TRUE);.vc_positive(x_tol,"x_tol",TRUE)
  .vc_positive(jacobian_tol,"jacobian_tol",TRUE)
  .vc_matrix_tolerances(atol,rtol,FALSE);.vc_flag(backtrack,"backtrack")
  initial <- .vc_vector_field_value(f,x0)
  if(length(initial)!=n) stop("A square system requires one output per input coordinate.")
  evaluate <- function(p) {
    z <- .vc_vector_field_value(f,p)
    if(length(z)!=n||!identical(names(z),names(initial))) stop("System output dimension/names changed.")
    z
  }
  point <- x0;value <- initial;residual <- max(abs(value))
  points <- matrix(point,nrow=1,dimnames=list(NULL,names(point)))
  values <- matrix(value,nrow=1,dimnames=list(NULL,names(value)))
  iterations <- data.frame(n=0L,residual=residual,step_norm=NA_real_,relative_step=NA_real_,
    damping=NA_real_,jacobian_rcond=NA_real_)
  diagnostics <- trials <- list()
  converged <- residual<=f_tol
  status <- if(converged)"converged" else "max_iter"
  if(!converged)for(k in seq_len(max_iter)) {
    d <- tryCatch(suppressWarnings(if(is.null(jacobian))
      vc_jacobian(evaluate,point,method,h,atol,rtol,FALSE) else {
        A <- jacobian(point)
        if(!is.matrix(A)||!is.numeric(A)||is.complex(A)||!identical(dim(A),c(n,n))||any(!is.finite(A)))
          stop("Supplied Jacobian must be a finite real n-by-n matrix.")
        list(estimate=A,all_reliable=NA)
      }),error=function(e)list(error=conditionMessage(e)))
    diagnostics[[k]] <- list(iteration=k-1L,source=if(is.null(jacobian))"numerical" else "supplied",result=d)
    if(!is.null(d$error)) {status <- "jacobian_unavailable";break}
    A <- d$estimate
    if(is.null(jacobian)&&!d$all_reliable) {status <- "jacobian_unreliable";break}
    condition <- tryCatch(rcond(A),error=function(e)0)
    diagnostics[[k]]$rcond <- condition
    if(!is.finite(condition)||condition<=jacobian_tol) {status <- "jacobian_singular";break}
    step <- tryCatch(as.double(solve(A,-value,tol=0)),error=function(e)numeric())
    if(length(step)!=n||any(!is.finite(step))) {status <- "linear_solve_failed";break}
    linear_residual <- max(abs(as.double(A%*%step)+value))
    diagnostics[[k]]$linear_residual <- linear_residual
    if(!is.finite(linear_residual)) {status <- "linear_solve_failed";break}
    accepted <- FALSE
    for(j in seq.int(0L,if(backtrack)max_backtrack else 0L)) {
      damping <- 2^(-j)
      next_point <- point+damping*step
      error <- "";next_value <- NULL;next_residual <- Inf
      if(any(!is.finite(next_point))) error <- "Nonfinite trial coordinates."
      else if(identical(as.double(next_point),as.double(point))) error <- "Trial step collapsed."
      else {
        next_value <- tryCatch(suppressWarnings(evaluate(next_point)),error=function(e) {
          error <<- conditionMessage(e);NULL
        })
        if(!is.null(next_value))next_residual <- max(abs(next_value))
      }
      accepted <- !nzchar(error)&&(!backtrack||next_residual<=f_tol||
        next_residual<=(1-1e-4*damping)*residual)
      trials[[length(trials)+1L]] <- list(iteration=k,halvings=j,damping=damping,
        point=next_point,value=next_value,residual=next_residual,accepted=accepted,error=error)
      if(accepted)break
    }
    if(!accepted) {
      status <- if(any(!is.finite(next_point)))"nonfinite_step" else
        if(identical(as.double(next_point),as.double(point)))"stagnation" else
        if(backtrack)"line_search_failed" else "invalid_trial"
      break
    }
    step_norm <- max(abs(next_point-point))
    relative_step <- max(abs(next_point/max(1,max(abs(next_point)))-point/max(1,max(abs(next_point)))))
    if(!is.finite(step_norm))step_norm <- Inf
    points <- rbind(points,next_point);values <- rbind(values,next_value)
    iterations <- rbind(iterations,data.frame(n=k,residual=next_residual,step_norm=step_norm,
      relative_step=relative_step,damping=damping,jacobian_rcond=condition))
    point <- next_point;value <- next_value;residual <- next_residual
    if(residual==0||(residual<=f_tol&&relative_step<=x_tol)) {
      converged <- TRUE;status <- "converged";break
    }
    if(relative_step<=x_tol) {status <- "stagnation";break}
    scale <- max(1,max(abs(point)))
    previous <- points[seq_len(k),,drop=FALSE]
    if(any(apply(previous,1,function(p)max(abs(p/scale-point/scale)))<=x_tol)) {
      status <- "cycle_detected";break
    }
  }
  rownames(points) <- rownames(values) <- NULL
  list(type="newton_system",x0=x0,estimate=point,value=value,converged=converged,
    status=status,points=points,values=values,iterations=iterations,
    jacobian_diagnostics=diagnostics,trials=trials,
    control=list(max_iter=max_iter,f_tol=f_tol,x_tol=x_tol,jacobian_tol=jacobian_tol,
      backtrack=backtrack,max_backtrack=max_backtrack,h=h,method=method,atol=atol,rtol=rtol),
    note="Finite stopping tests do not prove a root, smoothness, uniqueness or convergence from arbitrary starts.")
}
