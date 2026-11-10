.vc_flag <- function(x,name) {
  if(!is.logical(x)||length(x)!=1L||is.na(x)) stop("Invalid '",name,"'.")
}

.vc_norm <- function(x) {
  s <- max(abs(x))
  if(s==0)0 else s*sqrt(sum((x/s)^2))
}

#' Assess a multivariable stationary candidate by the second derivative test
#' @inheritParams vc_gradient
#' @param assume_c2 Logical declaration by the caller that f is twice
#'   continuously differentiable in a neighborhood. Default FALSE prevents
#'   an extremum label; numerical stencils cannot verify this hypothesis.
#' @param hessian_h Positive scalar/vector second-difference step, or NULL
#'   for vc_hessian defaults. h controls first differences independently.
#' @param tol_grad Positive absolute tolerance for the gradient norm and the
#'   finest left/right coordinate quotients. This is not an error bound.
#' @param tol_eigen Positive absolute eigenvalue margin, enlarged by an
#'   empirical matrix refinement/roundoff allowance before classification.
#' @return Plain list with type, point, value, classification, reason,
#'   stationary_evidence, gradient (unchanged complete engine result),
#'   side_max, hessian (complete result or NULL), eigenvalues, eigen_margin,
#'   matrix_allowance, assume_c2 and note. Labels are minimum, maximum,
#'   saddle, inconclusive or not_stationary; all are numerical assessments.
#' @details Central coordinate quotients must be stable and small on BOTH
#'   finest sides. This additional stationarity check can pass at a smooth
#'   quadratic minimum when vc_gradient's stricter relative side-agreement
#'   flag fails; that original flag is retained, never rewritten. A cusp such
#'   as abs(x)+y^2 is not accepted just because central cancellation gives zero.
#'   Finite checks can still miss nonsmoothness or discontinuity.
#'
#'   With declared C2 regularity, a stationary candidate and reliable Hessian,
#'   positive/negative definite estimates suggest a strict minimum/maximum;
#'   eigenvalues of both signs suggest a saddle. Semidefinite, degenerate,
#'   unresolved or unavailable Hessians are inconclusive, never flat extrema.
#'   The margin is tol_eigen plus dimension times the maximum entry change
#'   over consecutive refinements plus dimension times the maximum finest
#'   roundoff indicator. It is a conservative diagnostic, NOT a certified
#'   perturbation bound. Near-stationary numerical points need not be exact
#'   critical points; all labels remain conditional. No boundary or constrained
#'   extremum test, higher-order test or global classification is performed.
#' @examples
#' f <- function(v)v[1]^2+2*v[2]^2
#' f |> vc_classify_stationary(c(0,0),assume_c2=TRUE)
#' @export
vc_classify_stationary <- function(f,point,assume_c2=FALSE,h=NULL,hessian_h=NULL,
    tol_grad=1e-5,tol_eigen=1e-6,atol=1e-8,rtol=1e-4,warn=TRUE) {
  point <- .vc_point(point)
  if(length(point)<2) stop("Use vc_classify_critical for a one-variable function.")
  .vc_flag(assume_c2,"assume_c2")
  .vc_matrix_tolerances(atol,rtol,warn)
  .vc_positive(tol_grad,"tol_grad");.vc_positive(tol_eigen,"tol_eigen")
  .vc_matrix_step(point,hessian_h,TRUE)
  g <- vc_gradient(f,point,h=h,atol=atol,rtol=rtol,warn=FALSE)
  side_max <- max(vapply(g$results,function(z)max(abs(c(
    z$derivative$diagnostics$forward[3],z$derivative$diagnostics$backward[3]))),numeric(1)))
  stable <- all(vapply(g$results,function(z)z$derivative$stable,logical(1)))
  stationary <- .vc_norm(g$estimate)<=tol_grad && side_max<=tol_grad && stable
  label <- "inconclusive"
  reason <- "C2 regularity has not been declared."
  H <- NULL;ev <- numeric();allowance <- margin <- NA_real_
  if(!stationary) {
    label <- if(.vc_norm(g$estimate)>tol_grad)"not_stationary" else "inconclusive"
    reason <- "Gradient residual or coordinate-side/stability checks failed."
  } else {
    H <- tryCatch(vc_hessian(f,point,h=hessian_h,atol=atol,rtol=rtol,warn=FALSE),
      error=function(e)list(error=conditionMessage(e)))
    if(!is.null(H$error)) reason <- paste("Hessian unavailable:",H$error)
    else {
      ev <- eigen(H$estimate,symmetric=TRUE,only.values=TRUE)$values
      changes <- max(abs(H$refinements[[1]]-H$refinements[[2]]),
                     abs(H$refinements[[2]]-H$refinements[[3]]))
      rounding <- max(vapply(H$results,function(z)z$derivative$diagnostics$roundoff_indicator[3],numeric(1)))
      allowance <- length(point)*(changes+rounding)
      margin <- tol_eigen+allowance
      if(!H$all_reliable||!is.finite(margin)||any(!is.finite(ev)))
        reason <- "Hessian sensitivity/roundoff or eigenvalue checks failed."
      else if(assume_c2) {
        label <- if(all(ev>margin))"minimum" else if(all(ev< -margin))"maximum" else
          if(any(ev>margin)&&any(ev< -margin))"saddle" else "inconclusive"
        reason <- if(label=="inconclusive")"Second derivative test is inconclusive." else
          "Conditional second derivative test at a numerical stationary candidate."
      }
    }
  }
  if(warn&&label=="inconclusive") warning(reason,call.=FALSE)
  list(type="stationary_classification",point=point,value=g$value,
    classification=label,reason=reason,stationary_evidence=stationary,
    gradient=g,side_max=side_max,hessian=H,eigenvalues=ev,eigen_margin=margin,
    matrix_allowance=allowance,assume_c2=assume_c2,
    note="Conditional numerical evidence, not a smoothness, existence or extremum certificate.")
}

#' Search for multivariable stationary candidates from explicit starts
#' @inheritParams vc_classify_stationary
#' @param starts One finite coordinate vector or a finite numeric matrix with
#'   one starting point per row, at least two columns. Optional column names
#'   are passed to f and must be unique/nonempty. No random starts are added.
#' @param bounds Optional finite increasing n-by-2 numeric search-window matrix.
#'   This limits optimization, not the domain of f: central stencils may extend
#'   outside the window and must still be evaluable. Boundary hits are rejected.
#' @param tol_merge Positive Euclidean distance below which accepted candidates
#'   are merged. Every run, including duplicates, remains in runs.
#' @param max_iter Positive integer optimization limit passed to stats::optim.
#' @return Plain list with type, starts, bounds, points (rows of merged
#'   candidates), table (id, value, residual, classification), assessments,
#'   runs (start, terminal point, optimizer code/message, accepted, candidate_id,
#'   status, assessment and every objective evaluation), and note. Empty
#'   searches return a zero-row matrix/table. No optimizer code alone accepts
#'   a candidate; final stationarity/side checks are mandatory.
#' @details Recovers the original critical_points_2d/critical_points_nd strategy
#'   of minimizing squared numerical gradient norm and merging duplicates.
#'   Reuses vc_gradient and vc_classify_stationary. Unbounded searches use
#'   Nelder-Mead; finite windows use L-BFGS-B. Hessians are never obtained by
#'   repeated differentiation of gradients. No random state is read or changed.
#'   Failed function/stencil evaluations receive a finite penalty and remain
#'   visible in the objective history; errors in a run do not erase other runs.
#'   History rows are objective evaluations, NOT Newton iterations or a claimed
#'   descent trajectory. Residual convergence is assessed independently of the
#'   optimizer's stopping code. Local minimization can stop at a nonzero
#'   residual, miss roots, leave unbounded windows or fail on poor scales.
#'   There is no completeness, convergence or global-extremum guarantee.
#'   Nondifferentiable critical points and constrained/boundary extrema are
#'   outside this stationary-point search. Supply scientifically chosen starts.
#' @examples
#' f <- function(v)v[1]^2+2*v[2]^2
#' f |> vc_stationary_points(rbind(c(.1,.2),c(-.2,.1)),assume_c2=TRUE)
#' @export
vc_stationary_points <- function(f,starts,bounds=NULL,assume_c2=FALSE,h=NULL,
    hessian_h=NULL,tol_grad=1e-5,tol_eigen=1e-6,tol_merge=1e-4,max_iter=500L,
    atol=1e-8,rtol=1e-4,warn=TRUE) {
  if(!is.function(f)) stop("'f' must be a scalar field.")
  if(is.null(dim(starts))) starts <- matrix(.vc_point(starts,"starts"),nrow=1,
    dimnames=list(NULL,names(starts)))
  if(!is.matrix(starts)||!is.numeric(starts)||is.complex(starts)||!nrow(starts)||
     ncol(starts)<2||any(!is.finite(starts))) stop("Invalid 'starts' matrix.")
  .vc_point(stats::setNames(starts[1,],colnames(starts)),"starts")
  .vc_flag(assume_c2,"assume_c2");.vc_matrix_tolerances(atol,rtol,warn)
  .vc_positive(tol_grad,"tol_grad");.vc_positive(tol_eigen,"tol_eigen")
  .vc_positive(tol_merge,"tol_merge");.vc_integer(max_iter,"max_iter",1)
  n <- ncol(starts)
  labels <- if(is.null(colnames(starts)))paste0("x",seq_len(n)) else colnames(starts)
  .vc_matrix_step(starts[1,],h);.vc_matrix_step(starts[1,],hessian_h,TRUE)
  if(!is.null(bounds)) {
    if(!is.matrix(bounds)||!is.numeric(bounds)||is.complex(bounds)||
       !identical(dim(bounds),c(n,2L))||any(!is.finite(bounds))) stop("Invalid 'bounds'.")
    for(i in seq_len(n)) .vc_interval(bounds[i,])
    if(any(sweep(starts,2,bounds[,1],"<"))||any(sweep(starts,2,bounds[,2],">")))
      stop("Starts must lie within bounds.")
  }
  points <- matrix(numeric(),0,n,dimnames=list(NULL,labels))
  table <- data.frame(id=integer(),value=numeric(),residual=numeric(),classification=character())
  assessments <- runs <- list()
  for(k in seq_len(nrow(starts))) {
    history <- list()
    objective <- function(p) {
      names(p) <- colnames(starts)
      error <- ""
      residual <- tryCatch(.vc_norm(vc_gradient(f,p,h=h,atol=atol,rtol=rtol,warn=FALSE)$estimate),
        error=function(e) {error <<- conditionMessage(e);Inf})
      value <- if(is.finite(residual)&&residual<sqrt(.Machine$double.xmax))residual^2 else .Machine$double.xmax/100
      history[[length(history)+1L]] <<- list(point=p,residual=residual,objective=value,
        status=if(nzchar(error))error else if(!is.finite(residual))"Residual overflow" else "ok")
      value
    }
    start <- starts[k,];names(start) <- colnames(starts)
    fit <- tryCatch(if(is.null(bounds))stats::optim(start,objective,method="Nelder-Mead",
      control=list(maxit=max_iter,reltol=1e-12,abstol=tol_grad^2/100)) else
      stats::optim(start,objective,method="L-BFGS-B",lower=bounds[,1],upper=bounds[,2],
        control=list(maxit=max_iter,pgtol=tol_grad^2/100,factr=10)),error=function(e)e)
    a <- NULL;id <- NA_integer_;accepted <- FALSE
    if(inherits(fit,"error")) {
      terminal <- start;code <- NA_integer_;message <- conditionMessage(fit);status <- "optimizer_error"
    } else {
      terminal <- fit$par;code <- fit$convergence
      message <- if(is.null(fit$message))"" else fit$message
      boundary <- !is.null(bounds)&&any(terminal<=bounds[,1]+tol_merge|terminal>=bounds[,2]-tol_merge)
      a <- tryCatch(vc_classify_stationary(f,terminal,assume_c2,h,hessian_h,tol_grad,
        tol_eigen,atol,rtol,FALSE),error=function(e)list(error=conditionMessage(e)))
      accepted <- !boundary && is.null(a$error) && isTRUE(a$stationary_evidence)
      status <- if(boundary)"boundary_rejected" else if(!is.null(a$error))"assessment_error" else
        if(accepted)"candidate" else "stationarity_rejected"
      if(accepted) {
        distances <- if(nrow(points))apply(points,1,function(p).vc_norm(p-terminal)) else numeric()
        match_id <- which(distances<=tol_merge)
        if(length(match_id)) {id <- match_id[1];status <- "duplicate"} else {
          points <- rbind(points,as.double(terminal));id <- nrow(points)
          assessments[[id]] <- a
          table <- rbind(table,data.frame(id=id,value=a$value,residual=.vc_norm(a$gradient$estimate),
            classification=a$classification))
        }
      }
    }
    runs[[k]] <- list(start=start,point=terminal,optimizer_code=code,message=message,
      accepted=accepted,candidate_id=id,status=status,assessment=a,history=history)
  }
  if(warn&&!nrow(points)) warning("No stationary candidates accepted; inspect runs. This is not an absence proof.",call.=FALSE)
  list(type="stationary_points",starts=starts,bounds=bounds,points=points,table=table,
    assessments=assessments,runs=runs,assume_c2=assume_c2,
    note="Local numerical candidate search only; no completeness or convergence guarantee.")
}
