.vc_constraints <- function(g,point) {
  if(is.list(g)) {
    if(!length(g)||!all(vapply(g,is.function,logical(1))))stop("'g' must contain scalar functions.")
    fun <- function(p)vapply(g,function(gi).vc_field_value(gi,p),numeric(1))
  } else if(is.function(g))fun <- g else stop("'g' must be a vector function or nonempty list of scalar functions.")
  value <- .vc_vector_field_value(fun,point)
  evaluate <- function(p) {
    z <- .vc_vector_field_value(fun,p)
    if(length(z)!=length(value)||!identical(names(z),names(value)))stop("Constraint dimension/names changed.")
    z
  }
  list(value=value,evaluate=evaluate)
}

#' Assess an equality-constrained candidate using Lagrange multipliers
#' @inheritParams vc_classify_stationary
#' @param g Equality constraints g=0, either a deterministic vector-valued
#'   function accepting one coordinate vector, or a nonempty list of scalar
#'   functions. Optional constraint names must be fixed, unique and nonempty.
#' @param assume_c2 Caller declaration that BOTH f and all constraints are C2
#'   in a neighborhood. FALSE prevents extremum labels. Numerical evidence
#'   cannot establish regularity, exact feasibility or exact stationarity.
#' @param tol_feasible Positive absolute infinity-norm constraint tolerance.
#' @param tol_stationary Positive absolute infinity-norm tolerance for the
#'   Lagrange residual and its finest forward/backward coordinate quotients.
#' @param rank_tol Positive relative/absolute SVD rank threshold factor.
#'   Singular values must exceed rank_tol*max(1,largest singular value) plus
#'   an empirical Jacobian refinement allowance. Rank depends on scaling.
#' @param warn Logical; warn when the assessment is inconclusive.
#' @return Plain list with type, point, value, constraints, classification,
#'   reason, feasible, stationary_evidence, regular, rank, rank_threshold,
#'   singular_values, lambda, stationarity_residual, side_max, gradient,
#'   jacobian, hessians (objective/constraints), lagrangian_hessian,
#'   bordered_hessian, tangent_basis, reduced_hessian, eigenvalues,
#'   eigen_margin, control and note. Unavailable quantities are NULL/NA and
#'   engine failures remain in the corresponding diagnostics. Labels are
#'   minimum, maximum, saddle, inconclusive, infeasible or not_stationary.
#'   Multipliers at deficient rank are minimum-norm diagnostics, not unique.
#' @details Evaluates a supplied candidate; it does not search for candidates.
#'   Uses L=f-sum(lambda_i*g_i), so grad f=transpose(J)*lambda. vc_gradient,
#'   vc_jacobian and direct vc_hessian engines are reused unchanged. The SVD
#'   of J solves for multipliers without forming the normal equations J*J^T.
#'   Stored raw sensitivity flags remain unchanged. First-order evidence uses
#'   stable central estimates and BOTH small finest coordinate-side Lagrange
#'   residuals; individual zero partials can fail stricter relative side flags
#'   even for smooth fields. This criterion is finite evidence, not proof.
#'
#'   A full-row-rank Jacobian identifies an estimated tangent basis Z. With
#'   feasible/stationary evidence, declared C2 and reliable Hessian stencils,
#'   eigenvalues of transpose(Z)*H_L*Z suggest strict constrained minimum or
#'   maximum (definite), or saddle (both signs). Semidefinite, singular-rank,
#'   unresolved and zero-dimensional tangent cases remain inconclusive. The
#'   original bordered matrix is retained for teaching, but its arbitrary
#'   leading-minor sign pattern is not used as a universal classifier.
#'
#'   Eigenvalue margin is tol_eigen plus dimension times the maximum entry
#'   refinement change of H_L and weighted roundoff indicators. Jacobian rank
#'   allowance uses sqrt(m*n) times its maximum consecutive refinement change.
#'   These are heuristic sensitivity diagnostics, not certified matrix or
#'   tangent-space error bounds. Individual Hessian flags are required (zero
#'   multiplier constraint Hessians do not affect H_L reliability). Scaling,
#'   offsets and near-zero curvature may require informed h/hessian_h choices.
#'   No inequalities, boundary/global optimization, singular-constraint theory,
#'   higher-order classification or existence/completeness proof is supplied.
#' @examples
#' f <- function(v)sum(v^2)
#' g <- function(v)sum(v)-1
#' f |> vc_lagrange(g=g,point=c(.5,.5),assume_c2=TRUE,hessian_h=.01)
#' @export
vc_lagrange <- function(f,g,point,assume_c2=FALSE,h=NULL,hessian_h=NULL,
    tol_feasible=1e-8,tol_stationary=1e-5,rank_tol=1e-8,tol_eigen=1e-6,
    atol=1e-8,rtol=1e-4,warn=TRUE) {
  point <- .vc_point(point);n <- length(point)
  if(n<2)stop("Equality-constrained assessment requires at least two coordinates.")
  .vc_flag(assume_c2,"assume_c2");.vc_matrix_tolerances(atol,rtol,warn)
  for(nm in c("tol_feasible","tol_stationary","rank_tol","tol_eigen")) .vc_positive(get(nm),nm)
  .vc_matrix_step(point,h);.vc_matrix_step(point,hessian_h,TRUE)
  value <- .vc_field_value(f,point);constraints <- .vc_constraints(g,point)
  m <- length(constraints$value)
  out <- list(type="lagrange",point=point,value=value,constraints=constraints$value,
    classification="inconclusive",reason="C2 regularity has not been declared.",
    feasible=max(abs(constraints$value))<=tol_feasible,stationary_evidence=FALSE,
    regular=FALSE,rank=NA_integer_,rank_threshold=NA_real_,singular_values=NULL,
    lambda=NULL,stationarity_residual=NULL,side_max=NA_real_,gradient=NULL,jacobian=NULL,
    hessians=NULL,lagrangian_hessian=NULL,bordered_hessian=NULL,tangent_basis=NULL,
    reduced_hessian=NULL,eigenvalues=NULL,eigen_margin=NA_real_,
    control=list(assume_c2=assume_c2,h=h,hessian_h=hessian_h,tol_feasible=tol_feasible,
      tol_stationary=tol_stationary,rank_tol=rank_tol,tol_eigen=tol_eigen,atol=atol,rtol=rtol),
    note="Conditional numerical assessment only; no regularity, root or constrained-extremum certificate.")
  finish <- function(z) {
    if(warn&&z$classification=="inconclusive")warning(z$reason,call.=FALSE)
    z
  }
  if(!out$feasible) {
    out$classification <- "infeasible";out$reason <- "Constraint residual exceeds tol_feasible."
    return(finish(out))
  }
  safe <- function(expr)tryCatch(suppressWarnings(expr),error=function(e)list(error=conditionMessage(e)))
  grad <- safe(vc_gradient(f,point,h=h,atol=atol,rtol=rtol,warn=FALSE))
  jac <- safe(vc_jacobian(constraints$evaluate,point,h=h,atol=atol,rtol=rtol,warn=FALSE))
  out$gradient <- grad;out$jacobian <- jac
  if(!is.null(grad$error)||!is.null(jac$error)) {
    out$reason <- "First derivative stencils unavailable; inspect gradient/jacobian."
    return(finish(out))
  }
  J <- jac$estimate
  s <- safe(svd(J,nu=min(m,n),nv=n))
  if(!is.null(s$error)||any(!is.finite(s$d))) {
    out$reason <- "Constraint SVD unavailable or nonfinite."
    return(finish(out))
  }
  variation <- sqrt(m*n)*max(abs(jac$refinements[[1]]-jac$refinements[[2]]),
    abs(jac$refinements[[2]]-jac$refinements[[3]]))
  threshold <- rank_tol*max(1,s$d[1])+variation
  indices <- which(s$d>threshold)
  r <- length(indices)
  lambda <- if(r)as.double(s$u[,indices,drop=FALSE]%*%
    (as.double(crossprod(s$v[,indices,drop=FALSE],grad$estimate))/s$d[indices])) else numeric(m)
  out$lambda <- stats::setNames(lambda,rownames(J));out$rank <- r
  out$rank_threshold <- threshold;out$singular_values <- s$d;out$regular <- r==m
  out$tangent_basis <- if(r<n)s$v[,seq.int(r+1L,n),drop=FALSE] else matrix(numeric(),n,0)
  residual <- as.double(grad$estimate-crossprod(J,lambda))
  sides <- vapply(grad$results,function(z)c(z$derivative$diagnostics$forward[3],
    z$derivative$diagnostics$backward[3]),numeric(2))
  side_j <- replicate(2,matrix(0,m,n),simplify=FALSE)
  for(entry in jac$results)for(k in 1:2)side_j[[k]][entry$row,entry$column] <-
    entry$derivative$diagnostics[[if(k==1)"forward" else "backward"]][3]
  side_res <- lapply(1:2,function(k)sides[k,]-as.double(crossprod(side_j[[k]],lambda)))
  side_max <- max(abs(unlist(side_res)))
  stable <- all(jac$table$stable)&&all(vapply(grad$results,function(z)z$derivative$stable,logical(1)))
  out$stationarity_residual <- residual;out$side_max <- side_max
  out$stationary_evidence <- all(is.finite(c(lambda,residual,side_max)))&&
    max(abs(residual))<=tol_stationary&&side_max<=tol_stationary&&stable
  if(!out$regular) {
    out$reason <- "Constraint Jacobian lacks numerical full row rank; regular test is inconclusive."
    return(finish(out))
  }
  if(!out$stationary_evidence) {
    out$classification <- if(all(is.finite(residual))&&max(abs(residual))>tol_stationary)
      "not_stationary" else "inconclusive"
    out$reason <- "Lagrange residual or coordinate-side/stability checks failed."
    return(finish(out))
  }
  Hf <- safe(vc_hessian(f,point,h=hessian_h,atol=atol,rtol=rtol,warn=FALSE))
  Hg <- lapply(seq_len(m),function(i)safe(vc_hessian(function(p)constraints$evaluate(p)[i],
    point,h=hessian_h,atol=atol,rtol=rtol,warn=FALSE)))
  out$hessians <- list(objective=Hf,constraints=Hg)
  if(!is.null(Hf$error)||any(vapply(Hg,function(z)!is.null(z$error),logical(1)))) {
    out$reason <- "Hessian stencils unavailable; inspect hessians."
    return(finish(out))
  }
  combine <- function(k) {
    A <- Hf$refinements[[k]]
    for(i in seq_len(m))A <- A-lambda[i]*Hg[[i]]$refinements[[k]]
    A
  }
  refinements <- lapply(1:3,combine);HL <- refinements[[3]]
  out$lagrangian_hessian <- HL
  out$bordered_hessian <- rbind(cbind(matrix(0,m,m),J),cbind(t(J),HL))
  if(any(!is.finite(HL))) {
    out$reason <- "Lagrangian Hessian overflowed; rescale the problem."
    return(finish(out))
  }
  Z <- out$tangent_basis
  if(!ncol(Z)) {
    out$reason <- "Zero-dimensional tangent space; directional second-order classification is not supplied."
    return(finish(out))
  }
  reduced <- crossprod(Z,HL%*%Z);reduced <- reduced/2+t(reduced)/2
  out$reduced_hessian <- reduced
  ev <- safe(eigen(reduced,symmetric=TRUE,only.values=TRUE))
  if(!is.null(ev$error)||any(!is.finite(ev$values))) {
    out$reason <- "Reduced eigenvalue calculation unavailable or nonfinite."
    return(finish(out))
  }
  out$eigenvalues <- ev$values
  rounding <- function(H)max(vapply(H$results,function(z)z$derivative$diagnostics$roundoff_indicator[3],numeric(1)))
  allowance <- rounding(Hf)
  for(i in seq_len(m))if(lambda[i]!=0)allowance <- allowance+abs(lambda[i])*rounding(Hg[[i]])
  change <- max(abs(refinements[[1]]-refinements[[2]]),abs(refinements[[2]]-refinements[[3]]))
  margin <- tol_eigen+n*(change+allowance);out$eigen_margin <- margin
  reliable <- Hf$all_reliable&&all(vapply(seq_len(m),function(i)lambda[i]==0||Hg[[i]]$all_reliable,logical(1)))
  if(!reliable||!is.finite(margin))out$reason <- "Hessian sensitivity/roundoff checks failed."
  else if(assume_c2) {
    e <- ev$values
    out$classification <- if(all(e>margin))"minimum" else if(all(e< -margin))"maximum" else
      if(any(e>margin)&&any(e< -margin))"saddle" else "inconclusive"
    out$reason <- if(out$classification=="inconclusive")"Reduced second derivative test is inconclusive." else
      "Conditional regular equality-constrained second derivative test."
  }
  finish(out)
}
