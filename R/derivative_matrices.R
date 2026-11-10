.vc_matrix_step <- function(point,h,second=FALSE) {
  if(is.null(h)) return(if(second)pmax(1,abs(point))*.Machine$double.eps^(1/4) else NULL)
  if(!is.numeric(h)||is.complex(h)||!is.null(dim(h))||!length(h)%in%c(1L,length(point))||
     any(!is.finite(h))||any(h<=0)) stop("'h' must be NULL, positive scalar or one positive step per coordinate.")
  rep(as.double(h),length.out=length(point))
}

.vc_matrix_tolerances <- function(atol,rtol,warn) {
  .vc_positive(atol,"atol",TRUE)
  .vc_positive(rtol,"rtol",TRUE)
  if(!is.logical(warn)||length(warn)!=1||is.na(warn)) stop("Invalid 'warn'.")
}

.vc_vector_field_value <- function(f,p) {
  if(!is.function(f)) stop("'f' must be a function of one numeric vector.")
  .vc_point(f(p),"field output")
}

.vc_cached_field <- function(f,point,vector=FALSE) {
  storage <- new.env(parent=emptyenv())
  value <- if(vector).vc_vector_field_value(f,point) else .vc_field_value(f,point)
  key <- function(p)paste(sprintf("%a",as.double(p)),collapse=";")
  assign(key(point),value,storage)
  count <- 1L
  get_value <- function(p) {
    k <- key(p)
    if(exists(k,storage,inherits=FALSE)) return(get(k,storage,inherits=FALSE))
    z <- if(vector).vc_vector_field_value(f,p) else .vc_field_value(f,p)
    if(length(z)!=length(value)||!identical(names(z),names(value)))
      stop("Field output dimension/names changed across stencil points.")
    count <<- count+1L
    assign(k,z,storage)
    z
  }
  list(value=value,get=get_value,count=function()count)
}

#' Approximate the Jacobian of a vector-valued function
#' @inheritParams vc_gradient
#' @param f Deterministic function accepting one finite coordinate vector and
#'   returning a nonempty finite real vector of fixed dimension and names.
#'   Optional output names must be unique and nonempty.
#'   A scalar output is allowed and produces a one-row Jacobian.
#' @return Plain list with type, point, value, estimate (outputs by inputs),
#'   reliable (logical matrix), all_reliable, table, results (per-entry full
#'   derivative diagnostics), refinements (three estimate matrices),
#'   evaluation_count and note. Matrix names follow output/point names, or
#'   default to f1,... and x1,... . Table rows identify output and input index,
#'   labels, method, initial step, estimate, stable and reliable.
#' @details Reuses vc_derivative on each output/coordinate restriction, with
#'   shared cached vector evaluations. It does not guess f(x,y) signatures or
#'   silently recycle output dimensions. All three relative refinement levels
#'   are retained. Entry checks do not prove total differentiability, bound
#'   matrix error or establish a chain-rule hypothesis. One-sided differences
#'   are labeled and provide only their selected coordinate-side evidence.
#' @examples
#' f <- function(v)c(first=v[1]^2+v[2],second=v[2]*v[3])
#' f |> vc_jacobian(c(1,2,3))
#' @export
vc_jacobian <- function(f,point,method="central",h=NULL,atol=1e-8,rtol=1e-4,warn=TRUE) {
  point <- .vc_point(point)
  n <- length(point)
  .vc_matrix_tolerances(atol,rtol,warn)
  if(!is.character(method)||!length(method)||anyNA(method)||!length(method)%in%c(1L,n)||
     any(!method%in%c("central","forward","backward"))) stop("Invalid Jacobian 'method'.")
  method <- rep(method,length.out=n)
  h <- .vc_matrix_step(point,h)
  field <- .vc_cached_field(f,point,TRUE)
  m <- length(field$value)
  inputs <- if(is.null(names(point)))paste0("x",seq_len(n)) else names(point)
  outputs <- if(is.null(names(field$value)))paste0("f",seq_len(m)) else names(field$value)
  estimate <- matrix(NA_real_,m,n,dimnames=list(outputs,inputs))
  reliable <- matrix(FALSE,m,n,dimnames=dimnames(estimate))
  refinements <- replicate(3,estimate,simplify=FALSE)
  results <- rows <- list()
  for(j in seq_len(n)) for(i in seq_len(m)) {
    slice <- function(t) {p <- point;p[j] <- t;field$get(p)[i]}
    d <- vc_derivative(slice,point[j],method[j],if(is.null(h))NULL else h[j],atol,rtol,FALSE)
    estimate[i,j] <- d$estimate
    reliable[i,j] <- d$reliable
    for(k in 1:3) refinements[[k]][i,j] <- d$diagnostics$estimate[k]
    results[[length(results)+1L]] <- list(row=i,column=j,derivative=d)
    rows[[length(rows)+1L]] <- data.frame(row=i,column=j,output=outputs[i],input=inputs[j],
      method=method[j],h=d$h,estimate=d$estimate,stable=d$stable,reliable=d$reliable)
  }
  if(warn&&!all(reliable)) warning("Some Jacobian entry sensitivity checks failed; no differentiability or accuracy guarantee.",call.=FALSE)
  list(type="jacobian",point=point,value=field$value,estimate=estimate,reliable=reliable,
    all_reliable=all(reliable),table=do.call(rbind,rows),results=results,refinements=refinements,
    evaluation_count=field$count(),note="Finite entry checks do not prove total differentiability or bound Jacobian error.")
}

.vc_mixed_second <- function(field,point,i,j,hi,hj,atol,rtol) {
  steps_i <- hi/c(1,2,4)
  steps_j <- hj/c(1,2,4)
  estimates <- indicators <- widths_i <- widths_j <- numeric(3)
  for(k in 1:3) {
    li <- point[i]-steps_i[k];ri <- point[i]+steps_i[k]
    lj <- point[j]-steps_j[k];rj <- point[j]+steps_j[k]
    if(any(!is.finite(c(li,ri,lj,rj)))||li>=point[i]||ri<=point[i]||lj>=point[j]||rj<=point[j])
      stop("Mixed second-difference stencil collapsed or overflowed.")
    mi <- .vc_finite((point[i]-li)/2+(ri-point[i])/2,"First half-width")
    mj <- .vc_finite((point[j]-lj)/2+(rj-point[j])/2,"Second half-width")
    at <- function(a,b) {p <- point;p[i] <- a;p[j] <- b;field$get(p)}
    v <- c(at(ri,rj),at(ri,lj),at(li,rj),at(li,lj))
    numerator <- .vc_finite((v[1]-v[2])-(v[3]-v[4]),"Mixed difference numerator")
    estimates[k] <- .vc_finite(((numerator/mi)/2/mj)/2,"Mixed second difference")
    indicators[k] <- (4*.Machine$double.eps*max(abs(v))/mi)/mj
    if(!is.finite(indicators[k])) indicators[k] <- Inf
    widths_i[k] <- mi
    widths_j[k] <- mj
  }
  close <- function(a,b) {
    s <- max(abs(a),abs(b),atol,.Machine$double.xmin)
    abs(a/s-b/s)<=atol/s+rtol*max(abs(a/s),abs(b/s))
  }
  stable <- close(estimates[1],estimates[2])&&close(estimates[2],estimates[3])
  s <- max(abs(estimates[3]),atol,.Machine$double.xmin)
  resolved <- is.finite(indicators[3])&&indicators[3]/s<=atol/s+rtol*abs(estimates[3]/s)
  list(estimate=estimates[3],h=c(hi,hj),stable=stable,resolved=resolved,reliable=stable&&resolved,
    diagnostics=data.frame(step_i=steps_i,step_j=steps_j,half_width_i=widths_i,half_width_j=widths_j,
      estimate=estimates,roundoff_indicator=indicators),changes=abs(diff(estimates)))
}

#' Approximate a scalar-field Hessian with direct differences
#' @inheritParams vc_gradient
#' @param h Initial positive step(s) for second differences. NULL uses
#'   `pmax(1,abs(point))` times machine epsilon to the power one quarter.
#'   One step or one per coordinate is allowed. Central interior stencils
#'   are always required; no automatic boundary fallback is applied.
#' @return Plain list with type, point, value, estimate (symmetric matrix),
#'   reliable (logical matrix), all_reliable, table (upper-triangle entries,
#'   steps, estimate, stable, resolved, reliable), results (full direct stencil
#'   diagnostics), refinements (three matrices), evaluation_count and note.
#' @details Diagonal entries reuse the accepted direct .vc_second_difference
#'   engine used by vc_taylor. Mixed entries use four rectangle corners,
#'   never numerical differentiation of approximate gradients. Actual rectangle
#'   widths are used. Initial, half and quarter steps are retained; mixed
#'   roundoff indicators are `16*epsilon*max(abs(corners))/(width_i*width_j)`.
#'   These and diagonal indicators are heuristics, not certified error bounds.
#'   Near-zero curvature or large offsets may give conservative failed flags.
#'   Upper-triangle estimates are mirrored by construction. Symmetry is not a
#'   test of equality/existence of ordered mixed partials or of C2 regularity.
#'   All finite checks can miss nonsmoothness or sampling aliasing. No eigenvalue
#'   classification or extremum claim is made: a Hessian alone does not locate
#'   a critical point. The rectangle formula recovers the Hessian primitive
#'   embedded in original critical_points_nd/critical_points_2d, independently
#'   of their search/classification code.
#' @examples
#' f <- function(v)v[1]^2+v[1]*v[2]+2*v[2]^2
#' f |> vc_hessian(c(1,1),h=.01)
#' @export
vc_hessian <- function(f,point,h=NULL,atol=1e-8,rtol=1e-4,warn=TRUE) {
  point <- .vc_point(point)
  .vc_matrix_tolerances(atol,rtol,warn)
  h <- .vc_matrix_step(point,h,TRUE)
  field <- .vc_cached_field(f,point)
  n <- length(point)
  labels <- if(is.null(names(point)))paste0("x",seq_len(n)) else names(point)
  estimate <- matrix(NA_real_,n,n,dimnames=list(labels,labels))
  reliable <- matrix(FALSE,n,n,dimnames=dimnames(estimate))
  refinements <- replicate(3,estimate,simplify=FALSE)
  results <- rows <- list()
  for(i in seq_len(n)) for(j in i:n) {
    if(i==j) {
      slice <- function(t) {p <- point;p[i] <- t;field$get(p)}
      d <- .vc_second_difference(slice,point[i],h[i],atol,rtol)
    } else d <- .vc_mixed_second(field,point,i,j,h[i],h[j],atol,rtol)
    estimate[i,j] <- estimate[j,i] <- d$estimate
    reliable[i,j] <- reliable[j,i] <- d$reliable
    for(k in 1:3) refinements[[k]][i,j] <- refinements[[k]][j,i] <- d$diagnostics$estimate[k]
    results[[length(results)+1L]] <- list(row=i,column=j,derivative=d)
    rows[[length(rows)+1L]] <- data.frame(row=i,column=j,coordinate_i=labels[i],coordinate_j=labels[j],
      h_i=h[i],h_j=h[j],estimate=d$estimate,stable=d$stable,resolved=d$resolved,reliable=d$reliable)
  }
  if(warn&&!all(reliable)) warning("Some Hessian sensitivity/roundoff checks failed; no smoothness or accuracy guarantee.",call.=FALSE)
  list(type="hessian",point=point,value=field$value,estimate=estimate,reliable=reliable,
    all_reliable=all(reliable),table=do.call(rbind,rows),results=results,refinements=refinements,
    evaluation_count=field$count(),note="Symmetry is imposed by the stencil, not proof of C2 regularity or an extremum.")
}
