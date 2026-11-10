.vc_point <- function(point,name="point") {
  if(!is.numeric(point)||is.complex(point)||!is.null(dim(point))||!length(point)||any(!is.finite(point)))
    stop("'",name,"' must be a nonempty finite real vector.")
  nm <- names(point)
  if(!is.null(nm)&&(anyNA(nm)||any(!nzchar(nm))||anyDuplicated(nm))) stop("Point names must be nonempty and unique.")
  point <- as.double(point)
  names(point) <- nm
  point
}

.vc_field_value <- function(f,point) {
  if(!is.function(f)) stop("'f' must be a function of one numeric vector.")
  .vc_finite(f(point),"Scalar field value")
}

.vc_field_safe <- function(f,point)tryCatch(suppressWarnings(.vc_field_value(f,point)),error=function(e)NA_real_)

.vc_unit <- function(direction) {
  scale <- max(abs(direction))
  if(scale==0) stop("'direction' must be nonzero.")
  scaled <- direction/scale
  unit <- scaled/sqrt(sum(scaled^2))
  if(any(direction!=0 & unit==0)) stop("Direction normalization underflowed; rescale the coordinates.")
  unit
}

#' Approximate a coordinate partial derivative
#' @param f Deterministic scalar field accepting ONE numeric vector and returning
#'   one finite real scalar. For f(x,y), supply an explicit wrapper
#'   `function(v)f(v[1],v[2])`; no automatic signature guessing is performed.
#' @param point Nonempty finite real coordinate vector. Optional names must be
#'   unique and nonempty and are passed to f.
#' @param coordinate Integer coordinate index from one to length(point).
#' @inheritParams vc_derivative
#' @return Plain list with type, point, coordinate, value, estimate, derivative
#'   (complete vc_derivative diagnostics), reliable and note.
#' @details Holds every other coordinate fixed and reuses vc_derivative().
#'   Its actual-displacement and three-scale checks remain unchanged. Existence
#'   of partial derivatives is not proof of differentiability. Finite samples
#'   cannot establish the existence of the mathematical partial derivative.
#' @examples
#' (function(v)v[1]^2+2*v[2]^2) |> vc_partial(c(1,2),coordinate=2)
#' @export
vc_partial <- function(f,point,coordinate,method=c("central","forward","backward"),
                       h=NULL,atol=1e-8,rtol=1e-4,warn=TRUE) {
  point <- .vc_point(point)
  .vc_integer(coordinate,"coordinate",1)
  if(coordinate>length(point)) stop("'coordinate' exceeds the point dimension.")
  method <- match.arg(method)
  slice <- function(t) {
    p <- point
    p[coordinate] <- t
    .vc_field_value(f,p)
  }
  d <- vc_derivative(slice,point[coordinate],method,h,atol,rtol,warn)
  list(type="partial",point=point,coordinate=coordinate,value=.vc_field_value(f,point),
    estimate=unname(d$estimate),derivative=d,reliable=d$reliable,
    note="Coordinate differences do not prove existence of partials or total differentiability.")
}

#' Approximate the gradient of a scalar field
#' @inheritParams vc_partial
#' @param h NULL for per-coordinate vc_derivative heuristics, one positive
#'   finite step for all coordinates, or a positive finite vector of that
#'   dimension. No recycling of shorter vectors is allowed.
#' @param method One difference method, or one per coordinate. Explicit
#'   one-sided methods can be used at known coordinate-domain boundaries.
#' @return Plain list with type, point, value, estimate (named vector), results
#'   (full partial results), table (coordinate, label, method, h, estimate,
#'   reliable), reliable (all coordinate checks), and note.
#' @details Reuses vc_partial(), which reuses the accepted one-variable engine.
#'   A reliable flag means only that its finite coordinate checks passed.
#'   Even all existing coordinate partials do not imply total differentiability
#'   or justify computing every directional derivative by a dot product.
#'   The engine's conservative side checks may fail at smooth stationary points.
#' @examples
#' (function(v)v[1]^2+2*v[2]^2) |> vc_gradient(c(x=1,y=2))
#' @export
vc_gradient <- function(f,point,method="central",h=NULL,atol=1e-8,rtol=1e-4,warn=TRUE) {
  point <- .vc_point(point)
  n <- length(point)
  if(!is.character(method)||!length(method)||anyNA(method)||!length(method)%in%c(1L,n)||
     any(!method%in%c("central","forward","backward"))) stop("Invalid gradient 'method'.")
  methods <- rep(method,length.out=n)
  if(!is.null(h)) {
    if(!is.numeric(h)||is.complex(h)||!is.null(dim(h))||!length(h)%in%c(1L,n)||
       any(!is.finite(h))||any(h<=0)) stop("'h' must be NULL, a positive scalar or one step per coordinate.")
    h <- rep(as.double(h),length.out=n)
  }
  if(!is.logical(warn)||length(warn)!=1||is.na(warn)) stop("Invalid 'warn'.")
  results <- lapply(seq_len(n),function(i)vc_partial(f,point,i,methods[i],if(is.null(h))NULL else h[i],atol,rtol,FALSE))
  labels <- if(is.null(names(point)))paste0("x",seq_len(n)) else names(point)
  estimate <- vapply(results,function(z)z$estimate,numeric(1))
  reliable <- vapply(results,function(z)z$reliable,logical(1))
  names(estimate) <- labels
  if(warn&&!all(reliable)) warning("Some coordinate sensitivity checks failed; differentiability is not established.",call.=FALSE)
  list(type="gradient",point=point,value=.vc_field_value(f,point),estimate=estimate,
    results=results,table=data.frame(coordinate=seq_len(n),label=labels,method=methods,
      h=vapply(results,function(z)z$derivative$h,numeric(1)),estimate=unname(estimate),reliable=reliable),
    reliable=all(reliable),note="Coordinate sensitivity checks do not establish total differentiability.")
}

#' Approximate a directional derivative on its actual line restriction
#' @inheritParams vc_partial
#' @param direction Nonzero finite real vector matching point. It is normalized
#'   by scaled arithmetic to avoid squared-norm overflow/underflow. No change
#'   in direction magnitude changes the mathematical unit-direction derivative.
#' @param gradient Optional vc_gradient result at this same point/value for a
#'   diagnostic dot-product comparison. NULL does not compute any gradient.
#' @param h Initial step in the line parameter, or NULL for max(1,abs(point))
#'   times the one-variable epsilon heuristic. Every active coordinate must
#'   move representably at each nonzero stencil parameter; otherwise an error
#'   asks for rescaling or a larger step, rather than returning a false zero.
#' @return Plain list with type, point, value, direction (unit vector), method,
#'   estimate, derivative (full one-variable diagnostics), reliable, optional
#'   gradient_dot and difference (NA without gradient), gradient_comparison
#'   (not_requested, compatible, mismatch or inconclusive), and note.
#' @details Differentiates `t -> f(point+t*direction)` at zero with
#'   vc_derivative(), rather than assuming the gradient formula. With a supplied
#'   gradient, agreement is checked against atol+rtol times the larger magnitude,
#'   and both diagnostic flags must pass for compatible/mismatch labels.
#'   The identity with the dot product requires appropriate differentiability;
#'   agreement of finitely many probes is not a proof. Discontinuous or nonsmooth
#'   fields can have misleading symmetric quotients. One-sided methods yield
#'   only one-sided evidence. Normalized components lost to underflow are errors.
#'   For mathematically axis-aligned directions, use exact zeros rather than
#'   roundoff-sized components from trigonometric evaluations at cardinal angles.
#' @examples
#' f <- function(v)v[1]^2+2*v[2]^2
#' g <- f |> vc_gradient(c(1,2))
#' f |> vc_directional(c(1,2),direction=c(3,4),gradient=g)
#' @export
vc_directional <- function(f,point,direction,method=c("central","forward","backward"),
                           h=NULL,atol=1e-8,rtol=1e-4,warn=TRUE,gradient=NULL) {
  point <- .vc_point(point)
  direction <- .vc_point(direction,"direction")
  if(length(direction)!=length(point)) stop("Direction and point dimensions differ.")
  u <- .vc_unit(direction)
  method <- match.arg(method)
  if(is.null(h)) h <- max(1,abs(point))*.Machine$double.eps^(if(method=="central")1/3 else 1/2)
  slice <- function(t) {
    p <- point+t*u
    if(any(!is.finite(p))) stop("Directional stencil overflowed.")
    if(t!=0&&any(u!=0 & p==point)) stop("Active directional coordinates collapsed; enlarge h or rescale coordinates.")
    .vc_field_value(f,p)
  }
  d <- vc_derivative(slice,0,method,h,atol,rtol,warn)
  value <- .vc_field_value(f,point)
  dot <- delta <- NA_real_
  status <- "not_requested"
  if(!is.null(gradient)) {
    if(!is.list(gradient)||!identical(gradient$type,"gradient")||
       !identical(as.double(gradient$point),as.double(point))||!identical(gradient$value,value)||
       length(gradient$estimate)!=length(point)||any(!is.finite(gradient$estimate))||
       !is.logical(gradient$reliable)||length(gradient$reliable)!=1||is.na(gradient$reliable))
      stop("'gradient' must be a valid result at the same point and field value.")
    dot <- .vc_finite(sum(gradient$estimate*u),"Gradient dot product")
    delta <- .vc_finite(d$estimate-dot,"Directional comparison")
    threshold <- .vc_finite(atol+rtol*max(abs(dot),abs(d$estimate)),"Comparison threshold")
    status <- if(!d$reliable||!gradient$reliable)"inconclusive" else
      if(abs(delta)<=threshold)"compatible" else "mismatch"
  }
  list(type="directional",point=point,value=value,direction=u,method=method,
    estimate=d$estimate,derivative=d,reliable=d$reliable,gradient_dot=dot,difference=delta,
    gradient_comparison=status,note="Line-restriction diagnostics and gradient agreement do not prove total differentiability.")
}
