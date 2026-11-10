#' Construct a regular or custom partition
#' @param interval Finite increasing pair of bounds. NULL when breaks supplied.
#' @param n Positive integer number of regular cells; ignored with breaks.
#' @param breaks Optional finite strictly increasing vector of partition points.
#'   When interval is supplied, endpoints must agree exactly.
#' @return Plain list with type, interval, breaks, n, mesh and cells (left,
#'   right, width, midpoint). All cells have positive representable width.
#' @examples
#' vc_partition(c(0,2),n=4)
#' vc_partition(breaks=c(0,.2,1,2))
#' @export
vc_partition <- function(interval=NULL,n=10L,breaks=NULL) {
  if(is.null(breaks)) {
    interval <- .vc_interval(interval)
    .vc_integer(n,"n",1)
    if(n>=.Machine$integer.max) stop("'n' is too large.")
    breaks <- seq(interval[1],interval[2],length.out=n+1)
  } else {
    if(!is.numeric(breaks)||is.complex(breaks)||!is.null(dim(breaks))||length(breaks)<2||any(!is.finite(breaks)))
      stop("'breaks' must be a finite real vector with at least two points.")
    breaks <- as.double(breaks)
    if(!is.null(interval) && !identical(.vc_interval(interval),breaks[c(1,length(breaks))]))
      stop("Partition endpoints must match 'interval'.")
  }
  widths <- diff(breaks)
  if(any(!is.finite(widths))||any(widths<=0)) stop("Partition widths must be positive and finite.")
  interval <- .vc_interval(breaks[c(1,length(breaks))])
  left <- breaks[-length(breaks)]
  right <- breaks[-1L]
  midpoint <- left+widths/2
  if(any(!is.finite(midpoint))) stop("Nonfinite midpoint.")
  list(type="partition",interval=interval,breaks=breaks,n=length(widths),mesh=max(widths),
       cells=data.frame(left=left,right=right,width=widths,midpoint=midpoint))
}

#' Compute signed tagged Riemann sums
#' @inheritParams vc_secant
#' @param interval Finite increasing bounds, or NULL with partition supplied.
#' @param n Positive regular partition size; ignored with partition supplied.
#' @param method One of midpoint, left or right. These are sample rules,
#'   not lower/upper Darboux sums.
#' @param partition Optional vc_partition result. Its cells are reconstructed
#'   from breaks to validate consistency; interval, if given, must match.
#' @return Plain list with type, interval, partition, method, estimate, table
#'   (cell limits, width, midpoint, sample, value, contribution), and note.
#' @details Computes sum f(sample_i)*width_i. Negative heights contribute
#'   negatively: this is not geometric area. All sampled values/contributions
#'   and the total must be finite. A lower Darboux sum uses the infimum on
#'   every cell; upper uses its supremum. Endpoint heights (or minima/maxima
#'   of finitely many samples) do not establish these extrema for general f.
#'   Riemann convergence requires appropriate integrability and vanishing mesh;
#'   a finite table or an increasing cell count is not a proof of convergence.
#' @examples
#' function_to_sum <- function(x)x^2
#' function_to_sum |> vc_riemann(interval=c(0,2),n=100)
#' @export
vc_riemann <- function(f,interval=NULL,n=10L,method=c("midpoint","left","right"),partition=NULL) {
  method <- match.arg(method)
  if(!is.function(f)) stop("'f' must be a function.")
  if(is.null(partition)) partition <- vc_partition(interval,n) else {
    if(!is.list(partition)||!identical(partition$type,"partition")) stop("Invalid 'partition'.")
    partition <- vc_partition(interval,breaks=partition$breaks)
  }
  d <- partition$cells
  d$sample <- switch(method,left=d$left,right=d$right,midpoint=d$midpoint)
  d$value <- vapply(d$sample,function(x).vc_value(f,x),numeric(1))
  d$contribution <- d$value*d$width
  if(any(!is.finite(d$contribution))) stop("Nonfinite rectangle contribution.")
  estimate <- .vc_finite(sum(d$contribution),"Riemann sum")
  list(type="riemann",interval=partition$interval,partition=partition,method=method,
       estimate=estimate,table=d,note="Signed tagged sum; not a Darboux bound or geometric area.")
}

#' Tabulate Riemann refinement without claiming a limit
#' @inheritParams vc_riemann
#' @param n_values Strictly increasing positive integer regular cell counts.
#' @param reference Optional finite scalar analytical/reference value. It is
#'   supplied by the caller, not inferred from the finest sum.
#' @return Plain list with type, interval, method, results (full sum lists),
#'   table (n, mesh, estimate, change, error, abs_error) and note. First change
#'   and unavailable reference errors are NA; errors are estimate-reference.
#' @examples
#' f <- function(x)x^2
#' vc_riemann_convergence(f,c(0,2),reference=8/3)$table
#' @export
vc_riemann_convergence <- function(f,interval,n_values=c(4L,8L,16L,32L,64L),
                                   method=c("midpoint","left","right"),reference=NULL) {
  method <- match.arg(method)
  interval <- .vc_interval(interval)
  if(!is.numeric(n_values)||is.complex(n_values)||!length(n_values)||
     any(!is.finite(n_values))||any(n_values<1)||any(n_values!=floor(n_values))||any(diff(n_values)<=0))
    stop("'n_values' must be strictly increasing positive integers.")
  if(!is.null(reference)&&!.vc_scalar(reference)) stop("'reference' must be a finite real scalar.")
  results <- lapply(n_values,function(n)vc_riemann(f,interval,n,method))
  estimates <- vapply(results,function(x)x$estimate,numeric(1))
  errors <- if(is.null(reference)) rep(NA_real_,length(estimates)) else estimates-reference
  if(any(!is.finite(errors[!is.na(errors)]))) stop("Reference error overflowed.")
  list(type="riemann_convergence",interval=interval,method=method,results=results,
       table=data.frame(n=n_values,mesh=vapply(results,function(x)x$partition$mesh,numeric(1)),
                        estimate=estimates,change=c(NA_real_,diff(estimates)),error=errors,abs_error=abs(errors)),
       note="Observed refinement differences are not rigorous error bounds or a convergence proof.")
}
