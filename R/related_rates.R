#' Assess and solve a missing rate on an implicit level constraint
#' @param g Scalar constraint accepting one coordinate vector; level g=0.
#' @param point Finite coordinate vector, optionally named.
#' @param rates Numeric vector of the same length, in coordinate order. NA
#'   marks unknown rates (an all-NA logical vector is also accepted); all
#'   supplied rates must be finite. If both vectors
#'   are named, names must match exactly in order.
#' @param assume_differentiable Caller declaration of total differentiability
#'   near the point. Finite coordinate diagnostics cannot establish it.
#' @inheritParams vc_gradient
#' @param tol_feasible Positive absolute tolerance for g(point).
#' @param tol_rate Positive absolute tolerance for the first-order rate residual.
#' @param denom_tol Positive denominator threshold factor, multiplied by
#'   max(1,max(abs(gradient))) and enlarged by empirical derivative changes.
#' @return Plain list with status, point, value, feasible, supplied_rates,
#'   rates (unknowns filled only when solved), estimate, candidate, unknown,
#'   gradient (full unchanged engine), denominator, denominator_margin,
#'   known_contribution, residual, rate_indicator, conditional and note.
#' @details Uses grad(g) dot rates=0 conditionally on a differentiable path
#'   lying on g=0 and the caller's total-differentiability hypothesis. Solves
#'   one missing component only after feasibility, derivative reliability and
#'   denominator gates pass. Changes between derivative levels and finest side
#'   disagreement form heuristic denominator/error indicators, not bounds.
#'   Multiple missing rates remain underdetermined; a numerically zero missing
#'   coefficient may give first_order_inconsistent or underdetermined evidence,
#'   never a trajectory existence/nonexistence proof. Small nonzero coefficients
#'   remain denominator_unresolved, not zero by assertion. Fully supplied rates
#'   can be assessed without division. No hidden signature/name remapping.
#'   Explicit time dependence can be modeled by including time as a coordinate
#'   and supplying its rate as one. No higher-order or path integration solver.
#' @examples
#' g <- function(v)sum(v^2)-25
#' g |> vc_related_rates(c(x=4,y=3),c(-1,NA),assume_differentiable=TRUE)
#' @export
vc_related_rates <- function(g,point,rates,assume_differentiable=FALSE,h=NULL,
    tol_feasible=1e-8,tol_rate=1e-8,denom_tol=1e-8,atol=1e-8,rtol=1e-4,warn=TRUE) {
  point <- .vc_point(point)
  .vc_flag(assume_differentiable,"assume_differentiable")
  .vc_matrix_tolerances(atol,rtol,warn)
  .vc_positive(tol_feasible,"tol_feasible");.vc_positive(tol_rate,"tol_rate")
  .vc_positive(denom_tol,"denom_tol")
  if(is.logical(rates)&&is.null(dim(rates))&&all(is.na(rates))) {
    rate_names <- names(rates)
    rates <- as.double(rates)
    names(rates) <- rate_names
  }
  if(!is.numeric(rates)||is.complex(rates)||!is.null(dim(rates))||length(rates)!=length(point)||
     any(is.nan(rates))||any(!is.finite(rates[!is.na(rates)]))) stop("Invalid rates vector.")
  if(!is.null(names(rates))&&(anyNA(names(rates))||any(!nzchar(names(rates)))||anyDuplicated(names(rates))))
    stop("Rate names must be nonempty and unique.")
  if(!is.null(names(rates))&&!is.null(names(point))&&!identical(names(rates),names(point)))
    stop("Rate names must follow point coordinate order.")
  supplied <- rates
  rates <- as.double(rates);names(rates) <- names(point)
  unknown <- which(is.na(rates));known <- which(!is.na(rates))
  grad <- vc_gradient(g,point,h=h,atol=atol,rtol=rtol,warn=FALSE)
  a <- grad$estimate
  margins <- vapply(grad$results,function(z) {
    d <- z$derivative
    max(c(d$changes,abs(d$diagnostics$forward[3]-d$diagnostics$backward[3])))
  },numeric(1))
  b <- .vc_finite(sum(a[known]*rates[known]),"Known rate contribution")
  feasible <- abs(grad$value)<=tol_feasible
  denominator <- if(length(unknown)==1)unname(a[unknown]) else NA_real_
  margin <- if(length(unknown)==1)
    .vc_finite(denom_tol*max(1,abs(a))+margins[unknown],"Denominator margin") else NA_real_
  safe <- length(unknown)==1 && abs(denominator)>margin
  candidate <- if(safe).vc_finite(-b/denominator,"Candidate rate") else NA_real_
  status <- if(!feasible)"infeasible" else if(!assume_differentiable)"assumption_required" else
    if(!grad$reliable)"derivative_unreliable" else if(!length(unknown)) {
      if(abs(b)<=tol_rate)"consistent" else "first_order_inconsistent"
    } else if(length(unknown)>1)"underdetermined" else if(safe)"solved" else
      if(denominator!=0)"denominator_unresolved" else
        if(abs(b)>tol_rate)"first_order_inconsistent" else "underdetermined"
  estimate <- NA_real_
  indicator <- NA_real_
  if(status=="solved") {
    estimate <- candidate
    rates[unknown] <- estimate
    indicator <- .vc_finite((sum(abs(rates[known])*margins[known])+
      abs(estimate)*margins[unknown])/abs(denominator),"Rate sensitivity indicator")
  }
  residual <- if(anyNA(rates))NA_real_ else .vc_finite(sum(a*rates),"Rate residual")
  if(status=="solved" && abs(residual)>tol_rate) {
    status <- "residual_unresolved"
    rates[unknown] <- NA_real_
    estimate <- NA_real_
  }
  if(warn&&!status%in%c("solved","consistent"))warning("Related rates status: ",status,call.=FALSE)
  list(type="related_rates",status=status,point=point,value=grad$value,feasible=feasible,
    supplied_rates=supplied,rates=rates,estimate=estimate,candidate=candidate,unknown=unknown,
    gradient=grad,denominator=denominator,denominator_margin=margin,
    known_contribution=b,residual=residual,rate_indicator=indicator,
    conditional=assume_differentiable,note="Conditional first-order evidence; no path existence or certified error bound.")
}

#' Display the gradient and a stored related-rates velocity
#' @param result Result of vc_related_rates() in two dimensions.
#' @param g Same scalar constraint used for calculation.
#' @param xlim,ylim Optional contour bounds covering the point and velocity tip.
#' @param n Number of samples per contour coordinate, at least two.
#' @param vector_scale Positive display multiplier for both vectors.
#' @param title Optional character title.
#' @return Plotly widget; unresolved velocities are omitted. Sampled contours
#'   show the scalar field, not a certified constraint curve or path.
#' @details Reuses vc_plot_gradient without recalculating derivatives or rates.
#'   The velocity segment is a local first-order vector, not an integrated path.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   g <- function(v)sum(v^2)-25
#'   r <- vc_related_rates(g,c(4,3),c(-1,NA),TRUE)
#'   vc_plot_related_rates(r,g,vector_scale=.2)
#' }
#' @export
vc_plot_related_rates <- function(result,g,xlim=NULL,ylim=NULL,n=61L,vector_scale=1,title=NULL) {
  if(!requireNamespace("plotly",quietly=TRUE))stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"related_rates")||length(result$point)!=2)
    stop("This view requires a two-dimensional related rates result.")
  .vc_positive(vector_scale,"vector_scale")
  point <- result$point
  show <- result$status%in%c("solved","consistent") && all(is.finite(result$rates))
  end <- if(show)point+vector_scale*result$rates else point
  if(any(!is.finite(end)))stop("Velocity display overflowed.")
  extent <- point+vector_scale*result$gradient$estimate
  bounds <- function(i) {
    r <- range(c(point[i]+c(-1,1),end[i],extent[i]))
    r+c(-.2,.2)*diff(r)
  }
  if(is.null(xlim))xlim <- bounds(1)
  if(is.null(ylim))ylim <- bounds(2)
  if(is.null(title))title <- paste("Related rates:",result$status)
  p <- vc_plot_gradient(result$gradient,g,xlim=xlim,ylim=ylim,n=n,
    vector_scale=vector_scale,title=title)
  if(show)p <- plotly::add_trace(p,x=c(point[1],end[1]),y=c(point[2],end[2]),
    type="scatter",mode="lines+markers",name="Conditional velocity",inherit=FALSE)
  p
}
