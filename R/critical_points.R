#' Classify a possible critical point of a real function
#'
#' Evaluate stationarity, neighborhood signs and second-difference evidence.
#' All classifications are conditional numerical evidence, never proofs.
#' @inheritParams vc_tangent
#' @param radius Positive neighborhood radius. NULL uses 0.01*max(1,abs(x0)).
#'   Three radii (radius, radius/2, radius/4) are examined.
#' @param h2 Optional positive second-difference step, as in vc_taylor().
#' @param slope_tol Nonnegative absolute tolerance for a stationary candidate.
#' @return Plain list with type, x0, value, kind, classification, evidence,
#'   derivative, second_derivative, neighborhood, continuity_evidence and note.
#'   kind distinguishes stationary_candidate, nonsmooth_candidate, regular,
#'   flat_candidate and unresolved. classification is minimum, maximum,
#'   not_extremum, not_critical or inconclusive.
#' @details A finite central quotient is not sufficient for stationarity.
#'   Persistent disagreement of neighborhood one-sided slopes flags possible
#'   nonsmoothness. A shrinking gap can occur at a smooth stationary point
#'   even when vc_derivative's conservative side check fails.
#'   Stable near-zero central slopes plus reliable nonzero second differences
#'   support the second-derivative criterion under local smoothness. Otherwise
#'   consistent stable first-derivative signs at neighborhood points may
#'   support a first-sign criterion
#'   under continuity and side differentiability assumptions. Near-zero second
#'   derivative alone is inconclusive. Missing, noncontracting, conflicting or
#'   negligible neighborhood evidence is not classified automatically.
#'   The theoretical second-derivative test assumes differentiability near
#'   a stationary point and a nonzero second derivative there. The first-sign
#'   test requires signs on whole punctured neighborhoods, not merely samples.
#'   Three radii cannot prove these hypotheses, continuity or a local extremum.
#'   Warnings during derivative probes are suppressed; invalid probes become
#'   missing diagnostics, never stationary values.
#' @examples
#' vc_classify_critical(function(x) x^2, 0)
#' vc_classify_critical(function(x) x^3, 0)$classification
#' abs |> vc_classify_critical(x0=0)
#' @export
vc_classify_critical <- function(f,x0,radius=NULL,h=NULL,h2=NULL,
                                 slope_tol=1e-6,atol=1e-8,rtol=1e-4) {
  y <- .vc_value(f,x0)
  .vc_positive(slope_tol,"slope_tol",TRUE)
  .vc_positive(atol,"atol",TRUE)
  .vc_positive(rtol,"rtol",TRUE)
  .vc_optional_step(h,"h")
  .vc_optional_step(h2,"h2")
  if (is.null(radius)) radius <- .01*max(1,abs(x0))
  .vc_positive(radius,"radius")
  neighborhood <- .vc_continuity_evidence(f,x0,radius)
  d <- .vc_first(f,x0,h,atol,rtol)
  d2 <- .vc_second(f,x0,h2,atol,rtol)
  kind <- "unresolved"
  classification <- "inconclusive"
  evidence <- "insufficient_numerical_evidence"
  a <- neighborhood$data$backward
  b <- neighborhood$data$forward
  finite <- all(is.finite(c(a,b)))
  if (!is.null(d) && finite) {
    gaps <- abs(a-b)
    nonsmooth <- gaps[3]>atol && gaps[3]>=.999*gaps[1]
    stationary <- d$stable && abs(d$estimate)<=slope_tol && !nonsmooth
    flat <- all(abs(c(a,b))<=atol) && abs(d$estimate)<=slope_tol
    kind <- if (nonsmooth) "nonsmooth_candidate" else if (flat) "flat_candidate" else
      if (stationary) "stationary_candidate" else if (d$stable) "regular" else "unresolved"
    if (kind=="regular") classification <- "not_critical"
    radii <- neighborhood$data$radius
    ld <- lapply(x0-radii,function(x) .vc_first(f,x,h,atol,rtol))
    rd <- lapply(x0+radii,function(x) .vc_first(f,x,h,atol,rtol))
    estimate <- function(d) if(is.null(d)) NA_real_ else d$estimate
    stable <- function(d) !is.null(d) && d$stable
    neighborhood$data$left_derivative <- vapply(ld,estimate,numeric(1))
    neighborhood$data$right_derivative <- vapply(rd,estimate,numeric(1))
    neighborhood$data$left_stable <- vapply(ld,stable,logical(1))
    neighborhood$data$right_stable <- vapply(rd,stable,logical(1))
    eligible <- kind %in% c("stationary_candidate","nonsmooth_candidate") && neighborhood$continuity
    if (eligible && stationary && !is.null(d2) && d2$reliable && abs(d2$estimate)>atol) {
      classification <- if (d2$estimate>0) "minimum" else "maximum"
      evidence <- "conditional_second_derivative_test"
    } else if (eligible) {
      lp <- neighborhood$data$left_derivative
      rp <- neighborhood$data$right_derivative
      supported <- all(is.finite(c(lp,rp))) &&
        all(neighborhood$data$left_stable & neighborhood$data$right_stable)
      negative_left <- supported && all(lp < -atol)
      positive_left <- supported && all(lp > atol)
      negative_right <- supported && all(rp < -atol)
      positive_right <- supported && all(rp > atol)
      classification <- if (negative_left && positive_right) "minimum" else
        if (positive_left && negative_right) "maximum" else
          if ((positive_left && positive_right) || (negative_left && negative_right)) "not_extremum" else "inconclusive"
      if (classification!="inconclusive") evidence <- "conditional_neighborhood_sign_test"
    }
  }
  list(type="critical_classification",x0=x0,value=y,kind=kind,
       classification=classification,evidence=evidence,derivative=d,second_derivative=d2,
       neighborhood=neighborhood$data,continuity_evidence=neighborhood$continuity,
       note="Conditional numerical evidence only; smoothness and local extrema are not proved.")
}

#' Search for critical-point candidates on a finite interval
#' @inheritParams vc_classify_critical
#' @param interval Increasing pair of finite real search bounds.
#' @param n Integer number of grid samples, at least 5.
#' @param x_tol Positive root refinement tolerance; candidates within 10*x_tol
#'   are merged, keeping the smallest sampled derivative residual.
#' @param max_iter Positive integer limit for bracketed root refinement.
#'   Bounded local minimization uses stats::optimize() with its own stopping rule.
#' @return Plain list with type, interval, points, classifications (full lists),
#'   scan (x, value, derivative), flat_intervals, boundary_samples, complete
#'   (always FALSE), status and note. status distinguishes unavailable
#'   derivatives, found candidates and no detected candidates. refinements
#'   records method, bounds, location, residual, status and error message for
#'   every attempted refinement; acceptance here means residual acceptance,
#'   not a proof of criticality. points includes x, value, derivative, kind,
#'   classification and evidence. Boundaries are separate from interior points.
#' @details Uses vc_derivative(), sign-change root refinement and bounded
#'   minimization of absolute derivative near grid minima. All refined residuals
#'   are checked again. Low-derivative sampled runs are reported as possible
#'   flat intervals, not proved plateaus; isolated refined roots inside them
#'   are retained. Continuous stationary sets are not enumerated
#'   as infinitely many critical points. Neighborhood classification shares
#'   vc_classify_critical(). Known invalid stencils split derivative evidence;
#'   sign-change refinement assumes continuity and can miss discontinuities,
#'   even-multiplicity roots, narrow features or nonsmooth points. Completeness
#'   is never guaranteed. Candidates without derivatives are not stationary.
#' @examples
#' f <- function(x) x^3-3*x
#' (f |> vc_critical_points(interval=c(-2,2)))$points
#' @export
vc_critical_points <- function(f,interval,n=201L,h=NULL,x_tol=1e-8,
                               slope_tol=1e-6,max_iter=100L) {
  if (!is.function(f)) stop("'f' must be a function.")
  interval <- .vc_interval(interval)
  .vc_integer(n,"n",5)
  .vc_positive(x_tol,"x_tol")
  .vc_positive(slope_tol,"slope_tol",TRUE)
  .vc_integer(max_iter,"max_iter",1)
  .vc_optional_step(h,"h")
  grid <- seq(interval[1],interval[2],length.out=n)
  evaluate <- function(x) {
    d <- .vc_first(f,x,h)
    if (is.null(d)) NA_real_ else d$estimate
  }
  scan <- .vc_root_scan(evaluate,grid,slope_tol,x_tol,max_iter)
  labels <- ifelse(is.finite(scan$values) & abs(scan$values)<=slope_tol,"flat_candidate","other")
  flat <- .vc_runs(grid,labels)
  flat <- flat[flat$classification=="flat_candidate" & flat$samples>=3L,,drop=FALSE]
  roots <- scan$roots[scan$roots>interval[1] & scan$roots<interval[2]]
  details <- lapply(roots,function(x) vc_classify_critical(f,x,
    radius=min(diff(interval)/100,(x-interval[1])/2,(interval[2]-x)/2),h=h,slope_tol=slope_tol))
  points <- data.frame(x=numeric(),value=numeric(),derivative=numeric(),kind=character(),
                       classification=character(),evidence=character())
  if (length(details)) points <- do.call(rbind,lapply(details,function(a)
    data.frame(x=a$x0,value=a$value,derivative=if(is.null(a$derivative)) NA_real_ else a$derivative$estimate,
               kind=a$kind,classification=a$classification,evidence=a$evidence)))
  list(type="critical_points",interval=interval,points=points,classifications=details,
       scan=data.frame(x=grid,value=vapply(grid,function(x) .vc_safe_value(f,x),numeric(1)),derivative=scan$values),
       refinements=scan$refinements,flat_intervals=flat,boundary_samples=data.frame(x=interval,
         value=vapply(interval,function(x) .vc_safe_value(f,x),numeric(1))),complete=FALSE,
       status=if(!any(is.finite(scan$values))) "no_evaluable_derivatives" else
         if(nrow(points)) "candidates_found" else "no_candidates_detected",
       note="Finite-grid candidate search, not a complete list or proof of local/global extrema.")
}
