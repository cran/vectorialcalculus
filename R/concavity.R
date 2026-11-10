#' Sample concavity and search for smooth inflection candidates
#' @inheritParams vc_critical_points
#' @inheritParams vc_classify_critical
#' @return Plain list with type, interval, scan (x, value, second_derivative,
#'   reliable, classification), regions (sample spans), inflections (x, value,
#'   left_sign, right_sign, status), singular_candidates (same columns),
#'   rejected_candidates, complete=FALSE, status
#'   (sampled or no_evaluable_curvature), refinements (attempted root
#'   refinements with residual/status/messages), note.
#'   Scan labels are concave_up, concave_down, near_zero or unresolved.
#' @details Reuses the direct three-point second-difference engine from Taylor,
#'   with three-scale sensitivity and roundoff diagnostics. Only reliable
#'   nonzero signs label concavity; region bounds cover sample spans, not proved
#'   continuous intervals. Inflection candidates require opposite consistent
#'   reliable curvature signs at three radii, contracting function-value
#'   neighborhood evidence and a stable numerical first derivative. A zero
#'   second difference alone is not an inflection. Singular and unresolved
#'   candidates remain separate. Continuous curvature sign changes with
#'   unresolved first-derivative evidence are reported in singular_candidates,
#'   not asserted to be smooth inflections. Default neighborhood curvature
#'   steps are capped at one percent of distance from the candidate to avoid
#'   straddling its possible singularity. This finite search is incomplete;
#'   it need not detect inflections with nonexistent or singular derivatives.
#'   Positive/negative f'' on an interval gives the usual concavity criterion;
#'   this function labels only sampled evidence. Second-difference O(h^2)
#'   truncation presumes sufficient fourth-order local smoothness and spacing.
#'   Numerical evidence does not prove continuity or a change of concavity.
#'   Invalid/warning-producing numerical probes are recorded as unavailable
#'   when their values fail validation; probing warnings are suppressed.
#' @examples
#' vc_concavity(function(x) x^3, c(-1,1))$inflections
#' f <- function(x) x^4
#' (f |> vc_concavity(interval=c(-1,1)))$inflections
#' @export
vc_concavity <- function(f,interval,n=201L,h=NULL,h2=NULL,x_tol=1e-8,
                         atol=1e-8,rtol=1e-4,max_iter=100L) {
  if (!is.function(f)) stop("'f' must be a function.")
  interval <- .vc_interval(interval)
  .vc_integer(n,"n",5)
  .vc_positive(x_tol,"x_tol")
  .vc_positive(atol,"atol",TRUE)
  .vc_positive(rtol,"rtol",TRUE)
  .vc_optional_step(h,"h")
  .vc_optional_step(h2,"h2")
  .vc_integer(max_iter,"max_iter",1)
  grid <- seq(interval[1],interval[2],length.out=n)
  evaluate <- function(x) {
    d <- .vc_second(f,x,h2,atol,rtol)
    if (is.null(d)) NA_real_ else d$estimate
  }
  scan <- .vc_root_scan(evaluate,grid,atol,x_tol,max_iter)
  diagnostics <- lapply(grid,function(x) .vc_second(f,x,h2,atol,rtol))
  reliable <- vapply(diagnostics,function(d) !is.null(d) && d$reliable,logical(1))
  labels <- ifelse(!reliable,"unresolved",ifelse(scan$values>atol,"concave_up",
                   ifelse(scan$values < -atol,"concave_down","near_zero")))
  regions <- .vc_runs(grid,labels)
  inflections <- data.frame(x=numeric(),value=numeric(),left_sign=integer(),right_sign=integer(),status=character())
  singular <- inflections
  rejected <- data.frame(x=numeric(),reason=character())
  roots <- scan$roots[scan$roots>interval[1] & scan$roots<interval[2]]
  for (x in roots) {
    radius <- min(diff(interval)/100,(x-interval[1])/2,(interval[2]-x)/2)
    neighborhood <- .vc_continuity_evidence(f,x,radius)
    radii <- radius/c(1,2,4)
    sign_at <- function(z) {
      probe_step <- if(is.null(h2)) min(max(1,abs(z))*.Machine$double.eps^(1/4),abs(z-x)/100) else h2
      d <- .vc_second(f,z,probe_step,atol,rtol)
      if (is.null(d) || !d$reliable || abs(d$estimate)<=atol) 0L else as.integer(sign(d$estimate))
    }
    left <- vapply(x-radii,sign_at,integer(1))
    right <- vapply(x+radii,sign_at,integer(1))
    d1 <- .vc_first(f,x,h,atol,rtol)
    c0 <- vc_classify_critical(f,x,radius=radius,h=h,h2=h2,atol=atol,rtol=rtol)
    smooth <- !is.null(d1) && d1$stable && c0$kind!="nonsmooth_candidate"
    opposite <- all(left!=0L) && all(right!=0L) && length(unique(left))==1L &&
      length(unique(right))==1L && left[1]!=right[1]
    if (neighborhood$continuity && smooth && opposite) {
      inflections <- rbind(inflections,data.frame(x=x,value=.vc_value(f,x),left_sign=left[1],
        right_sign=right[1],status="conditional_smooth_inflection_candidate"))
    } else if (neighborhood$continuity && opposite) {
      singular <- rbind(singular,data.frame(x=x,value=.vc_value(f,x),left_sign=left[1],right_sign=right[1],
        status="concavity_change_with_unresolved_differentiability"))
    } else rejected <- rbind(rejected,data.frame(x=x,reason="no_consistent_continuous_concavity_change"))
  }
  list(type="concavity",interval=interval,
       scan=data.frame(x=grid,value=vapply(grid,function(x) .vc_safe_value(f,x),numeric(1)),
                       second_derivative=scan$values,reliable=reliable,classification=labels),
       diagnostics=diagnostics,regions=regions,inflections=inflections,singular_candidates=singular,rejected_candidates=rejected,
       refinements=scan$refinements,complete=FALSE,status=if(!any(is.finite(scan$values))) "no_evaluable_curvature" else "sampled",
       note="Sampled concavity and conditional smooth inflection candidates, not proofs.")
}
