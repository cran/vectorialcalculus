.vc_application_bounds <- function(interval,breaks) {
  interval <- sort(.vc_integral_interval(interval))
  if(is.null(breaks)) breaks <- numeric()
  if(!is.numeric(breaks)||is.complex(breaks)||!is.null(dim(breaks))||any(!is.finite(breaks))||
     any(diff(breaks)<=0)||any(breaks<=interval[1]|breaks>=interval[2]))
    stop("'breaks' must be finite, strictly increasing interior points.")
  list(interval=interval,points=c(interval[1],as.double(breaks),interval[2]))
}

.vc_application_integrate <- function(f,bounds,method,n,abs_tol,rel_tol,subdivisions) {
  # Validate controls even for degenerate intervals; divide absolute budget.
  vc_integrate(f,rep(bounds$interval[1],2),method,n,abs_tol,rel_tol,subdivisions)
  points <- bounds$points
  k <- length(points)-1L
  results <- lapply(seq_len(k),function(i)
    vc_integrate(f,points[c(i,i+1L)],method,n,abs_tol/k,rel_tol,subdivisions))
  estimate <- .vc_finite(sum(vapply(results,function(z)z$estimate,numeric(1))),"Application estimate")
  error <- .vc_finite(sum(vapply(results,function(z)z$error_estimate,numeric(1))),"Error indicator")
  tolerance <- .vc_finite(max(abs_tol,rel_tol*abs(estimate)),"Tolerance")
  list(estimate=estimate,error_estimate=error,tolerance=tolerance,
       tolerance_met=all(vapply(results,function(z)z$tolerance_met,logical(1))) && error<=tolerance,
       results=results,note="Quadrature diagnostics are not certified error bounds.")
}

.vc_boundary <- function(g) {
  if(is.null(g)) return(function(t)0)
  if(!is.function(g)) stop("'g' must be NULL or a function.")
  g
}

#' Compute geometric area between two graphs
#' @inheritParams vc_integrate
#' @param g Second scalar boundary function; NULL means zero.
#' @param breaks Optional finite strictly increasing interior integration split
#'   points, for example known intersections. These are not found automatically.
#' @return Plain list with type, sorted interval, breaks, method, estimate
#'   (geometric area), signed_estimate (integral of f-g on increasing bounds),
#'   area and signed (complete segmented quadrature diagnostics), and note.
#' @details Integrates `abs(f-g)`, not `abs(integral(f-g))`. Curve ordering may
#'   change. Bounds are sorted because area is unoriented; equal bounds give
#'   zero. Sampled values and totals must be finite. Known crossings can be
#'   supplied as breaks to help resolve the nonsmooth absolute difference.
#'   Finite sampling does not prove all crossings or features were resolved.
#' @examples
#' identity |> vc_area(c(-1,1),breaks=0)
#' @export
vc_area <- function(f,interval,g=NULL,breaks=NULL,method=c("adaptive","trapezoid","simpson"),
                    n=32L,abs_tol=1e-8,rel_tol=1e-8,subdivisions=100L) {
  if(!is.function(f)) stop("'f' must be a function.")
  g <- .vc_boundary(g)
  bounds <- .vc_application_bounds(interval,breaks)
  method <- match.arg(method)
  difference <- function(t).vc_finite(.vc_value(f,t)-.vc_value(g,t),"Boundary difference")
  area <- .vc_application_integrate(function(t)abs(difference(t)),bounds,method,n,abs_tol,rel_tol,subdivisions)
  if(area$estimate<0) stop("Quadrature returned an invalid negative geometric area.")
  signed <- .vc_application_integrate(difference,bounds,method,n,abs_tol,rel_tol,subdivisions)
  list(type="area",interval=bounds$interval,breaks=breaks,method=method,estimate=area$estimate,
       signed_estimate=signed$estimate,area=area,signed=signed,
       note="Geometric area is integral(abs(f-g)); signed_estimate uses increasing bounds.")
}

.vc_revolution_section <- function(f,g,t,axis_at,method) {
  a <- .vc_value(f,t)
  b <- .vc_value(g,t)
  lo <- min(a,b)
  hi <- max(a,b)
  if(method=="washers") {
    u <- .vc_finite(lo-axis_at,"Radius displacement")
    v <- .vc_finite(hi-axis_at,"Radius displacement")
    outer <- max(abs(u),abs(v))
    inner <- if(u<=0&&v>=0)0 else min(abs(u),abs(v))
    density <- if(outer==inner)0 else pi*((outer-inner)*outer)*(1+inner/outer)
  } else {
    outer <- .vc_finite(abs(t-axis_at),"Shell radius")
    inner <- NA_real_
    density <- 2*pi*outer*.vc_finite(hi-lo,"Shell height")
  }
  list(lower=lo,upper=hi,outer=outer,inner=inner,density=.vc_finite(density,"Volume density"))
}

#' Compute a volume of revolution from a planar strip region
#' @inheritParams vc_area
#' @param variable Independent coordinate of f and g: x means y=f(x), y=g(x);
#'   y means x=f(y), x=g(y).
#' @param axis Physical rotation axis: horizontal means y=axis_at; vertical
#'   means x=axis_at.
#' @param axis_at Finite scalar position of the axis.
#' @param quadrature Integration rule: adaptive, trapezoid or simpson.
#' @return Plain list with type, sorted interval, breaks, variable, axis,
#'   axis_at, method (washers or shells), quadrature (rule name), estimate
#'   (nonnegative volume), integral (full segmented quadrature diagnostics),
#'   and note. Disks are washers with zero inner radius.
#' @details The region consists of segments between f and g at each independent
#'   coordinate, with g=0 by default. Washers apply for x/horizontal or
#'   y/vertical; shells apply for x/vertical or y/horizontal. Washer radii
#'   account for a region crossing its axis: the inner radius is zero and the
#'   outer radius is the larger distance. Otherwise use the two endpoint
#'   distances. No subtraction of squared signed heights is used.
#'   Shells integrate `2*pi*abs(t-axis_at)*abs(f(t)-g(t))`; the interval must be
#'   entirely on one side of the axis, including an endpoint on it. Crossing
#'   the axis is rejected because shells from opposite sides can overlap;
#'   simply splitting and summing does not generally fix double counting.
#'   Sorted/equal bounds follow geometric-volume conventions. Quadrature may
#'   miss features; neither integration nor plots certify the region hypotheses.
#'   To recover rotation of the region between f and a shifted horizontal
#'   axis, explicitly supply g=function(t)axis_at (default g is zero).
#' @examples
#' identity |> vc_revolution(c(0,1)) # cone: pi/3
#' (function(x)sqrt(1-x^2)) |> vc_revolution(c(-1,1)) # sphere: 4*pi/3
#' @export
vc_revolution <- function(f,interval,g=NULL,variable=c("x","y"),axis=c("horizontal","vertical"),
                          axis_at=0,breaks=NULL,quadrature=c("adaptive","trapezoid","simpson"),
                          n=32L,abs_tol=1e-8,rel_tol=1e-8,subdivisions=100L) {
  if(!is.function(f)) stop("'f' must be a function.")
  g <- .vc_boundary(g)
  variable <- match.arg(variable)
  axis <- match.arg(axis)
  .vc_finite(axis_at,"axis_at")
  bounds <- .vc_application_bounds(interval,breaks)
  quadrature <- match.arg(quadrature)
  method <- if((variable=="x"&&axis=="horizontal")||(variable=="y"&&axis=="vertical"))"washers" else "shells"
  if(method=="shells"&&bounds$interval[1]<axis_at&&bounds$interval[2]>axis_at)
    stop("Shell interval crosses the axis; overlapping shells could double-count volume.")
  density <- function(t).vc_revolution_section(f,g,t,axis_at,method)$density
  integral <- .vc_application_integrate(density,bounds,quadrature,n,abs_tol,rel_tol,subdivisions)
  if(integral$estimate<0) stop("Quadrature returned an invalid negative volume.")
  list(type="revolution",interval=bounds$interval,breaks=breaks,variable=variable,axis=axis,
    axis_at=as.double(axis_at),method=method,quadrature=quadrature,estimate=integral$estimate,
    integral=integral,note="Unoriented volume; shell intervals must stay on one side of the axis.")
}

.vc_curve_application <- function(f,interval,df,h,axis_at,breaks,method,n,abs_tol,rel_tol,subdivisions) {
  if(!is.function(f)||(!is.null(df)&&!is.function(df))) stop("'f' and optional 'df' must be functions.")
  .vc_optional_step(h,"h")
  bounds <- .vc_application_bounds(interval,breaks)
  locations <- numeric()
  derivatives <- list()
  slope <- function(t) {
    i <- match(t,locations)
    if(!is.na(i)) return(derivatives[[i]]$estimate)
    .vc_value(f,t)
    d <- if(is.null(df)) vc_derivative(f,t,h=h,warn=FALSE) else
      list(estimate=.vc_value(df,t),reliable=NA)
    locations <<- c(locations,t)
    derivatives[[length(derivatives)+1L]] <<- d
    d$estimate
  }
  integrand <- function(t) {
    d <- slope(t)
    scale <- max(1,abs(d))
    speed <- scale*sqrt((1/scale)^2+(d/scale)^2)
    if(!is.null(axis_at)) speed <- speed*(2*pi*abs(.vc_finite(.vc_value(f,t)-axis_at,"Radius displacement")))
    .vc_finite(speed,"Curve density")
  }
  integral <- .vc_application_integrate(integrand,bounds,method,n,abs_tol,rel_tol,subdivisions)
  if(integral$estimate<0) stop("Quadrature returned an invalid negative curve measure.")
  sensitivity <- if(is.null(df)&&length(derivatives))
    all(vapply(derivatives,function(d)d$reliable,logical(1))) else NA
  list(type=if(is.null(axis_at))"arc_length" else "surface_area",interval=bounds$interval,
    breaks=breaks,method=method,axis_at=axis_at,estimate=integral$estimate,integral=integral,
    derivative_source=if(is.null(df))"numerical" else "supplied",
    derivative_sensitivity_ok=sensitivity,locations=locations,derivatives=derivatives,
    note="Requires a rectifiable graph and appropriate derivative; numerical diagnostics do not prove these hypotheses.")
}

#' Approximate the length of a graph
#' @inheritParams vc_area
#' @param df Optional scalar derivative function, assumed mathematically
#'   correct. NULL reuses vc_derivative() with central differences.
#' @param h Optional positive finite numerical derivative step.
#' @return Plain list with type, sorted interval, breaks, method, estimate,
#'   integral (full quadrature diagnostics), derivative_source,
#'   derivative_sensitivity_ok (NA for supplied derivatives or no samples),
#'   locations, derivatives (full diagnostics per sampled point), and note.
#' @details Integrates `sqrt(1+f'(t)^2)` for a rectifiable graph. Known corners
#'   can be supplied as breaks. Central numerical stencils may leave the domain
#'   near endpoints or cross corners; supply df where appropriate. Nonfinite
#'   sampled derivatives are errors. Failed sensitivity flags remain visible
#'   even if quadrature meets its tolerance. A supplied df is not verified.
#'   Quadrature error estimates do not include numerical derivative error;
#'   derivative_sensitivity_ok is a separate heuristic, not a bound.
#' @examples
#' identity |> vc_arc_length(c(0,1),df=function(x)1)
#' @export
vc_arc_length <- function(f,interval,df=NULL,h=NULL,breaks=NULL,
                          method=c("adaptive","trapezoid","simpson"),n=32L,
                          abs_tol=1e-8,rel_tol=1e-8,subdivisions=100L) {
  .vc_curve_application(f,interval,df,h,NULL,breaks,match.arg(method),n,abs_tol,rel_tol,subdivisions)
}

#' Approximate the lateral area generated by revolving a graph
#' @inheritParams vc_arc_length
#' @param axis_at Finite position of the parallel rotation axis in the dependent
#'   coordinate: y=axis_at for y=f(x), or x=axis_at for x=f(y).
#' @param end_caps Logical; add two circular end disks with endpoint radii
#'   abs(f(t)-axis_at). This assumes the filled region between this graph and
#'   its axis. It is not the cap rule for a two-graph annular region.
#' @return A surface_area list with the same diagnostics as vc_arc_length()
#'   and axis_at, lateral_area, cap_area, end_caps. Estimate is lateral area
#'   plus requested caps; integral diagnostics describe only the lateral part.
#' @details Integrates `2*pi*abs(f(t)-axis_at)*sqrt(1+f'(t)^2)`. This is the
#'   surface generated by one graph, not automatically the entire boundary of
#'   a washer/shell solid between two graphs. Parameterization must cover the
#'   surface without multiplicity for geometric-area interpretation. Caps are
#'   excluded by default. Equal bounds give zero lateral and cap areas without
#'   evaluating endpoints; coincident end disks are not counted twice.
#' @examples
#' identity |> vc_surface_area(c(0,1),df=function(x)1) # pi*sqrt(2)
#' @export
vc_surface_area <- function(f,interval,axis_at=0,df=NULL,h=NULL,breaks=NULL,end_caps=FALSE,
                            method=c("adaptive","trapezoid","simpson"),n=32L,
                            abs_tol=1e-8,rel_tol=1e-8,subdivisions=100L) {
  .vc_finite(axis_at,"axis_at")
  if(!is.logical(end_caps)||length(end_caps)!=1||is.na(end_caps)) stop("Invalid 'end_caps'.")
  out <- .vc_curve_application(f,interval,df,h,axis_at,breaks,match.arg(method),n,abs_tol,rel_tol,subdivisions)
  out$lateral_area <- out$estimate
  out$cap_area <- 0
  if(end_caps&&out$interval[1]<out$interval[2]) {
    radii <- vapply(out$interval,function(t)abs(.vc_finite(.vc_value(f,t)-axis_at,"Cap radius")),numeric(1))
    out$cap_area <- .vc_finite(sum(pi*radii*radii),"Cap area")
  }
  out$end_caps <- end_caps
  out$estimate <- .vc_finite(out$lateral_area+out$cap_area,"Total surface area")
  out
}
