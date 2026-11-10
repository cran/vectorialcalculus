#' Compare the two integrals in Green's circulation theorem
#' @param field Function of one Cartesian vector c(x,y) returning c(P,Q).
#' @param lower,upper Scalar functions defining the lower and upper y graphs.
#' @param interval Two finite, strictly increasing x endpoints.
#' @param lower_derivative,upper_derivative Optional derivatives of the y graphs.
#' @param rotation Optional scalar function of c(x,y) giving dQ/dx - dP/dy.
#'   Defaults to numerical partial derivatives from vc_jacobian().
#' @param orientation Boundary traversal. Clockwise negates the boundary
#'   circulation and the oriented interior integral.
#' @param boundary_method Rule for vc_line_integral().
#' @param n_boundary Initial fixed-rule cell count per boundary edge.
#' @inheritParams vc_double_integrate
#' @return List with circulation, signed interior integral, signed difference
#'   (circulation minus interior), comparison tolerance, both quadrature
#'   results, a summed error indicator, and separate derivative diagnostics.
#' @details The region is a single strip between two graphs. Counterclockwise
#'   traversal follows the lower graph to the right, the right edge upward,
#'   the upper graph to the left and the left edge downward. Collapsed
#'   vertical edges contribute zero. The caller supplies a suitable simple
#'   piecewise smooth boundary and a C1 field on a neighborhood of the region.
#'   Finite sampled boundary checks cannot establish those hypotheses.
#'   Central stencils may leave the domain. Numerical agreement and error
#'   indicators are diagnostic only. Regions with holes need decomposition.
#' @examples
#' vc_green(function(p)c(-p[2],p[1]),function(x)0,function(x)1,c(0,1),
#'   lower_derivative=function(x)0,upper_derivative=function(x)0,
#'   rotation=function(p)2,n_boundary=4,n_outer=2,n_inner=2)$circulation
#' @export
vc_green <- function(field,lower,upper,interval,lower_derivative=NULL,upper_derivative=NULL,
    rotation=NULL,orientation=c("counterclockwise","clockwise"),
    boundary_method=c("simpson","trapezoid","adaptive"),n_boundary=16L,
    method=c("simpson","trapezoid"),n_outer=8L,n_inner=8L,abs_tol=1e-8,rel_tol=1e-8) {
  if (!is.function(field) || !is.function(lower) || !is.function(upper))
    stop("'field', 'lower' and 'upper' must be functions.")
  for (f in list(lower_derivative,upper_derivative,rotation)) {
    if (!is.null(f) && !is.function(f)) stop("Supplied derivatives and rotation must be functions or NULL.")
  }
  interval <- .vc_integral_interval(interval)
  if (interval[1]>=interval[2]) stop("'interval' must be strictly increasing.")
  orientation <- match.arg(orientation)
  method <- match.arg(method)
  boundary_method <- match.arg(boundary_method)
  clockwise <- orientation=="clockwise"
  direction <- if(clockwise)-1 else 1
  evaluate <- function(p) {
    value <- .vc_point(field(p),"planar field output")
    if (length(value)!=2L) stop("Planar field must return two components.")
    unname(value)
  }
  heights <- function(x) {
    a <- .vc_value(lower,x)
    b <- .vc_value(upper,x)
    if (a>b) stop("Lower boundary exceeds upper boundary at a sampled section.")
    c(a,b)
  }
  a <- interval[1]
  b <- interval[2]
  left <- heights(a)
  right <- heights(b)
  edge <- function(curve,range,velocity=NULL) {
    if (clockwise) range <- rev(range)
    vc_line_integral(evaluate,curve,range,"work",velocity=velocity,
      method=boundary_method,n=n_boundary,abs_tol=abs_tol,rel_tol=rel_tol)
  }
  graph_velocity <- function(derivative) {
    if (is.null(derivative)) return(NULL)
    function(x)c(1,.vc_value(derivative,x))
  }
  boundary <- list(
    lower=edge(function(x)c(x,.vc_value(lower,x)),interval,graph_velocity(lower_derivative)),
    right=edge(function(y)c(b,y),right,function(y)c(0,1)),
    upper=edge(function(x)c(x,.vc_value(upper,x)),rev(interval),graph_velocity(upper_derivative)),
    left=edge(function(y)c(a,y),rev(left),function(y)c(0,1)))
  if (clockwise) boundary <- rev(boundary)
  circulation <- .vc_finite(sum(vapply(boundary,function(r)r$estimate,numeric(1))),"Circulation")
  boundary_error <- .vc_finite(sum(vapply(boundary,function(r)r$integral$error_estimate,numeric(1))),
    "Boundary error indicator")
  flags <- new.env(parent=emptyenv())
  curl_at <- function(p) {
    if (!is.null(rotation)) return(direction*.vc_finite(rotation(p),"Planar rotation"))
    key <- paste(sprintf("%a",p),collapse=";")
    if (exists(key,envir=flags,inherits=FALSE)) return(get(key,envir=flags)$value)
    j <- vc_jacobian(evaluate,p,warn=FALSE)
    value <- .vc_finite(direction*(j$estimate[2,1]-j$estimate[1,2]),"Planar rotation")
    assign(key,list(value=value,reliable=j$reliable[2,1]&&j$reliable[1,2]),envir=flags)
    value
  }
  interior <- vc_double_integrate(curl_at,interval,lower,upper,method=method,
    n_outer=n_outer,n_inner=n_inner,abs_tol=abs_tol,rel_tol=rel_tol)
  entries <- as.list(flags,all.names=TRUE)
  sensitivity <- if(length(entries)) all(vapply(entries,function(r)r$reliable,logical(1))) else NA
  difference <- .vc_finite(circulation-interior$estimate,"Green discrepancy")
  tolerance <- .vc_finite(max(abs_tol,rel_tol*max(abs(circulation),abs(interior$estimate))),
    "Comparison tolerance")
  error <- .vc_finite(boundary_error+interior$error_estimate,"Combined error indicator")
  list(type="green",interval=interval,orientation=orientation,circulation=circulation,
       interior_integral=interior$estimate,difference=difference,tolerance=tolerance,
       within_tolerance=abs(difference)<=tolerance,error_indicator=error,
       boundary=boundary,interior=interior,
       derivative_sensitivity=list(boundary=vapply(boundary,function(r)r$derivative_sensitivity_ok,logical(1)),
         interior=sensitivity,interior_evaluations=length(entries)),
       note="Agreement is a numerical comparison under caller-supplied Green theorem hypotheses.")
}

#' Display the stored region and oriented boundary for a Green comparison
#' @param result Result of vc_green().
#' @param region Result of vc_region() sampled from the same graphs and interval.
#' @param title Optional character title.
#' @return Plotly widget showing the region, sampled orientation arrows and
#'   both integral estimates. No field or boundary is reevaluated.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   lo <- function(x)0
#'   hi <- function(x)1
#'   result <- vc_green(function(p)c(-p[2],p[1]),lo,hi,c(0,1),
#'     rotation=function(p)2,n_outer=2,n_inner=2)
#'   vc_plot_green(result,vc_region(lo,hi,c(0,1)))
#' }
#' @export
vc_plot_green <- function(result,region,title=NULL) {
  if (!is.list(result) || !identical(result$type,"green")) stop("Invalid Green comparison.")
  if (!is.list(region) || !identical(region$type,"region")) stop("Invalid region sample.")
  if (!all(range(region$sections$outer)==result$interval))
    stop("Region and Green intervals must match.")
  p <- vc_plot_region(region,view="2d",show_grid=FALSE,
    title=if(is.null(title))"Green circulation comparison" else title)
  arrows <- list()
  for (edge in result$boundary) {
    d <- edge$samples
    if (nrow(d)<2L) next
    i <- max(1L,floor(.42*nrow(d)))
    j <- ceiling(.58*nrow(d))
    if (d$x[i]==d$x[j] && d$y[i]==d$y[j]) next
    arrows[[length(arrows)+1L]] <- list(x=d$x[j],y=d$y[j],
      ax=d$x[i],ay=d$y[i],xref="x",yref="y",axref="x",ayref="y",
      text="",showarrow=TRUE,arrowhead=3,arrowwidth=2,arrowcolor="#bb4728")
  }
  caption <- paste("Circulation =",format(result$circulation,digits=6),
    "| Interior =",format(result$interior_integral,digits=6),
    "| Difference =",format(result$difference,digits=3))
  arrows[[length(arrows)+1L]] <- list(text=caption,showarrow=FALSE,
    x=.5,y=0,xref="paper",yref="paper",yanchor="top",yshift=-30)
  plotly::layout(p,annotations=arrows,margin=list(b=85))
}
