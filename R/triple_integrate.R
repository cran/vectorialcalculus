#' Integrate over a nested Cartesian solid
#' @param f Scalar field accepting one vector in physical x,y,z order.
#' @param interval Finite outer bounds, reversed or equal allowed.
#' @param middle_lower,middle_upper Functions of the outer coordinate giving
#'   ordered middle bounds.
#' @param inner_lower,inner_upper Functions accepting a two-vector in
#'   outer,middle order, giving ordered inner bounds.
#' @param order Integration order from inner to outer; a permutation of xyz.
#' @param method Simpson or trapezoid fixed quadrature.
#' @param n_outer,n_middle,n_inner Initial cell counts, each refined internally.
#' @param abs_tol,rel_tol Nonnegative diagnostic tolerances, one positive.
#' @return List containing estimate, error_estimate, tolerance_met, outer (full
#'   vc_integrate result), slices (cached full vc_double_integrate results),
#'   outer_nodes, projection, order and interval.
#' @details Reuses unchanged scalar and double quadrature engines. Absolute
#'   finest outer weights propagate each double-slice error indicator; this
#'   sum plus outer refinement error is heuristic, not a certified bound.
#'   All slice flags also gate tolerance_met. Finite sampling does not prove
#'   Fubini hypotheses, integrability or ordering between sampled sections.
#'   Caller supplies a justified nested region and splits known singularities.
#'   Reversed outer bounds change sign; reversed sampled middle/inner bounds
#'   are rejected. Equal widths skip deeper functions. Signed integrals need
#'   not be geometric volumes. Storage and evaluations grow cubically with
#'   resolution. No cylindrical/spherical transformation or adaptive rule.
#' @examples
#' (function(v)1) |> vc_triple_integrate(c(0,1),function(x)0,function(x)1-x,
#'   function(v)0,function(v)1-sum(v),n_outer=2,n_middle=2,n_inner=2)
#' @export
vc_triple_integrate <- function(f,interval,middle_lower,middle_upper,
                                inner_lower,inner_upper,
                                order=c("zyx","zxy","yzx","yxz","xzy","xyz"),
                                method=c("simpson","trapezoid"),
                                n_outer=8L,n_middle=8L,n_inner=8L,
                                abs_tol=1e-8,rel_tol=1e-8) {
  if(!all(vapply(list(f,middle_lower,middle_upper,inner_lower,inner_upper),
    is.function,logical(1)))) stop("Field and boundaries must be functions.")
  order <- match.arg(order)
  method <- match.arg(method)
  interval <- .vc_integral_interval(interval)
  for(n in list(n_outer,n_middle,n_inner)) {
    .vc_integer(n,"cell count",if(method=="simpson")2 else 1)
    if(n>.Machine$integer.max/2 || (method=="simpson" && n%%2!=0))
      stop("Invalid cell count for refinement/rule.")
  }
  indices <- match(rev(strsplit(order,"")[[1]]),c("x","y","z"))
  nodes <- numeric()
  slices <- bounds <- list()
  slice <- function(t) {
    i <- match(t,nodes)
    if(!is.na(i)) return(slices[[i]]$estimate)
    a <- .vc_value(middle_lower,t)
    b <- .vc_value(middle_upper,t)
    if(a>b) stop("Middle lower bound exceeds upper bound at sampled section.")
    field <- function(v) {
      p <- numeric(3)
      p[indices] <- c(t,v[1],v[2])
      .vc_field_value(f,p)
    }
    z <- vc_double_integrate(field,c(a,b),
      function(s).vc_field_value(inner_lower,c(t,s)),
      function(s).vc_field_value(inner_upper,c(t,s)),
      method=method,n_outer=n_middle,n_inner=n_inner,abs_tol=abs_tol,rel_tol=rel_tol)
    nodes <<- c(nodes,t)
    slices[[length(slices)+1L]] <<- z
    bounds[[length(bounds)+1L]] <<- c(a,b)
    z$estimate
  }
  outer <- vc_integrate(slice,interval,method,n_outer,abs_tol,rel_tol)
  propagated <- 0
  if(!is.null(outer$table)) {
    ids <- match(outer$table$x,nodes)
    propagated <- .vc_finite(sum(abs(outer$table$weight)*
      vapply(slices[ids],function(z)z$error_estimate,numeric(1))),"Slice error indicator")
  }
  error <- .vc_finite(outer$error_estimate+propagated,"Combined error indicator")
  projection <- list(type="double_integral",order="yx",sections=data.frame(outer=nodes,
    lower=vapply(bounds,function(z)z[1],numeric(1)),
    upper=vapply(bounds,function(z)z[2],numeric(1))))
  list(type="triple_integral",estimate=outer$estimate,error_estimate=error,
    outer_error=outer$error_estimate,slice_error=propagated,tolerance=outer$tolerance,
    tolerance_met=outer$tolerance_met && error<=outer$tolerance &&
      all(vapply(slices,function(z)z$tolerance_met,logical(1))),
    outer=outer,slices=slices,outer_nodes=nodes,projection=projection,
    order=order,interval=interval,method=method,
    status=if(length(nodes))"nested_fixed_diagnostic" else "zero_interval",
    note="Signed integral; triple error indicators are not certified bounds.")
}

#' Display stored sections of a Cartesian triple-integration region
#' @param result Result of vc_triple_integrate().
#' @param view Stored 3D inner-coordinate sections or outer-middle projection.
#' @param title Optional character title.
#' @return Plotly widget without new field/boundary evaluations or quadrature.
#' @details Sections are stored fibers, not a watertight surface or quadrature
#'   cells. The projection joins sampled middle boundaries linearly, including
#'   degenerate fibers; it is a base-domain envelope, not a certified projection
#'   of the positive-volume support. Axis labels preserve physical coordinates.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   r <- vc_triple_integrate(sum,c(0,1),function(x)0,function(x)1,
#'     function(v)0,function(v)1,n_outer=2,n_middle=2,n_inner=2)
#'   vc_plot_triple_integral(r)
#' }
#' @export
vc_plot_triple_integral <- function(result,view=c("sections","projection"),title=NULL) {
  if(!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"triple_integral")) stop("Invalid triple result.")
  view <- match.arg(view)
  if(!is.null(title)&&(!is.character(title)||length(title)!=1||is.na(title))) stop("Invalid title.")
  labels <- rev(strsplit(result$order,"")[[1]])
  if(!length(result$outer_nodes)) stop("No stored sections to display.")
  if(view=="projection") {
    p <- vc_plot_double_integral(result$projection,title)
    return(plotly::layout(p,xaxis=list(title=labels[1]),yaxis=list(title=labels[2])))
  }
  indices <- match(labels,c("x","y","z"))
  xyz <- lapply(seq_along(result$slices),function(i) {
    d <- result$slices[[i]]$sections
    if(!nrow(d)) return(matrix(numeric(),0,3))
    a <- as.vector(rbind(rep(result$outer_nodes[i],nrow(d)),rep(result$outer_nodes[i],nrow(d)),NA_real_))
    b <- as.vector(rbind(d$outer,d$outer,NA_real_))
    c <- as.vector(rbind(d$lower,d$upper,NA_real_))
    m <- matrix(NA_real_,length(a),3)
    m[,indices] <- cbind(a,b,c)
    m
  })
  m <- do.call(rbind,xyz)
  if(!nrow(m)) stop("No inner fibers to display.")
  p <- plotly::plot_ly(x=m[,1],y=m[,2],z=m[,3],type="scatter3d",mode="lines",
    name="Stored inner fibers")
  plotly::layout(p,title=if(is.null(title))"Cartesian solid: sampled fibers" else title,
    scene=list(xaxis=list(title="x"),yaxis=list(title="y"),zaxis=list(title="z"),aspectmode="data"))
}
