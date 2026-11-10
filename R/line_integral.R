#' Integrate a scalar or vector field along a parametric curve
#' @inheritParams vc_curve_length
#' @param f Function of one Cartesian vector, or of position and parameter
#'   when time_dependent is TRUE. Returns one finite scalar or a vector
#'   matching the curve dimension, according to type.
#' @param curve Function of one scalar parameter returning a finite real
#'   Cartesian vector of length two or three.
#' @param interval Two finite parameter endpoints in traversal order.
#' @param velocity Optional derivative of curve with the same dimension.
#' @param type Scalar line integral with respect to arc length, or vector
#'   work using the dot product with the curve derivative.
#' @param time_dependent Pass the curve parameter explicitly as f(position,t).
#' @return List with estimate, interval, dimension, orientation, quadrature
#'   diagnostics on increasing bounds, samples in traversal order, derivative
#'   source, and separate numerical Jacobian diagnostics in sample order.
#' @details Scalar integrals are invariant under interval reversal, including
#'   for signed scalar fields. Work changes sign. Reversal traverses the same
#'   parameterized field values, including when time_dependent is TRUE. Equal
#'   endpoints give zero without evaluating either function. Central stencils
#'   can leave the domain or cross breaks. Supply velocity at known boundaries
#'   and split nonsmooth pieces with breaks when appropriate. Quadrature
#'   diagnostics cover the sampled integrand and omit derivative error.
#' @examples
#' vc_line_integral(function(p)c(-p[2],p[1]),
#'   function(t)c(cos(t),sin(t)),c(0,2*pi),type="work",
#'   velocity=function(t)c(-sin(t),cos(t)))$estimate
#' @export
vc_line_integral <- function(f,curve,interval,type=c("scalar","work"),
    velocity=NULL,h=NULL,breaks=NULL,method=c("adaptive","trapezoid","simpson"),
    n=32L,abs_tol=1e-8,rel_tol=1e-8,subdivisions=100L,time_dependent=FALSE) {
  if (!is.function(f) || !is.function(curve)) stop("'f' and 'curve' must be functions.")
  if (!is.null(velocity) && !is.function(velocity)) stop("'velocity' must be a function or NULL.")
  if (!is.null(h)) .vc_positive(h,"h")
  .vc_flag(time_dependent,"time_dependent")
  type <- match.arg(type)
  method <- match.arg(method)
  interval <- .vc_integral_interval(interval)
  bounds <- .vc_application_bounds(interval,breaks)
  orientation <- if(type=="work" && interval[2]<interval[1]) -1 else 1
  cache <- new.env(parent=emptyenv())
  keys <- character()
  dimension <- NA_integer_
  position <- function(t) {
    p <- .vc_point(curve(t),"curve output")
    if (!length(p)%in%c(2L,3L)) stop("Curve must have two or three coordinates.")
    if (is.na(dimension)) dimension <<- length(p)
    if (length(p)!=dimension) stop("Curve dimension changed.")
    unname(p)
  }
  integrand <- function(t) {
    key <- sprintf("%a",as.double(t))
    if (exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache)$integrand)
    derivative <- NULL
    if (is.null(velocity)) {
      derivative <- vc_jacobian(function(p)position(p[1]),t,h=h,warn=FALSE)
      p <- derivative$value
      v <- derivative$estimate[,1]
    } else {
      p <- position(t)
      v <- .vc_point(velocity(t),"velocity output")
      if (length(v)!=dimension) stop("Velocity dimension must match curve dimension.")
    }
    field <- if (time_dependent) f(p,t) else f(p)
    speed <- .vc_finite(.vc_norm(v),"Curve speed")
    density <- NA_real_
    if (type=="scalar") {
      density <- .vc_finite(field,"Scalar field value")
      value <- density*speed
    } else {
      field <- .vc_point(field,"vector field output")
      if (length(field)!=dimension) stop("Field dimension must match curve dimension.")
      value <- sum(field*v)
    }
    value <- .vc_finite(value,"Line integrand")
    assign(key,list(t=t,point=c(p,rep(0,3L-dimension)),speed=speed,
      density=density,integrand=value,derivative=derivative),envir=cache)
    keys <<- c(keys,key)
    value
  }
  integral <- .vc_application_integrate(integrand,bounds,method,n,abs_tol,rel_tol,subdivisions)
  records <- mget(keys,envir=cache,inherits=FALSE)
  parameters <- vapply(records,function(r)r$t,numeric(1))
  order <- order(parameters,decreasing=interval[2]<interval[1])
  records <- unname(records[order])
  p <- t(vapply(records,function(r)r$point,numeric(3)))
  samples <- data.frame(t=unname(parameters[order]),x=p[,1],y=p[,2],z=p[,3],
    speed=vapply(records,function(r)r$speed,numeric(1)),
    density=vapply(records,function(r)r$density,numeric(1)),
    integrand=vapply(records,function(r)r$integrand,numeric(1)))
  derivatives <- if(is.null(velocity))lapply(records,function(r)r$derivative) else list()
  sensitivity <- if(length(derivatives))
    all(vapply(derivatives,function(d)d$all_reliable,logical(1))) else NA
  list(type="line_integral",kind=type,interval=interval,dimension=dimension,
       orientation=orientation,estimate=orientation*integral$estimate,
       integral=integral,samples=samples,derivative_source=if(is.null(velocity))"numerical" else "supplied",
       derivatives=derivatives,derivative_sensitivity_ok=sensitivity)
}

#' Display stored samples from a line integral
#' @param result Result of vc_line_integral().
#' @param title Optional character title.
#' @return Plotly widget. For planar scalar integrals, shows the sampled
#'   curve and the signed curtain up to the field height. Other integrals
#'   show the sampled path with local integrand values in hover text.
#' @details Samples come from quadrature. Adaptive quadrature may omit
#'   endpoints, and connecting samples approximates the curve. Equal bounds
#'   leave no samples to display. No functions are reevaluated.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   r <- vc_line_integral(function(p)p[1],function(t)c(t,t^2),
#'     c(-1,1),method="simpson")
#'   vc_plot_line_integral(r)
#' }
#' @export
vc_plot_line_integral <- function(result,title=NULL) {
  if (!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if (!is.list(result) || !identical(result$type,"line_integral")) stop("Invalid line integral result.")
  if (!is.null(title) && (!is.character(title) || length(title)!=1L || is.na(title))) stop("Invalid title.")
  d <- result$samples
  if (!nrow(d)) stop("No quadrature samples to display.")
  label <- paste("t =",format(d$t),"/ integrand =",format(result$orientation*d$integrand))
  p <- plotly::plot_ly(x=d$x,y=d$y,z=d$z,type="scatter3d",mode="lines+markers",
    name="Sampled path",text=label,hoverinfo="x+y+z+text")
  if (result$kind=="scalar" && result$dimension==2L && nrow(d)>1L) {
    p <- plotly::add_trace(p,x=rbind(d$x,d$x),y=rbind(d$y,d$y),z=rbind(d$z,d$density),
      type="surface",name="Signed scalar curtain",showscale=FALSE,opacity=.5,inherit=FALSE)
    p <- plotly::add_trace(p,x=d$x,y=d$y,z=d$density,type="scatter3d",
      mode="lines",name="Field height",inherit=FALSE)
  }
  plotly::layout(p,title=if(is.null(title))paste("Line integral =",format(result$estimate,digits=6)) else title,
    scene=list(aspectmode="data",xaxis=list(title="x"),yaxis=list(title="y"),zaxis=list(title="z")))
}
