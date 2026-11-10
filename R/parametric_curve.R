.vc_curve_vector <- function(f, t) {
  value <- f(t)
  if (!is.numeric(value) || is.complex(value) || !is.null(dim(value)) ||
      length(value) != 3L || any(!is.finite(value))) {
    stop("Curve and velocity functions must return a finite real three-vector.")
  }
  as.double(value)
}

#' Sample an oriented spatial parametric curve
#' @param curve Function of one scalar parameter returning a finite real
#'   three-vector in Cartesian x,y,z order.
#' @param interval Two finite parameter endpoints. Reversed endpoints reverse
#'   the sample order. Equal endpoints produce one sample.
#' @param n Integer number of samples, at least two. Equal endpoints override
#'   this count and produce a single sample.
#' @return List with type curve3d, interval and samples (t,x,y,z).
#' @details Samples preserve the supplied parameter order. A plotted polyline
#'   approximates the sampled curve. Coincident spatial points are allowed.
#' @examples
#' vc_curve_sample(function(t)c(cos(t),sin(t),t),c(0,2*pi),n=20)
#' @export
vc_curve_sample <- function(curve, interval, n=101L) {
  if (!is.function(curve)) stop("'curve' must be a function.")
  interval <- .vc_integral_interval(interval)
  .vc_integer(n, "n", 2)
  parameters <- if (interval[1]==interval[2]) interval[1] else
    seq(interval[1], interval[2], length.out=n)
  if (any(!is.finite(parameters)) || anyDuplicated(parameters))
    stop("Parameter grid is not representable.")
  points <- t(vapply(parameters, function(value).vc_curve_vector(curve,value), numeric(3)))
  list(type="curve3d", interval=interval,
       samples=data.frame(t=parameters,x=points[,1],y=points[,2],z=points[,3]))
}

#' Display stored spatial curve samples
#' @param result Result of vc_curve_sample().
#' @param mode Plotly mode for lines, markers or both.
#' @param line,marker Plotly style lists.
#' @param title Optional character title.
#' @return A Plotly widget in physical Cartesian coordinates. Parameter values
#'   are available in hover text. A one-point curve uses a marker.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   vc_curve_sample(function(t)c(t,t^2,t^3),c(0,1)) |> vc_plot_curve()
#' }
#' @export
vc_plot_curve <- function(result, mode=c("lines","markers","lines+markers"),
                          line=list(color="#24476b",width=4), marker=list(), title=NULL) {
  if (!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if (!is.list(result) || !identical(result$type,"curve3d")) stop("Invalid curve result.")
  mode <- match.arg(mode)
  if (!is.list(line) || !is.list(marker)) stop("Line and marker styles must be lists.")
  if (!is.null(title) && (!is.character(title) || length(title)!=1L || is.na(title)))
    stop("Invalid title.")
  d <- result$samples
  if (nrow(d)==1L) mode <- "markers"
  p <- plotly::plot_ly(x=d$x,y=d$y,z=d$z,type="scatter3d",
    mode=mode,line=if(mode!="markers")line else NULL,
    marker=if(mode!="lines")marker else NULL,
    text=paste("t =",format(d$t)),hoverinfo="x+y+z+text",name="Sampled curve")
  plotly::layout(p,title=if(is.null(title))"Sampled parametric curve" else title,
    scene=list(xaxis=list(title="x(t)"),yaxis=list(title="y(t)"),
      zaxis=list(title="z(t)"),aspectmode="data"))
}

#' Integrate speed along a spatial parametric curve
#' @inheritParams vc_curve_sample
#' @inheritParams vc_arc_length
#' @param interval Two finite parameter endpoints. Length uses their increasing
#'   order, regardless of traversal direction.
#' @param n Number of initial subintervals for fixed quadrature, as in
#'   vc_integrate(). Adaptive quadrature ignores this count.
#' @param velocity Optional function returning the Cartesian derivative as a
#'   finite real three-vector. The caller supplies its mathematical validity.
#' @return List with type curve_length, sorted interval, estimate, integral
#'   diagnostics, speed samples, derivative source, numerical Jacobian results
#'   (in sample order) and a separate derivative sensitivity flag.
#' @details Length uses increasing parameter bounds and is nonnegative.
#'   A zero interval gives zero without evaluations. With no supplied velocity,
#'   central Jacobian stencils can leave the interval or cross a break. Supply
#'   velocity when the curve is undefined at needed stencil points. Split known
#'   nonsmooth pieces using breaks and prefer adaptive quadrature at corners.
#'   The speed formula requires a suitable piecewise smooth parameterization.
#'   Quadrature diagnostics cover the sampled speed and omit derivative error.
#' @examples
#' vc_curve_length(function(t)c(cos(t),sin(t),t),c(0,2*pi),
#'   velocity=function(t)c(-sin(t),cos(t),1))
#' @export
vc_curve_length <- function(curve, interval, velocity=NULL, h=NULL, breaks=NULL,
    method=c("adaptive","trapezoid","simpson"), n=32L,
    abs_tol=1e-8, rel_tol=1e-8, subdivisions=100L) {
  if (!is.function(curve)) stop("'curve' must be a function.")
  if (!is.null(velocity) && !is.function(velocity)) stop("'velocity' must be a function or NULL.")
  if (!is.null(h)) .vc_positive(h,"h")
  bounds <- .vc_application_bounds(interval,breaks)
  method <- match.arg(method)
  cache <- new.env(parent=emptyenv())
  locations <- speeds <- numeric()
  derivatives <- list()
  speed <- function(t) {
    key <- sprintf("%a",as.double(t))
    if (exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache))
    if (is.null(velocity)) {
      derivative <- vc_jacobian(function(p).vc_curve_vector(curve,p[1]),t,h=h,warn=FALSE)
      v <- derivative$estimate[,1]
      derivatives[[length(derivatives)+1L]] <<- derivative
    } else {
      .vc_curve_vector(curve,t)
      v <- .vc_curve_vector(velocity,t)
    }
    value <- .vc_finite(.vc_norm(v),"Curve speed")
    locations <<- c(locations,t)
    speeds <<- c(speeds,value)
    assign(key,value,envir=cache)
    value
  }
  integral <- .vc_application_integrate(speed,bounds,method,n,abs_tol,rel_tol,subdivisions)
  if (integral$estimate<0) stop("Quadrature returned an invalid negative curve length.")
  sensitivity <- if(length(derivatives))
    all(vapply(derivatives,function(d)d$all_reliable,logical(1))) else NA
  list(type="curve_length",interval=bounds$interval,estimate=integral$estimate,
       integral=integral,samples=data.frame(t=locations,speed=speeds),
       derivative_source=if(is.null(velocity))"numerical" else "supplied",
       derivatives=derivatives,derivative_sensitivity_ok=sensitivity,
       note="Quadrature and derivative sensitivity provide separate numerical diagnostics.")
}
