#' Integrate area, scalar density or flux on a parameterized surface
#' @param surface Function of one vector c(u,v) returning a finite Cartesian
#'   three-vector.
#' @param u_range,v_range Two finite endpoints for each parameter. Integration
#'   uses increasing bounds. Equal bounds give zero without surface evaluations.
#' @param f NULL for area, a scalar density for scalar integration, or a
#'   three-component vector field for flux. Functions accept one Cartesian vector.
#' @param type Scalar integration or oriented flux.
#' @param jacobian Optional function of c(u,v) returning a finite real 3 by 2
#'   matrix, with surface derivatives as columns. Defaults to vc_jacobian().
#' @param orientation Forward uses the u derivative crossed with the v
#'   derivative. Reverse negates flux. It has no effect on scalar integration.
#' @param h Optional numerical Jacobian step or two coordinate steps.
#' @inheritParams vc_double_integrate
#' @return List with estimate, full nested quadrature diagnostics, sample
#'   parameters, positions, area factors, signed integrands and numerical
#'   derivative sensitivity. Supplied derivatives have NA sensitivity.
#' @details Scalar integration uses the norm of the cross product and counts
#'   parameterization multiplicity. Flux uses the oriented, unnormalized
#'   cross product. A zero cross product contributes zero without requiring
#'   a unit normal. Numerical stencils can leave the parameter domain. The
#'   quadrature diagnostics omit derivative error. Samples do not establish
#'   smoothness, regularity, injectivity or a globally consistent orientation.
#' @examples
#' vc_surface_integral(function(p)c(p[1],p[2],p[1]+p[2]),c(0,1),c(0,1),
#'   n_outer=2,n_inner=2)$estimate
#' @export
vc_surface_integral <- function(surface,u_range,v_range,f=NULL,type=c("scalar","flux"),
    jacobian=NULL,orientation=c("forward","reverse"),h=NULL,
    method=c("simpson","trapezoid"),n_outer=8L,n_inner=8L,abs_tol=1e-8,rel_tol=1e-8) {
  if (!is.function(surface)) stop("'surface' must be a function.")
  if (!is.null(f) && !is.function(f)) stop("'f' must be a function or NULL.")
  if (!is.null(jacobian) && !is.function(jacobian)) stop("'jacobian' must be a function or NULL.")
  type <- match.arg(type)
  orientation <- match.arg(orientation)
  method <- match.arg(method)
  if (type=="flux" && is.null(f)) stop("Flux requires a vector field.")
  u_range <- sort(.vc_integral_interval(u_range))
  v_range <- sort(.vc_integral_interval(v_range))
  h <- .vc_matrix_step(c(0,0),h)
  sign <- if (orientation=="forward") 1 else -1
  position <- function(p).vc_spatial_field(surface(p))
  cache <- new.env(parent=emptyenv())
  keys <- character()
  integrand <- function(p) {
    key <- paste(sprintf("%a",p),collapse=";")
    if (exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache)$integrand)
    reliable <- NA
    if (is.null(jacobian)) {
      derivative <- vc_jacobian(position,p,h=h,warn=FALSE)
      xyz <- derivative$value
      j <- derivative$estimate
      reliable <- derivative$all_reliable
    } else {
      xyz <- position(p)
      j <- jacobian(p)
      if (!is.numeric(j) || is.complex(j) || !identical(dim(j),c(3L,2L)) ||
          any(!is.finite(j))) stop("Jacobian must be a finite real 3 by 2 matrix.")
    }
    normal <- .vc_curve_cross(j[,1],j[,2])
    if (any(!is.finite(normal))) stop("Surface cross product is not finite.")
    area <- .vc_finite(.vc_norm(normal),"Surface area factor")
    value <- if (type=="flux") sign*sum(.vc_spatial_field(f(xyz))*normal) else
      area*(if(is.null(f))1 else .vc_finite(f(xyz),"Scalar density"))
    value <- .vc_finite(value,"Surface integrand")
    assign(key,data.frame(u=p[1],v=p[2],x=xyz[1],y=xyz[2],z=xyz[3],
      area_factor=area,integrand=value,reliable=reliable),envir=cache)
    keys <<- c(keys,key)
    value
  }
  integral <- vc_double_integrate(integrand,u_range,function(u)v_range[1],
    function(u)v_range[2],method=method,n_outer=n_outer,n_inner=n_inner,
    abs_tol=abs_tol,rel_tol=rel_tol)
  samples <- data.frame(u=numeric(),v=numeric(),x=numeric(),y=numeric(),z=numeric(),
    area_factor=numeric(),integrand=numeric(),reliable=logical())
  if (length(keys)) samples <- do.call(rbind,unname(mget(keys,envir=cache,inherits=FALSE)))
  sensitivity <- if (is.null(jacobian) && nrow(samples)) all(samples$reliable) else NA
  list(type="surface_integral",kind=type,orientation=orientation,u_range=u_range,
       v_range=v_range,estimate=integral$estimate,integral=integral,samples=samples,
       derivative_source=if(is.null(jacobian))"numerical" else "supplied",
       derivative_sensitivity_ok=sensitivity)
}

#' Integrate a scalar density or oriented flux on a graph
#' @param graph Function of one vector c(x,y) returning the scalar height.
#' @param x_range,y_range Finite endpoint pairs for the Cartesian base rectangle.
#' @param gradient Optional function of c(x,y) returning the two height partials.
#' @param orientation Up uses the normal with positive z component. Down
#'   negates flux. Scalar integration is unchanged.
#' @inheritParams vc_surface_integral
#' @return The vc_surface_integral() result with up/down orientation.
#' @examples
#' vc_surface_graph_integral(function(p)sum(p),c(0,1),c(0,1),
#'   f=function(p)c(0,0,1),type="flux",n_outer=2,n_inner=2)$estimate
#' @export
vc_surface_graph_integral <- function(graph,x_range,y_range,f=NULL,type=c("scalar","flux"),
    gradient=NULL,orientation=c("up","down"),h=NULL,
    method=c("simpson","trapezoid"),n_outer=8L,n_inner=8L,abs_tol=1e-8,rel_tol=1e-8) {
  if (!is.function(graph)) stop("'graph' must be a function.")
  if (!is.null(gradient) && !is.function(gradient)) stop("'gradient' must be a function or NULL.")
  orientation <- match.arg(orientation)
  j <- if (!is.null(gradient)) function(p) {
    g <- .vc_point(gradient(p),"gradient output")
    if (length(g)!=2L) stop("Gradient output must have two components.")
    rbind(diag(2),g)
  }
  result <- vc_surface_integral(function(p)c(p,.vc_finite(graph(p),"Graph height")),
    x_range,y_range,f,match.arg(type),j,if(orientation=="up")"forward" else "reverse",
    h,match.arg(method),n_outer,n_inner,abs_tol,rel_tol)
  result$orientation <- orientation
  result
}
