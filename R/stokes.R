#' Compare circulation and curl flux in Stokes' theorem
#' @inheritParams vc_surface_integral
#' @param field Spatial vector field accepting one Cartesian vector.
#' @param curl Optional function of position giving the three Cartesian curl
#'   components. Defaults to vc_curl().
#' @param u_range,v_range Finite strictly increasing parameter endpoint pairs.
#' @param orientation Forward uses the u derivative crossed with the v
#'   derivative. Reverse negates the normal and reverses the boundary.
#' @param curl_h Optional numerical step or three coordinate steps for curl.
#' @inheritParams vc_green
#' @return List with boundary circulation, curl flux, signed difference,
#'   comparison tolerance, both quadrature results and derivative diagnostics.
#' @details The boundary follows the four edges of the parameter rectangle.
#'   Forward traversal starts along increasing u at minimum v. A supplied
#'   surface Jacobian is also used for edge velocities. The caller supplies
#'   the smoothness, orientation and boundary hypotheses of Stokes' theorem.
#'   Multiple charts, trimming and topology inference are outside this API.
#'   Central stencils can leave the supplied domains. Agreement is diagnostic.
#' @examples
#' patch <- function(p)c(p[1],p[2],p[1]+p[2])
#' vc_stokes(function(p)c(-p[2],p[1],0),patch,c(0,1),c(0,1),
#'   n_boundary=4,n_outer=2,n_inner=2)$curl_flux
#' @export
vc_stokes <- function(field,surface,u_range,v_range,jacobian=NULL,curl=NULL,
    orientation=c("forward","reverse"),h=NULL,curl_h=NULL,
    boundary_method=c("simpson","trapezoid","adaptive"),n_boundary=16L,
    method=c("simpson","trapezoid"),n_outer=8L,n_inner=8L,abs_tol=1e-8,rel_tol=1e-8) {
  if (!is.function(field) || !is.function(surface)) stop("'field' and 'surface' must be functions.")
  if (!is.null(jacobian) && !is.function(jacobian)) stop("'jacobian' must be a function or NULL.")
  if (!is.null(curl) && !is.function(curl)) stop("'curl' must be a function or NULL.")
  u <- .vc_integral_interval(u_range)
  v <- .vc_integral_interval(v_range)
  if (u[1]>=u[2] || v[1]>=v[2]) stop("Parameter ranges must be strictly increasing.")
  orientation <- match.arg(orientation)
  method <- match.arg(method)
  boundary_method <- match.arg(boundary_method)
  h <- .vc_matrix_step(c(0,0),h)
  curl_h <- .vc_matrix_step(c(0,0,0),curl_h)
  reverse <- orientation=="reverse"
  checked_jacobian <- if (!is.null(jacobian)) function(p) {
    j <- jacobian(p)
    if (!is.numeric(j) || is.complex(j) || !identical(dim(j),c(3L,2L)) ||
        any(!is.finite(j))) stop("Jacobian must be a finite real 3 by 2 matrix.")
    j
  }
  edge <- function(path,range,column) {
    if (reverse) range <- rev(range)
    velocity <- if (!is.null(checked_jacobian)) function(t)checked_jacobian(path(t))[,column]
    step <- if (!is.null(h)) h[column]
    vc_line_integral(function(p).vc_spatial_field(field(p)),
      function(t).vc_spatial_field(surface(path(t))),range,"work",velocity=velocity,h=step,
      method=boundary_method,n=n_boundary,abs_tol=abs_tol,rel_tol=rel_tol)
  }
  boundary <- list(
    v_min=edge(function(t)c(t,v[1]),u,1L),
    u_max=edge(function(t)c(u[2],t),v,2L),
    v_max=edge(function(t)c(t,v[2]),rev(u),1L),
    u_min=edge(function(t)c(u[1],t),rev(v),2L))
  if (reverse) boundary <- rev(boundary)
  circulation <- .vc_finite(sum(vapply(boundary,function(r)r$estimate,numeric(1))),"Circulation")
  error <- .vc_finite(sum(vapply(boundary,function(r)r$integral$error_estimate,numeric(1))),
    "Boundary error indicator")
  curl_flags <- logical()
  rotation <- if (!is.null(curl)) curl else function(p) {
    r <- vc_curl(field,p,h=curl_h,warn=FALSE)
    curl_flags <<- c(curl_flags,all(r$reliable))
    r$estimate
  }
  surface_result <- vc_surface_integral(surface,u,v,rotation,"flux",checked_jacobian,
    orientation,h,method,n_outer,n_inner,abs_tol,rel_tol)
  flux <- surface_result$estimate
  difference <- .vc_finite(circulation-flux,"Stokes discrepancy")
  tolerance <- .vc_finite(max(abs_tol,rel_tol*max(abs(circulation),abs(flux))),
    "Comparison tolerance")
  list(type="stokes",u_range=u,v_range=v,orientation=orientation,circulation=circulation,
       curl_flux=flux,difference=difference,tolerance=tolerance,within_tolerance=abs(difference)<=tolerance,
       error_indicator=.vc_finite(error+surface_result$integral$error_estimate,"Combined error indicator"),
       boundary=boundary,surface=surface_result,
       derivative_sensitivity=list(boundary=vapply(boundary,function(r)r$derivative_sensitivity_ok,logical(1)),
         surface=surface_result$derivative_sensitivity_ok,curl=if(length(curl_flags))all(curl_flags) else NA,
         curl_evaluations=length(curl_flags)),
       note="Agreement is a numerical comparison under caller-supplied Stokes theorem hypotheses.")
}

#' Display stored Stokes circulation and curl flux
#' @param result Result of vc_stokes().
#' @param patch Result of vc_surface_sample() for the same surface and ranges.
#' @param title Optional character title.
#' @return Plotly widget with the patch, oriented boundary and both estimates.
#'   No field or surface is reevaluated.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   patch <- function(p)c(p[1],p[2],p[1]+p[2])
#'   result <- vc_stokes(function(p)c(-p[2],p[1],0),patch,c(0,1),c(0,1),
#'     n_boundary=4,n_outer=2,n_inner=2)
#'   vc_plot_stokes(result,vc_surface_sample(patch,c(0,1),c(0,1)))
#' }
#' @export
vc_plot_stokes <- function(result,patch,title=NULL) {
  if (!is.list(result) || !identical(result$type,"stokes")) stop("Invalid Stokes comparison.")
  if (!is.list(patch) || !identical(patch$type,"surface3d")) stop("Invalid surface patch.")
  if (!all(range(patch$u)==result$u_range) || !all(range(patch$v)==result$v_range))
    stop("Surface and Stokes parameter ranges must match.")
  p <- vc_plot_surface(patch,show_curve=FALSE,title=if(is.null(title))"Stokes comparison" else title)
  for (name in names(result$boundary)) {
    d <- result$boundary[[name]]$samples
    p <- plotly::add_trace(p,x=d$x,y=d$y,z=d$z,type="scatter3d",mode="lines",
      name=name,line=list(width=6,color="#ad442d"),inherit=FALSE)
    i <- max(1L,floor(.42*nrow(d)))
    j <- ceiling(.58*nrow(d))
    xyz <- as.matrix(d[c("x","y","z")])
    vector <- xyz[j,]-xyz[i,]
    if (all(vector==0)) next
    p <- plotly::add_trace(p,x=xyz[j,1],y=xyz[j,2],z=xyz[j,3],
      u=vector[1],v=vector[2],w=vector[3],type="cone",anchor="tip",
      sizemode="absolute",sizeref=.5,showscale=FALSE,showlegend=FALSE,
      colorscale=list(list(0,"#ad442d"),list(1,"#ad442d")),inherit=FALSE)
  }
  caption <- paste("Circulation =",format(result$circulation,digits=6),
    "| Curl flux =",format(result$curl_flux,digits=6),
    "| Difference =",format(result$difference,digits=3),"/",result$orientation)
  plotly::layout(p,annotations=list(list(text=caption,showarrow=FALSE,
    x=.5,y=0,xref="paper",yref="paper",yanchor="top")),margin=list(b=55))
}
