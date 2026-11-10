.vc_spatial_field <- function(value) {
  value <- .vc_point(value,"field output")
  if (length(value)!=3L) stop("Field output must have three Cartesian components.")
  unname(value)
}

#' Estimate the divergence of a square vector field
#' @inheritParams vc_jacobian
#' @param f Function of a Cartesian vector returning the same number of
#'   finite real components as the input.
#' @return List with estimate, refinements, reliable and the full Jacobian.
#' @details The estimate is the trace of the Jacobian. Reliability summarizes
#'   only the diagonal entries. Numerical sensitivity does not certify
#'   differentiability or bound the divergence error.
#' @examples
#' vc_divergence(function(p)c(p[1],2*p[2]),c(1,2))$estimate
#' @export
vc_divergence <- function(f,point,method="central",h=NULL,atol=1e-8,rtol=1e-4,warn=TRUE) {
  .vc_flag(warn,"warn")
  if (!is.function(f)) stop("'f' must be a function.")
  square <- function(p) {
    value <- .vc_vector_field_value(f,p)
    if (length(value)!=length(p)) stop("Field output length must match point dimension.")
    value
  }
  jacobian <- vc_jacobian(square,point,method,h,atol,rtol,warn=FALSE)
  reliable <- all(diag(jacobian$reliable))
  trace <- function(m).vc_finite(sum(diag(m)),"Divergence")
  if (warn && !reliable) warning("Some divergence derivative sensitivity checks failed.",call.=FALSE)
  list(type="divergence",estimate=trace(jacobian$estimate),
       refinements=vapply(jacobian$refinements,trace,numeric(1)),
       reliable=reliable,jacobian=jacobian)
}

#' Estimate the Cartesian curl of a spatial vector field
#' @inheritParams vc_jacobian
#' @param f Function of a Cartesian vector returning three finite real components.
#' @param point Finite real Cartesian vector of length three.
#' @return List with estimate, refinements, component reliability and full Jacobian.
#' @details Curl uses the right-handed Cartesian convention. Reliability
#'   summarizes the two partial derivatives used per component. To evaluate
#'   a time-dependent field at fixed time, capture that time in f.
#' @examples
#' vc_curl(function(p)c(-p[2],p[1],0),c(1,2,0))$estimate
#' @export
vc_curl <- function(f,point,method="central",h=NULL,atol=1e-8,rtol=1e-4,warn=TRUE) {
  .vc_flag(warn,"warn")
  point <- .vc_point(point)
  if (length(point)!=3L) stop("'point' must have three Cartesian coordinates.")
  if (!is.function(f)) stop("'f' must be a function.")
  jacobian <- vc_jacobian(function(p).vc_spatial_field(f(p)),point,method,h,atol,rtol,warn=FALSE)
  positive <- cbind(c(3,1,2),c(2,3,1))
  negative <- cbind(c(2,3,1),c(3,1,2))
  curl <- function(m) {
    value <- m[positive]-m[negative]
    if (any(!is.finite(value))) stop("Curl is not finite.")
    stats::setNames(value,c("x","y","z"))
  }
  reliable <- stats::setNames(jacobian$reliable[positive]&jacobian$reliable[negative],c("x","y","z"))
  if (warn && !all(reliable)) warning("Some curl derivative sensitivity checks failed.",call.=FALSE)
  list(type="curl",estimate=curl(jacobian$estimate),
       refinements=vapply(jacobian$refinements,curl,numeric(3)),
       reliable=reliable,jacobian=jacobian)
}

#' Sample a spatial vector field in a nested Cartesian solid
#' @inheritParams vc_solid_cartesian
#' @param f Function of one Cartesian vector returning three finite components.
#' @return List with type vector_field3d and samples x,y,z,u,v,w,magnitude.
#' @details Subdivisions include both endpoints and z varies fastest. The
#'   sampled solid validates boundary order at the nodes. Degenerate fibers
#'   can repeat positions.
#' @examples
#' vc_vector_field_sample(function(p)c(-p[2],p[1],0),function(x)-1,
#'   function(x)1,function(p)-1,function(p)1,c(-1,1),n=2,n_u=2,n_v=1)
#' @export
vc_vector_field_sample <- function(f,lower,upper,bottom,top,interval=NULL,
                                   n=8L,n_u=6L,n_v=6L,partition=NULL) {
  if (!is.function(f)) stop("'f' must be a function.")
  solid <- vc_solid_cartesian(lower,upper,bottom,top,interval,n,n_u,n_v,partition)
  grid <- solid$base$grid
  fractions <- solid$v
  lo <- as.vector(solid$bottom)
  hi <- as.vector(solid$top)
  if (as.double(length(fractions))*length(lo)>.Machine$integer.max) stop("Field grid too large.")
  z <- outer(fractions,hi-lo)
  z <- sweep(z,2,lo,"+")
  z[length(fractions),] <- hi
  points <- cbind(rep(as.vector(grid$x),each=length(fractions)),
                  rep(as.vector(grid$y),each=length(fractions)),as.vector(z))
  values <- t(vapply(seq_len(nrow(points)),function(i).vc_spatial_field(f(points[i,])),numeric(3)))
  magnitude <- apply(values,1,function(v).vc_finite(.vc_norm(v),"Field magnitude"))
  list(type="vector_field3d",samples=data.frame(x=points[,1],y=points[,2],z=points[,3],
    u=values[,1],v=values[,2],w=values[,3],magnitude=magnitude))
}

#' Trace a spatial field trajectory with fixed-step RK4
#' @param f Function of one Cartesian vector, or of a vector and time when
#'   time_dependent is TRUE. Returns three finite Cartesian components.
#' @param seed Initial Cartesian position at the first endpoint.
#' @param interval Two finite times in traversal order.
#' @param n Number of steps, at least one.
#' @param time_dependent Use the signature f(position,time) for a pathline.
#' @return List with type streamline3d, samples t,x,y,z,speed and field mode.
#' @details Autonomous fields give streamlines. Reversed times integrate
#'   backward. Equal endpoints give one sample. The method uses the supplied
#'   field as velocity and has no adaptive error estimate or boundary events.
#'   The visualization sampling domain does not clip the trajectory.
#' @examples
#' vc_streamline(function(p)c(-p[2],p[1],0),c(1,0,0),c(0,pi/2))$samples
#' @export
vc_streamline <- function(f,seed,interval,n=100L,time_dependent=FALSE) {
  if (!is.function(f)) stop("'f' must be a function.")
  seed <- .vc_point(seed,"seed")
  if (length(seed)!=3L) stop("'seed' must have three Cartesian coordinates.")
  interval <- .vc_integral_interval(interval)
  .vc_integer(n,"n",1)
  .vc_flag(time_dependent,"time_dependent")
  times <- if (interval[1]==interval[2]) interval[1] else
    seq(interval[1],interval[2],length.out=as.double(n)+1)
  steps <- diff(times)
  if (any(!is.finite(times)) || any(!is.finite(steps)) || any(steps==0))
    stop("Time grid is not representable.")
  field <- function(p,t) {
    if (any(!is.finite(p))) stop("Trajectory left the representable coordinate range.")
    .vc_spatial_field(if (time_dependent) f(p,t) else f(p))
  }
  positions <- matrix(NA_real_,length(times),3)
  speeds <- numeric(length(times))
  positions[1,] <- seed
  for (i in seq_along(times)) {
    p <- positions[i,]
    t <- times[i]
    k1 <- field(p,t)
    speeds[i] <- .vc_finite(.vc_norm(k1),"Trajectory speed")
    if (i==length(times)) break
    dt <- steps[i]
    mid <- t+dt/2
    if (mid==t || mid==times[i+1]) stop("RK4 midpoint is not representable.")
    k2 <- field(p+dt/2*k1,mid)
    k3 <- field(p+dt/2*k2,mid)
    k4 <- field(p+dt*k3,times[i+1])
    positions[i+1,] <- p+dt*(k1/6+k2/3+k3/3+k4/6)
  }
  list(type="streamline3d",time_dependent=time_dependent,
       samples=data.frame(t=times,x=positions[,1],y=positions[,2],z=positions[,3],speed=speeds))
}

#' Display stored field arrows and an optional trajectory
#' @param result Result of vc_vector_field_sample().
#' @param streamline Optional vc_streamline() result.
#' @param scale Positive maximum arrow length in Cartesian units. By default
#'   uses eight percent of the largest sampled axis span.
#' @param title Optional plot title.
#' @return Plotly widget using stored positions and vectors. Arrow lengths
#'   are proportional to magnitude and zero vectors are shown as points.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(p)c(-p[2],p[1],0)
#'   s <- vc_vector_field_sample(f,function(x)-1,function(x)1,
#'     function(p)0,function(p)0,c(-1,1),n=4,n_u=4,n_v=1)
#'   vc_plot_vector_field(s,vc_streamline(f,c(1,0,0),c(0,pi)))
#' }
#' @export
vc_plot_vector_field <- function(result,streamline=NULL,scale=NULL,title=NULL) {
  if (!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if (!is.list(result) || !identical(result$type,"vector_field3d")) stop("Invalid vector field samples.")
  if (!is.null(streamline) && (!is.list(streamline) || !identical(streamline$type,"streamline3d")))
    stop("Invalid streamline result.")
  if (!is.null(title) && (!is.character(title) || length(title)!=1L || is.na(title))) stop("Invalid title.")
  d <- result$samples
  p <- as.matrix(d[c("x","y","z")])
  if (is.null(scale)) {
    span <- max(apply(p,2,function(v)diff(range(v))))
    scale <- if (span==0) 1 else .08*span
  }
  .vc_positive(scale,"scale")
  fig <- plotly::plot_ly()
  moving <- d$magnitude>0
  if (any(moving)) {
    starts <- p[moving,,drop=FALSE]
    vectors <- as.matrix(d[moving,c("u","v","w")])
    magnitudes <- d$magnitude[moving]
    peak <- max(magnitudes)
    offsets <- scale*(vectors/peak)
    tips <- starts+offsets
    if (any(!is.finite(tips)) || any(rowSums(tips!=starts)==0))
      stop("Arrow endpoints are not representable; change 'scale'.")
    side <- t(vapply(seq_len(nrow(vectors)),function(i) {
      unit <- vectors[i,]/magnitudes[i]
      axis <- diag(3)[which.min(abs(unit)),]
      normal <- .vc_curve_cross(unit,axis)
      normal/.vc_norm(normal)
    },numeric(3)))
    wings <- side*(.09*(scale*(magnitudes/peak)))
    back <- tips-.23*offsets
    left <- back+wings
    right <- back-wings
    if (any(!is.finite(left)) || any(!is.finite(right)))
      stop("Arrowheads are not representable; change 'scale'.")
    lines <- lapply(seq_len(3),function(i)as.vector(rbind(
      starts[,i],tips[,i],NA_real_,left[,i],tips[,i],right[,i],NA_real_)))
    fig <- plotly::add_trace(fig,x=lines[[1]],y=lines[[2]],z=lines[[3]],
      type="scatter3d",mode="lines",name="Field arrows",connectgaps=FALSE,inherit=FALSE)
  }
  if (any(!moving)) fig <- plotly::add_trace(fig,x=d$x[!moving],y=d$y[!moving],
    z=d$z[!moving],type="scatter3d",mode="markers",name="Zero vectors",inherit=FALSE)
  if (!is.null(streamline)) {
    s <- streamline$samples
    fig <- plotly::add_trace(fig,x=s$x,y=s$y,z=s$z,type="scatter3d",mode="lines+markers",
      name=if(streamline$time_dependent)"Pathline" else "Streamline",inherit=FALSE,
      text=paste("t =",format(s$t),"/ speed =",format(s$speed)),hoverinfo="x+y+z+text")
  }
  plotly::layout(fig,title=if(is.null(title))"Spatial vector field" else title,
    scene=list(aspectmode="data",xaxis=list(title="x"),yaxis=list(title="y"),zaxis=list(title="z")))
}
