.vc_osculating_section <- function(frame) {
  bending <- .vc_curve_bending(frame)
  radius <- NA_real_
  center <- rep(NA_real_,3)
  status <- bending$status
  if (status=="ok") {
    radius <- 1/bending$curvature
    candidate <- frame$point+radius*frame$normal
    if (!is.finite(radius) || any(!is.finite(candidate)) || all(candidate==frame$point)) {
      status <- "unrepresentable"
      if (!is.finite(radius)) radius <- NA_real_
    } else center <- candidate
  }
  list(frame=frame,curvature=bending$curvature,radius=radius,center=center,status=status)
}

#' Sample osculating circles of a spatial curve
#' @inheritParams vc_frenet_frame
#' @param t Nonempty finite real vector of parameter values in the desired order.
#' @param n Number of circle segments, at least three.
#' @return List with type osculating_circle, t and circles in parameter order.
#'   Each circle stores the frame, curvature, radius, center, status and a
#'   three-column ring matrix when the circle is available.
#' @details A circle at a resolved point has radius 1/curvature and center
#'   point plus radius times the principal normal. The sampled ring starts
#'   and ends at the curve point and initially follows the tangent. A ring
#'   is unavailable for unresolved curvature or speed, or when its geometry
#'   cannot be represented. The frame diagnostics are retained. The
#'   construction needs only the first two derivatives.
#' @examples
#' curve <- function(t)c(cos(t),sin(t),t/2)
#' vc_osculating_circle(curve,c(0,pi/2))$circles[[1]]$radius
#' @export
vc_osculating_circle <- function(curve,t,velocity=NULL,acceleration=NULL,
                                  h=NULL,h2=NULL,n=72L) {
  t <- .vc_point(t,"t")
  .vc_integer(n,"n",3)
  theta <- seq(0,2*pi,length.out=as.double(n)+1)
  circles <- lapply(t,function(s) {
    frame <- vc_frenet_frame(curve,s,velocity,acceleration,h,h2)
    section <- .vc_osculating_section(frame)
    section["ring"] <- list(NULL)
    if (section$status=="ok") {
      directions <- outer(-cos(theta),frame$normal)+outer(sin(theta),frame$tangent)
      ring <- sweep(section$radius*directions,2,section$center,"+")
      ring[1,] <- ring[nrow(ring),] <- frame$point
      next_ring <- ring[-1,,drop=FALSE]
      previous_ring <- ring[-nrow(ring),,drop=FALSE]
      if (any(!is.finite(ring)) || any(rowSums(next_ring!=previous_ring)==0)) {
        section$status <- "unrepresentable"
      } else section$ring <- ring
    }
    section
  })
  list(type="osculating_circle",t=t,circles=circles)
}

#' Sample the osculating ribbon of a spatial curve
#' @inheritParams vc_frenet_frame
#' @inheritParams vc_curve_sample
#' @param n_u Number of samples across each normal segment, at least two.
#' @param u_max Positive fraction of the curvature radius at most one.
#' @return List with type osculating_ribbon, stored curve, normal fractions u,
#'   sections and grid matrices x, y, z with rows in u and columns in parameter
#'   order. Unavailable sections are entire columns of NA.
#' @details At each available parameter the segment is the curve point plus
#'   u times the radius times the principal normal. At u=1 it reaches the
#'   center of curvature. The stored curve retains its samples at unavailable
#'   sections. A reversed interval reverses columns and an equal interval
#'   samples one section. Surface interpolation is a sampled visualization.
#' @examples
#' curve <- function(t)c(cos(t),sin(t),t/2)
#' vc_osculating_ribbon(curve,c(0,pi),n=9)$sections[[1]]$center
#' @export
vc_osculating_ribbon <- function(curve,interval,velocity=NULL,acceleration=NULL,
                                  h=NULL,h2=NULL,n=61L,n_u=9L,u_max=1) {
  .vc_integer(n_u,"n_u",2)
  .vc_positive(u_max,"u_max")
  if (u_max>1) stop("'u_max' must be at most one.")
  sampled <- vc_curve_sample(curve,interval,n)
  u <- seq(0,u_max,length.out=n_u)
  sections <- lapply(sampled$samples$t,function(t)
    .vc_osculating_section(vc_frenet_frame(curve,t,velocity,acceleration,h,h2)))
  fibers <- vector("list",length(sections))
  for (i in seq_along(sections)) {
    section <- sections[[i]]
    fibers[[i]] <- matrix(NA_real_,n_u,3)
    if (section$status=="ok") {
      p <- section$frame$point
      fiber <- sweep(outer(u,section$radius*section$frame$normal),2,p,"+")
      if (all(is.finite(fiber)) && any(fiber[n_u,]!=p)) {
        fibers[[i]] <- fiber
      } else sections[[i]]$status <- "unrepresentable"
    }
  }
  grid <- lapply(seq_len(3),function(axis)
    vapply(fibers,function(fiber)fiber[,axis],numeric(n_u)))
  names(grid) <- c("x","y","z")
  list(type="osculating_ribbon",curve=sampled,u=u,sections=sections,grid=grid)
}

#' Display stored osculating circles
#' @inheritParams vc_plot_curve
#' @param result Result of vc_osculating_circle().
#' @param curve_samples Optional result of vc_curve_sample() for the same curve
#'   covering the circle parameters.
#' @param fill Display either rings or triangulated disks.
#' @param show_radius Logical flag to draw a segment from each center to its
#'   curve point.
#' @return Plotly widget with stored anchors and available circles.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   curve <- function(t)c(cos(t),sin(t),t/2)
#'   vc_plot_osculating_circle(vc_osculating_circle(curve,pi/2),
#'     vc_curve_sample(curve,c(0,pi)),show_radius=TRUE)
#' }
#' @export
vc_plot_osculating_circle <- function(result,curve_samples=NULL,fill=c("ring","disk"),
                                       show_radius=FALSE,title=NULL) {
  if (!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if (!is.list(result) || !identical(result$type,"osculating_circle")) stop("Invalid osculating circles.")
  fill <- match.arg(fill)
  .vc_flag(show_radius,"show_radius")
  if (is.null(title)) title <- "Osculating circles"
  if (is.null(curve_samples)) {
    if (!is.character(title) || length(title)!=1L || is.na(title)) stop("Invalid title.")
    p <- plotly::layout(plotly::plot_ly(),title=title,scene=list(aspectmode="data"))
  } else {
    p <- vc_plot_curve(curve_samples,title=title)
    if (any(result$t<min(curve_samples$interval) | result$t>max(curve_samples$interval)))
      stop("Circle parameters are outside the sampled interval.")
  }
  anchors <- vapply(result$circles,function(circle)circle$frame$point,numeric(3))
  labels <- vapply(result$circles,function(circle)
    paste("t =",format(circle$frame$t),"/",circle$status),character(1))
  p <- plotly::add_trace(p,x=anchors[1,],y=anchors[2,],z=anchors[3,],
    type="scatter3d",mode="markers",name="Circle points",text=labels,
    hoverinfo="x+y+z+text",inherit=FALSE)
  for (circle in result$circles) {
    ring <- circle$ring
    if (is.null(ring)) next
    name <- paste("t =",format(circle$frame$t))
    if (fill=="ring") {
      p <- plotly::add_trace(p,x=ring[,1],y=ring[,2],z=ring[,3],
        type="scatter3d",mode="lines",name=name,inherit=FALSE)
    } else {
      vertices <- rbind(circle$center,ring)
      edge <- seq_len(nrow(ring)-1L)
      p <- plotly::add_trace(p,x=vertices[,1],y=vertices[,2],z=vertices[,3],
        i=rep(0L,length(edge)),j=edge,k=edge+1L,type="mesh3d",
        name=name,opacity=.5,color=I("#cf553e"),inherit=FALSE)
    }
    if (show_radius) {
      segment <- rbind(circle$center,circle$frame$point)
      p <- plotly::add_trace(p,x=segment[,1],y=segment[,2],z=segment[,3],
        type="scatter3d",mode="lines",name="Radius",showlegend=FALSE,
        line=list(color="#555555",dash="dash"),inherit=FALSE)
    }
  }
  p
}

#' Display a stored osculating ribbon
#' @inheritParams vc_plot_curve
#' @param result Result of vc_osculating_ribbon().
#' @param opacity Ribbon opacity between zero and one.
#' @param show_centers Logical flag to display the sampled centers of curvature.
#' @return Plotly widget with the stored curve and available ribbon sections.
#'   Surface gaps remain open and isolated sections are shown as line segments.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   curve <- function(t)c(cos(t),sin(t),t/2)
#'   vc_plot_osculating_ribbon(vc_osculating_ribbon(curve,c(0,pi),n=19))
#' }
#' @export
vc_plot_osculating_ribbon <- function(result,opacity=.4,show_centers=TRUE,title=NULL) {
  if (!is.list(result) || !identical(result$type,"osculating_ribbon")) stop("Invalid osculating ribbon.")
  if (!.vc_scalar(opacity) || opacity<0 || opacity>1) stop("Invalid opacity.")
  .vc_flag(show_centers,"show_centers")
  p <- vc_plot_curve(result$curve,title=if(is.null(title))"Osculating ribbon" else title)
  grid <- result$grid
  available <- is.finite(grid$x[1,])
  left <- c(FALSE,available[-length(available)])
  right <- c(available[-1],FALSE)
  if (any(available & (left | right))) {
    p <- plotly::add_trace(p,x=grid$x,y=grid$y,z=grid$z,type="surface",
      name="Ribbon",opacity=opacity,showscale=FALSE,connectgaps=FALSE,inherit=FALSE)
  }
  isolated <- which(available & !left & !right)
  if (length(isolated)) {
    segments <- lapply(grid,function(axis)as.vector(rbind(axis[,isolated,drop=FALSE],NA_real_)))
    p <- plotly::add_trace(p,x=segments$x,y=segments$y,z=segments$z,
      type="scatter3d",mode="lines",connectgaps=FALSE,name="Isolated sections",inherit=FALSE)
  }
  if (show_centers && any(available)) {
    centers <- vapply(result$sections,function(section)
      if (section$status=="ok") section$center else rep(NA_real_,3),numeric(3))
    p <- plotly::add_trace(p,x=centers[1,],y=centers[2,],z=centers[3,],
      type="scatter3d",mode="lines+markers",connectgaps=FALSE,name="Curvature centers",
      marker=list(size=3),inherit=FALSE)
  }
  p
}
