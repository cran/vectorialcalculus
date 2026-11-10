#' Sample a planar region independently of integration
#' @param lower,upper Scalar boundary functions of x, ordered at sampled nodes.
#' @param interval Finite increasing x bounds, or NULL with partition supplied.
#' @param n Positive number of regular x strips, ignored with partition.
#' @param n_u Positive number of subdivisions across each vertical section.
#' @param partition Optional vc_partition result; reconstructed from breaks.
#' @return Plain list with type region, partition, sections (outer,lower,upper),
#'   u blending fractions, grid (x,y matrices with u rows and x columns),
#'   strips (the scalar partition cell table) and note. No volume/area computed.
#' @details Reuses vc_partition unchanged. Validates sampled ordering without
#'   swapping boundaries. Linear interpolation between sampled boundary values
#'   is a representation, not proof of topology, ordering between samples or
#'   exact curved boundaries. A single interval between two graphs is supported;
#'   disconnected vertical cross sections require justified decomposition.
#'   Equal boundary heights are allowed as degenerate sections. Every stored
#'   coordinate must be finite; fine representable grids are not guaranteed.
#' @examples
#' vc_region(function(x)x^2,sqrt,c(0,1),n=10,n_u=4)
#' @export
vc_region <- function(lower,upper,interval=NULL,n=40L,n_u=20L,partition=NULL) {
  if(!is.function(lower)||!is.function(upper))stop("Both boundaries must be functions.")
  .vc_integer(n_u,"n_u",1)
  if(n_u>=.Machine$integer.max)stop("Blending resolution too large.")
  if(is.null(partition))p <- vc_partition(interval,n) else {
    if(!is.list(partition)||!identical(partition$type,"partition"))stop("Invalid partition.")
    p <- vc_partition(interval,breaks=partition$breaks)
  }
  if((as.double(n_u)+1)*length(p$breaks)>.Machine$integer.max)stop("Too many region grid nodes.")
  x <- p$breaks
  a <- vapply(x,function(t).vc_value(lower,t),numeric(1))
  b <- vapply(x,function(t).vc_value(upper,t),numeric(1))
  if(any(a>b))stop("Sampled lower boundary exceeds upper boundary.")
  width <- b-a
  if(any(!is.finite(width)))stop("Section width overflowed.")
  u <- seq(0,1,length.out=n_u+1)
  y <- vapply(seq_along(x),function(i)a[i]+u*width[i],numeric(length(u)))
  if(any(!is.finite(y)))stop("Region grid overflowed.")
  # Preserve actual sampled endpoints rather than reconstructed rounded sums.
  y[1,] <- a;y[length(u),] <- b
  list(type="region",partition=p,sections=data.frame(outer=x,lower=a,upper=b),u=u,
    grid=list(x=matrix(rep(x,each=length(u)),nrow=length(u)),y=y),strips=p$cells,
    note="Sampled region; interpolated graphs are not exact curved boundaries or topology certificates.")
}

#' Display a stored planar region in 2D or as a horizontal plate
#' @param result Result of vc_region().
#' @param view Filled 2D region or sampled 3D plate.
#' @param plane_z Finite display height for the 3D plate, not an integration value.
#' @param show_grid Logical; show stored vertical sections in the 3D view.
#' @param title Optional character title.
#' @return Plotly widget based only on stored region samples; no boundary
#'   reevaluation, integration or area calculation is performed.
#' @details In 2D reuses vc_plot_double_integral's sampled-domain renderer.
#'   In 3D the grid has u rows/x columns and every height equals plane_z.
#'   The plate has zero thickness; display height is not scalar-field height.
#'   Interpolation may miss crossings/oscillations between samples. Degenerate
#'   sections may collapse surface triangles. Human rendering review is separate.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   r <- vc_region(function(x)x^2,sqrt,c(0,1),n=8,n_u=4)
#'   vc_plot_region(r,view="3d",plane_z=2)
#' }
#' @export
vc_plot_region <- function(result,view=c("2d","3d"),plane_z=0,show_grid=TRUE,title=NULL) {
  if(!requireNamespace("plotly",quietly=TRUE))stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"region"))stop("Invalid region result.")
  view <- match.arg(view)
  .vc_finite(plane_z,"Plane height");.vc_flag(show_grid,"show_grid")
  if(!is.null(title)&&(!is.character(title)||length(title)!=1||is.na(title)))stop("Invalid title.")
  if(view=="2d")return(vc_plot_double_integral(list(type="double_integral",order="yx",
    sections=result$sections),title=if(is.null(title))"Sampled planar region" else title))
  x <- result$grid$x;y <- result$grid$y
  p <- plotly::plot_ly(x=x,y=y,z=matrix(plane_z,nrow(x),ncol(x)),type="surface",
    colorscale=list(list(0,"#188367"),list(1,"#188367")),showscale=FALSE,opacity=.55,
    name="Sampled plate")
  if(show_grid) {
    d <- result$sections
    a <- as.vector(rbind(d$outer,d$outer,NA_real_))
    b <- as.vector(rbind(d$lower,d$upper,NA_real_))
    z <- rep(plane_z,length(a));z[is.na(a)] <- NA_real_
    p <- plotly::add_trace(p,x=a,y=b,z=z,type="scatter3d",mode="lines",inherit=FALSE,
      name="Stored sections",line=list(color="#24476b"))
  }
  plotly::layout(p,title=if(is.null(title))paste("Sampled plate at z =",plane_z) else title,
    scene=list(xaxis=list(title="x"),yaxis=list(title="y"),zaxis=list(title="z"),aspectmode="data"))
}
