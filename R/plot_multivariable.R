.vc_multivariable_plot <- function(result,f,type) {
  if(!is.list(result)||!identical(result$type,type)) stop("Invalid multivariable result.")
  point <- .vc_point(result$point)
  if(!length(point)%in%c(2L,3L)) stop("This view requires dimension two or three.")
  if(!identical(.vc_field_value(f,point),result$value)) stop("'f' does not match the stored point value.")
  point
}

#' Plot a gradient with optional stored directional playback
#' @param result Result from vc_gradient(), in dimension two or three.
#' @inheritParams vc_partial
#' @param directions Optional list of vc_directional results at the same point.
#'   Playback uses those stored unit directions/estimates; none are computed.
#' @param xlim,ylim Optional finite increasing contour bounds covering the
#'   point. NULL covers the displayed gradient endpoint and unit directions.
#'   Used only in dimension two.
#' @param n Integer number of contour samples per coordinate, at least two.
#' @param vector_scale Positive finite display multiplier on the gradient vector.
#'   Hover retains the unscaled estimate. Unit directions have unit display length.
#' @param colors Named field, gradient, direction and point Plotly/CSS colors.
#' @param alpha Finite contour opacity between zero and one.
#' @param line_width Positive finite vector line width.
#' @param animate Logical; multiple stored directions get playback and slider.
#'   FALSE shows the last direction; a single direction remains static.
#' @param frame_ms Positive finite frame duration in milliseconds.
#' @param title Optional scalar title; default reports numerical sensitivity.
#' @param scene Plotly scene list for the three-dimensional vector view.
#' @return Plotly widget: sampled contours plus gradient in dimension two,
#'   or gradient/point in three-dimensional domain coordinates. In 3D the
#'   scalar field itself is not drawn as an asserted four-dimensional surface.
#' @details Recovers the original gradient_scalar vector and gradient_direction2d
#'   directional playback representations. Vectors are drawn in the domain;
#'   endpoints/hover identify direction and magnitude. A zero numerical gradient
#'   is not assigned an arbitrary preferred direction. No maximum/minimum or
#'   differentiability assertion follows from the arrows. Invalid contour
#'   evaluations are gaps, and axis aspect ratios preserve data proportions.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(v)v[1]^2+2*v[2]^2
#'   f |> vc_gradient(c(1,1)) |> vc_plot_gradient(f=f)
#' }
#' @export
vc_plot_gradient <- function(result,f,directions=NULL,xlim=NULL,ylim=NULL,n=61L,vector_scale=1,
    colors=c(field="#24476b",gradient="#c44736",direction="#188367",point="#222222"),
    alpha=.6,line_width=3,animate=TRUE,frame_ms=400,title=NULL,scene=list(aspectmode="data")) {
  .vc_application_style(colors,c("field","gradient","direction","point"),alpha,line_width,title)
  point <- .vc_multivariable_plot(result,f,"gradient")
  .vc_integer(n,"n",2)
  .vc_positive(vector_scale,"vector_scale")
  .vc_positive(frame_ms,"frame_ms")
  if(!is.logical(animate)||length(animate)!=1||is.na(animate)) stop("Invalid 'animate'.")
  if(!is.list(scene)) stop("Invalid 'scene'.")
  end <- point+vector_scale*result$estimate
  if(length(end)!=length(point)||any(!is.finite(end))) stop("Gradient display endpoint is invalid.")
  p <- plotly::plot_ly()
  if(length(point)==2) {
    bounds <- function(i) {
      r <- range(c(point[i]+c(-1,1),end[i]))
      r+c(-.2,.2)*diff(r)
    }
    if(is.null(xlim)) xlim <- bounds(1)
    if(is.null(ylim)) ylim <- bounds(2)
    xlim <- .vc_interval(xlim)
    ylim <- .vc_interval(ylim)
    if(point[1]<xlim[1]||point[1]>xlim[2]||point[2]<ylim[1]||point[2]>ylim[2]) stop("Contour bounds must contain the point.")
    x <- seq(xlim[1],xlim[2],length.out=n)
    y <- seq(ylim[1],ylim[2],length.out=n)
    z <- vapply(x,function(a)vapply(y,function(b) {
      v <- point;v[] <- c(a,b);.vc_field_safe(f,v)
    },numeric(1)),numeric(n))
    if(!any(is.finite(z))) stop("No finite contour samples.")
    p <- plotly::add_trace(p,x=x,y=y,z=z,type="contour",name="Scalar field (sampled)",
      opacity=alpha,showscale=FALSE,colorscale=list(list(0,"#ffffff"),list(1,colors[["field"]])))
  }
  add_vector <- function(a,b,name,color,text) {
    if(length(point)==2) p <<- plotly::add_trace(p,x=c(a[1],b[1]),y=c(a[2],b[2]),
      type="scatter",mode="lines+markers",inherit=FALSE,name=name,text=text,
      line=list(color=color,width=line_width),marker=list(color=color,size=c(4,9))) else
      p <<- plotly::add_trace(p,x=c(a[1],b[1]),y=c(a[2],b[2]),z=c(a[3],b[3]),
      type="scatter3d",mode="lines+markers",inherit=FALSE,name=name,text=text,
      line=list(color=color,width=line_width),marker=list(color=color,size=c(3,6)))
  }
  add_vector(point,end,"Numerical gradient",colors[["gradient"]],
    paste("gradient:",paste(format(result$estimate,digits=7),collapse=", "),"; reliable:",result$reliable))
  if(length(directions)) {
    if(!is.list(directions)) stop("'directions' must be a list of stored results.")
    if(!animate) directions <- directions[length(directions)]
    data <- do.call(rbind,lapply(seq_along(directions),function(i) {
      d <- directions[[i]]
      .vc_multivariable_plot(d,f,"directional")
      if(!identical(as.double(d$point),as.double(point))) stop("Stored directions must use the same point.")
      v <- point+d$direction
      if(any(!is.finite(v))) stop("Direction display overflowed.")
      data.frame(x=c(point[1],v[1]),y=c(point[2],v[2]),
        z=if(length(point)==3)c(point[3],v[3]) else c(0,0),frame=i,
        text=paste0("unit direction; derivative=",format(d$estimate,digits=7),
                    "; sensitivity=",d$reliable,"; gradient comparison=",d$gradient_comparison))
    }))
    args <- list(p=p,data=data,x=~x,y=~y,type=if(length(point)==2)"scatter" else "scatter3d",
      mode="lines+markers",inherit=FALSE,name="Stored unit direction",text=~text,
      line=list(color=colors[["direction"]],width=line_width),marker=list(color=colors[["direction"]],size=6))
    if(length(point)==3) args$z <- ~z
    if(length(directions)>1) args$frame <- ~frame
    p <- do.call(plotly::add_trace,args)
    if(length(directions)>1) {
      p <- plotly::animation_opts(p,frame=frame_ms,transition=0,redraw=TRUE)
      p <- plotly::animation_button(p,x=0,y=-.15)
      p <- plotly::animation_slider(p,x=.18,len=.8,y=-.15,currentvalue=list(prefix="Direction index: "))
    }
  }
  if(is.null(title)) title <- paste("Numerical gradient; coordinate sensitivity:",result$reliable)
  if(length(point)==2) p <- plotly::layout(p,title=title,xaxis=list(title="x",range=xlim),
    yaxis=list(title="y",range=ylim,scaleanchor="x",scaleratio=1),margin=list(b=120)) else
    p <- plotly::layout(p,title=title,scene=scene,margin=list(b=120))
  plotly::config(p,displaylogo=FALSE,scrollZoom=TRUE)
}

#' Plot a directional line restriction on a sampled scalar surface
#' @param result Two-dimensional vc_directional() result.
#' @inheritParams vc_plot_gradient
#' @param t_range Finite increasing line-parameter bounds containing zero.
#' @param xlim,ylim Optional finite increasing surface bounds covering the
#'   point. NULL uses point plus/minus two; t_range controls the separate curve.
#' @param colors Named surface, curve, tangent and point Plotly/CSS colors.
#' @return Plotly 3D surface with the actual line-restriction curve,
#'   numerical tangent through the stored point, and stored base point.
#' @details Recovers the directional curve/surface idea of legacy
#'   directional_derivative3d. Curve height is sampled f(point+t*u), whereas
#'   the tangent height uses only the stored derivative estimate. The direction
#'   is a unit domain vector, so the parameter is domain distance. Failed
#'   diagnostics and one-sided methods are explicitly labeled. Missing samples
#'   are gaps, not proof that no intervening singularity exists.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(v)v[1]^2+2*v[2]^2
#'   f |> vc_directional(c(1,1),c(1,1)) |> vc_plot_directional(f=f)
#' }
#' @export
vc_plot_directional <- function(result,f,xlim=NULL,ylim=NULL,t_range=c(-1,1),n=61L,
    colors=c(surface="#24476b",curve="#188367",tangent="#c44736",point="#222222"),
    alpha=.45,line_width=3,title=NULL,scene=list(aspectmode="data")) {
  .vc_application_style(colors,c("surface","curve","tangent","point"),alpha,line_width,title)
  point <- .vc_multivariable_plot(result,f,"directional")
  if(length(point)!=2) stop("The scalar-surface view requires a two-dimensional domain.")
  .vc_integer(n,"n",2)
  if(!is.list(scene)) stop("Invalid 'scene'.")
  t_range <- .vc_interval(t_range)
  if(t_range[1]>0||t_range[2]<0) stop("'t_range' must contain zero.")
  if(is.null(xlim)) xlim <- point[1]+c(-2,2)
  if(is.null(ylim)) ylim <- point[2]+c(-2,2)
  xlim <- .vc_interval(xlim)
  ylim <- .vc_interval(ylim)
  if(point[1]<xlim[1]||point[1]>xlim[2]||point[2]<ylim[1]||point[2]>ylim[2]) stop("Surface bounds must contain the point.")
  x <- seq(xlim[1],xlim[2],length.out=n)
  y <- seq(ylim[1],ylim[2],length.out=n)
  z <- vapply(x,function(a)vapply(y,function(b) {
    v <- point;v[] <- c(a,b);.vc_field_safe(f,v)
  },numeric(1)),numeric(n))
  if(!any(is.finite(z))) stop("No finite surface samples.")
  t <- sort(unique(c(0,seq(t_range[1],t_range[2],length.out=n))))
  path <- lapply(t,function(s)point+s*result$direction)
  if(any(!is.finite(unlist(path)))) stop("Plot path overflowed.")
  height <- vapply(path,function(p).vc_field_safe(f,p),numeric(1))
  tangent <- result$value+result$estimate*t
  if(any(!is.finite(tangent))) stop("Tangent display overflowed.")
  p <- plotly::add_surface(plotly::plot_ly(),x=x,y=y,z=z,opacity=alpha,showscale=FALSE,
    colorscale=list(list(0,"#ffffff"),list(1,colors[["surface"]])),name="Sampled scalar surface")
  for(k in c("curve","tangent")) p <- plotly::add_trace(p,
    x=vapply(path,function(v)v[1],numeric(1)),y=vapply(path,function(v)v[2],numeric(1)),
    z=if(k=="curve")height else tangent,type="scatter3d",mode="lines",inherit=FALSE,
    connectgaps=FALSE,name=if(k=="curve")"Actual line restriction (sampled)" else "Numerical directional tangent",
    line=list(color=colors[[k]],width=line_width))
  p <- plotly::add_trace(p,x=point[1],y=point[2],z=result$value,type="scatter3d",mode="markers",
    inherit=FALSE,name="Stored base point",marker=list(color=colors[["point"]],size=5))
  if(is.null(title)) title <- paste("Directional derivative:",format(result$estimate,digits=7),
    ";",result$method,"; sensitivity:",result$reliable)
  plotly::config(plotly::layout(p,title=title,scene=scene),displaylogo=FALSE,scrollZoom=TRUE)
}
