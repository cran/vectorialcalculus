#' Plot a candidate tangent plane, coordinate tangents and linearization error
#' @param result Two-dimensional vc_linearization() result.
#' @inheritParams vc_plot_gradient
#' @param xlim,ylim Optional finite increasing bounds containing the point;
#'   NULL covers the point and successfully evaluated probe endpoints.
#' @param n Integer number of surface/plane samples per coordinate, at least two.
#' @param alpha Finite original-surface opacity between zero and one.
#' @param scene Plotly 3D scene list for camera, axes and aspect settings.
#' @param colors Named surface, plane, x_tangent, y_tangent, normal, increment
#'   and point Plotly/CSS colors.
#' @param plane_alpha Finite plane opacity between zero and one.
#' @param t_len Positive finite half-length of coordinate tangent segments.
#' @param show_surface,show_plane,show_tangents,show_sections,show_normal Logical
#'   display switches. Sections are sampled true coordinate curves; tangents
#'   are linear candidates. Normal uses the stored unit vector (-fx,-fy,1).
#' @param normal_scale Positive finite display length for the normal vector.
#' @param animate Logical; multiple valid probe rows get playback/slider;
#'   FALSE displays the last valid row. One valid row is static.
#' @param title Optional scalar title; default labels a candidate affine plane
#'   and the separate coordinate-sensitivity flag. It never certifies a tangent.
#' @return Plotly widget with sampled scalar surface, candidate affine plane,
#'   coordinate tangent models and optional normal. Each displayed probe joins
#'   its predicted and actual heights above the same actual domain endpoint.
#' @details Calculation is not rerun. Surface and affine values are sampled
#'   separately with rows corresponding to y and columns to x. Invalid samples
#'   are gaps. Probe rows with unavailable values remain in result diagnostics
#'   but are omitted from geometry and counted in the default title. Hover
#'   retains scale, differential, finite change, remainder and normalized
#'   remainder. The plane is a mathematical tangent only if appropriate total
#'   differentiability holds; finite probes and partial checks do not prove it.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(v)v[1]^2+2*v[2]^2
#'   f |> vc_linearization(c(1,1),dx=c(.1,.2)) |> vc_plot_linearization(f=f)
#' }
#' @export
vc_plot_linearization <- function(result,f,xlim=NULL,ylim=NULL,n=61L,
    colors=c(surface="#24476b",plane="#188367",x_tangent="#c44736",y_tangent="#8250a3",
             normal="#555555",increment="#dd8800",point="#222222"),
    alpha=.5,plane_alpha=.35,line_width=3,t_len=.5,show_surface=TRUE,show_plane=TRUE,
    show_tangents=TRUE,show_sections=TRUE,show_normal=FALSE,normal_scale=1,animate=TRUE,frame_ms=500,
    title=NULL,scene=list(aspectmode="data")) {
  .vc_application_style(colors,c("surface","plane","x_tangent","y_tangent","normal","increment","point"),
                         alpha,line_width,title)
  point <- .vc_multivariable_plot(result,f,"linearization")
  if(length(point)!=2) stop("The plane view requires a two-dimensional domain.")
  .vc_integer(n,"n",2)
  .vc_positive(t_len,"t_len")
  .vc_positive(normal_scale,"normal_scale")
  .vc_positive(frame_ms,"frame_ms")
  flag <- function(z)is.logical(z)&&length(z)==1&&!is.na(z)
  if(!all(vapply(list(show_surface,show_plane,show_tangents,show_sections,show_normal,animate),flag,logical(1))))
    stop("Invalid logical display switch.")
  if(!.vc_scalar(plane_alpha)||plane_alpha<0||plane_alpha>1) stop("Invalid 'plane_alpha'.")
  if(!is.list(scene)) stop("Invalid 'scene'.")
  d <- result$comparison
  valid <- if(is.null(d))integer() else which(d$status%in%c("ok","zero_increment") &
    is.finite(d$approximation)&is.finite(d$value))
  bounds <- function(i) {
    v <- c(point[i]+c(-1,1),vapply(result$probes[valid],function(z)z$endpoint[i],numeric(1)))
    r <- range(v)
    r+c(-.1,.1)*diff(r)
  }
  if(is.null(xlim)) xlim <- bounds(1)
  if(is.null(ylim)) ylim <- bounds(2)
  xlim <- .vc_interval(xlim)
  ylim <- .vc_interval(ylim)
  if(point[1]<xlim[1]||point[1]>xlim[2]||point[2]<ylim[1]||point[2]>ylim[2])
    stop("Display bounds must contain the base point.")
  x <- seq(xlim[1],xlim[2],length.out=n)
  y <- seq(ylim[1],ylim[2],length.out=n)
  sample <- function(fun)vapply(x,function(a)vapply(y,function(b) {
    v <- point;v[] <- c(a,b);.vc_field_safe(fun,v)
  },numeric(1)),numeric(n))
  p <- plotly::plot_ly()
  if(show_surface) {
    z <- sample(f)
    if(!any(is.finite(z))) stop("No finite surface samples.")
    p <- plotly::add_surface(p,x=x,y=y,z=z,inherit=FALSE,opacity=alpha,showscale=FALSE,
      connectgaps=FALSE,name="Sampled scalar surface",
      colorscale=list(list(0,"#ffffff"),list(1,colors[["surface"]])))
  }
  if(show_plane) {
    z <- sample(result$model)
    if(!any(is.finite(z))) stop("No finite affine plane samples; reduce the display range.")
    p <- plotly::add_surface(p,x=x,y=y,z=z,inherit=FALSE,opacity=plane_alpha,showscale=FALSE,
      connectgaps=FALSE,name="Candidate affine plane",
      colorscale=list(list(0,colors[["plane"]]),list(1,colors[["plane"]])))
  }
  slopes <- as.double(result$gradient$estimate)
  if(length(slopes)!=2||any(!is.finite(slopes))) stop("Invalid stored gradient.")
  if(show_sections) for(i in 1:2) {
    a <- if(i==1)x else rep(point[1],n)
    b <- if(i==2)y else rep(point[2],n)
    z <- vapply(seq_len(n),function(k) {
      v <- point;v[] <- c(a[k],b[k]);.vc_field_safe(f,v)
    },numeric(1))
    p <- plotly::add_trace(p,x=a,y=b,z=z,type="scatter3d",mode="lines",inherit=FALSE,
      connectgaps=FALSE,name=if(i==1)"Actual x section" else "Actual y section",
      line=list(color=colors[[if(i==1)"x_tangent" else "y_tangent"]],width=line_width,dash="dot"))
  }
  if(show_tangents) for(i in 1:2) {
    t <- c(-t_len,t_len)
    a <- rep(point[1],2);b <- rep(point[2],2)
    if(i==1) a <- a+t else b <- b+t
    z <- result$value+slopes[i]*t
    if(any(!is.finite(c(a,b,z)))) stop("Coordinate tangent display overflowed.")
    p <- plotly::add_trace(p,x=a,y=b,z=z,type="scatter3d",mode="lines",inherit=FALSE,
      name=if(i==1)"Candidate x tangent" else "Candidate y tangent",
      line=list(color=colors[[if(i==1)"x_tangent" else "y_tangent"]],width=line_width))
  }
  if(show_normal) {
    base <- c(point,result$value)
    if(is.null(result$unit_normal)) stop("Stored unit normal is unavailable; rescale coordinates.")
    end <- base+normal_scale*result$unit_normal
    if(any(!is.finite(end))) stop("Normal display overflowed.")
    p <- plotly::add_trace(p,x=c(base[1],end[1]),y=c(base[2],end[2]),z=c(base[3],end[3]),
      type="scatter3d",mode="lines+markers",inherit=FALSE,name="Candidate unit normal",
      line=list(color=colors[["normal"]],width=line_width),marker=list(size=4))
  }
  p <- plotly::add_trace(p,x=point[1],y=point[2],z=result$value,type="scatter3d",mode="markers",
    inherit=FALSE,name="Stored base point",marker=list(color=colors[["point"]],size=5))
  if(length(valid)) {
    selected <- if(animate)valid else valid[length(valid)]
    data <- do.call(rbind,lapply(seq_along(selected),function(k) {
      i <- selected[k]
      endpoint <- result$probes[[i]]$endpoint
      data.frame(x=rep(endpoint[1],2),y=rep(endpoint[2],2),z=c(d$approximation[i],d$value[i]),frame=k,
        text=paste0(c("Affine model","Actual function"),"; probe=",d$probe[i],"; scale=",d$scale[i],
          "; differential=",format(d$differential[i],digits=7),"; actual change=",format(d$change[i],digits=7),
          "; remainder=",format(d$remainder[i],digits=7),"; normalized=",format(d$normalized_remainder[i],digits=7)))
    }))
    if(length(selected)>1) {
      p <- plotly::add_trace(p,data=data,x=~x,y=~y,z=~z,frame=~frame,text=~text,
        type="scatter3d",mode="lines+markers",inherit=FALSE,name="Linear versus actual increment",
        hoverinfo="x+y+z+text",line=list(color=colors[["increment"]],width=line_width),marker=list(size=6))
      p <- plotly::animation_opts(p,frame=frame_ms,transition=0,redraw=TRUE)
      p <- plotly::animation_button(p,x=0,y=-.15)
      p <- plotly::animation_slider(p,x=.18,len=.8,y=-.15,currentvalue=list(prefix="Probe index: "))
    } else p <- plotly::add_trace(p,x=data$x,y=data$y,z=data$z,text=data$text,
      type="scatter3d",mode="lines+markers",inherit=FALSE,name="Linear versus actual increment",
      hoverinfo="x+y+z+text",line=list(color=colors[["increment"]],width=line_width),marker=list(size=6))
  }
  omitted <- if(is.null(d))0 else nrow(d)-length(valid)
  if(is.null(title)) title <- paste("Candidate affine plane; coordinate sensitivity:",result$gradient_reliable,
    "; omitted invalid probes:",omitted)
  plotly::config(plotly::layout(p,title=title,scene=scene,margin=list(t=90,b=120)),displaylogo=FALSE,scrollZoom=TRUE)
}
