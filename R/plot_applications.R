.vc_application_style <- function(colors,required,alpha,line_width,title) {
  if(!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if(!is.character(colors)||!all(required %in% names(colors))||anyNA(colors)||
     anyDuplicated(names(colors))||any(!nzchar(colors))) stop("Invalid 'colors'.")
  if(!.vc_scalar(alpha)||alpha<0||alpha>1) stop("'alpha' must be between 0 and 1.")
  .vc_positive(line_width,"line_width")
  if(!is.null(title)&&(!is.character(title)||length(title)!=1||is.na(title))) stop("Invalid 'title'.")
}

#' Display geometric area between sampled graphs
#' @param result Result from vc_area().
#' @param f,g The same scalar boundary functions used in calculation; NULL g
#'   means zero. Functions are sampled, but quadrature is not recomputed.
#' @param n Integer curve sample count, at least two. Known stored breaks
#'   are also sampled to show supplied crossings.
#' @param colors Named f, g and region Plotly/CSS colors.
#' @param alpha Finite fill opacity between zero and one.
#' @param line_width Positive finite curve width.
#' @param line_dash Curve line style: solid, dash, dot or dashdot.
#' @param title Optional scalar title; default reports stored geometric and
#'   signed estimates and the area tolerance diagnostic.
#' @param xlab,ylab Scalar axis labels.
#' @param xlim,ylim Optional finite increasing display ranges; NULL uses data.
#' @return Plotly widget shading the sampled lower/upper envelopes, with both
#'   original graphs. Fill is a geometric region, not a new area calculation.
#' @details Crossings between samples can be missed. Invalid samples are gaps;
#'   neither quadrature nor rendering proves all features were resolved. Axes
#'   are linear to retain geometric meaning. Match f/g to the stored analysis.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   identity |> vc_area(c(-1,1),breaks=0) |> vc_plot_area(f=identity)
#' }
#' @export
vc_plot_area <- function(result,f,g=NULL,n=401L,
    colors=c(f="#24476b",g="#c44736",region="#188367"),alpha=.35,line_width=2,
    line_dash=c("solid","dash","dot","dashdot"),title=NULL,xlab="x",ylab="y",xlim=NULL,ylim=NULL) {
  .vc_application_style(colors,c("f","g","region"),alpha,line_width,title)
  if(!is.list(result)||!identical(result$type,"area")) stop("Invalid area result.")
  if(!is.function(f)) stop("'f' must be a function.")
  g <- .vc_boundary(g)
  .vc_integer(n,"n",2)
  line_dash <- match.arg(line_dash)
  label <- function(z)is.character(z)&&length(z)==1&&!is.na(z)
  if(!label(xlab)||!label(ylab)) stop("Invalid axis labels.")
  if(!is.null(xlim)) xlim <- .vc_interval(xlim)
  if(!is.null(ylim)) ylim <- .vc_interval(ylim)
  x <- sort(unique(c(seq(result$interval[1],result$interval[2],length.out=n),result$breaks)))
  a <- vapply(x,function(t).vc_safe_value(f,t),numeric(1))
  b <- vapply(x,function(t).vc_safe_value(g,t),numeric(1))
  if(!any(is.finite(a)&is.finite(b))) stop("No finite boundary samples.")
  p <- plotly::add_trace(plotly::plot_ly(),x=x,y=pmin(a,b),type="scatter",mode="lines",
    line=list(width=0),showlegend=FALSE,name="Sampled lower envelope",connectgaps=FALSE)
  p <- plotly::add_trace(p,x=x,y=pmax(a,b),type="scatter",mode="lines",inherit=FALSE,
    fill="tonexty",fillcolor=colors[["region"]],opacity=alpha,line=list(width=0),
    name="Geometric region (sampled)",connectgaps=FALSE)
  for(k in c("f","g")) p <- plotly::add_trace(p,x=x,y=if(k=="f")a else b,
    type="scatter",mode="lines",inherit=FALSE,name=k,connectgaps=FALSE,
    line=list(color=colors[[k]],width=line_width,dash=line_dash))
  if(is.null(title)) title <- paste0("Geometric area = ",format(result$estimate,digits=8),
    "; signed integral(f-g) = ",format(result$signed_estimate,digits=8),
    "; area diagnostic tolerance met: ",result$area$tolerance_met)
  plotly::layout(p,title=title,xaxis=list(title=xlab,range=xlim),yaxis=list(title=ylab,range=ylim,zeroline=TRUE))
}

.vc_revolution_xyz <- function(axial,radius,theta,axis,axis_at) {
  a <- matrix(rep(axial,each=length(theta)),nrow=length(theta))
  c <- axis_at+outer(cos(theta),radius)
  z <- outer(sin(theta),radius)
  list(x=if(axis=="horizontal")a else c,y=if(axis=="horizontal")c else a,z=z)
}

#' Display a solid of revolution and optional section playback
#' @param result Result from vc_revolution().
#' @inheritParams vc_plot_area
#' @param nt Integer angular sample count, at least four, including both 0
#'   and 2*pi to close circles.
#' @param colors Named outer, inner, axis, profile and section Plotly/CSS colors.
#' @param show_axis,show_profiles,show_caps Logical switches for rotation axis,
#'   generating boundary profiles and end faces. Caps affect graphics only.
#' @param curve_thetas Finite real vector of profile angles in radians.
#' @param curve_opacity Finite profile opacity between zero and one.
#' @param sections Optional finite vector of independent coordinates within
#'   interval for section rings/connectors. NULL draws none.
#' @param animate Logical; multiple sections get playback and a slider;
#'   FALSE displays the last section, one section is static.
#' @param frame_ms Positive finite frame duration in milliseconds.
#' @param scene Plotly scene list for axis labels, aspect ratio, camera and
#'   display scales. Default uses geometric data proportions and linear axes.
#' @param title Optional scalar title; default reports the stored volume,
#'   method and quadrature tolerance diagnostic.
#' @return Plotly 3D widget with washer outer/inner surfaces or shell upper/lower
#'   boundaries, optional end faces, generating profiles and rotation axis.
#'   Sections are rings/connectors of the stored procedure, not new integrals.
#' @details Coordinates reuse the same radii/height rules as the calculation.
#'   Sampled surfaces can miss features and are not a geometric proof. Curves
#'   are generating profiles; when a region straddles its axis, an interior
#'   profile need not lie on the outer boundary. Cap display does not compute
#'   cap area or change the reported volume/lateral area. Use the same f/g as
#'   the calculation; default g remains zero even for a shifted rotation axis.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   identity |> vc_revolution(c(0,1)) |> vc_plot_revolution(f=identity,sections=c(0,.5,1))
#' }
#' @export
vc_plot_revolution <- function(result,f,g=NULL,n=81L,nt=61L,
    colors=c(outer="#24476b",inner="#188367",axis="#333333",profile="#c44736",section="#dd8800"),
    alpha=.6,line_width=3,title=NULL,show_axis=TRUE,show_profiles=TRUE,show_caps=TRUE,
    curve_thetas=0,curve_opacity=1,sections=NULL,animate=TRUE,frame_ms=600,
    scene=list(aspectmode="data",xaxis=list(title="x"),yaxis=list(title="y"),zaxis=list(title="z"))) {
  .vc_application_style(colors,c("outer","inner","axis","profile","section"),alpha,line_width,title)
  if(!is.list(result)||!identical(result$type,"revolution")) stop("Invalid revolution result.")
  if(!is.function(f)) stop("'f' must be a function.")
  g <- .vc_boundary(g)
  .vc_integer(n,"n",2)
  .vc_integer(nt,"nt",4)
  .vc_positive(frame_ms,"frame_ms")
  logical_scalar <- function(z)is.logical(z)&&length(z)==1&&!is.na(z)
  if(!all(vapply(list(show_axis,show_profiles,show_caps,animate),logical_scalar,logical(1)))) stop("Invalid logical switch.")
  if(!.vc_scalar(curve_opacity)||curve_opacity<0||curve_opacity>1) stop("Invalid 'curve_opacity'.")
  if(!is.numeric(curve_thetas)||is.complex(curve_thetas)||any(!is.finite(curve_thetas))) stop("Invalid 'curve_thetas'.")
  if(!is.list(scene)) stop("'scene' must be a list.")
  if(!is.null(sections)&&(!is.numeric(sections)||is.complex(sections)||!is.null(dim(sections))||
      !length(sections)||any(!is.finite(sections))||any(sections<result$interval[1]|sections>result$interval[2])))
    stop("Invalid 'sections'.")
  t <- sort(unique(c(seq(result$interval[1],result$interval[2],length.out=n),result$breaks)))
  theta <- seq(0,2*pi,length.out=nt)
  sample <- function(t).vc_revolution_section(f,g,t,result$axis_at,result$method)
  d <- lapply(t,sample)
  get <- function(name)vapply(d,function(z)z[[name]],numeric(1))
  xyz <- function(a,r,th=theta).vc_revolution_xyz(a,r,th,result$axis,result$axis_at)
  p <- plotly::plot_ly()
  surface <- function(a,r,color,name) {
    z <- xyz(a,r)
    p <<- plotly::add_surface(p,x=z$x,y=z$y,z=z$z,inherit=FALSE,showscale=FALSE,
      colorscale=list(list(0,color),list(1,color)),opacity=alpha,name=name,showlegend=TRUE)
  }
  if(result$method=="washers") {
    surface(t,get("outer"),colors[["outer"]],"Outer radial boundary")
    if(any(get("inner")>0)) surface(t,get("inner"),colors[["inner"]],"Inner radial boundary")
  } else {
    surface(get("upper"),get("outer"),colors[["outer"]],"Upper axial boundary")
    surface(get("lower"),get("outer"),colors[["inner"]],"Lower axial boundary")
  }
  if(show_caps) for(i in c(1,length(t))) {
    if(result$method=="washers") surface(rep(t[i],2),c(d[[i]]$inner,d[[i]]$outer),colors[["inner"]],"End annulus") else
      surface(c(d[[i]]$lower,d[[i]]$upper),rep(d[[i]]$outer,2),colors[["inner"]],"Radial end cylinder")
  }
  if(show_axis) {
    axial <- if(result$method=="washers")range(t) else range(c(get("lower"),get("upper")))
    z <- xyz(axial,c(0,0),0)
    p <- plotly::add_trace(p,x=as.numeric(z$x),y=as.numeric(z$y),z=as.numeric(z$z),
      type="scatter3d",mode="lines",inherit=FALSE,name="Rotation axis",line=list(color=colors[["axis"]],width=line_width))
  }
  if(show_profiles) for(th in curve_thetas) for(fun in list(f,g)) {
    v <- vapply(t,function(u).vc_value(fun,u),numeric(1))
    z <- if(result$method=="washers") xyz(t,v-result$axis_at,th) else xyz(v,t-result$axis_at,th)
    p <- plotly::add_trace(p,x=as.numeric(z$x),y=as.numeric(z$y),z=as.numeric(z$z),
      type="scatter3d",mode="lines",inherit=FALSE,name=paste("Generating profile at",format(th),"rad"),
      opacity=curve_opacity,line=list(color=colors[["profile"]],width=line_width))
  }
  if(length(sections)) {
    if(!animate) sections <- sections[length(sections)]
    data <- do.call(rbind,lapply(seq_along(sections),function(i) {
      u <- sections[i]
      q <- sample(u)
      axial <- if(result$method=="washers")rep(u,2) else c(q$lower,q$upper)
      radial <- if(result$method=="washers")c(q$outer,q$inner) else rep(q$outer,2)
      rings <- xyz(axial,radial)
      connector <- xyz(axial,radial,0)
      flatten <- function(r,c)c(r[,1],NA_real_,r[,2],NA_real_,as.numeric(c))
      data.frame(x=flatten(rings$x,connector$x),y=flatten(rings$y,connector$y),
        z=flatten(rings$z,connector$z),frame=i,
        text=paste0("t=",u,"; ",result$method,"; volume density=",q$density))
    }))
    if(length(sections)>1) {
      p <- plotly::add_trace(p,data=data,x=~x,y=~y,z=~z,frame=~frame,text=~text,
        type="scatter3d",mode="lines",inherit=FALSE,name="Selected section",hoverinfo="text",
        line=list(color=colors[["section"]],width=line_width))
      p <- plotly::animation_opts(p,frame=frame_ms,transition=0,redraw=TRUE)
      p <- plotly::animation_button(p,x=0,y=-.1)
      p <- plotly::animation_slider(p,x=.18,len=.8,y=-.1,currentvalue=list(prefix="Section index: "))
    } else p <- plotly::add_trace(p,x=data$x,y=data$y,z=data$z,text=data$text,
      type="scatter3d",mode="lines",inherit=FALSE,name="Selected section",hoverinfo="text",
      line=list(color=colors[["section"]],width=line_width))
  }
  if(is.null(title)) title <- paste0("Revolution (",result$method,"): V = ",format(result$estimate,digits=8),
    "; diagnostic tolerance met: ",result$integral$tolerance_met)
  p <- plotly::layout(p,title=title,scene=scene,margin=list(b=100))
  plotly::config(p,displaylogo=FALSE,scrollZoom=TRUE)
}
