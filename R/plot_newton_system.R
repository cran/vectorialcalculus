#' Animate two-dimensional Newton system iterates and zero contours
#' @param result Two-dimensional vc_newton_system result, including failed runs.
#' @param f The same vector-valued function used in the calculation.
#' @param xlim,ylim Optional finite increasing display bounds. NULL covers
#'   accepted iterates plus a margin; rejected trials do not expand the view.
#' @param n Integer samples per coordinate axis, at least two. Zero contours
#'   are sampled independently; missed/tangential roots are possible.
#' @param colors Named first, second, path, current and rejected colors.
#' @param line_width Positive finite curve/path width.
#' @param marker_size Positive finite iterate marker size.
#' @param alpha Opacity of optional rejected-trial markers, between zero and one.
#' @param show_rejected Logical; display finite rejected trial coordinates,
#'   labeled as rejected and not connected to the accepted iterate path.
#' @param animate Logical; TRUE plays stored accepted iterates, FALSE shows
#'   the complete path. A run with zero updates has a static initial marker.
#' @param frame_ms Positive animation frame duration in milliseconds.
#' @param title Optional scalar title; default reports stored termination status.
#' @param xlab,ylab Scalar coordinate-axis labels.
#' @return Plotly widget with two zero level sets and accepted Newton path.
#'   Each frame includes the current iterate, residual and accepted damping.
#'   Contour matrices have y rows and x columns. The terminal frame is included
#'   even after failure; no updates or Jacobians are recomputed for plotting.
#' @details Invalid samples are gaps. f is checked against stored system
#'   values at every accepted point. Contour intersections illustrate possible
#'   roots, not proofs of existence or uniqueness. Rejected-trial markers,
#'   when requested, are a static overlay of all rejections and do not imply
#'   they occurred at each playback time. Only 2D systems are drawn; the engine
#'   supports higher dimensions, whose Jacobians can use vc_plot_matrix.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(v)c(v[1]^2+v[2]^2-1,v[1]-v[2])
#'   f |> vc_newton_system(c(.8,.2)) |> vc_plot_newton_system(f=f)
#' }
#' @export
vc_plot_newton_system <- function(result,f,xlim=NULL,ylim=NULL,n=81L,
    colors=c(first="#188367",second="#8250a3",path="#24476b",current="#c44736",rejected="#dd8800"),
    line_width=2,marker_size=9,alpha=.45,show_rejected=FALSE,animate=TRUE,
    frame_ms=600,title=NULL,xlab="x",ylab="y") {
  .vc_application_style(colors,c("first","second","path","current","rejected"),alpha,line_width,title)
  if(!is.list(result)||!identical(result$type,"newton_system")||length(result$x0)!=2)
    stop("A two-dimensional vc_newton_system result is required.")
  .vc_integer(n,"n",2);.vc_positive(marker_size,"marker_size");.vc_positive(frame_ms,"frame_ms")
  .vc_flag(show_rejected,"show_rejected");.vc_flag(animate,"animate")
  if(!is.character(xlab)||length(xlab)!=1||is.na(xlab)||!is.character(ylab)||length(ylab)!=1||is.na(ylab))
    stop("Invalid axis labels.")
  named <- function(p)stats::setNames(p,names(result$x0))
  for(i in seq_len(nrow(result$points))) {
    z <- .vc_vector_field_value(f,named(result$points[i,]))
    if(length(z)!=2||!identical(names(z),colnames(result$values))||
       any(abs(z-result$values[i,])>1e-8*pmax(1,abs(z),abs(result$values[i,]))))
      stop("'f' does not match the stored system values.")
  }
  bounds <- function(i) {
    r <- range(result$points[,i]);span <- max(1,diff(r))
    r+c(-.2,.2)*span
  }
  if(is.null(xlim))xlim <- bounds(1)
  if(is.null(ylim))ylim <- bounds(2)
  xlim <- .vc_interval(xlim);ylim <- .vc_interval(ylim)
  x <- seq(xlim[1],xlim[2],length.out=n);y <- seq(ylim[1],ylim[2],length.out=n)
  samples <- lapply(x,function(a)lapply(y,function(b)tryCatch({
    z <- .vc_vector_field_value(f,named(c(a,b)))
    if(length(z)!=2||!identical(names(z),colnames(result$values)))stop("Invalid system shape.")
    z
  },error=function(e)c(NA_real_,NA_real_))))
  p <- plotly::plot_ly()
  for(i in 1:2) {
    z <- vapply(samples,function(col)vapply(col,function(v)v[i],numeric(1)),numeric(n))
    if(!any(is.finite(z)))stop("No finite contour samples.")
    color <- colors[[if(i==1)"first" else "second"]]
    p <- plotly::add_trace(p,x=x,y=y,z=z,type="contour",inherit=FALSE,
      name=paste0("F",i," = 0"),showscale=FALSE,connectgaps=FALSE,
      autocontour=FALSE,contours=list(start=0,end=0,size=1,coloring="lines",showlabels=TRUE),
      colorscale=list(list(0,color),list(1,color)),line=list(width=line_width),
      hovertemplate=paste0("x=%{x}<br>y=%{y}<br>F",i,"=%{z}<extra></extra>"))
  }
  m <- nrow(result$points)
  labels <- paste0("Iteration ",result$iterations$n,"<br>residual=",format(result$iterations$residual),
    "<br>damping=",format(result$iterations$damping),"<br>Status: ",result$status)
  if(animate&&m>1) {
    d <- do.call(rbind,lapply(seq_len(m),function(i)data.frame(
      x=result$points[seq_len(i),1],y=result$points[seq_len(i),2],text=labels[seq_len(i)],frame=i-1L)))
    p <- plotly::add_trace(p,data=d,x=~x,y=~y,text=~text,frame=~frame,
      type="scatter",mode="lines+markers",name="Accepted Newton iterates",inherit=FALSE,
      line=list(color=colors[["path"]],width=line_width),marker=list(size=marker_size),hoverinfo="text")
    current <- data.frame(x=result$points[,1],y=result$points[,2],text=labels,frame=seq_len(m)-1L)
    p <- plotly::add_trace(p,data=current,x=~x,y=~y,text=~text,frame=~frame,
      type="scatter",mode="markers",name="Current iterate",inherit=FALSE,
      marker=list(size=marker_size+3,color=colors[["current"]]),hoverinfo="text")
    p <- plotly::animation_opts(p,frame=frame_ms,transition=0,redraw=TRUE)
    p <- plotly::animation_slider(p,currentvalue=list(prefix="Newton iteration: "))
  } else p <- plotly::add_trace(p,x=result$points[,1],y=result$points[,2],text=labels,
    type="scatter",mode="lines+markers",inherit=FALSE,name="Accepted Newton iterates",
    line=list(color=colors[["path"]],width=line_width),marker=list(size=marker_size,color=colors[["current"]]),hoverinfo="text")
  rejected <- Filter(function(t)!t$accepted&&all(is.finite(t$point)),result$trials)
  if(show_rejected&&length(rejected)) {
    points <- do.call(rbind,lapply(rejected,function(t)t$point))
    p <- plotly::add_trace(p,x=points[,1],y=points[,2],type="scatter",mode="markers",inherit=FALSE,
      name="All rejected trials (static overlay)",opacity=alpha,
      marker=list(size=marker_size,symbol="x",color=colors[["rejected"]]),
      text=vapply(rejected,function(t)paste0("Rejected trial, update ",t$iteration,
        "<br>damping=",t$damping,"<br>residual=",t$residual,"<br>",t$error),character(1)),hoverinfo="text")
  }
  if(is.null(title))title <- paste("Newton system:",result$status,"- finite stopping evidence only")
  plotly::layout(p,title=title,xaxis=list(title=xlab,range=xlim),yaxis=list(title=ylab,range=ylim))
}
