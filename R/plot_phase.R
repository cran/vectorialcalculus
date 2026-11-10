#' Plot critical-point candidates and derivative evidence
#' @param result Result from vc_critical_points(), vc_concavity() or
#'   vc_newton(), matching the selected plot function.
#' @inheritParams vc_plot_secant
#' @param xlim Finite increasing plot bounds. NULL uses the stored search
#'   interval, or the accepted Newton iterates with padding.
#' @param colors Named vector with curve, line, point and diagnostic colors.
#' @return Interactive plotly htmlwidget. Analysis plots have function and
#'   derivative-evidence panels; Newton adds actual iteration steps and a
#'   residual history panel. Calculations are not rerun.
#' @details Markers are conditional numerical candidates, not proven extrema.
#'   Sampled curves may miss singularities; invalid samples are gaps. Existing
#'   scan values and classifications are used without silently recomputing them.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(x) x^3-3*x
#'   f |> vc_critical_points(c(-2,2)) |> vc_plot_critical_points(f=f)
#' }
#' @export
vc_plot_critical_points <- function(result,f,xlim=NULL,n=401L,
                                    colors=c(curve="#24476b",line="#c44736",point="#188367",diagnostic="#8250a3"),
                                    title=NULL,xlab="x",ylab="f(x)") {
  .vc_plot_analysis(result,f,"critical_points",xlim,n,colors,title,xlab,ylab)
}

#' Plot sampled concavity and smooth inflection candidates
#' @inheritParams vc_plot_critical_points
#' @return Interactive plotly widget with sampled curvature labels in hover
#'   and inflection candidate markers; no proof of inflection is implied.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(x) x^3
#'   f |> vc_concavity(c(-1,1)) |> vc_plot_concavity(f=f)
#' }
#' @export
vc_plot_concavity <- function(result,f,xlim=NULL,n=401L,
                              colors=c(curve="#24476b",line="#c44736",point="#188367",diagnostic="#8250a3"),
                              title=NULL,xlab="x",ylab="f(x)") {
  .vc_plot_analysis(result,f,"concavity",xlim,n,colors,title,xlab,ylab)
}

#' Plot Newton iterations, tangents and numerical residuals
#' @inheritParams vc_plot_critical_points
#' @param animate Logical; animate computed tangent updates when at least
#'   two computed outgoing steps exist. Initial-root, single-step and
#'   failed-before-step results remain static; no next iterate is fabricated.
#' @param frame_ms Positive finite frame duration in milliseconds.
#' @return Interactive plotly widget with iteration geometry and residual
#'   history. Play/pause and a frame slider are available for multiple steps.
#' @details Tangent segments join the current function point to its actual
#'   computed Newton x-intercept. Final points without a computed outgoing
#'   step have no tangent-update frame. Hover and title include residual/status;
#'   numerical convergence does not guarantee a mathematical root.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(x) x^2-2
#'   f |> vc_newton(1) |> vc_plot_newton(f=f)
#' }
#' @export
vc_plot_newton <- function(result,f,xlim=NULL,n=401L,animate=TRUE,frame_ms=500,
                           colors=c(curve="#24476b",line="#c44736",point="#188367",diagnostic="#8250a3"),
                           title=NULL,xlab="x",ylab="f(x)") {
  if (!is.logical(animate) || length(animate)!=1L || is.na(animate)) stop("'animate' must be TRUE or FALSE.")
  .vc_positive(frame_ms,"frame_ms")
  .vc_plot_analysis(result,f,"newton",xlim,n,colors,title,xlab,ylab,animate,frame_ms)
}

.vc_phase_plot_validate <- function(result,f,type,n,colors,title,xlab,ylab) {
  if (!requireNamespace("plotly",quietly=TRUE)) stop("Package 'plotly' is required only for visualization.")
  if (!is.list(result) || !identical(result$type,type)) stop("'result' has the wrong analysis type.")
  if (!is.function(f)) stop("'f' must be a function.")
  .vc_integer(n,"n",2)
  if (!is.character(colors) || !all(c("curve","line","point","diagnostic") %in% names(colors)) ||
      anyNA(colors) || anyDuplicated(names(colors)) || any(!nzchar(colors))) stop("Invalid named 'colors'.")
  label <- function(x) is.character(x) && length(x)==1L && !is.na(x)
  if ((!is.null(title) && !label(title)) || !label(xlab) || !label(ylab)) stop("Invalid plot labels.")
  if (type=="newton") {
    xs <- result$iterations$x
    ys <- result$iterations$fx
  } else {
    xs <- result$scan$x
    ys <- result$scan$value
  }
  finite <- is.finite(ys)
  actual <- vapply(xs[finite],function(x) .vc_safe_value(f,x),numeric(1))
  if (!length(actual) || !isTRUE(all.equal(as.double(actual),as.double(ys[finite]),tolerance=0)))
    stop("'f' does not match stored values.")
}

.vc_plot_analysis <- function(result,f,type,xlim,n,colors,title,xlab,ylab,
                              animate=FALSE,frame_ms=500) {
  .vc_phase_plot_validate(result,f,type,n,colors,title,xlab,ylab)
  if (is.null(xlim)) {
    if (type!="newton") xlim <- result$interval else {
      bounds <- range(result$iterations$x)
      pad <- max(diff(bounds)/4,.5)
      xlim <- bounds+c(-pad,pad)
    }
  }
  xlim <- .vc_interval(xlim)
  xs <- seq(xlim[1],xlim[2],length.out=n)
  ys <- vapply(xs,function(x) .vc_safe_value(f,x),numeric(1))
  if (!any(is.finite(ys))) stop("No finite function samples.")
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p,x=xs,y=ys,type="scatter",mode="lines",name="f(x) (sampled)",
                         connectgaps=FALSE,line=list(color=colors[["curve"]]))
  e <- plotly::plot_ly()
  if (type=="newton") {
    d <- result$iterations
    p <- plotly::add_trace(p,x=d$x,y=d$fx,text=paste0("n=",d$n,"; residual=",format(d$residual,digits=6)),
      type="scatter",mode="markers",name="Accepted iterates",marker=list(color=colors[["point"]],size=8))
    p <- plotly::add_trace(p,x=xlim,y=c(0,0),type="scatter",mode="lines",name="Zero axis",
                           line=list(color=colors[["diagnostic"]],dash="dot"))
    e <- plotly::add_trace(e,x=d$n,y=d$residual,type="scatter",mode="lines+markers",name="Numerical residual",
                           line=list(color=colors[["diagnostic"]]))
  } else {
    points <- if (type=="concavity") rbind(result$inflections,result$singular_candidates) else result$points
    if (nrow(points)) {
      labels <- if (type=="concavity") points$status else paste(points$kind,points$classification,points$evidence)
      marker_colors <- if(type=="concavity") rep(colors[["point"]],nrow(points)) else
        ifelse(points$classification=="maximum",colors[["line"]],
          ifelse(points$classification=="minimum",colors[["point"]],colors[["diagnostic"]]))
      symbols <- if(type=="concavity") ifelse(grepl("unresolved",points$status),"diamond-open","circle") else
        ifelse(points$kind=="nonsmooth_candidate","diamond","circle")
      p <- plotly::add_trace(p,x=points$x,y=points$value,text=labels,type="scatter",mode="markers",
        name="Conditional candidates",marker=list(color=marker_colors,size=10,symbol=symbols))
    }
    values <- if (type=="concavity") result$scan$second_derivative else result$scan$derivative
    labels <- if (type=="concavity") result$scan$classification else rep("Central quotient; not a differentiability test",length(values))
    e <- plotly::add_trace(e,x=result$scan$x,y=values,text=labels,type="scatter",mode="lines+markers",
                           name=if(type=="concavity") "Second differences" else "First differences",
                           connectgaps=FALSE,line=list(color=colors[["diagnostic"]]),marker=list(size=3))
    if(type=="concavity") for(label in c("concave_up","concave_down")) {
      use <- result$scan$classification==label & is.finite(values)
      if(any(use)) e <- plotly::add_trace(e,x=result$scan$x,y=ifelse(use,values,NA_real_),
        type="scatter",mode="lines+markers",connectgaps=FALSE,
        name=if(label=="concave_up") "Concave up (sampled)" else "Concave down (sampled)",
        line=list(color=colors[[if(label=="concave_up") "point" else "line"]],width=3),marker=list(size=5))
    }
  }
  p <- plotly::subplot(p,e,nrows=2,shareX=type!="newton",titleY=TRUE,heights=c(.65,.35),margin=.08)
  if (type=="newton" && nrow(result$iterations)>1L) {
    d <- result$iterations
    steps <- do.call(rbind,lapply(seq_len(nrow(d)-1L),function(i)
      data.frame(x=c(d$x[i],d$x[i+1L]),y=c(d$fx[i],0),frame=i,
                 text=paste0("n=",d$n[i],"; slope=",format(d$dfx[i],digits=6)))))
    if (animate && nrow(d)>2L) {
      p <- plotly::add_trace(p,data=steps,x=~x,y=~y,frame=~frame,text=~text,
        xaxis="x",yaxis="y",type="scatter",mode="lines+markers",name="Computed tangent update",
        line=list(color=colors[["line"]]),marker=list(size=7))
      p <- plotly::animation_opts(p,frame=frame_ms,transition=0,redraw=TRUE)
      p <- plotly::animation_button(p,x=0,y=-.16)
      p <- plotly::animation_slider(p,x=.18,len=.8,y=-.16,currentvalue=list(prefix="Computed step: "))
    } else {
      segments <- do.call(rbind,lapply(split(steps,steps$frame),function(d)
        data.frame(x=c(d$x,NA_real_),y=c(d$y,NA_real_))))
      p <- plotly::add_trace(p,x=segments$x,y=segments$y,xaxis="x",yaxis="y",type="scatter",mode="lines",
        name="Computed tangent updates",connectgaps=FALSE,line=list(color=colors[["line"]]))
    }
  }
  if (is.null(title)) title <- switch(type,newton=paste("Newton-Raphson:",result$status),
    critical_points="Critical-point candidates",concavity="Sampled concavity and inflection candidates")
  p <- plotly::layout(p,title=list(text=paste0(title,"<br><sup>Numerical evidence, not a mathematical guarantee</sup>")),
    xaxis=list(title=xlab,range=xlim),yaxis=list(title=ylab),
    yaxis2=list(title=if(type=="newton") "Residual" else "Derivative evidence"),
    paper_bgcolor="white",plot_bgcolor="white",margin=list(t=100,b=130),
    legend=list(orientation="h",y=1.04,yanchor="bottom"))
  if (type=="newton") p <- plotly::layout(p,xaxis2=list(title="Iteration"))
  p$height <- 800
  plotly::config(p,displaylogo=FALSE,scrollZoom=TRUE)
}
