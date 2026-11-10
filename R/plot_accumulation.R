#' Explore accumulation and signed regions interactively
#' @param result Result from vc_accumulation().
#' @param f The same scalar function used for accumulation. The sampled
#'   function is drawn without recomputing any integrals.
#' @param fundamental Optional vc_fundamental() result with exactly the same
#'   lower limit and x vector. Its derivative estimates are overlaid on f;
#'   hover and open markers retain inconclusive/mismatch diagnostics.
#' @param n Integer number of curve/fill samples, at least two per frame.
#' @param colors Named curve, positive, negative, accumulation and derivative
#'   Plotly/CSS colors. Positive/negative identify heights above/below zero,
#'   not the sign of a reversed integral.
#' @param alpha Finite fill opacity from zero to one.
#' @param line_width Positive finite curve width.
#' @param line_dash One of solid, dash, dot or dashdot for sampled curves.
#' @param animate Logical. Multiple stored upper limits get playback and a
#'   slider; FALSE draws the last stored limit. A single limit is static.
#' @param frame_ms Positive finite animation duration in milliseconds.
#' @param xlim Optional finite increasing display bounds covering lower and
#'   every stored x. NULL selects their range with padding for a single point.
#' @param title Optional scalar title.
#' @param xlab,ylab,alab Scalar labels for x, f(x) and A(x).
#' @return Two-panel Plotly widget: sampled f with zero-based positive/negative
#'   fills between lower and selected x, and sampled accumulation with a
#'   synchronized zero-to-A segment. Orientation, integral estimate, estimated error and
#'   status appear in hover. All stored accumulation nodes remain visible.
#' @details Only stored integrals and comparisons are displayed. Function
#'   sampling and linear joins between accumulation nodes are illustrative,
#'   not new quadrature. Fill is not an unsigned-area computation; for x below
#'   lower the integral reverses sign even though heights remain unchanged.
#'   Invalid f samples are gaps. Theorems, continuity and integrability are
#'   not visually certified. Linear axes preserve signed-area interpretation.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(x)x^2
#'   f |> vc_accumulation(seq(0,2,length.out=9)) |> vc_plot_accumulation(f=f)
#' }
#' @export
vc_plot_accumulation <- function(result,f,fundamental=NULL,n=101L,
    colors=c(curve="#24476b",positive="#188367",negative="#c44736",
             accumulation="#8250a3",derivative="#dd8800"),
    alpha=.35,line_width=2,line_dash=c("solid","dash","dot","dashdot"),
    animate=TRUE,frame_ms=600,xlim=NULL,title=NULL,xlab="x",ylab="f(x)",alab="A(x)") {
  if(!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"accumulation")) stop("Invalid accumulation result.")
  if(!is.function(f)) stop("'f' must be a function.")
  .vc_integer(n,"n",2)
  .vc_positive(line_width,"line_width")
  .vc_positive(frame_ms,"frame_ms")
  line_dash <- match.arg(line_dash)
  if(!.vc_scalar(alpha)||alpha<0||alpha>1) stop("'alpha' must be between 0 and 1.")
  if(!is.logical(animate)||length(animate)!=1||is.na(animate)) stop("Invalid 'animate'.")
  if(!is.character(colors)||!all(c("curve","positive","negative","accumulation","derivative") %in% names(colors))||
     anyNA(colors)||anyDuplicated(names(colors))||any(!nzchar(colors))) stop("Invalid 'colors'.")
  label <- function(z)is.character(z)&&length(z)==1&&!is.na(z)
  if((!is.null(title)&&!label(title))||!all(vapply(list(xlab,ylab,alab),label,logical(1)))) stop("Invalid labels.")
  d <- result$table
  xs <- .vc_accumulation_x(d$x)
  .vc_finite(result$lower,"lower")
  if(length(d$estimate)!=length(xs)||any(!is.finite(d$estimate))) stop("Invalid stored estimates.")
  bounds <- range(c(result$lower,xs))
  if(is.null(xlim)) {
    xlim <- bounds
    if(xlim[1]==xlim[2]) xlim <- xlim+c(-1,1)*max(1,abs(xlim[1]))*.05
  }
  xlim <- .vc_interval(xlim)
  if(xlim[1]>bounds[1]||xlim[2]<bounds[2]) stop("'xlim' must cover the stored limits.")
  grid <- seq(xlim[1],xlim[2],length.out=n)
  values <- vapply(grid,function(z).vc_safe_value(f,z),numeric(1))
  if(!any(is.finite(values))) stop("No finite function samples.")
  p <- plotly::add_trace(plotly::plot_ly(),x=grid,y=values,type="scatter",mode="lines",
    name="f(x) (sampled)",connectgaps=FALSE,
    line=list(color=colors[["curve"]],width=line_width,dash=line_dash))
  if(!is.null(fundamental)) {
    if(!is.list(fundamental)||!identical(fundamental$type,"fundamental")||
       !identical(fundamental$lower,result$lower)||!identical(fundamental$table$x,d$x))
      stop("'fundamental' must use the same lower limit and x vector.")
    q <- fundamental$table
    actual <- vapply(q$x,function(z).vc_value(f,z),numeric(1))
    if(!isTRUE(all.equal(actual,q$value,tolerance=0))) stop("'f' does not match the comparison values.")
    p <- plotly::add_trace(p,x=q$x,y=q$estimate,type="scatter",mode="markers",inherit=FALSE,
      name="Numerical A'(x)",text=paste(q$status,"; difference=",q$difference,
        "; quadrature indicator=",q$quadrature_indicator),hoverinfo="x+y+text",
      marker=list(color=colors[["derivative"]],size=8,
        symbol=ifelse(q$status=="consistent","circle","diamond-open")))
  }
  text <- paste0("lower=",result$lower,"; x=",d$x,"; signed A(x)=",format(d$estimate,digits=8),
    "; estimated error=",format(d$error_estimate,digits=4),"; ",d$status,
    "; diagnostic tolerance met=",d$tolerance_met)
  profile <- data.frame(x=c(result$lower,d$x),y=c(0,d$estimate))
  profile <- profile[order(profile$x),]
  profile <- profile[!duplicated(profile$x),]
  a <- plotly::add_trace(plotly::plot_ly(),x=profile$x,y=profile$y,type="scatter",mode="lines",
    name="A(x) (sampled; linear joins)",line=list(color=colors[["accumulation"]],width=line_width,dash=line_dash))
  a <- plotly::add_trace(a,x=d$x,y=d$estimate,type="scatter",mode="markers",inherit=FALSE,name="Stored accumulation",
    text=text,hoverinfo="x+y+text",marker=list(color=colors[["accumulation"]],size=5))
  p <- plotly::subplot(p,a,nrows=2,shareX=TRUE,titleY=TRUE,heights=c(.6,.4),margin=.08)
  selected <- if(animate)seq_along(xs) else length(xs)
  animated <- length(selected)>1
  fills <- do.call(rbind,lapply(selected,function(i) {
    limits <- sort(c(result$lower,xs[i]))
    xx <- seq(limits[1],limits[2],length.out=n)
    yy <- vapply(xx,function(z).vc_safe_value(f,z),numeric(1))
    data.frame(x=xx,positive=pmax(yy,0),negative=pmin(yy,0),frame=i,text=text[i])
  }))
  for(s in c("positive","negative")) {
    z <- data.frame(x=fills$x,y=fills[[s]],frame=fills$frame,text=fills$text)
    if(animated) p <- plotly::add_trace(p,data=z,x=~x,y=~y,frame=~frame,text=~text,
      xaxis="x",yaxis="y",type="scatter",mode="lines",inherit=FALSE,fill="tozeroy",name=paste(s,"heights"),
      connectgaps=FALSE,opacity=alpha,line=list(color=colors[[s]],width=1),hoverinfo="x+y+text") else
      p <- plotly::add_trace(p,x=z$x,y=z$y,text=z$text,xaxis="x",yaxis="y",type="scatter",mode="lines",
        fill="tozeroy",name=paste(s,"heights"),inherit=FALSE,connectgaps=FALSE,opacity=alpha,
        line=list(color=colors[[s]],width=1),hoverinfo="x+y+text")
  }
  marker <- do.call(rbind,lapply(selected,function(i)
    data.frame(x=rep(xs[i],2),y=c(0,d$estimate[i]),frame=i,text=text[i])))
  if(animated) {
    p <- plotly::add_trace(p,data=marker,x=~x,y=~y,frame=~frame,text=~text,xaxis="x",yaxis="y2",
      type="scatter",mode="lines+markers",inherit=FALSE,name="Selected upper limit",hoverinfo="x+y+text",
      line=list(color=colors[["negative"]],width=line_width),marker=list(color=colors[["negative"]],size=9))
    p <- plotly::animation_opts(p,frame=frame_ms,transition=0,redraw=TRUE)
    p <- plotly::animation_button(p,x=0,y=-.16)
    p <- plotly::animation_slider(p,x=.18,len=.8,y=-.16,currentvalue=list(prefix="Limit index: "))
  } else p <- plotly::add_trace(p,x=marker$x,y=marker$y,text=marker$text,xaxis="x",yaxis="y2",
    type="scatter",mode="lines+markers",inherit=FALSE,name="Selected upper limit",hoverinfo="x+y+text",
    line=list(color=colors[["negative"]],width=line_width),marker=list(color=colors[["negative"]],size=9))
  if(is.null(title)) title <- paste0("Signed accumulation: A(x) = integral from ",result$lower," to x")
  p <- plotly::layout(p,title=title,xaxis=list(title=xlab,range=xlim),
    yaxis=list(title=ylab,zeroline=TRUE,rangemode="tozero"),
    yaxis2=list(title=alab,zeroline=TRUE,rangemode="tozero"),margin=list(t=90,b=130))
  p$height <- 800
  plotly::config(p,displaylogo=FALSE,scrollZoom=TRUE)
}
