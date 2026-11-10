#' Draw signed Riemann rectangles with optional refinement animation
#' @param result Result from vc_riemann() or vc_riemann_convergence().
#' @inheritParams vc_plot_secant
#' @param colors Named curve, positive, negative and edge Plotly/CSS colors.
#' @param alpha Finite fill opacity between zero and one.
#' @param line_width Nonnegative finite rectangle edge width.
#' @param animate Logical; animate multiple stored refinements. A single
#'   stored sum is static. No new sums are computed by this function.
#' @param frame_ms Positive finite animation frame duration in milliseconds.
#' @return Plotly widget with sampled curve and zero-based signed rectangles.
#' @details Rectangles extend from zero to f(sample), including below zero.
#'   This baseline is deliberately fixed to preserve the integral interpretation.
#'   Invalid curve samples are gaps. The curve is sampled, not an integrability
#'   proof. Hover shows n, rule, signed sum and rectangle contribution.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(x)x^2
#'   f |> vc_riemann_convergence(c(0,2),reference=8/3) |> vc_plot_riemann(f=f)
#' }
#' @export
vc_plot_riemann <- function(result,f,n=401L,
                            colors=c(curve="#24476b",positive="#188367",negative="#c44736",edge="#333333"),
                            alpha=.45,line_width=1,animate=TRUE,frame_ms=700,
                            title=NULL,xlab="x",ylab="f(x)") {
  if(!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if(!is.list(result)||!result$type %in% c("riemann","riemann_convergence")) stop("Invalid Riemann result.")
  .vc_integer(n,"n",2)
  .vc_positive(line_width,"line_width",TRUE)
  .vc_positive(frame_ms,"frame_ms")
  if(!.vc_scalar(alpha)||alpha<0||alpha>1) stop("'alpha' must be between 0 and 1.")
  if(!is.logical(animate)||length(animate)!=1||is.na(animate)) stop("Invalid 'animate'.")
  if(!is.character(colors)||!all(c("curve","positive","negative","edge") %in% names(colors))||
     anyNA(colors)||anyDuplicated(names(colors))||any(!nzchar(colors))) stop("Invalid 'colors'.")
  label <- function(z)is.character(z)&&length(z)==1&&!is.na(z)
  if((!is.null(title)&&!label(title))||!label(xlab)||!label(ylab)) stop("Invalid labels.")
  sums <- if(result$type=="riemann") list(result) else result$results
  if(!length(sums)) stop("No stored sums.")
  for(a in sums) {
    values <- vapply(a$table$sample,function(x).vc_value(f,x),numeric(1))
    if(!isTRUE(all.equal(values,a$table$value,tolerance=0))) stop("'f' does not match stored samples.")
  }
  if(!animate) sums <- sums[length(sums)]
  rows <- list()
  for(k in seq_along(sums)) {
    a <- sums[[k]]
    for(i in seq_len(nrow(a$table))) {
      d <- a$table[i,]
      rows[[length(rows)+1L]] <- data.frame(x=c(d$left,d$right,d$right,d$left,d$left,NA_real_),
        y=c(0,0,d$value,d$value,0,NA_real_),frame=k,sign=if(d$value<0)"negative" else "positive",
        text=paste0("n=",a$partition$n,"; ",a$method,"; signed sum=",format(a$estimate,digits=8),
                    "; contribution=",format(d$contribution,digits=8)))
    }
  }
  rectangles <- do.call(rbind,rows)
  xs <- seq(result$interval[1],result$interval[2],length.out=n)
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p,x=xs,y=vapply(xs,function(x).vc_safe_value(f,x),numeric(1)),
      type="scatter",mode="lines",name="Sampled function",connectgaps=FALSE,line=list(color=colors[["curve"]]))
  for(s in unique(rectangles$sign)) {
    d <- rectangles[rectangles$sign==s,,drop=FALSE]
    if(length(sums)>1) for(k in seq_along(sums)) {
      if(!any(d$frame==k)) d <- rbind(d,data.frame(x=NA_real_,y=NA_real_,frame=k,sign=s,text="No contributions of this sign"))
    }
    if(length(sums)>1) p <- plotly::add_trace(p,data=d,x=~x,y=~y,frame=~frame,text=~text,
      type="scatter",mode="lines",fill="toself",name=paste(s,"contributions"),
      fillcolor=colors[[s]],opacity=alpha,line=list(color=colors[["edge"]],width=line_width)) else
      p <- plotly::add_trace(p,x=d$x,y=d$y,text=d$text,type="scatter",mode="lines",fill="toself",
      name=paste(s,"contributions"),fillcolor=colors[[s]],opacity=alpha,
      line=list(color=colors[["edge"]],width=line_width))
  }
  if(length(sums)>1) {
    p <- plotly::animation_opts(p,frame=frame_ms,transition=0,redraw=TRUE)
    p <- plotly::animation_button(p,x=0,y=-.15)
    p <- plotly::animation_slider(p,x=.18,len=.8,y=-.15,currentvalue=list(prefix="Refinement: "))
  }
  if(is.null(title)) title <- "Signed tagged Riemann sums (not Darboux bounds)"
  plotly::layout(p,title=title,xaxis=list(title=xlab),yaxis=list(title=ylab,zeroline=TRUE),
                 margin=list(b=120),legend=list(orientation="h"))
}
