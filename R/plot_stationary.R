#' Inspect two-dimensional stationary candidates and search evaluations
#' @param result Result of vc_stationary_points with two input coordinates.
#' @param f Same scalar field used in the search, accepting one vector.
#' @param xlim,ylim Finite increasing display bounds; NULL uses the search
#'   window, or the range of starts/candidates enlarged by one coordinate unit.
#' @param n Integer samples per axis, at least two. Matrix rows represent y
#'   and columns x. Nonfinite or unavailable field samples remain gaps.
#' @param run Integer run index whose objective evaluation history is shown.
#' @param show_history,animate Logical switches. Playback shows stored objective
#'   evaluations, not optimizer iterations or a Newton/descent path.
#' @param max_frames Positive integer maximum displayed history samples.
#'   Deterministic thinning keeps first/last samples; full history stays in result.
#' @param frame_ms Positive playback duration in milliseconds.
#' @param colors Named minimum, maximum, saddle, inconclusive and history colors.
#' @param colorscale Plotly colorscale name or list for the sampled contours.
#' @param alpha Opacity between zero and one for search-evaluation markers.
#' @param marker_size Positive candidate marker size.
#' @param line_width Positive history-line width.
#' @param title Optional scalar title.
#' @param xlab,ylab Scalar axis labels.
#' @return Plotly contour widget with conditional candidate labels and optional
#'   evaluation playback. Hover retains residuals, stopping codes and statuses.
#' @details No derivative or optimization is recomputed. Only f is sampled for
#'   context and checked against stored candidate values. Search evaluations may
#'   move uphill and include optimizer probes; their connecting polyline is not
#'   a mathematical flow. Failed residual evaluations are omitted from geometry
#'   but remain in the calculation history. All candidate labels are conditional
#'   numerical evidence; this plot does not establish regularity or completeness.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(v)v[1]^2+2*v[2]^2
#'   r <- vc_stationary_points(f,c(.1,.2),assume_c2=TRUE)
#'   r |> vc_plot_stationary(f=f)
#' }
#' @export
vc_plot_stationary <- function(result,f,xlim=NULL,ylim=NULL,n=61L,run=1L,
    show_history=TRUE,animate=TRUE,max_frames=40L,frame_ms=300,
    colors=c(minimum="#188367",maximum="#c44736",saddle="#8250a3",
      inconclusive="#dd8800",history="#24476b"),colorscale="Viridis",alpha=.6,
    marker_size=10,line_width=2,title=NULL,xlab="x",ylab="y") {
  .vc_application_style(colors,c("minimum","maximum","saddle","inconclusive","history"),alpha,line_width,title)
  if(!is.list(result)||!identical(result$type,"stationary_points")||ncol(result$starts)!=2)
    stop("A two-dimensional vc_stationary_points result is required.")
  if(!is.function(f)) stop("'f' must be the same scalar field.")
  .vc_integer(n,"n",2);.vc_integer(run,"run",1);.vc_integer(max_frames,"max_frames",2)
  if(run>length(result$runs)) stop("Invalid 'run' index.")
  .vc_flag(show_history,"show_history");.vc_flag(animate,"animate")
  .vc_positive(frame_ms,"frame_ms");.vc_positive(marker_size,"marker_size")
  if(!is.character(xlab)||length(xlab)!=1||is.na(xlab)||!is.character(ylab)||length(ylab)!=1||is.na(ylab))
    stop("Invalid axis labels.")
  named <- function(v)stats::setNames(v,colnames(result$starts))
  if(nrow(result$points))for(i in seq_len(nrow(result$points))) {
    z <- .vc_field_value(f,named(result$points[i,]))
    if(abs(z-result$table$value[i])>1e-8*max(1,abs(z),abs(result$table$value[i])))
      stop("Field does not match the stored candidate values.")
  }
  bounds <- function(i)if(!is.null(result$bounds))result$bounds[i,] else
    range(c(result$starts[,i],result$points[,i]))+c(-1,1)
  if(is.null(xlim))xlim <- bounds(1)
  if(is.null(ylim))ylim <- bounds(2)
  xlim <- .vc_interval(xlim);ylim <- .vc_interval(ylim)
  x <- seq(xlim[1],xlim[2],length.out=n);y <- seq(ylim[1],ylim[2],length.out=n)
  z <- vapply(x,function(a)vapply(y,function(b).vc_field_safe(f,named(c(a,b))),numeric(1)),numeric(n))
  if(!any(is.finite(z))) stop("No finite contour samples.")
  p <- plotly::plot_ly(x=x,y=y,z=z,type="contour",colorscale=colorscale,
    connectgaps=FALSE,name="Sampled field",showscale=TRUE,
    hovertemplate="x=%{x}<br>y=%{y}<br>f=%{z}<extra>Sampled field</extra>")
  for(label in c("minimum","maximum","saddle","inconclusive")) {
    ids <- which(result$table$classification==label)
    if(length(ids))p <- plotly::add_trace(p,x=result$points[ids,1],y=result$points[ids,2],
      type="scatter",mode="markers",inherit=FALSE,name=paste("Conditional",label),
      marker=list(size=marker_size,color=colors[[label]],symbol="diamond"),
      text=paste0("Candidate ",ids,"<br>f=",format(result$table$value[ids]),
        "<br>gradient residual=",format(result$table$residual[ids]),"<br>",label," (conditional)"),
      hoverinfo="text")
  }
  h <- result$runs[[run]]$history
  valid <- which(vapply(h,function(a)is.finite(a$residual),logical(1)))
  if(show_history&&length(valid)) {
    selected <- valid[unique(round(seq(1,length(valid),length.out=min(max_frames,length(valid)))))]
    points <- do.call(rbind,lapply(h[selected],function(a)a$point))
    text <- paste0("Objective evaluation ",selected,"<br>residual=",
      vapply(h[selected],function(a)format(a$residual),character(1)),
      "<br>Optimizer stopping code: ",result$runs[[run]]$optimizer_code,
      "<br>Run status: ",result$runs[[run]]$status)
    if(animate&&length(selected)>1) {
      d <- do.call(rbind,lapply(seq_along(selected),function(i)data.frame(
        x=points[seq_len(i),1],y=points[seq_len(i),2],text=text[seq_len(i)],frame=i)))
      p <- plotly::add_trace(p,data=d,x=~x,y=~y,text=~text,frame=~frame,
        type="scatter",mode="lines+markers",inherit=FALSE,name="Objective evaluations",
        opacity=alpha,line=list(color=colors[["history"]],width=line_width),hoverinfo="text")
      p <- plotly::animation_opts(p,frame=frame_ms,transition=0,redraw=TRUE)
      p <- plotly::animation_slider(p,currentvalue=list(prefix="Stored evaluation sample: "))
    } else p <- plotly::add_trace(p,x=points[,1],y=points[,2],text=text,
      type="scatter",mode="lines+markers",inherit=FALSE,name="Objective evaluations",
      opacity=alpha,line=list(color=colors[["history"]],width=line_width),hoverinfo="text")
  }
  if(is.null(title))title <- "Stationary candidates: conditional labels; objective probes are not iterations"
  plotly::layout(p,title=title,xaxis=list(title=xlab,range=xlim),yaxis=list(title=ylab,range=ylim))
}
