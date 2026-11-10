#' Inspect stored stationary candidates on a sampled 3D graph
#' @param result Two-coordinate result of vc_stationary_points().
#' @param f Same scalar field, accepting one coordinate vector.
#' @param xlim,ylim Finite increasing display bounds; default search bounds or
#'   starts/candidate extent enlarged by one unit.
#' @param n Integer samples per coordinate, at least two. Matrix rows are y,
#'   columns x. Unavailable/nonfinite surface values are gaps.
#' @param colors Named minimum,maximum,saddle,inconclusive Plotly colors.
#' @param opacity Surface opacity between zero and one.
#' @param marker_size Positive candidate marker size.
#' @param title Optional character title.
#' @return Plotly surface and stored candidate markers with stored heights,
#'   conditional classification, residual and reason in hover text.
#' @details No optimization, derivative, Hessian or classification is repeated.
#'   Only f is sampled for the graph and checked at stored candidates. Surface
#'   interpolation is not exact geometry or a smoothness/completeness proof.
#'   Labels remain conditional on the calculation hypotheses and diagnostics;
#'   inconclusive candidates stay inconclusive. No optimizer history is shown
#'   as a path. Use vc_plot_stationary for objective-evaluation playback.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(v)v[1]^2+2*v[2]^2
#'   r <- vc_stationary_points(f,c(.1,.2),assume_c2=TRUE,warn=FALSE)
#'   vc_plot_stationary_surface(r,f)
#' }
#' @export
vc_plot_stationary_surface <- function(result,f,xlim=NULL,ylim=NULL,n=41L,
    colors=c(minimum="#188367",maximum="#c44736",saddle="#8250a3",inconclusive="#dd8800"),
    opacity=.65,marker_size=6,title=NULL) {
  if(!requireNamespace("plotly",quietly=TRUE))stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"stationary_points")||ncol(result$starts)!=2)
    stop("A two-dimensional stationary search result is required.")
  if(!is.function(f))stop("The same scalar field must be supplied.")
  .vc_integer(n,"n",2);.vc_positive(marker_size,"marker_size")
  if(!.vc_scalar(opacity)||opacity<0||opacity>1)stop("Invalid opacity.")
  labels <- c("minimum","maximum","saddle","inconclusive")
  if(!is.character(colors)||!all(labels%in%names(colors))||anyNA(colors)||
     any(!nzchar(colors))||anyDuplicated(names(colors)))stop("Invalid colors.")
  if(!is.null(title)&&(!is.character(title)||length(title)!=1||is.na(title)))stop("Invalid title.")
  named <- function(p)stats::setNames(p,colnames(result$starts))
  if(nrow(result$points))for(i in seq_len(nrow(result$points))) {
    value <- .vc_field_value(f,named(result$points[i,]))
    if(abs(value-result$table$value[i])>1e-8*max(1,abs(value),abs(result$table$value[i])))
      stop("Field does not match stored candidate values.")
  }
  bounds <- function(i)if(!is.null(result$bounds))result$bounds[i,] else
    range(c(result$starts[,i],result$points[,i]))+c(-1,1)
  if(is.null(xlim))xlim <- bounds(1)
  if(is.null(ylim))ylim <- bounds(2)
  xlim <- .vc_interval(xlim);ylim <- .vc_interval(ylim)
  x <- seq(xlim[1],xlim[2],length.out=n);y <- seq(ylim[1],ylim[2],length.out=n)
  z <- vapply(x,function(a)vapply(y,function(b).vc_field_safe(f,named(c(a,b))),numeric(1)),numeric(n))
  if(!any(is.finite(z)))stop("No finite surface samples.")
  p <- plotly::plot_ly(x=x,y=y,z=z,type="surface",opacity=opacity,connectgaps=FALSE,
    name="Sampled field",showscale=FALSE)
  for(label in labels) {
    ids <- which(result$table$classification==label)
    if(length(ids))p <- plotly::add_trace(p,x=result$points[ids,1],y=result$points[ids,2],
      z=result$table$value[ids],type="scatter3d",mode="markers",inherit=FALSE,
      name=paste("Conditional",label),marker=list(size=marker_size,color=colors[[label]]),
      text=paste0("Candidate ",ids,"<br>f=",format(result$table$value[ids]),
        "<br>gradient residual=",format(result$table$residual[ids]),
        "<br>",vapply(result$assessments[ids],function(a)a$reason,character(1))),hoverinfo="text")
  }
  plotly::layout(p,title=if(is.null(title))"Stationary candidates: conditional labels on sampled graph" else title,
    scene=list(xaxis=list(title="x",range=xlim),yaxis=list(title="y",range=ylim),
      zaxis=list(title="f(x,y)"),aspectmode="data"))
}
