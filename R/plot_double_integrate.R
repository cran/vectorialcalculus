#' Display stored Cartesian integration sections
#' @param result Result of vc_double_integrate().
#' @param title Optional character title.
#' @return Plotly widget of sampled boundary polygon and stored sections.
#' @details No function evaluations or quadrature occur here. Straight segments
#'   connect stored boundary samples, not an exact curved boundary. Physical
#'   x,y axes are preserved for either order. No signed height is represented.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   r <- vc_double_integrate(function(v)sum(v),c(0,1),function(x)0,function(x)x)
#'   vc_plot_double_integral(r)
#' }
#' @export
vc_plot_double_integral <- function(result,title=NULL) {
  if(!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"double_integral")||
     !result$order %in% c("yx","xy")) stop("Invalid double integral result.")
  if(!is.null(title)&&(!is.character(title)||length(title)!=1||is.na(title)))
    stop("Invalid title.")
  d <- result$sections
  if(!is.data.frame(d)||!nrow(d)) stop("No stored sections to display.")
  d <- d[order(d$outer),,drop=FALSE]
  a <- c(d$outer,rev(d$outer),d$outer[1])
  b <- c(d$lower,rev(d$upper),d$lower[1])
  xy <- function(a,b)if(result$order=="yx")list(x=a,y=b) else list(x=b,y=a)
  polygon <- xy(a,b)
  p <- plotly::plot_ly(x=polygon$x,y=polygon$y,type="scatter",mode="lines",
    fill="toself",name="Sampled domain")
  lines <- xy(as.vector(rbind(d$outer,d$outer,NA_real_)),
    as.vector(rbind(d$lower,d$upper,NA_real_)))
  p <- plotly::add_trace(p,x=lines$x,y=lines$y,type="scatter",mode="lines",
    name="Stored sections",inherit=FALSE)
  plotly::layout(p,title=if(is.null(title))"Cartesian domain: sampled sections" else title,
    xaxis=list(title="x"),yaxis=list(title="y",scaleanchor="x",scaleratio=1))
}

