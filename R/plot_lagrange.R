#' View equality constraints, objective contours and Lagrange gradients
#' @param result Two-dimensional vc_lagrange result, including inconclusive or
#'   infeasible assessments. This is a supplied-candidate view, not an optimizer.
#' @param f,g Same objective and constraints used in calculation.
#' @param xlim,ylim Optional finite increasing bounds; NULL uses point +/- 1.
#' @param n Integer samples per coordinate, at least two. Stored contour
#'   matrices use y rows and x columns. Unavailable samples are gaps.
#' @param vector_scale Positive common display multiplier for gradient vectors.
#'   Magnitudes are not normalized; hover retains the unscaled estimates.
#' @param show_gradients Logical; show objective and each constraint gradient
#'   when the stored first-derivative data are available.
#' @param show_weighted Logical; additionally show lambda_i times grad(g_i).
#'   Their SUM, not each individual vector, matches grad(f) at a regular
#'   stationary candidate. They share the candidate as their display origin.
#' @param colors Named point, gradient and weighted colors.
#' @param constraint_colors One color or one per constraint; the same color is
#'   used for its zero set and unweighted gradient.
#' @param colorscale Plotly name/list for the objective contours.
#' @param alpha Objective contour opacity between zero and one.
#' @param line_width Positive line/arrow width.
#' @param marker_size Positive candidate marker size.
#' @param title Optional scalar title; default reports conditional classification,
#'   reason and numerical constraint-rank diagnostic.
#' @param xlab,ylab Scalar coordinate-axis labels.
#' @return Plotly contour widget with candidate marker, zero constraints and
#'   stored gradient vectors, plus arrow annotations. Calculation and derivative
#'   assessment are not rerun. Legend controls toggle traces; arrowhead
#'   annotations are static overlays controlled by show_gradients/show_weighted.
#' @details Functions are checked against stored candidate values, then sampled
#'   for context. A grid may miss singular/tangential/narrow zero sets; absence
#'   of a rendered zero contour does not mean no feasible point exists. Finite
#'   arrow endpoints that overflow are omitted. Labels remain conditional on
#'   the mathematical hypotheses and numerical evidence in result; this plot
#'   does not certify regularity, stationarity or an extremum.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(v)sum(v^2)
#'   g <- function(v)sum(v)-1
#'   f |> vc_lagrange(g,c(.5,.5),assume_c2=TRUE,hessian_h=.01) |>
#'     vc_plot_lagrange(f=f,g=g)
#' }
#' @export
vc_plot_lagrange <- function(result,f,g,xlim=NULL,ylim=NULL,n=81L,vector_scale=.3,
    show_gradients=TRUE,show_weighted=FALSE,
    colors=c(point="#222222",gradient="#c44736",weighted="#dd8800"),
    constraint_colors="#8250a3",colorscale="Viridis",alpha=.6,line_width=2,marker_size=10,
    title=NULL,xlab="x",ylab="y") {
  .vc_application_style(colors,c("point","gradient","weighted"),alpha,line_width,title)
  if(!is.list(result)||!identical(result$type,"lagrange")||length(result$point)!=2)
    stop("A two-dimensional vc_lagrange result is required.")
  .vc_integer(n,"n",2);.vc_positive(vector_scale,"vector_scale");.vc_positive(marker_size,"marker_size")
  .vc_flag(show_gradients,"show_gradients");.vc_flag(show_weighted,"show_weighted")
  if(!is.character(xlab)||length(xlab)!=1||is.na(xlab)||!is.character(ylab)||length(ylab)!=1||is.na(ylab))
    stop("Invalid axis labels.")
  point <- result$point
  constraints <- .vc_constraints(g,point);m <- length(constraints$value)
  if(!isTRUE(all.equal(.vc_field_value(f,point),result$value,tolerance=1e-12))||
     !isTRUE(all.equal(constraints$value,result$constraints,tolerance=1e-12)))
    stop("Functions do not match stored candidate values.")
  if(!is.character(constraint_colors)||!length(constraint_colors)%in%c(1L,m)||
     anyNA(constraint_colors)||any(!nzchar(constraint_colors)))stop("Invalid 'constraint_colors'.")
  constraint_colors <- rep(constraint_colors,length.out=m)
  if(is.null(xlim))xlim <- point[1]+c(-1,1)
  if(is.null(ylim))ylim <- point[2]+c(-1,1)
  xlim <- .vc_interval(xlim);ylim <- .vc_interval(ylim)
  x <- seq(xlim[1],xlim[2],length.out=n);y <- seq(ylim[1],ylim[2],length.out=n)
  named <- function(v)stats::setNames(v,names(point))
  sample <- function(fun)vapply(x,function(a)vapply(y,function(b)
    .vc_field_safe(fun,named(c(a,b))),numeric(1)),numeric(n))
  z <- sample(f)
  if(!any(is.finite(z)))stop("No finite objective contour samples.")
  p <- plotly::plot_ly(x=x,y=y,z=z,type="contour",connectgaps=FALSE,colorscale=colorscale,
    opacity=alpha,name="Objective contours",showscale=TRUE,
    hovertemplate="x=%{x}<br>y=%{y}<br>f=%{z}<extra></extra>")
  for(i in seq_len(m)) {
    z <- sample(function(v)constraints$evaluate(v)[i])
    p <- plotly::add_trace(p,x=x,y=y,z=z,type="contour",inherit=FALSE,
      name=paste0("g",i," = 0"),showscale=FALSE,connectgaps=FALSE,autocontour=FALSE,
      contours=list(start=0,end=0,size=1,coloring="lines",showlabels=TRUE),
      colorscale=list(list(0,constraint_colors[i]),list(1,constraint_colors[i])),
      line=list(width=line_width),hovertemplate=paste0("g",i,"=%{z}<extra></extra>"))
  }
  p <- plotly::add_trace(p,x=point[1],y=point[2],type="scatter",mode="markers",inherit=FALSE,
    name="Assessed candidate",marker=list(size=marker_size,color=colors[["point"]]),
    text=paste0("f=",result$value,"<br>constraint residual=",max(abs(result$constraints)),
      "<br>classification=",result$classification,"<br>",result$reason),hoverinfo="text")
  arrows <- list()
  add_vector <- function(p,v,label,color) {
    endpoint <- point+vector_scale*as.double(v)
    if(any(!is.finite(endpoint)))return(p)
    arrows[[length(arrows)+1L]] <<- list(x=endpoint[1],y=endpoint[2],ax=point[1],ay=point[2],
      xref="x",yref="y",axref="x",ayref="y",text="",showarrow=TRUE,arrowhead=2,
      arrowcolor=color,arrowwidth=line_width)
    plotly::add_trace(p,x=c(point[1],endpoint[1]),y=c(point[2],endpoint[2]),type="scatter",
      mode="lines",inherit=FALSE,name=label,line=list(color=color,width=line_width),
      text=paste0(label,"<br>unscaled vector=(",paste(format(v),collapse=", "),")"),hoverinfo="text")
  }
  if(show_gradients&&!is.null(result$gradient$estimate))
    p <- add_vector(p,result$gradient$estimate,"Objective gradient",colors[["gradient"]])
  if(!is.null(result$jacobian$estimate))for(i in seq_len(m)) {
    if(show_gradients)p <- add_vector(p,result$jacobian$estimate[i,],paste0("Constraint gradient ",i),constraint_colors[i])
    if(show_weighted&&!is.null(result$lambda))p <- add_vector(p,
      result$lambda[i]*result$jacobian$estimate[i,],paste0("lambda",i," * gradient g",i),colors[["weighted"]])
  }
  if(is.null(title))title <- paste("Lagrange:",result$classification,"(conditional); rank",result$rank,
    "of",m,"<br>",result$reason)
  plotly::layout(p,title=title,annotations=arrows,xaxis=list(title=xlab,range=xlim),yaxis=list(title=ylab,range=ylim))
}
