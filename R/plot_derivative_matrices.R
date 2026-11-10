#' Inspect a numerical derivative matrix and its stored refinements
#' @param result A vc_jacobian() or vc_hessian() result.
#' @param colors Named negative, zero, positive and unresolved Plotly/CSS colors.
#' @param alpha Finite heatmap opacity between zero and one.
#' @param digits Integer significant digits in displayed cell labels, 1 to 15;
#'   hover retains a more precise estimate.
#' @param show_values Logical; display rounded values as cell annotations.
#' @param show_diagnostics Logical; mark entries whose full sensitivity checks
#'   failed. Flags refer to the three-scale diagnostic, not only the visible level.
#' @param slider Logical; allow selection among the three stored refinement
#'   matrices. FALSE draws only the finest matrix.
#' @param title Optional scalar title. Default identifies matrix type, the
#'   selected relative refinement level and the numerical-only interpretation.
#' @return Plotly matrix widget with correctly labeled rows and columns,
#'   numeric hover, optional values, unresolved-entry markers and refinement
#'   selector. Slider changes both the matrix and its value annotations.
#' @details Jacobian rows are outputs and columns are input coordinates;
#'   Hessian rows/columns are input coordinates. The first row is at the top.
#'   No derivatives are recomputed. Relative levels mean initial steps, halves
#'   and quarters, not a common absolute step for differently scaled coordinates.
#'   Hessian symmetry is imposed by its rectangle stencil, not established as
#'   equality/existence of ordered mixed partials. No extremum classification
#'   or total differentiability certificate is inferred from colors or flags.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(v)v[1]^2+v[1]*v[2]+2*v[2]^2
#'   f |> vc_hessian(c(1,1),h=.01) |> vc_plot_matrix()
#' }
#' @export
vc_plot_matrix <- function(result,colors=c(negative="#24476b",zero="#ffffff",positive="#c44736",unresolved="#dd8800"),
                           alpha=.95,digits=5L,show_values=TRUE,show_diagnostics=TRUE,slider=TRUE,title=NULL) {
  .vc_application_style(colors,c("negative","zero","positive","unresolved"),alpha,1,title)
  if(!is.list(result)||!result$type%in%c("jacobian","hessian")||!is.matrix(result$estimate)||
     !is.numeric(result$estimate)||is.complex(result$estimate)||!length(result$estimate)||any(!is.finite(result$estimate)))
    stop("Invalid derivative matrix result.")
  .vc_integer(digits,"digits",1)
  if(digits>15) stop("'digits' must not exceed 15.")
  flag <- function(z)is.logical(z)&&length(z)==1&&!is.na(z)
  if(!all(vapply(list(show_values,show_diagnostics,slider),flag,logical(1)))) stop("Invalid display switch.")
  if(!is.matrix(result$reliable)||!is.logical(result$reliable)||anyNA(result$reliable)||
     !identical(dim(result$reliable),dim(result$estimate))) stop("Invalid reliability matrix.")
  matrices <- result$refinements
  if(length(matrices)!=3||!all(vapply(matrices,function(z)is.matrix(z)&&is.numeric(z)&&!is.complex(z)&&
      identical(dim(z),dim(result$estimate))&&all(is.finite(z)),logical(1)))) stop("Invalid stored refinements.")
  nr <- nrow(result$estimate);nc <- ncol(result$estimate)
  x <- seq_len(nc);y <- seq_len(nr)
  xlabels <- colnames(result$estimate);ylabels <- rownames(result$estimate)
  if(is.null(xlabels)) xlabels <- paste0("x",x)
  if(is.null(ylabels)) ylabels <- paste0("row",y)
  limit <- max(abs(unlist(matrices)))
  if(limit==0) limit <- 1
  annotations <- function(k) {
    if(!show_values) return(list())
    unlist(lapply(y,function(i)lapply(x,function(j)list(x=j,y=i,xref="x",yref="y",showarrow=FALSE,
      text=format(matrices[[k]][i,j],digits=digits,trim=TRUE),
      font=list(color=if(abs(matrices[[k]][i,j])>.55*limit)"#ffffff" else "#222222")))),recursive=FALSE)
  }
  caption <- function(k) {
    base <- if(is.null(title))paste("Numerical",result$type) else title
    paste0(base,"; relative refinement ",k,"/3<br><sup>Entry checks are not a differentiability or smoothness proof</sup>")
  }
  levels <- if(slider)1:3 else 3L
  p <- plotly::plot_ly()
  for(k in levels) {
    detail <- matrix("",nr,nc)
    for(entry in result$results) {
      i <- entry$row;j <- entry$column;d <- entry$derivative
      if(result$type=="jacobian") extra <- paste0("; step=",format(d$diagnostics$step[k],digits=6),"; stable=",d$stable) else {
        steps <- if(i==j)rep(d$diagnostics$step[k],2) else c(d$diagnostics$step_i[k],d$diagnostics$step_j[k])
        extra <- paste0("; stencil coordinates=",xlabels[i],",",xlabels[j],"; steps=",paste(format(steps,digits=6),collapse=","),
          "; stable=",d$stable,"; resolved=",d$resolved,"; roundoff indicator=",format(d$diagnostics$roundoff_indicator[k],digits=6))
      }
      detail[i,j] <- extra
      if(result$type=="hessian"&&i!=j) detail[j,i] <- paste0(extra,"; mirrored upper stencil")
    }
    text <- matrix("",nr,nc)
    for(i in y) for(j in x) text[i,j] <- paste0("row ",ylabels[i],"; column ",xlabels[j],
      "; estimate=",format(matrices[[k]][i,j],digits=12),"; three-scale checks=",result$reliable[i,j],detail[i,j])
    p <- plotly::add_trace(p,x=x,y=y,z=matrices[[k]],text=text,type="heatmap",inherit=FALSE,
      name=paste("Refinement",k),visible=k==3,opacity=alpha,zmid=0,zmin=-limit,zmax=limit,
      colorscale=list(list(0,colors[["negative"]]),list(.5,colors[["zero"]]),list(1,colors[["positive"]])),
      hoverinfo="text",showscale=TRUE)
  }
  unresolved <- which(!result$reliable,arr.ind=TRUE)
  marked <- show_diagnostics&&nrow(unresolved)>0
  if(marked) p <- plotly::add_trace(p,x=unresolved[,2],y=unresolved[,1],type="scatter",mode="markers",
    inherit=FALSE,name="Unresolved entry checks",marker=list(symbol="square-open",size=32,color=colors[["unresolved"]],line=list(width=3)),
    text="Sensitivity or roundoff check failed; inspect stored diagnostics",hoverinfo="text")
  steps <- lapply(1:3,function(k)list(method="update",label=c("initial","half","quarter")[k],
    args=list(list(visible=c(seq_along(levels)==k,if(marked)TRUE else logical())),
              list(annotations=annotations(k),title=caption(k)))))
  p <- plotly::layout(p,title=caption(3),annotations=annotations(3),
    xaxis=list(title="Input coordinate",tickmode="array",tickvals=x,ticktext=xlabels,side="top",range=c(.5,nc+.5)),
    yaxis=list(title=if(result$type=="jacobian")"Output component" else "Input coordinate",
      tickmode="array",tickvals=y,ticktext=ylabels,range=c(nr+.5,.5)),margin=list(t=120,b=100))
  if(slider) p <- plotly::layout(p,sliders=list(list(active=2,currentvalue=list(prefix="Step level: "),steps=steps)))
  plotly::config(p,displaylogo=FALSE)
}
