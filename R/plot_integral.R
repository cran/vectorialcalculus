#' Visualize signed numerical quadrature
#' @param result Result of vc_integrate().
#' @inheritParams vc_plot_riemann
#' @param colors Named curve, approximation and nodes Plotly/CSS colors.
#' @param line_width Positive finite curve/approximation line width.
#' @param xscale,yscale Plotly axis scale: linear or log. Log scales may hide
#'   nonpositive data and cannot faithfully show signed areas across zero.
#' @return Plotly widget. Fixed rules show finest stored nodes and piecewise
#'   linear/parabolic interpolation with fill to zero. Adaptive results show
#'   only a sampled function with fill to zero, not QUADPACK's internal mesh.
#' @details Plotting never recomputes quadrature. Stored node values are checked
#'   against f. Simpson arcs interpolate each consecutive pair of cells.
#'   Fill represents a signed-integral region, not an unsigned area calculation.
#'   A reversed orientation changes the integral sign without changing f.
#'   Hover reports the stored estimate, diagnostic and orientation. Sampled
#'   curves/fills can miss features and do not establish integrability.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(x)x^2
#'   f |> vc_integrate(c(0,2),method="simpson",n=4) |> vc_plot_integral(f=f)
#' }
#' @export
vc_plot_integral <- function(result,f,n=401L,
                             colors=c(curve="#24476b",approximation="#188367",nodes="#c44736"),
                             alpha=.3,line_width=2,title=NULL,xlab="x",ylab="f(x)",
                             xscale=c("linear","log"),yscale=c("linear","log")) {
  if(!requireNamespace("plotly",quietly=TRUE)) stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"integral")||
     !result$method %in% c("trapezoid","simpson","adaptive")) stop("Invalid integral result.")
  if(!is.function(f)) stop("'f' must be a function.")
  .vc_integer(n,"n",2)
  .vc_positive(line_width,"line_width")
  if(!.vc_scalar(alpha)||alpha<0||alpha>1) stop("'alpha' must be between 0 and 1.")
  if(!is.character(colors)||!all(c("curve","approximation","nodes") %in% names(colors))||
     anyNA(colors)||anyDuplicated(names(colors))||any(!nzchar(colors))) stop("Invalid 'colors'.")
  label <- function(z)is.character(z)&&length(z)==1&&!is.na(z)
  if((!is.null(title)&&!label(title))||!label(xlab)||!label(ylab)) stop("Invalid labels.")
  xscale <- match.arg(xscale)
  yscale <- match.arg(yscale)
  bounds <- sort(.vc_integral_interval(result$interval))
  xs <- seq(bounds[1],bounds[2],length.out=n)
  ys <- vapply(xs,function(x).vc_safe_value(f,x),numeric(1))
  p <- plotly::plot_ly()
  text <- paste0("Rule: ",result$method,"; orientation: ",paste(result$interval,collapse=" to "),
    "; estimate: ",format(result$estimate,digits=9),"; estimated error: ",
    format(result$error_estimate,digits=4),"; ",result$status)
  d <- result$table
  if(!is.null(d)) {
    values <- vapply(d$x,function(x).vc_value(f,x),numeric(1))
    if(!isTRUE(all.equal(values,d$value,tolerance=0))) stop("'f' does not match stored nodes.")
    if(result$method=="trapezoid") {
      ax <- d$x
      ay <- d$value
    } else {
      ax <- ay <- numeric()
      for(i in seq(1,nrow(d)-2,by=2)) {
        xx <- d$x[i+(0:2)]
        yy <- d$value[i+(0:2)]
        t <- seq(0,1,length.out=25)
        u <- (xx[2]-xx[1])/(xx[3]-xx[1])
        ax <- c(ax,xx[1]+t*(xx[3]-xx[1]))
        ay <- c(ay,yy[1]*(t-u)*(t-1)/u+
          yy[2]*t*(t-1)/(u*(u-1))+yy[3]*t*(t-u)/(1-u))
      }
    }
    p <- plotly::add_trace(p,x=ax,y=ay,type="scatter",mode="lines",fill="tozeroy",
      name=paste(result$method,"interpolant"),text=text,hoverinfo="x+y+text",
      opacity=alpha,line=list(color=colors[["approximation"]],width=line_width))
    p <- plotly::add_trace(p,x=d$x,y=d$value,type="scatter",mode="markers",name="Stored nodes",
      marker=list(color=colors[["nodes"]]),text=paste0("signed weight=",d$weight,
        "; contribution=",d$contribution),hoverinfo="x+y+text")
  }
  p <- plotly::add_trace(p,x=xs,y=ys,type="scatter",mode="lines",connectgaps=FALSE,
    name=if(result$method=="adaptive")"Sampled function (not adaptive mesh)" else "Sampled function",
    fill=if(is.null(d)&&result$status!="zero_interval")"tozeroy" else "none",
    text=text,hoverinfo="x+y+text",line=list(color=colors[["curve"]],width=line_width))
  if(is.null(title)) title <- paste("Signed integration:",result$method)
  plotly::layout(p,title=title,xaxis=list(title=xlab,type=xscale),
                 yaxis=list(title=ylab,type=yscale,zeroline=TRUE))
}
