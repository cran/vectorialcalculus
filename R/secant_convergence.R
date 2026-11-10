#' Tabulate secants at decreasing steps
#' @inheritParams vc_tangent
#' @param steps Nonempty strictly decreasing positive finite real step vector.
#'   A single step is allowed. Steps are not silently sorted or deduplicated.
#' @return Plain list with type, x0, y0, method, steps, results (vc_secant
#'   objects), tangent (vc_tangent with full derivative diagnostics), table
#'   (step, left, right, slope, slope_difference), and note. slope_difference
#'   compares with a numerical tangent slope, not an exact derivative error.
#' @details Recovers the decreasing-step procedure of legacy secant_tangent.
#'   Forward/backward secants join x0 to x0 plus/minus each step. Central
#'   secants join x0-step and x0+step: their lines pass through both endpoints,
#'   not generally through (x0,f(x0)). Reuses vc_secant and vc_tangent without
#'   changing their algorithms. Collapsed/nonfinite endpoints are errors.
#'   Numerical sequences and apparent agreement do not prove differentiability
#'   or convergence. Failed tangent diagnostics are retained without warnings.
#' @examples
#' (function(x)x^2) |> vc_secant_convergence(x0=1)
#' abs |> vc_secant_convergence(x0=0,method="central")
#' @export
vc_secant_convergence <- function(f,x0,steps=2^(-(1:8)),method=c("forward","backward","central"),
                                   h=NULL,atol=1e-8,rtol=1e-4) {
  method <- match.arg(method)
  x0 <- .vc_finite(x0,"x0")
  if(!is.numeric(steps)||is.complex(steps)||!is.null(dim(steps))||!length(steps)||
     any(!is.finite(steps))||any(steps<=0)||any(diff(steps)>=0)) stop("'steps' must be strictly decreasing positive finite values.")
  steps <- as.double(steps)
  tangent <- vc_tangent(f,x0,method=method,h=h,atol=atol,rtol=rtol,warn=FALSE)
  results <- lapply(steps,function(s) {
    left <- if(method=="forward")x0 else .vc_finite(x0-s,"Left endpoint")
    right <- if(method=="backward")x0 else .vc_finite(x0+s,"Right endpoint")
    if((method!="forward"&&left>=x0)||(method!="backward"&&right<=x0)) stop("Secant endpoint collapsed; enlarge the step.")
    vc_secant(f,left,right)
  })
  slopes <- vapply(results,function(z)z$slope,numeric(1))
  delta <- slopes-tangent$slope
  if(any(!is.finite(delta))) stop("Slope comparison overflowed.")
  list(type="secant_convergence",x0=x0,y0=tangent$y0,method=method,steps=steps,
    results=results,tangent=tangent,table=data.frame(step=steps,
      left=vapply(results,function(z)z$x0,numeric(1)),right=vapply(results,function(z)z$x1,numeric(1)),
      slope=slopes,slope_difference=delta),note="Finite secants and tangent diagnostics do not prove differentiability or convergence.")
}

#' Animate stored secants approaching a numerical tangent
#' @param result Result from vc_secant_convergence().
#' @inheritParams vc_plot_secant
#' @param animate Logical; multiple stored steps get playback/slider. FALSE
#'   displays the last step; one step is always static.
#' @param frame_ms Positive finite frame duration in milliseconds. Transition
#'   is zero to avoid interpolating artificial intermediate secants.
#' @param line_width Positive finite animated secant width.
#' @return Plotly widget with sampled function, numerical tangent and its
#'   diagnostic label, actual secant lines and actual endpoint pairs.
#' @details No secants or derivatives are recomputed. Stored endpoints are
#'   checked against f. In central mode the secant is the chord of both
#'   endpoints, correcting the legacy line-through-x0 display. The tangent
#'   remains an approximation, including at cusps or domain boundaries.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(x)x^2
#'   f |> vc_secant_convergence(1) |> vc_plot_secant_convergence(f=f)
#' }
#' @export
vc_plot_secant_convergence <- function(result,f,xlim=NULL,n=401L,
    colors=c(curve="#24476b",line="#c44736",point="#222222",differential="#188367",change="#8250a3"),
    title=NULL,xlab="x",ylab="y",animate=TRUE,frame_ms=400,line_width=3) {
  if(!is.list(result)||!identical(result$type,"secant_convergence")||!length(result$results)) stop("Invalid secant sequence.")
  if(!is.logical(animate)||length(animate)!=1||is.na(animate)) stop("Invalid 'animate'.")
  .vc_positive(frame_ms,"frame_ms")
  .vc_positive(line_width,"line_width")
  for(z in result$results) {
    if(!isTRUE(all.equal(c(.vc_value(f,z$x0),.vc_value(f,z$x1)),c(z$y0,z$y1),tolerance=0)))
      stop("'f' does not match stored secant endpoints.")
  }
  bounds <- range(c(result$table$left,result$table$right,result$x0))
  if(is.null(xlim)) xlim <- bounds+c(-.5,.5)*.vc_finite(diff(bounds),"Endpoint span")
  xlim <- .vc_interval(xlim)
  if(xlim[1]>bounds[1]||xlim[2]<bounds[2]) stop("'xlim' must cover all stored secant endpoints.")
  if(is.null(title)) title <- paste("Decreasing-step secants:",result$method)
  p <- .vc_plot_model(result$tangent,f,"tangent",xlim,n,colors,title,xlab,ylab)
  indices <- if(animate)seq_along(result$results) else length(result$results)
  lines <- points <- list()
  for(i in indices) {
    z <- result$results[[i]]
    y <- z$y0+z$slope*(xlim-z$x0)
    if(any(!is.finite(y))) stop("Secant line overflowed; narrow 'xlim'.")
    text <- paste0("step=",result$steps[i],"; slope=",format(z$slope,digits=9),"; ",result$method)
    lines[[length(lines)+1L]] <- data.frame(x=xlim,y=y,frame=i,text=text)
    points[[length(points)+1L]] <- data.frame(x=c(z$x0,z$x1),y=c(z$y0,z$y1),frame=i,text=text)
  }
  lines <- do.call(rbind,lines)
  points <- do.call(rbind,points)
  if(length(indices)>1) {
    p <- plotly::add_trace(p,data=lines,x=~x,y=~y,frame=~frame,text=~text,inherit=FALSE,
      type="scatter",mode="lines",name="Stored secant",line=list(color=colors[["differential"]],width=line_width))
    p <- plotly::add_trace(p,data=points,x=~x,y=~y,frame=~frame,text=~text,inherit=FALSE,
      type="scatter",mode="markers",name="Secant endpoints",marker=list(color=colors[["change"]],size=9))
    p <- plotly::animation_opts(p,frame=frame_ms,transition=0,redraw=TRUE)
    p <- plotly::animation_button(p,x=0,y=-.15)
    p <- plotly::animation_slider(p,x=.18,len=.8,y=-.15,currentvalue=list(prefix="Step index: "))
  } else {
    p <- plotly::add_trace(p,x=lines$x,y=lines$y,text=lines$text,inherit=FALSE,type="scatter",mode="lines",
      name="Stored secant",line=list(color=colors[["differential"]],width=line_width))
    p <- plotly::add_trace(p,x=points$x,y=points$y,text=points$text,inherit=FALSE,type="scatter",mode="markers",
      name="Secant endpoints",marker=list(color=colors[["change"]],size=9))
  }
  plotly::layout(p,margin=list(b=120))
}
