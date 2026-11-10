#' Probe a finite numerical limit from one or both sides
#' @inheritParams vc_secant
#' @param x0 Finite real approach point. f(x0) is never evaluated.
#' @param side One of both, left or right.
#' @param h Positive initial distance. NULL uses 0.1*max(1,abs(x0)).
#' @param steps Integer number of halving distances, between 4 and 1000.
#' @param atol,rtol Nonnegative absolute/relative comparison tolerances.
#' @return Plain list with type, x0, side, samples (x, distance, value, side),
#'   left, right (side summaries with estimate and stable), estimate (NA unless
#'   applicable checks pass), status, and note. status is compatible_finite,
#'   side_disagreement, unresolved or unavailable.
#' @details Each side's last three values must agree within
#'   atol+rtol*max(abs(a),abs(b)); both-side estimates must also agree.
#'   This finite Cauchy-like check is evidence only. Divergence, oscillations,
#'   floating-point cancellation and aliased sampling can evade finite tests.
#'   Unbounded or unresolved sequences are not replaced by a finite limit or
#'   an asserted infinity. Invalid values remain NA. Collapsed/underflowed or
#'   overflowing approach points are errors. Existence is never proved.
#' @examples
#' f <- function(x) sin(x)/x
#' f |> vc_limit(x0=0)
#' vc_limit(function(x) abs(x)/x,0)$status
#' @export
vc_limit <- function(f,x0,side=c("both","left","right"),h=NULL,steps=20L,
                     atol=1e-6,rtol=1e-4) {
  if (!is.function(f)) stop("'f' must be a function.")
  if (!.vc_scalar(x0)) stop("'x0' must be a finite real scalar.")
  side <- match.arg(side)
  .vc_integer(steps,"steps",4)
  if (steps>1000L) stop("'steps' must be <=1000.")
  .vc_positive(atol,"atol",TRUE)
  .vc_positive(rtol,"rtol",TRUE)
  if (is.null(h)) h <- .1*max(1,abs(x0))
  .vc_positive(h,"h")
  distances <- h/2^(seq_len(steps)-1L)
  if (any(distances==0)) stop("Approach distance underflow; reduce 'steps'.")
  close <- function(a,b) {
    scale <- max(abs(a),abs(b),atol,.Machine$double.xmin)
    abs(a/scale-b/scale)<=atol/scale+rtol*max(abs(a/scale),abs(b/scale))
  }
  samples <- data.frame(x=numeric(),distance=numeric(),value=numeric(),side=character())
  summaries <- list(left=NULL,right=NULL)
  for (s in if(side=="both") c("left","right") else side) {
    xs <- x0+if(s=="left") -distances else distances
    if (any(!is.finite(xs)) || any(xs==x0)) stop("Approach points collapsed or overflowed; change 'h' or 'steps'.")
    values <- vapply(xs,function(x) .vc_safe_value(f,x),numeric(1))
    samples <- rbind(samples,data.frame(x=xs,distance=abs(xs-x0),value=values,side=s))
    last <- values[(steps-2L):steps]
    available <- all(is.finite(last))
    stable <- available && close(last[1],last[2]) && close(last[2],last[3])
    summaries[[s]] <- list(estimate=if(stable) last[3] else NA_real_,stable=stable,available=available)
  }
  active <- summaries[if(side=="both") c("left","right") else side]
  available <- all(vapply(active,function(a)a$available,logical(1)))
  stable <- all(vapply(active,function(a)a$stable,logical(1)))
  status <- if(!available) "unavailable" else if(!stable) "unresolved" else "compatible_finite"
  if (side=="both" && stable && !close(summaries$left$estimate,summaries$right$estimate)) status <- "side_disagreement"
  estimate <- if (status=="compatible_finite") {
    if(side=="both") summaries$left$estimate/2+summaries$right$estimate/2 else active[[1]]$estimate
  } else NA_real_
  list(type="limit",x0=x0,side=side,samples=samples,left=summaries$left,right=summaries$right,
       estimate=estimate,status=status,note="Finite probes do not prove the existence or value of a limit.")
}

#' Plot finite limit probes without evaluating the approach point
#' @param result Result from vc_limit().
#' @inheritParams vc_plot_secant
#' @param colors Named vector with curve, left, right and limit colors.
#' @return Plotly widget with sampled function, side probes and a distance
#'   panel. A finite compatible estimate is an open marker, not f(x0).
#' @details The approach point is excluded from all function evaluations.
#'   The lower panel uses logarithmic distance with smaller distances to the
#'   right. No proof of limit existence or infinity is implied.
#' @examples
#' if (requireNamespace("plotly",quietly=TRUE)) {
#'   f <- function(x) sin(x)/x
#'   f |> vc_limit(0) |> vc_plot_limit(f=f)
#' }
#' @export
vc_plot_limit <- function(result,f,xlim=NULL,n=401L,
                          colors=c(curve="#24476b",left="#188367",right="#c44736",limit="#8250a3"),
                          title=NULL,xlab="x",ylab="f(x)") {
  if (!requireNamespace("plotly",quietly=TRUE)) stop("Package 'plotly' is required only for visualization.")
  if (!is.list(result) || !identical(result$type,"limit")) stop("'result' must come from vc_limit().")
  if (!is.function(f)) stop("'f' must be a function.")
  .vc_integer(n,"n",2)
  if (!is.character(colors) || !all(c("curve","left","right","limit") %in% names(colors)) ||
      anyNA(colors) || anyDuplicated(names(colors)) || any(!nzchar(colors))) stop("Invalid 'colors'.")
  label <- function(x) is.character(x) && length(x)==1L && !is.na(x)
  if ((!is.null(title) && !label(title)) || !label(xlab) || !label(ylab)) stop("Invalid labels.")
  d <- result$samples
  actual <- vapply(d$x,function(x) .vc_safe_value(f,x),numeric(1))
  if (!isTRUE(all.equal(actual,d$value,tolerance=0))) stop("'f' does not match stored probes.")
  if (is.null(xlim)) xlim <- result$x0+c(-1,1)*max(d$distance)*1.2
  xlim <- .vc_interval(xlim)
  if(result$x0<xlim[1] || result$x0>xlim[2]) stop("'xlim' must contain the approach point.")
  xs <- sort(unique(c(seq(xlim[1],xlim[2],length.out=n),result$x0)))
  ys <- vapply(xs,function(x) if(x==result$x0) NA_real_ else .vc_safe_value(f,x),numeric(1))
  if(!any(is.finite(ys)) && !any(is.finite(d$value))) stop("No finite function samples or probes.")
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p,x=xs,y=ys,type="scatter",mode="lines",name="Sampled function (x0 excluded)",
    connectgaps=FALSE,line=list(color=colors[["curve"]]))
  e <- plotly::plot_ly()
  for(s in unique(d$side)) {
    z <- d[d$side==s,,drop=FALSE]
    p <- plotly::add_trace(p,x=z$x,y=z$value,type="scatter",mode="markers",name=paste(s,"probes"),
      marker=list(color=colors[[s]],size=6))
    e <- plotly::add_trace(e,x=z$distance,y=z$value,type="scatter",mode="lines+markers",name=paste(s,"approach"),
      connectgaps=FALSE,line=list(color=colors[[s]]))
  }
  if (is.finite(result$estimate)) p <- plotly::add_trace(p,x=result$x0,y=result$estimate,
    type="scatter",mode="markers",name="Compatible limit candidate (not f(x0))",
    marker=list(color=colors[["limit"]],size=12,symbol="circle-open"))
  p <- plotly::subplot(p,e,nrows=2,shareX=FALSE,titleY=TRUE,heights=c(.65,.35),margin=.08)
  if (is.null(title)) title <- paste("Finite limit probes:",result$status)
  p <- plotly::layout(p,title=list(text=paste0(title,"<br><sup>Finite probes are not a proof</sup>")),
    xaxis=list(title=xlab,range=xlim),xaxis2=list(title="Distance to x0 (decreasing)",type="log",autorange="reversed"),
    yaxis=list(title=ylab),yaxis2=list(title="Probe values"),margin=list(t=100,b=80),
    legend=list(orientation="h",y=1.04,yanchor="bottom"))
  p$height <- 800
  plotly::config(p,displaylogo=FALSE,scrollZoom=TRUE)
}
