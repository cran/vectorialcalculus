#' Construct a rectangular product partition
#' @param xlim,ylim Finite increasing coordinate bounds; NULL with custom breaks.
#' @param nx,ny Positive regular cell counts, ignored with corresponding breaks.
#' @param xbreaks,ybreaks Optional finite strictly increasing coordinate breaks.
#' @return Plain list with two validated vc_partition results and cell table
#'   (ix,iy,x0,x1,y0,y1,dx,dy,xmid,ymid,area). x varies fastest.
#' @details Reuses vc_partition unchanged. Mesh is the maximum cell diameter.
#' @examples
#' vc_partition2d(c(0,1),c(0,2),nx=2,ny=3)
#' @export
vc_partition2d <- function(xlim=NULL,ylim=NULL,nx=8L,ny=8L,xbreaks=NULL,ybreaks=NULL) {
  x <- vc_partition(xlim,nx,xbreaks);y <- vc_partition(ylim,ny,ybreaks)
  if(as.double(x$n)*y$n>.Machine$integer.max)stop("Too many product cells.")
  d <- expand.grid(ix=seq_len(x$n),iy=seq_len(y$n))
  a <- x$cells[d$ix,];b <- y$cells[d$iy,]
  d$x0 <- a$left;d$x1 <- a$right;d$y0 <- b$left;d$y1 <- b$right
  d$dx <- a$width;d$dy <- b$width;d$xmid <- a$midpoint;d$ymid <- b$midpoint
  d$area <- d$dx*d$dy
  if(any(!is.finite(d$area)|d$area<=0))stop("Cell area overflowed or underflowed.")
  # Scaled hypotenuse avoids squaring large widths.
  s <- pmax(d$dx,d$dy)
  mesh <- max(s*sqrt((d$dx/s)^2+(d$dy/s)^2))
  .vc_finite(mesh,"Product mesh")
  list(type="partition2d",x=x,y=y,cells=d,mesh=mesh)
}

#' Compute signed double Riemann sums and sampled cell statistics
#' @param f Scalar field accepting one x,y vector.
#' @inheritParams vc_partition2d
#' @param partition Optional vc_partition2d result, reconstructed from breaks.
#' @param method Midpoint, lower-left or upper-right tag; or finite sample
#'   minimum, maximum or mean. Sample extrema are NOT certified Darboux bounds.
#' @param sample_n Positive subdivisions per coordinate for sample statistics;
#'   one samples the four corners. Ignored for single-tag rules.
#' @param lower,upper Optional pair of scalar y-boundary functions of x. The
#'   supplied bounding box must contain their sampled bounds, in order.
#' @return Plain list with estimate, partition, method, table including value,
#'   sample_x/sample_y, contribution, inside_fraction, and samples (per-cell
#'   tables of actual x,y,inside,value evaluations). Mean tags have NA coordinates.
#' @details For nonrectangular regions computes the Riemann sum of f times
#'   the sampled domain indicator on the supplied box: outside samples give
#'   zero without evaluating f. Full rectangular cell areas remain in use;
#'   boundary cells are not exact clipped pieces. inside_fraction counts samples,
#'   not the fraction of geometric area. Membership/ordering are
#'   checked only at sampled x. Caller supplies a valid bounding box and
#'   justifies integrability of the zero extension and boundary behavior.
#'   Negative heights/contributions are signed, not geometric volumes.
#'   Finite statistics do not establish extrema, integrability or convergence.
#' @examples
#' (function(v)sum(v)) |> vc_riemann2d(c(0,1),c(0,1),nx=4,ny=4)
#' @export
vc_riemann2d <- function(f,xlim=NULL,ylim=NULL,nx=8L,ny=8L,
    method=c("midpoint","lower_left","upper_right","sample_min","sample_max","sample_mean"),
    partition=NULL,sample_n=1L,lower=NULL,upper=NULL) {
  if(!is.function(f))stop("'f' must be a scalar field.")
  method <- match.arg(method)
  .vc_integer(sample_n,"sample_n",1)
  if(sample_n>=.Machine$integer.max)stop("Sample resolution too large.")
  if(xor(is.null(lower),is.null(upper))||(!is.null(lower)&&
    (!is.function(lower)||!is.function(upper))))stop("Supply both boundary functions or neither.")
  if(is.null(partition))p <- vc_partition2d(xlim,ylim,nx,ny) else {
    if(!is.list(partition)||!identical(partition$type,"partition2d"))stop("Invalid product partition.")
    p <- vc_partition2d(xlim,ylim,xbreaks=partition$x$breaks,ybreaks=partition$y$breaks)
  }
  d <- p$cells
  samples <- lapply(seq_len(nrow(d)),function(i) {
    z <- d[i,]
    q <- if(method=="midpoint")data.frame(x=z$xmid,y=z$ymid) else
      if(method=="lower_left")data.frame(x=z$x0,y=z$y0) else
      if(method=="upper_right")data.frame(x=z$x1,y=z$y1) else {
        if((as.double(sample_n)+1)^2>.Machine$integer.max)stop("Too many cell samples.")
        expand.grid(x=seq(z$x0,z$x1,length.out=sample_n+1),
                    y=seq(z$y0,z$y1,length.out=sample_n+1))
      }
    q$inside <- TRUE
    if(!is.null(lower)) {
      a <- vapply(q$x,function(x).vc_value(lower,x),numeric(1))
      b <- vapply(q$x,function(x).vc_value(upper,x),numeric(1))
      if(any(a>b))stop("Sampled lower boundary exceeds upper boundary.")
      if(any(a<p$y$interval[1]|b>p$y$interval[2]))stop("Sampled boundaries leave the bounding box.")
      q$inside <- q$y>=a & q$y<=b
    }
    q$value <- vapply(seq_len(nrow(q)),function(j)
      if(q$inside[j]).vc_field_value(f,c(q$x[j],q$y[j])) else 0,numeric(1))
    q
  })
  d$value <- d$sample_x <- d$sample_y <- d$inside_fraction <- numeric(nrow(d))
  for(i in seq_len(nrow(d))) {
    q <- samples[[i]]
    k <- switch(method,sample_min=which.min(q$value),sample_max=which.max(q$value),1L)
    d$value[i] <- if(method=="sample_mean").vc_finite(mean(q$value),"Sample mean") else q$value[k]
    d$sample_x[i] <- if(method=="sample_mean")NA_real_ else q$x[k]
    d$sample_y[i] <- if(method=="sample_mean")NA_real_ else q$y[k]
    d$inside_fraction[i] <- mean(q$inside)
  }
  d$contribution <- d$value*d$area
  if(any(!is.finite(d$contribution)))stop("Nonfinite prism contribution.")
  list(type="riemann2d",estimate=.vc_finite(sum(d$contribution),"Double Riemann sum"),
    partition=p,method=method,table=d,samples=samples,
    region=if(is.null(lower))"rectangle" else "sampled_indicator_on_box",
    note="Signed sampled sum; no certified Darboux bound, clipped geometry or convergence proof.")
}

#' Tabulate double Riemann refinement
#' @inheritParams vc_riemann2d
#' @param n_values Strictly increasing positive cell counts per coordinate.
#' @param reference Optional finite analytical value.
#' @return List with complete results and table of n, cell count, mesh,
#'   estimate, change, error and abs_error. No limit is inferred from the table.
#' @examples
#' vc_riemann2d_convergence(function(v)sum(v),c(0,1),c(0,1),reference=1)
#' @export
vc_riemann2d_convergence <- function(f,xlim,ylim,n_values=c(2L,4L,8L,16L),
    method="midpoint",sample_n=1L,lower=NULL,upper=NULL,reference=NULL) {
  if(!is.numeric(n_values)||is.complex(n_values)||!length(n_values)||
     any(!is.finite(n_values))||any(n_values<1|n_values!=floor(n_values))||any(diff(n_values)<=0))
    stop("Invalid refinement counts.")
  if(!is.null(reference)&&!.vc_scalar(reference))stop("Invalid reference.")
  results <- lapply(n_values,function(n)vc_riemann2d(f,xlim,ylim,n,n,method,
    sample_n=sample_n,lower=lower,upper=upper))
  e <- vapply(results,function(z)z$estimate,numeric(1))
  change <- c(NA_real_,diff(e))
  error <- if(is.null(reference))rep(NA_real_,length(e)) else e-reference
  if(any(!is.finite(change[-1]))||any(!is.finite(error[!is.na(error)])))stop("Refinement differences overflowed.")
  list(type="riemann2d_convergence",results=results,table=data.frame(n=n_values,
    cells=n_values^2,mesh=vapply(results,function(z)z$partition$mesh,numeric(1)),
    estimate=e,change=change,error=error,abs_error=abs(error)),reference=reference,
    note="Finite refinement is not a proof of Riemann integrability or convergence.")
}

#' Draw stored signed Riemann prisms independently of calculation
#' @param result Result of vc_riemann2d().
#' @param title Optional character title.
#' @return Plotly mesh widget. Every rectangular cell extends from zero to
#'   its stored signed height. No new field or boundary evaluation is made.
#' @details For restricted regions these are prisms of the zero-extension
#'   approximation on the box, not exact clipped boundary cells. Sampled extrema
#'   are labeled as such. A prism below zero contributes negatively.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   vc_riemann2d(function(v)v[1]-v[2],c(0,1),c(0,1),nx=2,ny=2) |>
#'     vc_plot_riemann2d()
#' }
#' @export
vc_plot_riemann2d <- function(result,title=NULL) {
  if(!requireNamespace("plotly",quietly=TRUE))stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"riemann2d"))stop("Invalid double Riemann result.")
  if(!is.null(title)&&(!is.character(title)||length(title)!=1||is.na(title)))stop("Invalid title.")
  d <- result$table
  p <- plotly::plot_ly()
  for(k in seq_len(nrow(d))) {
    z <- d[k,]
    x <- rep(c(z$x0,z$x1,z$x1,z$x0),2)
    y <- rep(c(z$y0,z$y0,z$y1,z$y1),2)
    heights <- c(rep(0,4),rep(z$value,4))
    p <- plotly::add_trace(p,x=x,y=y,z=heights,type="mesh3d",inherit=FALSE,
      i=c(0,0,4,4,0,0,1,1,2,2,3,3),j=c(1,2,5,6,1,5,2,6,3,7,0,4),
      k=c(2,3,6,7,5,4,6,5,7,6,4,7),
      color=if(z$value<0)"#c44736" else "#188367",opacity=.55,showlegend=FALSE,
      name=paste("Cell",k),text=paste(result$method,"signed contribution",format(z$contribution)))
  }
  plotly::layout(p,title=if(is.null(title))paste("Signed prisms:",result$method) else title,
    scene=list(xaxis=list(title="x"),yaxis=list(title="y"),zaxis=list(title="f sample"),aspectmode="data"))
}
