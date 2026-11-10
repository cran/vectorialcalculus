.vc_positive <- function(x, name, zero = FALSE) {
  if (!.vc_scalar(x) || if (zero) x < 0 else x <= 0)
    stop("'",name,"' must be ",if (zero) "nonnegative" else "positive", " and finite.")
}

.vc_integer <- function(x, name, minimum) {
  if (!.vc_scalar(x) || x < minimum || x != floor(x) || x > .Machine$integer.max)
    stop("'",name,"' must be an integer >= ",minimum,".")
}

.vc_interval <- function(interval) {
  if (!is.numeric(interval) || is.complex(interval) || length(interval)!=2L ||
      any(!is.finite(interval))) stop("'interval' must be finite increasing bounds with finite width.")
  interval <- as.double(interval)
  if(interval[1]>=interval[2] || !is.finite(diff(interval)))
    stop("'interval' must be finite increasing bounds with finite width.")
  interval
}

.vc_optional_step <- function(h,name) if (!is.null(h)) .vc_positive(h,name)

.vc_safe_value <- function(f,x) {
  tryCatch(suppressWarnings(.vc_value(f,x)),error=function(e) NA_real_)
}

.vc_first <- function(f,x,h=NULL,atol=1e-8,rtol=1e-4) {
  tryCatch(suppressWarnings(vc_derivative(f,x,h=h,atol=atol,rtol=rtol,warn=FALSE)),error=function(e) NULL)
}

.vc_second <- function(f,x,h2=NULL,atol=1e-8,rtol=1e-4) {
  if (is.null(h2)) h2 <- max(1,abs(x))*.Machine$double.eps^(1/4)
  tryCatch(suppressWarnings(.vc_second_difference(f,x,h2,atol,rtol)),error=function(e) NULL)
}

.vc_continuity_evidence <- function(f,x,radius) {
  radii <- radius/c(1,2,4)
  y <- .vc_safe_value(f,x)
  left <- vapply(x-radii,function(z) .vc_safe_value(f,z),numeric(1))
  right <- vapply(x+radii,function(z) .vc_safe_value(f,z),numeric(1))
  a <- x-(x-radii)
  b <- (x+radii)-x
  gap <- pmax(abs(left-y),abs(right-y))
  usable <- is.finite(y) && all(is.finite(c(left,right,a,b,gap))) && all(a>0 & b>0)
  continuity <- usable && (gap[3] <= .8*gap[1] || all(gap==0))
  data <- data.frame(radius=radii,left=left,right=right,
                     backward=(y-left)/a,forward=(right-y)/b)
  list(data=data,continuity=continuity)
}

.vc_root_scan <- function(evaluate,grid,residual_tol,x_tol,max_iter) {
  values <- vapply(grid,function(x) {
    v <- tryCatch(evaluate(x),error=function(e) NA_real_)
    if (.vc_scalar(v)) v else NA_real_
  },numeric(1))
  roots <- numeric()
  events <- list()
  record <- function(method,bounds,x=NA_real_,residual=NA_real_,status,message="") {
    events[[length(events)+1L]] <<- data.frame(method=method,lower=bounds[1],upper=bounds[2],
      x=x,residual=residual,status=status,message=message)
  }
  add <- function(x,method,bounds) {
    residual <- tryCatch(evaluate(x),error=function(e) NA_real_)
    finite <- .vc_scalar(residual)
    accepted <- finite && abs(residual)<=residual_tol
    record(method,bounds,x,if(finite) abs(residual) else NA_real_,
           if(accepted) "accepted" else if(finite) "residual_rejected" else "nonfinite_residual")
    if (accepted) roots <<- c(roots,x)
  }
  for(i in seq_len(length(grid)-1L)) {
    if (all(is.finite(values[c(i,i+1L)])) && sign(values[i])!=sign(values[i+1L])) {
      bounds <- grid[c(i,i+1L)]
      z <- tryCatch(stats::uniroot(evaluate,bounds,tol=x_tol,
                                   maxiter=max_iter,check.conv=TRUE)$root,error=function(e) e)
      if(inherits(z,"error")) record("bracket",bounds,status="refinement_error",message=conditionMessage(z))
      else if(is.finite(z)) add(z,"bracket",bounds)
    }
  }
  for(i in 2:(length(grid)-1L)) {
    v <- abs(values[c(i-1L,i,i+1L)])
    if (!all(is.finite(v))) next
    if (v[2]<=v[1] && v[2]<=v[3] && (v[2]<v[1] || v[2]<v[3])) {
      bounds <- grid[c(i-1L,i+1L)]
      if (v[2]==0) add(grid[i],"grid_minimum",rep(grid[i],2))
      z <- tryCatch(stats::optimize(function(x) {
        q <- evaluate(x)
        if (.vc_scalar(q)) abs(q) else .Machine$double.xmax
      },bounds,tol=x_tol)$minimum,error=function(e) e)
      if(inherits(z,"error")) record("local_absolute_minimum",bounds,status="refinement_error",message=conditionMessage(z))
      else if(is.finite(z)) add(z,"local_absolute_minimum",bounds)
    }
  }
  if (length(roots)) {
    roots <- sort(unique(roots))
    groups <- cumsum(c(TRUE,diff(roots)>10*x_tol))
    roots <- vapply(split(roots,groups),function(xs) {
      residual <- vapply(xs,function(x) abs(evaluate(x)),numeric(1))
      xs[which.min(residual)]
    },numeric(1))
  }
  refinements <- if(length(events)) do.call(rbind,events) else data.frame(method=character(),
    lower=numeric(),upper=numeric(),x=numeric(),residual=numeric(),status=character(),message=character())
  list(roots=unname(roots),values=values,refinements=refinements)
}

.vc_runs <- function(x,labels) {
  runs <- rle(labels)
  ends <- cumsum(runs$lengths)
  starts <- c(1L,ends[-length(ends)]+1L)
  data.frame(lower=x[starts],upper=x[ends],classification=runs$values,samples=runs$lengths)
}
