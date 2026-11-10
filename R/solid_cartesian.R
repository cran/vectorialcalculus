#' Sample the six faces of a nested Cartesian solid
#' @param lower,upper Scalar y-boundary functions of x, ordered at sampled nodes.
#' @param bottom,top Scalar z-boundary functions accepting one x,y vector.
#' @param interval Finite increasing x bounds, or NULL with partition supplied.
#' @param n,n_u,n_v Positive subdivisions in x, across y and across z.
#' @param partition Optional vc_partition result for custom x strips.
#' @return Plain list with type solid_cartesian, base (full vc_region result),
#'   bottom/top matrices, v fractions and faces: x_min, x_max, y_lower,
#'   y_upper, z_lower, z_upper. Each face has finite x,y,z matrices with
#'   matching dimensions. No area, volume, normals or quadrature computed.
#' @details Reuses vc_region unchanged. Samples z bounds at every base grid
#'   node and retains exact shared endpoints on all face seams. Ordered bounds
#'   are checked at samples, not between nodes. Interpolated parametric faces
#'   may miss crossings or folds; this is not a certified watertight manifold.
#'   Degenerate sections/heights are allowed and can collapse face triangles.
#'   Disconnected sections need justified decomposition. No unused slice
#'   arguments are copied from the original API; it did not implement slices.
#' @examples
#' vc_solid_cartesian(function(x)0,function(x)1-x,function(v)0,
#'   function(v)1-sum(v),c(0,1),n=4,n_u=4,n_v=4)
#' @export
vc_solid_cartesian <- function(lower,upper,bottom,top,interval=NULL,
                               n=20L,n_u=12L,n_v=12L,partition=NULL) {
  if(!is.function(bottom)||!is.function(top))stop("Both height boundaries must be functions.")
  .vc_integer(n_v,"n_v",1)
  if(n_v>=.Machine$integer.max)stop("Vertical blending resolution too large.")
  base <- vc_region(lower,upper,interval,n,n_u,partition)
  x <- base$grid$x;y <- base$grid$y
  evaluate <- function(g)matrix(vapply(seq_along(x),function(i)
    .vc_field_value(g,c(x[i],y[i])),numeric(1)),nrow=nrow(x))
  a <- evaluate(bottom);b <- evaluate(top)
  if(any(a>b))stop("Sampled bottom exceeds top.")
  if(any(!is.finite(b-a)))stop("Vertical width overflowed.")
  v <- seq(0,1,length.out=n_v+1)
  side <- function(xs,ys,lo,hi) {
    if(as.double(length(v))*length(xs)>.Machine$integer.max)stop("Face grid too large.")
    z <- vapply(seq_along(xs),function(i)lo[i]+v*(hi[i]-lo[i]),numeric(length(v)))
    if(any(!is.finite(z)))stop("Face coordinate overflowed.")
    z[1,] <- lo;z[length(v),] <- hi
    list(x=matrix(rep(xs,each=length(v)),nrow=length(v)),
         y=matrix(rep(ys,each=length(v)),nrow=length(v)),z=z)
  }
  left <- 1L;right <- ncol(x);low <- 1L;high <- nrow(x)
  faces <- list(x_min=side(x[,left],y[,left],a[,left],b[,left]),
    x_max=side(x[,right],y[,right],a[,right],b[,right]),
    y_lower=side(x[low,],y[low,],a[low,],b[low,]),
    y_upper=side(x[high,],y[high,],a[high,],b[high,]),
    z_lower=list(x=x,y=y,z=a),z_upper=list(x=x,y=y,z=b))
  list(type="solid_cartesian",base=base,bottom=a,top=b,v=v,faces=faces,
    note="Sampled parametric faces; shared seams do not certify topology, exact geometry or volume.")
}

#' Display stored Cartesian solid faces or wireframes
#' @param result Result of vc_solid_cartesian().
#' @param mode Faces, wireframe or both.
#' @param faces Nonempty unique character subset of stored face names.
#' @param opacity Surface opacity between zero and one.
#' @param wire_step Positive integer stride for wire lines, retaining endpoints.
#' @param title Optional character title.
#' @return Plotly widget consuming stored face data only. No boundary evaluation,
#'   integral calculation or automatic internal slices are performed.
#' @details Wire rows/columns use the stored parameter grids. Faces can be
#'   selected independently. Interpolation is a display approximation; normal
#'   orientation, manifold structure and exact closure are not certified.
#' @examples
#' if(requireNamespace("plotly",quietly=TRUE)) {
#'   s <- vc_solid_cartesian(function(x)0,function(x)1,function(v)0,
#'     function(v)1,c(0,1),n=2,n_u=2,n_v=2)
#'   vc_plot_solid_cartesian(s,mode="both")
#' }
#' @export
vc_plot_solid_cartesian <- function(result,mode=c("faces","wireframe","both"),
    faces=names(result$faces),opacity=.4,wire_step=2L,title=NULL) {
  if(!requireNamespace("plotly",quietly=TRUE))stop("Plotly is required only for graphics.")
  if(!is.list(result)||!identical(result$type,"solid_cartesian"))stop("Invalid Cartesian solid.")
  mode <- match.arg(mode);.vc_integer(wire_step,"wire_step",1)
  if(!.vc_scalar(opacity)||opacity<0||opacity>1)stop("Invalid opacity.")
  if(!is.character(faces)||!length(faces)||anyNA(faces)||anyDuplicated(faces)||
     any(!faces%in%names(result$faces)))stop("Invalid face selection.")
  if(!is.null(title)&&(!is.character(title)||length(title)!=1||is.na(title)))stop("Invalid title.")
  p <- plotly::plot_ly()
  for(name in faces) {
    f <- result$faces[[name]]
    if(mode!="wireframe")p <- plotly::add_trace(p,x=f$x,y=f$y,z=f$z,type="surface",
      inherit=FALSE,opacity=opacity,showscale=FALSE,name=name)
    if(mode!="faces") {
      indices <- function(n)unique(c(seq(1,n,by=wire_step),n))
      rows <- indices(nrow(f$x));cols <- indices(ncol(f$x))
      lines <- lapply(c("x","y","z"),function(axis) {
        m <- f[[axis]]
        c(unlist(lapply(rows,function(i)c(m[i,],NA_real_))),
          unlist(lapply(cols,function(j)c(m[,j],NA_real_))))
      })
      p <- plotly::add_trace(p,x=lines[[1]],y=lines[[2]],z=lines[[3]],
        type="scatter3d",mode="lines",inherit=FALSE,name=paste(name,"wire"),
        line=list(color="#24476b"),showlegend=FALSE)
    }
  }
  plotly::layout(p,title=if(is.null(title))"Cartesian solid: sampled faces" else title,
    scene=list(xaxis=list(title="x"),yaxis=list(title="y"),zaxis=list(title="z"),aspectmode="data"))
}
