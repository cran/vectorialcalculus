.vc_increment <- function(dx,point) {
  dx <- .vc_point(dx,"dx")
  if(length(dx)!=length(point)) stop("Increment and point dimensions differ.")
  if(!is.null(names(dx))&&!is.null(names(point))&&!identical(names(dx),names(point)))
    stop("Named increments must follow the point coordinate order.")
  dx
}

.vc_increment_rows <- function(dx,point) {
  if(is.null(dim(dx))) return(list(.vc_increment(dx,point)))
  if(!is.matrix(dx)||!is.numeric(dx)||is.complex(dx)||ncol(dx)!=length(point)||
     !nrow(dx)||any(!is.finite(dx))) stop("'dx' must be an increment vector or finite row-per-increment matrix.")
  lapply(seq_len(nrow(dx)),function(i).vc_increment(dx[i,],point))
}

.vc_model_norm <- function(v) {
  scale <- max(abs(v))
  if(scale==0) return(0)
  .vc_finite(scale*sqrt(sum((v/scale)^2)),"Displacement norm")
}

#' Construct a numerical total differential and affine local model
#' @inheritParams vc_gradient
#' @param dx Optional increment vector, or matrix with one increment per row.
#'   Coordinates are positional; when both point and increments are named,
#'   names must agree in order. NULL constructs evaluators without probes.
#' @param scales Strictly decreasing positive finite probe multipliers;
#'   default compares each increment at full, half and quarter scale.
#' @param gradient Optional vc_gradient result at the same point and field
#'   value. When supplied, its diagnostics are retained and method, h, atol,
#'   rtol and warn are ignored; the gradient is not recomputed.
#' @return Plain list with type, point, value, gradient (full result),
#'   gradient_reliable, differential (function of an increment), model
#'   (function of an absolute point), comparison (table or NULL), probes
#'   (requested/actual increments, endpoints), scales, normal, unit_normal,
#'   normal_status, plane_coefficients, coefficient_status and note. For the
#'   candidate graph hyperplane, normal is (-gradient,1); plane_coefficients
#'   appends the offset to that vector. An unrepresentable offset is NA and
#'   failed unit normalization returns NULL, with explicit geometry statuses;
#'   the centered model remains usable. Table columns
#'   are probe, scale, distance, requested_differential, differential,
#'   approximation, value, change, remainder, normalized_remainder,
#'   displacement_error, status and message. Unavailable diagnostics are NA.
#' @details Reuses vc_gradient without modifying the accepted derivative engine.
#'   The candidate differential is `sum(gradient*dx)`; the model is
#'   `f(point)+sum(gradient*(new_point-point))`. It is the mathematical total
#'   differential/tangent model only under appropriate total differentiability.
#'   Coordinate sensitivity checks do not establish this hypothesis.
#'   Probes compare the model with actual finite changes at point+scaled_dx.
#'   Comparison differential uses the representable displacement; the requested
#'   functional value is also retained. Absolute remainder divided by actual
#'   displacement norm is reported, not promoted to a differentiability test.
#'   Even decreasing ratios on finitely many directions are not a proof.
#'   Gradient approximation error and roundoff can produce a residual floor.
#'   Zero increments have no normalized remainder. Domain errors, collapsed
#'   active coordinates, scaling underflow/overflow and nonfinite diagnostics
#'   remain explicit probe statuses; they do not become apparently valid zeros.
#'   Callable evaluators reject invalid dimensions/nonfinite results. Probing
#'   may evaluate outside the coordinate stencils, so it can fail even when
#'   gradient diagnostics pass. A supplied gradient is checked at the point,
#'   not proven to belong to the same field everywhere.
#' @examples
#' f <- function(v)v[1]^2+2*v[2]^2
#' a <- f |> vc_linearization(c(1,1),dx=c(.1,.2))
#' a$differential(c(.1,.2))
#' a$model(c(1.1,1.2))
#' a$comparison
#' @export
vc_linearization <- function(f,point,dx=NULL,scales=c(1,.5,.25),gradient=NULL,
                             method="central",h=NULL,atol=1e-8,rtol=1e-4,warn=TRUE) {
  point <- .vc_point(point)
  value <- .vc_field_value(f,point)
  if(!is.numeric(scales)||is.complex(scales)||!is.null(dim(scales))||!length(scales)||
     any(!is.finite(scales))||any(scales<=0)||any(diff(scales)>=0)) stop("'scales' must be strictly decreasing positive finite values.")
  if(is.null(gradient)) gradient <- vc_gradient(f,point,method,h,atol,rtol,warn) else {
    if(!is.list(gradient)||!identical(gradient$type,"gradient")||
       !identical(as.double(gradient$point),as.double(point))||!identical(gradient$value,value)||
       !is.numeric(gradient$estimate)||is.complex(gradient$estimate)||
       length(gradient$estimate)!=length(point)||any(!is.finite(gradient$estimate))||
       !is.logical(gradient$reliable)||length(gradient$reliable)!=1||is.na(gradient$reliable))
      stop("'gradient' must be a valid result at the same point and field value.")
  }
  slopes <- as.double(gradient$estimate)
  differential <- function(dx) {
    dx <- .vc_increment(dx,point)
    .vc_finite(sum(slopes*dx),"Linear functional value")
  }
  model <- function(new_point) {
    new_point <- .vc_increment(new_point,point)
    delta <- new_point-point
    if(any(!is.finite(delta))) stop("Model displacement overflowed.")
    .vc_finite(value+differential(delta),"Affine model value")
  }
  probes <- list()
  rows <- list()
  if(!is.null(dx)) {
    increments <- .vc_increment_rows(dx,point)
    for(i in seq_along(increments)) for(s in scales) {
      requested <- increments[[i]]*s
      endpoint <- actual <- rep(NA_real_,length(point))
      row <- data.frame(probe=i,scale=s,distance=NA_real_,requested_differential=NA_real_,
        differential=NA_real_,approximation=NA_real_,value=NA_real_,change=NA_real_,
        remainder=NA_real_,normalized_remainder=NA_real_,displacement_error=NA_real_,
        status="not_evaluated",message="")
      evaluate <- function() {
        if(any(!is.finite(requested))) return("step_overflow")
        if(any(increments[[i]]!=0 & requested==0)) return("step_underflow")
        row$requested_differential <<- differential(requested)
        endpoint <<- point+requested
        if(any(!is.finite(endpoint))) return("endpoint_not_finite")
        actual <<- endpoint-point
        if(any(!is.finite(actual))) return("displacement_not_finite")
        if(any(requested!=0 & actual==0)) return("collapsed_increment")
        row$distance <<- .vc_model_norm(actual)
        row$displacement_error <<- .vc_model_norm(actual-requested)
        row$differential <<- differential(actual)
        row$approximation <<- model(endpoint)
        sampled <- tryCatch(suppressWarnings(.vc_field_value(f,endpoint)),error=function(e)e)
        if(inherits(sampled,"error")) {
          row$message <<- conditionMessage(sampled)
          return("value_not_finite")
        }
        row$value <<- sampled
        row$change <<- .vc_finite(sampled-value,"Actual finite change")
        row$remainder <<- .vc_finite(row$change-row$differential,"Linear remainder")
        if(row$distance==0) return("zero_increment")
        row$normalized_remainder <<- .vc_finite(abs(row$remainder)/row$distance,"Normalized remainder")
        "ok"
      }
      row$status <- tryCatch(evaluate(),error=function(e) {
        row$message <<- conditionMessage(e)
        "numerical_error"
      })
      rows[[length(rows)+1L]] <- row
      probes[[length(probes)+1L]] <- list(requested_dx=requested,actual_dx=actual,endpoint=endpoint)
    }
  }
  normal <- c(-slopes,1)
  normal_status <- "ok"
  unit_normal <- tryCatch(.vc_unit(normal),error=function(e) {
    normal_status <<- conditionMessage(e)
    NULL
  })
  offset <- sum(slopes*point)-value
  coefficient_status <- if(is.finite(offset))"ok" else "offset_not_finite"
  if(!is.finite(offset)) offset <- NA_real_
  list(type="linearization",point=point,value=value,gradient=gradient,
    gradient_reliable=gradient$reliable,differential=differential,model=model,
    comparison=if(length(rows))do.call(rbind,rows) else NULL,probes=probes,scales=as.double(scales),
    normal=normal,unit_normal=unit_normal,normal_status=normal_status,
    plane_coefficients=c(normal,offset=offset),coefficient_status=coefficient_status,
    note="Candidate affine model only: coordinate checks and finite remainder probes do not prove total differentiability.")
}
