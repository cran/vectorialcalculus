#' Cylindrical solid defined by radial and vertical bounds (with optional plot)
#'
#' @description
#' Builds and optionally plots, using \pkg{plotly}, a three-dimensional solid
#' described in cylindrical coordinates. The solid is defined by:
#' \itemize{
#'   \item an angular variable \code{theta} in the interval
#'         \code{[th_min, th_max]},
#'   \item a radial variable \code{r} between \code{R1(theta)} and
#'         \code{R2(theta)},
#'   \item and a vertical coordinate \code{z} between \code{Z1(r, theta)} and
#'         \code{Z2(r, theta)}.
#' }
#'
#' The surface is rendered by sampling a curvilinear grid in the parameters
#' \code{(theta, u, v)}, where \code{u} and \code{v} vary in \code{[0, 1]} and
#' are used as linear blending variables along the radial and vertical
#' directions, respectively.
#'
#' When volume computation is requested, the function numerically approximates
#' the triple integral of the form
#' \code{integral theta from th_min to th_max} of
#' \code{integral r from R1(theta) to R2(theta)} of
#' \code{integral z from Z1(r, theta) to Z2(r, theta)} of
#' \code{r dz dr dtheta}, which is the standard volume element in cylindrical
#' coordinates.
#'
#' @param R1,R2 Functions \code{function(theta)} giving the inner and outer
#'        radius bounds, respectively.
#' @param Z1,Z2 Functions \code{function(r, theta)} giving the lower and upper
#'        \code{z} bounds.
#' @param th_min,th_max Angular limits (numeric scalars) with
#'        \code{th_max > th_min}.
#' @param plot Logical. If \code{TRUE}, the solid is plotted using \pkg{plotly}.
#' @param n_theta,n_u,n_v Mesh resolution in \code{theta} (angle), \code{u}
#'        (radial blend) and \code{v} (vertical blend).
#' @param mode Character string, one of \code{"faces"}, \code{"wireframe"} or
#'        \code{"both"}, indicating whether to draw only the surface,
#'        only a wireframe or both.
#' @param colorscale Plotly colorscale for the surface. It can be a named scale,
#'        a single color, or a character vector of colors interpreted as a
#'        gradient.
#' @param opacity Surface opacity, a numeric value between 0 and 1.
#' @param show_surface_grid Logical. If \code{TRUE}, draws a grid over the
#'        surface.
#' @param surface_grid_color,surface_grid_width Color and line width used for
#'        the surface grid.
#' @param edge_line,wire_line Line style lists used for the edges and the
#'        wireframe lines when those are drawn.
#' @param scene,bg Plotly 3D scene configuration and background colors. The
#'        background list typically has entries \code{paper} and \code{plot}.
#' @param compute_volume Logical. If \code{TRUE}, the volume of the solid is
#'        approximated numerically.
#' @param vol_method Character string. Either \code{"adaptive"}, which uses
#'        nested \code{stats::integrate}, or \code{"grid"}, which uses the
#'        trapezoidal rule over a regular mesh.
#' @param ntheta_vol,nr_vol Mesh sizes in the angular and radial directions used
#'        when \code{vol_method = "grid"}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{theta_seq}, \code{u_seq}, \code{v_seq}: the parameter sequences
#'         used for sampling the surface,
#'   \item \code{fig}: a \pkg{plotly} object when \code{plot = TRUE}, or
#'         \code{NULL} otherwise,
#'   \item \code{volume}: \code{NULL} if \code{compute_volume = FALSE}, or a
#'         list containing the numeric volume estimate and metadata (method and
#'         grid parameters) when \code{compute_volume = TRUE}.
#' }
#'
#' @examples
#' \donttest{
#' # Example: a quarter-twisted cup
#' # R in [0, 1 + 0.2*cos(theta)],  z in [0, 1 + 0.5*r]
#' R1 <- function(theta) 0
#' R2 <- function(theta) 1 + 0.2*cos(theta)
#' Z1 <- function(r, theta) 0
#' Z2 <- function(r, theta) 1 + 0.5*r
#' solid_cylindrical3d(
#'   R1, R2, Z1, Z2, th_min = 0, th_max = pi/2,
#'   plot = FALSE, mode = "both",
#'   colorscale = c("white", "#2a9d8f"), opacity = 0.35,
#'   show_surface_grid = TRUE,
#'   compute_volume = TRUE, vol_method = "adaptive"
#' )$volume
#' }
#'
#' @export
#' @importFrom stats integrate
solid_cylindrical3d <- function(
    R1, R2, Z1, Z2,
    th_min, th_max,
    plot = TRUE,
    n_theta = 160, n_u = 70, n_v = 70,
    mode = c("faces","wireframe","both"),
    colorscale = "Blues",
    opacity    = 0.35,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    edge_line = list(color = "black", width = 2),
    wire_line = list(color = "rgba(0,0,0,0.35)", width = 1),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white"),
    compute_volume = FALSE,
    vol_method = c("adaptive","grid"),
    ntheta_vol = 400, nr_vol = 400
) {
  mode       <- match.arg(mode)
  vol_method <- match.arg(vol_method)

  if (!is.function(R1) || !is.function(R2)) {
    stop("'R1' and 'R2' must be functions of theta.", call. = FALSE)
  }
  if (!is.function(Z1) || !is.function(Z2)) {
    stop("'Z1' and 'Z2' must be functions of (r, theta).", call. = FALSE)
  }
  if (!is.numeric(th_min) || !is.numeric(th_max) || th_max <= th_min) {
    stop("'th_max' must be > 'th_min'.", call. = FALSE)
  }

  # Coerce mesh sizes to integers and validate
  n_theta    <- as.integer(n_theta)
  n_u        <- as.integer(n_u)
  n_v        <- as.integer(n_v)
  ntheta_vol <- as.integer(ntheta_vol)
  nr_vol     <- as.integer(nr_vol)

  for (nm in c("n_theta","n_u","n_v","ntheta_vol","nr_vol")) {
    val <- get(nm)
    if (!is.finite(val) || val < 2L) {
      stop(nm, " must be an integer >= 2.", call. = FALSE)
    }
  }

  # helpers (colorscale parsing)
  is_color <- function(x) {
    if (!is.character(x) || length(x) != 1L) return(FALSE)
    if (grepl("^rgba?\\(", x, ignore.case = TRUE)) return(TRUE)
    if (grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)) return(TRUE)
    ok <- tryCatch({
      grDevices::col2rgb(x)
      TRUE
    }, error = function(...) FALSE)
    ok
  }

  as_rgba <- function(col, alpha = NULL) {
    if (grepl("^rgba?\\(", col, ignore.case = TRUE)) return(col)
    rgb <- grDevices::col2rgb(col) / 255
    a   <- if (is.null(alpha)) 1 else max(0, min(1, alpha))
    sprintf("rgba(%g,%g,%g,%g)", 255 * rgb[1], 255 * rgb[2], 255 * rgb[3], a)
  }

  as_colorscale <- function(x, alpha = NULL) {
    if (is.list(x) && length(x) >= 2L && is.numeric(x[[1]][[1]])) return(x)
    if (is.character(x) && length(x) == 1L) {
      if (is_color(x)) {
        ccol <- as_rgba(x, alpha)
        return(list(list(0, ccol), list(1, ccol)))
      } else {
        return(x)
      }
    }
    if (is.character(x) && length(x) > 1L) {
      cols <- vapply(x, as_rgba, character(1), alpha = alpha)
      pos  <- seq(0, 1, length.out = length(cols))
      return(lapply(seq_along(cols), function(i) list(pos[i], cols[i])))
    }
    stop("Unrecognized 'colorscale' format.", call. = FALSE)
  }

  # parameter grids
  theta_seq <- seq(th_min, th_max, length.out = n_theta)
  u_seq     <- seq(0, 1, length.out = n_u)   # radius blend
  v_seq     <- seq(0, 1, length.out = n_v)   # z blend

  # 1) theta = const faces
  build_face_theta <- function(th) {
    rr_in  <- R1(th)
    rr_out <- R2(th)
    r_line <- (1 - u_seq) * rr_in + u_seq * rr_out
    R <- matrix(rep(r_line, times = n_v), nrow = n_u)
    Z <- matrix(NA_real_, n_u, n_v)
    for (j in seq_len(n_v)) {
      z_low  <- Z1(r_line, th)
      z_high <- Z2(r_line, th)
      Z[, j] <- (1 - v_seq[j]) * z_low + v_seq[j] * z_high
    }
    X <- R * cos(th)
    Y <- R * sin(th)
    list(X = X, Y = Y, Z = Z)
  }

  # 2) r = R1(theta), r = R2(theta)
  build_face_r <- function(which = c("inner","outer")) {
    which <- match.arg(which)
    Rfix  <- if (which == "inner") {
      vapply(theta_seq, R1, numeric(1))
    } else {
      vapply(theta_seq, R2, numeric(1))
    }
    R  <- matrix(rep(Rfix, each = n_v), nrow = n_v)
    TH <- matrix(rep(theta_seq, each = n_v), nrow = n_v)
    Z  <- matrix(NA_real_, n_v, n_theta)
    for (j in seq_len(n_theta)) {
      z_low  <- Z1(Rfix[j], theta_seq[j])
      z_high <- Z2(Rfix[j], theta_seq[j])
      Z[, j] <- (1 - v_seq) * z_low + v_seq * z_high
    }
    X <- R * cos(TH)
    Y <- R * sin(TH)
    list(X = X, Y = Y, Z = Z)
  }

  # 3) z = Z1(r,theta), z = Z2(r,theta)
  build_face_z <- function(which = c("low","high")) {
    which <- match.arg(which)
    TH    <- matrix(rep(theta_seq, each = n_u), nrow = n_u)
    rr_in  <- vapply(theta_seq, R1, numeric(1))
    rr_out <- vapply(theta_seq, R2, numeric(1))
    R <- matrix(NA_real_, n_u, n_theta)
    for (j in seq_len(n_theta)) {
      R[, j] <- (1 - u_seq) * rr_in[j] + u_seq * rr_out[j]
    }
    Z <- matrix(NA_real_, n_u, n_theta)
    for (j in seq_len(n_theta)) {
      z_low  <- Z1(R[, j], theta_seq[j])
      z_high <- Z2(R[, j], theta_seq[j])
      Z[, j] <- if (which == "low") z_low else z_high
    }
    X <- R * cos(TH)
    Y <- R * sin(TH)
    list(X = X, Y = Y, Z = Z)
  }

  add_edges <- function(fig, F, line = edge_line) {
    X <- F$X; Y <- F$Y; Z <- F$Z
    i1 <- 1L; i2 <- nrow(X); j1 <- 1L; j2 <- ncol(X)
    fig |>
      plotly::add_trace(x = X[i1,], y = Y[i1,], z = Z[i1,],
                        type = "scatter3d", mode = "lines",
                        line = line, hoverinfo = "none", showlegend = FALSE) |>
      plotly::add_trace(x = X[i2,], y = Y[i2,], z = Z[i2,],
                        type = "scatter3d", mode = "lines",
                        line = line, hoverinfo = "none", showlegend = FALSE) |>
      plotly::add_trace(x = X[, j1], y = Y[, j1], z = Z[, j1],
                        type = "scatter3d", mode = "lines",
                        line = line, hoverinfo = "none", showlegend = FALSE) |>
      plotly::add_trace(x = X[, j2], y = Y[, j2], z = Z[, j2],
                        type = "scatter3d", mode = "lines",
                        line = line, hoverinfo = "none", showlegend = FALSE)
  }

  add_wireframe <- function(fig, F, line = wire_line, step = 6L) {
    X <- F$X; Y <- F$Y; Z <- F$Z
    nr <- nrow(X); nc <- ncol(X)
    r_idx <- seq(1L, nr, by = max(1L, step))
    c_idx <- seq(1L, nc, by = max(1L, step))
    for (i in r_idx) {
      fig <- fig |>
        plotly::add_trace(
          x = X[i,], y = Y[i,], z = Z[i,],
          type = "scatter3d", mode = "lines",
          line = line, hoverinfo = "none", showlegend = FALSE
        )
    }
    for (j in c_idx) {
      fig <- fig |>
        plotly::add_trace(
          x = X[, j], y = Y[, j], z = Z[, j],
          type = "scatter3d", mode = "lines",
          line = line, hoverinfo = "none", showlegend = FALSE
        )
    }
    fig
  }

  add_surface <- function(fig, F) {
    contours_arg <- if (isTRUE(show_surface_grid)) list(
      x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
      y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
      z = list(show = FALSE)
    ) else NULL

    fig |>
      plotly::add_surface(
        x = F$X, y = F$Y, z = F$Z,
        colorscale = as_colorscale(colorscale),
        showscale  = FALSE,
        opacity    = opacity,
        contours   = contours_arg
      )
  }

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("You need 'plotly' installed to plot.", call. = FALSE)
    } else {
      fig <- plotly::plot_ly()

      # theta faces
      Fth_min <- build_face_theta(th_min)
      Fth_max <- build_face_theta(th_max)
      for (Fth in list(Fth_min, Fth_max)) {
        if (mode != "wireframe") fig <- add_surface(fig, Fth)
        if (mode != "faces")     fig <- add_wireframe(fig, Fth)
        fig <- add_edges(fig, Fth)
      }

      # inner/outer radius faces
      Fin  <- build_face_r("inner")
      Fout <- build_face_r("outer")
      for (Fr in list(Fin, Fout)) {
        if (mode != "wireframe") fig <- add_surface(fig, Fr)
        if (mode != "faces")     fig <- add_wireframe(fig, Fr)
        fig <- add_edges(fig, Fr)
      }

      # z = Z1, z = Z2 faces
      Fz1 <- build_face_z("low")
      Fz2 <- build_face_z("high")
      for (Fz in list(Fz1, Fz2)) {
        if (mode != "wireframe") fig <- add_surface(fig, Fz)
        if (mode != "faces")     fig <- add_wireframe(fig, Fz)
        fig <- add_edges(fig, Fz)
      }

      fig <- fig |>
        plotly::layout(
          title = "Cylindrical solid",
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )
      print(fig)
    }
  }

  # ----- Volume -------------------------------------------------------
  volume <- NULL
  if (isTRUE(compute_volume)) {
    R1v <- Vectorize(R1)
    R2v <- Vectorize(R2)
    Z1v <- Vectorize(Z1)
    Z2v <- Vectorize(Z2)

    if (vol_method == "adaptive") {

      # NOTE: stats::integrate passes a VECTOR of theta values to the integrand.
      # This wrapper makes the inner integral robust by evaluating it per scalar theta.
      inner_r <- function(theta) {
        vapply(theta, function(th) {
          rL <- R1v(th)
          rU <- R2v(th)

          if (rU < rL) {
            tmp <- rL; rL <- rU; rU <- tmp
          }
          if (!is.finite(rL) || !is.finite(rU) || rU <= rL) return(0)

          g <- function(r) (Z2v(r, th) - Z1v(r, th)) * r
          stats::integrate(g, lower = rL, upper = rU, rel.tol = 1e-6)$value
        }, numeric(1))
      }

      val <- stats::integrate(
        function(th) inner_r(th),
        lower = th_min, upper = th_max, rel.tol = 1e-6
      )$value

      volume <- list(estimate = val, method = "adaptive")

    } else {
      ths <- seq(th_min, th_max, length.out = ntheta_vol)
      dth <- (th_max - th_min) / (ntheta_vol - 1L)
      accum_th <- numeric(length(ths))

      for (i in seq_along(ths)) {
        th <- ths[i]
        rL <- R1v(th)
        rU <- R2v(th)
        if (rU < rL) {
          tmp <- rL; rL <- rU; rU <- tmp
        }
        if (!is.finite(rL) || !is.finite(rU) || rU <= rL) {
          accum_th[i] <- 0
          next
        }
        rs <- seq(rL, rU, length.out = nr_vol)
        if (length(rs) < 2L) {
          accum_th[i] <- 0
          next
        }
        dr   <- (rU - rL) / (nr_vol - 1L)
        vals <- (Z2v(rs, th) - Z1v(rs, th)) * rs
        Ir   <- dr * (sum(vals) - 0.5 * vals[1L] - 0.5 * vals[length(vals)])
        accum_th[i] <- Ir
      }

      val <- dth * (sum(accum_th) - 0.5 * accum_th[1L] - 0.5 * accum_th[length(accum_th)])
      volume <- list(
        estimate = val,
        method   = "grid",
        ntheta   = ntheta_vol,
        nr       = nr_vol
      )
    }
  }

  list(
    theta_seq = theta_seq,
    u_seq     = u_seq,
    v_seq     = v_seq,
    fig       = fig,
    volume    = volume
  )
}
