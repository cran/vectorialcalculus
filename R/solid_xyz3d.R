#' Solid defined by bounds in x, y and z
#'
#' @description
#' Constructs a three-dimensional solid defined by bounds in the variables
#' \code{x}, \code{y} y \code{z}, and optionally renders it using
#' \pkg{plotly}. The solid is described by:
#' \itemize{
#'   \item an interval for \code{x} between \code{a} y \code{b},
#'   \item lower and upper functions in the \code{y} direction,
#'   \item lower and upper functions in the \code{z} direction that may
#'         depend on both \code{x} y \code{y}.
#' }
#' The function uses a curvilinear-prism parametrization to build meshes for
#' the six faces of the solid. It supports different display modes (faces,
#' wireframe, or both), optional numerical volume computation, and internal
#' slices on coordinate planes (slice arguments are reserved for future
#' extensions and are currently ignored).
#'
#' @details
#' The solid is sampled on a three-parameter grid. Two of the parameters
#' describe the position on the base region in the \code{x}-\code{y} plane,
#' and the third parameter interpolates between the lower and upper
#' \code{z} bounds. From this parametrization the function constructs the
#' six bounding faces, corresponding to the two extreme values of \code{x},
#' the two extreme values of \code{y}, and the two extreme values of
#' \code{z}.
#'
#' Rendering options allow:
#' \itemize{
#'   \item drawing only the faces of the solid,
#'   \item drawing only a wireframe of the mesh,
#'   \item combining both faces and wireframe,
#'   \item assigning individual color scales and opacities to each face,
#'   \item showing or hiding surface grids and edges.
#' }
#'
#' When internal slices are requested, the intention is to intersect the solid
#' with planes of the form \code{x = constant}, \code{y = constant} o
#' \code{z = constant}. The corresponding slice arguments are reserved for
#' future versions of the function and are not yet implemented.
#'
#' If \code{compute_volume = TRUE}, the function also computes an
#' approximate volume of the solid using either:
#' \itemize{
#'   \item a nested adaptive integration based on \code{stats::integrate},
#'   \item or a trapezoidal rule on a regular grid in the \code{x} y
#'         \code{y} directions.
#' }
#'
#' @param H1,H2 Functions of one variable \code{x} giving the lower and
#' upper bounds in the \code{y} direction.
#' @param G1,G2 Functions of two variables \code{x} y \code{y} giving the
#' lower and upper bounds in the \code{z} direction.
#' @param a,b Numeric endpoints of the interval for \code{x}. It is assumed
#' that \code{b > a}.
#' @param plot Logical; if \code{TRUE}, the solid is rendered with
#' \pkg{plotly}.
#' @param n_x,n_u,n_v Integers giving the mesh resolution in the principal
#' parameter along \code{x} and in the two internal parameters of the face
#' meshes.
#' @param mode Character string; one of \code{"faces"}, \code{"wireframe"}
#' or \code{"both"}, indicating whether to draw surfaces, wireframe, or a
#' combination of both.
#' @param show_faces Logical vector indicating which of the six faces to
#' display. The order is \code{c("x=a","x=b","y=H1","y=H2","z=G1","z=G2")}.
#' A single logical value is also allowed and will be recycled.
#' @param colorscales Color specification for faces. It can be:
#' \itemize{
#'   \item a single \pkg{plotly} colorscale name applied to all faces,
#'   \item a single flat color (R color name, hexadecimal code, or
#'         \code{"rgba(...)"}),
#'   \item a vector of colors that define a gradient,
#'   \item or a list or vector of length six, assigning a scale or color
#'         specification to each face separately.
#' }
#' @param opacities Numeric values controlling face opacity; may be a single
#' value or a vector of length six.
#' @param show_surface_grid Logical; if \code{TRUE}, draws grid lines on the
#' faces.
#' @param surface_grid_color,surface_grid_width Color and width for surface
#' grid lines.
#' @param show_edges Logical; if \code{TRUE}, draws the edges of each face.
#' @param edge_line List with style options for edges (for example, color,
#' width and dash pattern).
#' @param wire_step Integer greater or equal to one; controls how many mesh
#' lines are skipped between wireframe lines.
#' @param wire_line List with style options for wireframe lines.
#' @param scene List with 3D scene options for \pkg{plotly}. By default, an
#' aspect ratio based on the data is used.
#' @param bg List specifying background colors for the figure, typically
#' with entries \code{paper} and \code{plot}.
#' @param compute_volume Logical; if \code{TRUE}, computes an approximate
#' volume of the solid.
#' @param vol_method Character string selecting the volume integration
#' method: \code{"adaptive"} for nested calls to \code{stats::integrate},
#' or \code{"grid"} for a trapezoidal rule on a regular grid.
#' @param nx_vol,ny_vol Integer grid sizes used when
#' \code{vol_method = "grid"}.
#' @param slice List describing slices to be drawn, with components
#' \code{x}, \code{y} y \code{z}. Each component can be \code{NULL}, a
#' single numeric value or a numeric vector of slice positions. (Reserved
#' for future use; currently ignored.)
#' @param slice_mode Character string indicating how to render slices:
#' \code{"surface"}, \code{"wireframe"} or \code{"both"}. (Reserved for
#' future use; currently ignored.)
#' @param slice_nx,slice_nu,slice_nv Mesh resolutions used to build the
#' slices. (Reserved for future use; currently ignored.)
#' @param slice_colorscales List with color scales for slices in the
#' \code{x}, \code{y} y \code{z} directions, in the same formats accepted
#' by \code{colorscales}. (Reserved for future use; currently ignored.)
#' @param slice_opacity Numeric opacity for slices, between 0 and 1.
#' (Reserved for future use; currently ignored.)
#' @param slice_show_grid Logical; if \code{TRUE}, draws grid lines on the
#' slices. (Reserved for future use; currently ignored.)
#' @param slice_grid_color,slice_grid_width Color and width for slice grid
#' lines. (Reserved for future use; currently ignored.)
#' @param slice_wire_step Integer controlling the spacing of wireframe
#' lines on slices. (Reserved for future use; currently ignored.)
#' @param slice_wire_line List with style options for slice wireframe
#' lines. (Reserved for future use; currently ignored.)
#'
#' @return
#' A list with:
#' \itemize{
#'   \item \code{x_seq}, \code{u_seq}, \code{v_seq}: the parameter sequences
#'         used to build the mesh,
#'   \item \code{fig}: a \pkg{plotly} object when \code{plot = TRUE},
#'         otherwise \code{NULL},
#'   \item \code{volume}: either \code{NULL} or a list with an approximate
#'         volume estimate and related metadata when
#'         \code{compute_volume = TRUE}.
#' }
#'
#' @examples
#' # Note: examples avoid plotting for CRAN checks
#' H1 <- function(x) -1 - x
#' H2 <- function(x)  1 - x^2
#' G1 <- function(x, y) y
#' G2 <- function(x, y) y + 1
#' s <- solid_xyz3d(
#'   H1, H2, G1, G2,
#'   a = -1, b = 1,
#'   plot = FALSE,
#'   compute_volume = TRUE,
#'   vol_method = "grid",
#'   nx_vol = 50, ny_vol = 50
#' )
#' s$volume$estimate
#'
#' @export
solid_xyz3d <- function(
    H1, H2, G1, G2,
    a, b,
    plot = TRUE,
    n_x = 120, n_u = 60, n_v = 60,
    mode = c("faces","wireframe","both"),
    show_faces = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    colorscales = c("Blues","Blues","Greens","Greens","Reds","Reds"),
    opacities   = 0.35,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    show_edges = TRUE,
    edge_line  = list(color = "black", width = 2),
    wire_step  = 6,
    wire_line  = list(color = "black", width = 1),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white"),
    # Volume
    compute_volume = FALSE,
    vol_method = c("adaptive","grid"),
    nx_vol = 300, ny_vol = 300,
    # Slices (reservados para versiones futuras)
    slice = list(x = NULL, y = NULL, z = NULL),
    slice_mode = c("surface","wireframe","both"),
    slice_nx = 200, slice_nu = 120, slice_nv = 120,
    slice_colorscales = list(x = "Oranges", y = "Purples", z = "Greens"),
    slice_opacity = 0.55,
    slice_show_grid  = TRUE,
    slice_grid_color = "rgba(80,80,80,0.25)",
    slice_grid_width = 1,
    slice_wire_step = 8,
    slice_wire_line = list(color = "black", width = 2, dash = "dot")
)
{
  # ---------- helpers ----------
  assert_fun <- function(f, name) {
    if (!is.function(f)) {
      stop(sprintf("'%s' must be a function.", name), call. = FALSE)
    }
  }
  assert_num <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      stop(sprintf("'%s' must be a finite numeric scalar.", name), call. = FALSE)
    }
  }
  assert_pos_int <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 1 || abs(x - round(x)) > .Machine$double.eps^0.5) {
      stop(sprintf("'%s' must be a positive integer.", name), call. = FALSE)
    }
  }
  as_vec6_exact <- function(x, name) {
    if (length(x) == 1L) {
      rep(x, 6)
    } else if (length(x) == 6L) {
      x
    } else {
      stop(sprintf("'%s' must have length 1 or 6.", name), call. = FALSE)
    }
  }

  is_color <- function(x) {
    if (!is.character(x) || length(x) != 1) return(FALSE)
    if (grepl("^rgba?\\(", x, ignore.case = TRUE)) return(TRUE)
    if (grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)) return(TRUE)
    ok <- TRUE
    tryCatch(grDevices::col2rgb(x), error = function(...) ok <<- FALSE)
    ok
  }
  as_rgba <- function(col, alpha = NULL) {
    if (grepl("^rgba?\\(", col, ignore.case = TRUE)) return(col)
    rgb <- grDevices::col2rgb(col) # 0..255 integers
    a <- if (is.null(alpha)) 1 else max(0, min(1, alpha))
    sprintf("rgba(%d,%d,%d,%g)", rgb[1], rgb[2], rgb[3], a)
  }
  is_plotly_scale_list <- function(x) {
    is.list(x) && length(x) >= 2 && all(vapply(
      x,
      function(it) {
        is.list(it) && length(it) == 2 &&
          is.numeric(it[[1]]) &&
          is.character(it[[2]]) && length(it[[2]]) == 1L
      },
      logical(1)
    ))
  }
  as_colorscale <- function(x, alpha = NULL) {
    if (is.null(x)) return(NULL)
    if (is_plotly_scale_list(x)) return(x)
    if (is.character(x) && length(x) == 1) {
      if (is_color(x)) {
        ccol <- as_rgba(x, alpha)
        return(list(list(0, ccol), list(1, ccol)))
      } else {
        return(x) # Plotly named scale
      }
    }
    if (is.character(x) && length(x) > 1) {
      cols <- vapply(x, as_rgba, character(1), alpha = alpha)
      pos  <- seq(0, 1, length.out = length(cols))
      return(Map(function(p, c) list(p, c), pos, cols))
    }
    stop("Unrecognized 'colorscale' format.", call. = FALSE)
  }

  # ---------- validations ----------
  assert_fun(H1,"H1"); assert_fun(H2,"H2")
  assert_fun(G1,"G1"); assert_fun(G2,"G2")
  assert_num(a,"a");  assert_num(b,"b")
  if (b <= a) stop("'b' must be > 'a'.", call. = FALSE)

  for (nm in c("n_x","n_u","n_v","wire_step",
               "nx_vol","ny_vol",
               "slice_nx","slice_nu","slice_nv","slice_wire_step")) {
    assert_pos_int(get(nm), nm)
  }

  mode        <- match.arg(mode)
  vol_method  <- match.arg(vol_method)
  slice_mode  <- match.arg(slice_mode)

  # Para el método de volumen por rejilla se requieren al menos 2 puntos
  if (vol_method == "grid") {
    if (nx_vol < 2) stop("'nx_vol' must be >= 2 for grid volume.", call. = FALSE)
    if (ny_vol < 2) stop("'ny_vol' must be >= 2 for grid volume.", call. = FALSE)
  }

  show_faces <- as_vec6_exact(show_faces, "show_faces")
  opacities  <- as_vec6_exact(opacities,  "opacities")

  # ---------- colors (6 faces + slices) ----------
  if (length(colorscales) == 6L) {
    colorscales <- lapply(colorscales, as_colorscale)
  } else {
    cs_global <- as_colorscale(colorscales)
    colorscales <- rep(list(cs_global), 6L)
  }
  slice_cs_x <- as_colorscale(slice_colorscales$x)
  slice_cs_y <- as_colorscale(slice_colorscales$y)
  slice_cs_z <- as_colorscale(slice_colorscales$z)
  # (slice_* se reservan para uso futuro)

  # ---------- parametrización ----------
  y_blend <- function(x, u) (1 - u) * H1(x) + u * H2(x)
  z_blend <- function(x, y, v) (1 - v) * G1(x, y) + v * G2(x, y)

  x_seq <- seq(a, b, length.out = n_x)
  u_seq <- seq(0, 1, length.out = n_u)
  v_seq <- seq(0, 1, length.out = n_v)

  # clamp steps to mesh sizes (defensivo)
  wire_step       <- max(1L, min(wire_step, max(n_u, n_v)))
  slice_wire_step <- max(1L, min(slice_wire_step, max(slice_nu, slice_nv)))

  # ---------- drawing helpers ----------
  add_edges <- function(fig, X, Y, Z, line = edge_line) {
    i1 <- 1; i2 <- nrow(X); j1 <- 1; j2 <- ncol(X)
    fig |>
      plotly::add_trace(
        x = X[i1,], y = Y[i1,], z = Z[i1,],
        type = "scatter3d", mode = "lines",
        line = line, hoverinfo = "none", showlegend = FALSE
      ) |>
      plotly::add_trace(
        x = X[i2,], y = Y[i2,], z = Z[i2,],
        type = "scatter3d", mode = "lines",
        line = line, hoverinfo = "none", showlegend = FALSE
      ) |>
      plotly::add_trace(
        x = X[,j1], y = Y[,j1], z = Z[,j1],
        type = "scatter3d", mode = "lines",
        line = line, hoverinfo = "none", showlegend = FALSE
      ) |>
      plotly::add_trace(
        x = X[,j2], y = Y[,j2], z = Z[,j2],
        type = "scatter3d", mode = "lines",
        line = line, hoverinfo = "none", showlegend = FALSE
      )
  }

  add_wireframe <- function(fig, X, Y, Z, step, line = wire_line) {
    nr <- nrow(X); nc <- ncol(X)
    r_idx <- seq(1, nr, by = max(1L, step))
    c_idx <- seq(1, nc, by = max(1L, step))
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
          x = X[,j], y = Y[,j], z = Z[,j],
          type = "scatter3d", mode = "lines",
          line = line, hoverinfo = "none", showlegend = FALSE
        )
    }
    fig
  }

  add_surface <- function(fig, X, Y, Z, cs, op, show_grid, grid_col, grid_w) {
    contours_arg <- if (isTRUE(show_grid)) list(
      x = list(show = TRUE, color = grid_col, width = grid_w),
      y = list(show = TRUE, color = grid_col, width = grid_w),
      z = list(show = FALSE)
    ) else NULL
    fig |>
      plotly::add_surface(
        x = X, y = Y, z = Z,
        colorscale = cs, showscale = FALSE,
        opacity = op, contours = contours_arg
      )
  }

  build_face_x <- function(xfix) {
    U <- seq(0, 1, length.out = n_u)
    V <- seq(0, 1, length.out = n_v)
    X <- matrix(xfix, nrow = n_u, ncol = n_v)
    y_line <- y_blend(xfix, U)
    Y <- matrix(rep(y_line, times = n_v), nrow = n_u)
    Z <- matrix(NA_real_, n_u, n_v)
    for (j in seq_len(n_v)) {
      Z[, j] <- z_blend(xfix, y_line, V[j])
    }
    list(X = X, Y = Y, Z = Z)
  }

  build_face_y <- function(y_fun) {
    V <- seq(0, 1, length.out = n_v)
    X <- matrix(rep(x_seq, each = n_v), nrow = n_v)
    y_line <- y_fun(x_seq)
    Y <- matrix(rep(y_line, each = n_v), nrow = n_v)
    Z <- matrix(NA_real_, n_v, n_x)
    for (j in seq_len(n_x)) {
      Z[, j] <- z_blend(x_seq[j], y_line[j], V)
    }
    list(X = X, Y = Y, Z = Z)
  }

  build_face_z <- function(z_fun) {
    U <- seq(0, 1, length.out = n_u)
    X <- matrix(rep(x_seq, each = n_u), nrow = n_u)
    Y <- matrix(NA_real_, n_u, n_x)
    Z <- matrix(NA_real_, n_u, n_x)
    for (j in seq_len(n_x)) {
      yline <- y_blend(x_seq[j], U)
      Y[, j] <- yline
      Z[, j] <- z_fun(x_seq[j], yline)
    }
    list(X = X, Y = Y, Z = Z)
  }

  # ---------- drawing ----------
  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("You need 'plotly' installed to draw.", call. = FALSE)
    } else {
      fig <- plotly::plot_ly()

      if (isTRUE(show_faces[1])) {
        F <- build_face_x(a)
        if (mode != "wireframe") {
          fig <- add_surface(fig, F$X, F$Y, F$Z,
                             colorscales[[1]], opacities[1],
                             show_surface_grid, surface_grid_color, surface_grid_width)
        }
        if (mode != "faces") {
          fig <- add_wireframe(fig, F$X, F$Y, F$Z, wire_step)
        }
        if (isTRUE(show_edges)) {
          fig <- add_edges(fig, F$X, F$Y, F$Z)
        }
      }

      if (isTRUE(show_faces[2])) {
        F <- build_face_x(b)
        if (mode != "wireframe") {
          fig <- add_surface(fig, F$X, F$Y, F$Z,
                             colorscales[[2]], opacities[2],
                             show_surface_grid, surface_grid_color, surface_grid_width)
        }
        if (mode != "faces") {
          fig <- add_wireframe(fig, F$X, F$Y, F$Z, wire_step)
        }
        if (isTRUE(show_edges)) {
          fig <- add_edges(fig, F$X, F$Y, F$Z)
        }
      }

      if (isTRUE(show_faces[3])) {
        F <- build_face_y(H1)
        if (mode != "wireframe") {
          fig <- add_surface(fig, F$X, F$Y, F$Z,
                             colorscales[[3]], opacities[3],
                             show_surface_grid, surface_grid_color, surface_grid_width)
        }
        if (mode != "faces") {
          fig <- add_wireframe(fig, F$X, F$Y, F$Z, wire_step)
        }
        if (isTRUE(show_edges)) {
          fig <- add_edges(fig, F$X, F$Y, F$Z)
        }
      }

      if (isTRUE(show_faces[4])) {
        F <- build_face_y(H2)
        if (mode != "wireframe") {
          fig <- add_surface(fig, F$X, F$Y, F$Z,
                             colorscales[[4]], opacities[4],
                             show_surface_grid, surface_grid_color, surface_grid_width)
        }
        if (mode != "faces") {
          fig <- add_wireframe(fig, F$X, F$Y, F$Z, wire_step)
        }
        if (isTRUE(show_edges)) {
          fig <- add_edges(fig, F$X, F$Y, F$Z)
        }
      }

      if (isTRUE(show_faces[5])) {
        F <- build_face_z(function(x, y) G1(x, y))
        if (mode != "wireframe") {
          fig <- add_surface(fig, F$X, F$Y, F$Z,
                             colorscales[[5]], opacities[5],
                             show_surface_grid, surface_grid_color, surface_grid_width)
        }
        if (mode != "faces") {
          fig <- add_wireframe(fig, F$X, F$Y, F$Z, wire_step)
        }
        if (isTRUE(show_edges)) {
          fig <- add_edges(fig, F$X, F$Y, F$Z)
        }
      }

      if (isTRUE(show_faces[6])) {
        F <- build_face_z(function(x, y) G2(x, y))
        if (mode != "wireframe") {
          fig <- add_surface(fig, F$X, F$Y, F$Z,
                             colorscales[[6]], opacities[6],
                             show_surface_grid, surface_grid_color, surface_grid_width)
        }
        if (mode != "faces") {
          fig <- add_wireframe(fig, F$X, F$Y, F$Z, wire_step)
        }
        if (isTRUE(show_edges)) {
          fig <- add_edges(fig, F$X, F$Y, F$Z)
        }
      }

      fig <- fig |>
        plotly::layout(
          title = sprintf("Solid (mode = %s)", mode),
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )
      if (interactive()) print(fig)
    }
  }

  # ---------- volume ----------
  volume <- NULL
  if (isTRUE(compute_volume)) {
    # thickness always nonnegative
    f_gap  <- function(x, y) abs(G2(x, y) - G1(x, y))
    f_gapv <- Vectorize(f_gap, SIMPLIFY = TRUE)
    H1v <- Vectorize(H1); H2v <- Vectorize(H2)

    if (vol_method == "adaptive") {
      inner <- function(x) {
        yl <- pmin(H1v(x), H2v(x))
        yh <- pmax(H1v(x), H2v(x))
        g  <- function(y) f_gapv(x, y)
        sapply(seq_along(x), function(i) {
          res <- try(
            stats::integrate(
              g,
              lower = yl[i], upper = yh[i],
              rel.tol = 1e-6, stop.on.error = FALSE
            ),
            silent = TRUE
          )
          if (inherits(res, "try-error") || !is.finite(res$value)) 0 else res$value
        })
      }
      res <- try(
        stats::integrate(
          function(x) inner(x),
          lower = a, upper = b,
          rel.tol = 1e-6, stop.on.error = FALSE
        ),
        silent = TRUE
      )
      vol_val <- if (inherits(res, "try-error") || !is.finite(res$value)) NA_real_ else res$value
      volume  <- list(estimate = vol_val, method = "adaptive")

    } else {
      xs <- seq(a, b, length.out = nx_vol)
      dx <- (b - a) / (nx_vol - 1)
      vol_x <- numeric(nx_vol)

      for (i in seq_along(xs)) {
        xl <- xs[i]
        yl <- H1(xl)
        yh <- H2(xl)
        if (yh < yl) {
          tmp <- yl; yl <- yh; yh <- tmp
        }
        ys <- seq(yl, yh, length.out = ny_vol)
        if (length(ys) < 2) {
          vol_x[i] <- 0
          next
        }
        dy <- (yh - yl) / (ny_vol - 1)
        vals <- f_gapv(xl, ys)
        Iy <- dy * (sum(vals) - 0.5 * vals[1] - 0.5 * vals[length(vals)])
        vol_x[i] <- Iy
      }
      vol_val <- dx * (sum(vol_x) - 0.5 * vol_x[1] - 0.5 * vol_x[length(vol_x)])
      volume  <- list(estimate = vol_val, method = "grid", nx = nx_vol, ny = ny_vol)
    }
  }

  list(
    x_seq = x_seq,
    u_seq = u_seq,
    v_seq = v_seq,
    fig   = fig,
    volume = volume
  )
}
