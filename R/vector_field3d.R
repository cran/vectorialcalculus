#' 3D vector field in a curvilinear prism
#'
#' @description
#' Displays a three-dimensional vector field inside a solid region whose
#' bounds in the variables \code{x}, \code{y} and \code{z} are defined by
#' user-supplied functions. The region is described by an interval for
#' \code{x}, lower and upper bounds in the \code{y} direction depending on
#' \code{x}, and lower and upper bounds in the \code{z} direction that may
#' depend on both \code{x} and \code{y}. The function samples this volume on
#' a regular grid and draws arrows representing the vector field using
#' \pkg{plotly}.
#'
#' @details
#' The domain is parameterized by three normalized parameters, one for each
#' direction. For each grid point, the corresponding physical coordinates
#' in \code{x}, \code{y} and \code{z} are obtained by linear interpolation
#' between the lower and upper bounds. The vector field is evaluated at
#' each of these points.
#'
#' Arrow lengths are scaled using a saturated version of the vector norm.
#' This avoids extremely long arrows when the magnitude of the field varies
#' strongly across the region. A bias parameter controls how quickly the
#' lengths approach saturation: small magnitudes produce short arrows and
#' large magnitudes are capped so that they remain visible without
#' dominating the picture.
#'
#' Depending on the selected mode, the function can:
#' \itemize{
#'   \item draw only line segments representing the arrow shafts,
#'   \item draw only arrow heads (cones or chevrons),
#'   \item or combine both shafts and heads.
#' }
#'
#' The plotted figure can be customized through colors, opacity settings,
#' line widths and standard \pkg{plotly} scene options. If plotting is
#' disabled, the function still returns the sampled data for further
#' processing.
#'
#' @param F A function \code{function(x, y, z)} returning a numeric vector
#'   of length three, interpreted as \code{c(Fx, Fy, Fz)}.
#' @param H1,H2 Functions of one variable \code{x} giving the lower and
#'   upper bounds in the \code{y} direction.
#' @param G1,G2 Functions of two variables \code{x} and \code{y} giving the
#'   lower and upper bounds in the \code{z} direction.
#' @param a,b Numeric endpoints of the interval for \code{x}. It is assumed
#'   that \code{b > a}.
#' @param NX,NY,NZ Integers greater than or equal to one giving the grid
#'   density in the three parameter directions. Each parameter is sampled
#'   using a regular sequence between zero and one with \code{N + 1} points.
#' @param plot Logical; if \code{TRUE}, the vector field is drawn using
#'   \pkg{plotly}.
#' @param arrows Character string indicating the arrow mode. Allowed values
#'   are \code{"both"}, \code{"line"}, \code{"cone"} and \code{"none"}.
#' @param arrow_scale Global length scale for arrows, expressed as a
#'   fraction of the largest span of the bounding box.
#' @param normalize_bias Numeric saturation bias used in the scaling of the
#'   vector norm. Larger values make the arrow lengths saturate earlier.
#' @param normal_color Color for the arrow shafts (line segments).
#' @param normal_width Numeric width of the arrow shafts.
#' @param arrow_color Color for the arrow heads (cones or chevrons).
#' @param arrow_opacity Numeric opacity for arrow heads when cones are
#'   available.
#' @param arrow_size Relative size of arrow heads with respect to
#'   \code{arrow_scale}.
#' @param scene List with 3D scene options for \pkg{plotly}.
#' @param bg List with background colors for the figure, typically with
#'   entries \code{paper} and \code{plot}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{points}: a data frame with base positions \code{x, y, z}
#'         and the magnitude of the field at each point,
#'   \item \code{segments}: a data frame with columns
#'         \code{x0, y0, z0, x1, y1, z1} describing the arrow shafts,
#'   \item \code{fig}: a \pkg{plotly} object when \code{plot = TRUE},
#'         otherwise \code{NULL}.
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' H1 <- function(x) -1 - x
#' H2 <- function(x)  1 - x^2
#' G1 <- function(x, y) y
#' G2 <- function(x, y) y + 1
#'
#' F <- function(x, y, z) c(-y, x, 1)
#'
#' vector_field3d(
#'   F, H1 = H1, H2 = H2, G1 = G1, G2 = G2,
#'   a = -1, b = 1, NX = 8, NY = 6, NZ = 6,
#'   plot = TRUE, arrows = "both",
#'   arrow_scale = 0.08, normalize_bias = 1,
#'   normal_color = "rgba(0,0,0,0.6)", normal_width = 2,
#'   arrow_color  = "#1d3557", arrow_opacity = 0.95, arrow_size = 0.35
#' )
#' \dontshow{\}}
#'
#' @export
vector_field3d <- function(
    F, H1, H2, G1, G2,
    a, b, NX = 8, NY = 6, NZ = 6,
    plot = TRUE,
    arrows = c("both","line","cone","none"),
    arrow_scale = 0.08,
    normalize_bias = 1,
    normal_color = "black",
    normal_width = 1.5,
    arrow_color  = "black",
    arrow_opacity = 0.9,
    arrow_size    = 0.35,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  arrows <- match.arg(arrows)

  # ---- Basic validations
  if (!is.function(F))  stop("'F' must be function(x,y,z) -> c(Fx,Fy,Fz).", call. = FALSE)
  if (!is.function(H1) || !is.function(H2)) stop("'H1' and 'H2' must be functions of x.", call. = FALSE)
  if (!is.function(G1) || !is.function(G2)) stop("'G1' and 'G2' must be functions of (x,y).", call. = FALSE)
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || b <= a) {
    stop("'a' and 'b' must be scalars with b > a.", call. = FALSE)
  }
  for (nm in c("NX","NY","NZ")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1 ||
        abs(v - round(v)) > .Machine$double.eps^0.5) {
      stop(sprintf("'%s' must be an integer >= 1.", nm), call. = FALSE)
    }
  }

  # ---- Color helpers for cones
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
    rgb <- grDevices::col2rgb(col) / 255
    a_  <- if (is.null(alpha)) 1 else max(0, min(1, alpha))
    sprintf("rgba(%g,%g,%g,%g)", 255*rgb[1], 255*rgb[2], 255*rgb[3], a_)
  }
  as_colorscale <- function(x, alpha = NULL) {
    if (is_color(x)) {
      ccol <- as_rgba(x, alpha)
      return(list(list(0, ccol), list(1, ccol)))
    }
    if (is.character(x) && length(x) > 1) {
      cols <- vapply(x, as_rgba, character(1), alpha = alpha)
      pos  <- seq(0, 1, length.out = length(cols))
      return(lapply(seq_along(cols), function(i) list(pos[i], cols[i])))
    }
    x
  }

  # ---- Grids in [0,1] for i,j,k and mapping to the volume
  Ii <- seq(0, 1, length.out = NX + 1)
  Jj <- seq(0, 1, length.out = NY + 1)
  Kk <- seq(0, 1, length.out = NZ + 1)

  n_tot <- length(Ii) * length(Jj) * length(Kk)
  P0 <- matrix(NA_real_, n_tot, 3)  # base points
  V  <- matrix(NA_real_, n_tot, 3)  # vectors F
  idx <- 0L

  for (i in Ii) {
    x <- (1 - i) * a + i * b
    y1 <- H1(x); y2 <- H2(x)
    for (j in Jj) {
      y <- (1 - j) * y1 + j * y2
      g1 <- G1(x, y); g2 <- G2(x, y)
      for (k in Kk) {
        z <- (1 - k) * g1 + k * g2
        idx <- idx + 1L
        P0[idx, ] <- c(x, y, z)
        vec <- F(x, y, z)
        if (!is.numeric(vec) || length(vec) != 3L) {
          stop("F(x,y,z) must return c(Fx,Fy,Fz).", call. = FALSE)
        }
        V[idx, ] <- vec
      }
    }
  }

  # ---- Bounding box and length scale
  rx <- range(P0[,1], finite = TRUE)
  ry <- range(P0[,2], finite = TRUE)
  rz <- range(P0[,3], finite = TRUE)
  span <- max(diff(rx), diff(ry), diff(rz))
  if (!is.finite(span) || span <= 0) span <- 1

  # ---- Saturated lengths and directions
  vmag <- sqrt(rowSums(V*V))
  Vsat <- V / sqrt(vmag^2 + normalize_bias)  # saturation
  lsat <- sqrt(rowSums(Vsat*Vsat))           # factor in [0,1)
  diru <- Vsat
  nzr  <- lsat > 0
  diru[nzr, ]  <- Vsat[nzr, , drop = FALSE] / lsat[nzr]
  diru[!nzr, ] <- 0

  L  <- arrow_scale * span * lsat
  P1 <- P0 + diru * L

  # ---- Output data (useful if plot = FALSE)
  points_df <- data.frame(
    x = P0[,1], y = P0[,2], z = P0[,3],
    mag = vmag
  )
  seg_df <- data.frame(
    x0 = P0[,1], y0 = P0[,2], z0 = P0[,3],
    x1 = P1[,1], y1 = P1[,2], z1 = P1[,3]
  )

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      has_cone <- "add_cone" %in% getNamespaceExports("plotly")
      plt <- plotly::plot_ly()

      # Shafts (lines)
      if (arrows %in% c("both","line")) {
        xs <- as.numeric(rbind(seg_df$x0, seg_df$x1, NA))
        ys <- as.numeric(rbind(seg_df$y0, seg_df$y1, NA))
        zs <- as.numeric(rbind(seg_df$z0, seg_df$z1, NA))
        plt <- plt |>
          plotly::add_trace(
            x = xs, y = ys, z = zs,
            type = "scatter3d", mode = "lines",
            line = list(color = normal_color, width = normal_width),
            hoverinfo = "none", showlegend = FALSE
          )
      }

      # Heads: cones if available; otherwise chevrons
      if (arrows %in% c("both","cone")) {
        if (isTRUE(has_cone)) {
          cone_len <- max(1e-8, arrow_scale * span * arrow_size)
          cs_cone  <- as_colorscale(arrow_color, alpha = arrow_opacity)
          plt <- plt |>
            plotly::add_trace(
              type = "cone",
              x = P1[,1], y = P1[,2], z = P1[,3],
              u = diru[,1], v = diru[,2], w = diru[,3],
              anchor   = "tip",
              sizemode = "absolute",
              sizeref  = cone_len,
              colorscale = cs_cone,
              showscale  = FALSE,
              opacity    = arrow_opacity
            )
        } else {
          # Fallback: chevron (two short lines) at the tip
          head_len <- max(1e-8, arrow_scale * span * arrow_size)
          head_w   <- head_len * 0.5

          perp_unit <- function(n) {
            if (abs(n[1]) <= abs(n[2]) && abs(n[1]) <= abs(n[3])) v <- c(0, -n[3], n[2])
            else if (abs(n[2]) <= abs(n[1]) && abs(n[2]) <= abs(n[3])) v <- c(-n[3], 0, n[1])
            else v <- c(-n[2], n[1], 0)
            v / sqrt(sum(v*v))
          }

          hx <- hy <- hz <- numeric(0)
          for (i in seq_len(nrow(P0))) {
            if (lsat[i] <= 0) next
            tip <- P1[i, ]; n <- diru[i, ]
            p   <- perp_unit(n)
            base  <- tip - head_len * n
            left  <- base + head_w * p
            right <- base - head_w * p
            hx <- c(hx, tip[1], left[1], NA, tip[1], right[1], NA)
            hy <- c(hy, tip[2], left[2], NA, tip[2], right[2], NA)
            hz <- c(hz, tip[3], left[3], NA, tip[3], right[3], NA)
          }
          plt <- plt |>
            plotly::add_trace(
              x = hx, y = hy, z = hz,
              type = "scatter3d", mode = "lines",
              line = list(color = arrow_color, width = max(1, normal_width)),
              hoverinfo = "none", showlegend = FALSE
            )
        }
      }

      plt <- plt |>
        plotly::layout(
          title = "3D vector field",
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      fig <- plt
      print(plt)
    }
  }

  list(points = points_df, segments = seg_df, fig = fig)
}
