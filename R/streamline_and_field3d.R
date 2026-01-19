#' Vector field and streamline in 3D (single combined figure)
#'
#' @description
#' Draws a three-dimensional vector field inside a curvilinear volume and
#' overlays a streamline that follows the field, all in a single
#' \pkg{plotly} figure. The streamline is obtained by integrating an
#' ordinary differential equation using a fixed-step Runge-Kutta method
#' of order four (RK4), starting from an initial point.
#'
#' @details
#' The volume is defined by:
#' \itemize{
#'   \item an interval for \code{x} between \code{a} and \code{b},
#'   \item lower and upper functions for \code{y} that depend on \code{x},
#'   \item lower and upper functions for \code{z} that may depend on
#'         both \code{x} and \code{y}.
#' }
#'
#' The function builds a regular grid in three parameters and maps each
#' grid point to the physical coordinates \code{(x, y, z)} using linear
#' blends between the corresponding bounds. The vector field is then
#' evaluated at each of these points to obtain the base positions and
#' vectors used to draw the arrows.
#'
#' Arrow lengths are scaled using a saturated version of the vector norm.
#' This means that small magnitudes produce short arrows, whereas very
#' large magnitudes are limited by a bias parameter so that they do not
#' dominate the entire plot. A global scale factor controls the typical
#' arrow length relative to the size of the domain.
#'
#' The streamline is defined as the trajectory of a particle whose
#' velocity at each point is given by the vector field evaluated along
#' the path. The field may optionally depend on time; if the function
#' \code{field} has a time argument, it is used during integration. The
#' integration runs from time zero up to a final time, with positive or
#' negative direction depending on the sign of the final time.
#'
#' The resulting figure combines:
#' \itemize{
#'   \item a set of arrows representing the vector field,
#'   \item a space curve representing the streamline,
#'   \item optional markers along the streamline and highlighted start
#'         and end points.
#' }
#'
#' @param field A function representing the vector field. It can be given as
#' \code{function(x, y, z)} or \code{function(x, y, z, t)}, and must
#' return a numeric vector of length three \code{c(Fx, Fy, Fz)}.
#' @param H1,H2 Functions of one variable \code{x} giving the lower and
#' upper bounds in the \code{y} direction.
#' @param G1,G2 Functions of two variables \code{x} and \code{y} giving
#' the lower and upper bounds in the \code{z} direction.
#' @param a,b Numeric endpoints of the interval for \code{x}. It is
#' assumed that \code{b > a}.
#' @param NX,NY,NZ Integers greater than or equal to one specifying the
#' sampling density of the field in the three parameter directions.
#' @param p Numeric vector of length three giving the initial point of
#' the streamline, in the form \code{c(x0, y0, z0)}.
#' @param t_final Final integration time for the streamline. A negative value
#' integrates backward in time.
#' @param step Step size for the fixed-step RK4 integration. Must be
#' strictly positive.
#'
#' @param arrows Character string indicating the arrow mode. Allowed
#' values are \code{"none"}, \code{"line"}, \code{"cone"} and
#' \code{"both"}.
#' @param arrow_scale Global arrow length scale, expressed as a fraction
#' of the largest span of the bounding box.
#' @param normalize_bias Numeric saturation bias used in the scaling of
#' the vector norm. Larger values make arrow lengths saturate earlier.
#' @param normal_color Color of the arrow shafts (line segments).
#' @param normal_width Numeric width of the arrow shafts.
#' @param arrow_color Color of the arrow heads (cones or chevrons).
#' @param arrow_opacity Opacity of the arrow heads.
#' @param arrow_size Relative size of the arrow heads with respect to
#' \code{arrow_scale}.
#'
#' @param traj_color Color of the streamline.
#' @param traj_width Width of the streamline.
#' @param traj_markers Logical; if \code{TRUE}, draws markers along the
#' streamline.
#' @param traj_marker_size Size of the markers drawn on the streamline.
#'
#' @param scene List with 3D scene settings for \pkg{plotly}, such as
#' aspect mode and axis titles.
#' @param bg List with background colors for the figure, typically with
#' entries \code{paper} and \code{plot}.
#' @param ... Reserved for backward compatibility. Do not use.

#'
#' @return A list with:
#' \itemize{
#'   \item \code{field_points}: a data frame with base positions
#'         \code{x, y, z} and the magnitude of the field at each point,
#'   \item \code{field_segments}: a data frame with columns
#'         \code{x0, y0, z0, x1, y1, z1} describing the arrow shafts,
#'   \item \code{traj}: a data frame with the streamline data, including
#'         time, coordinates and local speed,
#'   \item \code{fig}: a \pkg{plotly} object containing the combined
#'         field and streamline visualization.
#' }
#'
#' @examples
#' \donttest{
#' H1 <- function(x) -1
#' H2 <- function(x)  1
#' G1 <- function(x, y) -0.5
#' G2 <- function(x, y)  0.5
#' field <- function(x, y, z) c(-y, x, 0.6)
#'
#' out <- streamline_and_field3d(
#'   field, H1, H2, G1, G2,
#'   a = -2, b = 2, NX = 10, NY = 8, NZ = 5,
#'   p = c(1, 0, 0), t_final = 2, step = 0.05,
#'   arrows = "both", arrow_scale = 0.12, normalize_bias = 1,
#'   normal_color = "rgba(0,0,0,0.55)", normal_width = 2,
#'   arrow_color  = "#1d3557", arrow_opacity = 0.95, arrow_size = 0.4,
#'   traj_color   = "#e63946", traj_width = 5, traj_markers = TRUE
#' )
#' }
#'
#' @export
streamline_and_field3d <- function(
    field, H1, H2, G1, G2,
    a, b, NX = 8, NY = 6, NZ = 6,
    p, t_final, step,
    # field style
    arrows = c("both", "line", "cone", "none"),
    arrow_scale = 0.08,
    normalize_bias = 1,
    normal_color = "rgba(0,0,0,0.55)",
    normal_width = 2,
    arrow_color  = "#1d3557",
    arrow_opacity = 0.95,
    arrow_size    = 0.35,
    # streamline style
    traj_color = "#e63946",
    traj_width = 5,
    traj_markers = TRUE,
    traj_marker_size = 2,
    # scene / backgrounds
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white"),
    ...
) {
  # ---- Backward compatibility: accept legacy arg names F and T via ...
  dots <- list(...)
  if (is.null(field) && !is.null(dots$F)) {
    field <- dots$F
  }
  if (missing(t_final) && !is.null(dots$T)) {
    t_final <- dots$T
  }

  arrows <- match.arg(arrows)

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for plotting.", call. = FALSE)
  }

  # --- basic checks
  if (!is.function(field)) {
    stop("'field' must be function(x,y,z) or function(x,y,z,t).", call. = FALSE)
  }
  if (!is.function(H1) || !is.function(H2)) {
    stop("'H1' and 'H2' must be functions of x.", call. = FALSE)
  }
  if (!is.function(G1) || !is.function(G2)) {
    stop("'G1' and 'G2' must be functions of (x,y).", call. = FALSE)
  }
  if (!is.numeric(a) || !is.numeric(b) ||
      length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || b <= a) {
    stop("'a' and 'b' must be scalars with b > a.", call. = FALSE)
  }
  for (nm in c("NX", "NY", "NZ")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1) {
      stop(sprintf("'%s' must be integer >= 1.", nm), call. = FALSE)
    }
  }
  if (!is.numeric(p) || length(p) != 3L || any(!is.finite(p))) {
    stop("'p' must be c(x0,y0,z0).", call. = FALSE)
  }
  if (!is.numeric(t_final) || length(t_final) != 1L || !is.finite(t_final)) {
    stop("'t_final' must be a numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(step) || length(step) != 1L || !is.finite(step) || step <= 0) {
    stop("'step' must be > 0.", call. = FALSE)
  }

  # --- helpers
  has_cone <- "add_cone" %in% getNamespaceExports("plotly")
  perp_unit <- function(n) {
    if (abs(n[1]) <= abs(n[2]) && abs(n[1]) <= abs(n[3])) {
      v <- c(0, -n[3], n[2])
    } else if (abs(n[2]) <= abs(n[1]) && abs(n[2]) <= abs(n[3])) {
      v <- c(-n[3], 0, n[1])
    } else {
      v <- c(-n[2], n[1], 0)
    }
    v / sqrt(sum(v * v))
  }

  # --- grids & curvilinear mapping
  Ii <- seq(0, 1, length.out = NX + 1)
  Jj <- seq(0, 1, length.out = NY + 1)
  Kk <- seq(0, 1, length.out = NZ + 1)

  n_tot <- length(Ii) * length(Jj) * length(Kk)
  P0 <- matrix(NA_real_, n_tot, 3)  # bases
  V  <- matrix(NA_real_, n_tot, 3)  # vectors
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
        vec <- field(x, y, z)
        if (!is.numeric(vec) || length(vec) != 3L) {
          stop("field(x,y,z) must return numeric length-3 c(Fx,Fy,Fz).", call. = FALSE)
        }
        V[idx, ] <- vec
      }
    }
  }

  # --- bounding box & saturated arrow lengths
  rx <- range(P0[, 1], finite = TRUE)
  ry <- range(P0[, 2], finite = TRUE)
  rz <- range(P0[, 3], finite = TRUE)
  span <- max(diff(rx), diff(ry), diff(rz))
  if (!is.finite(span) || span <= 0) span <- 1

  vmag <- sqrt(rowSums(V * V))
  Vsat <- V / sqrt(vmag^2 + normalize_bias)
  lsat <- sqrt(rowSums(Vsat * Vsat))
  diru <- Vsat
  nzr  <- lsat > 0
  diru[nzr, ]  <- Vsat[nzr, , drop = FALSE] / lsat[nzr]
  diru[!nzr, ] <- 0

  L  <- arrow_scale * span * lsat
  P1 <- P0 + diru * L

  field_points   <- data.frame(x = P0[, 1], y = P0[, 2], z = P0[, 3], mag = vmag)
  field_segments <- data.frame(
    x0 = P0[, 1], y0 = P0[, 2], z0 = P0[, 3],
    x1 = P1[, 1], y1 = P1[, 2], z1 = P1[, 3]
  )

  # --- streamline (internal RK4)
  field_eval <- function(x, y, z, t) {
    if (length(formals(field)) >= 4L) field(x, y, z, t) else field(x, y, z)
  }

  if (abs(t_final) < .Machine$double.eps) {
    v0 <- field_eval(p[1], p[2], p[3], 0)
    speed0 <- sqrt(sum(v0 * v0))
    traj <- data.frame(
      t = 0,
      x = p[1], y = p[2], z = p[3],
      speed = speed0
    )
  } else {
    n_steps <- ceiling(abs(t_final) / step)
    signT   <- if (t_final >= 0) 1 else -1
    dt_def  <- signT * step
    t <- numeric(n_steps + 1L); t[1] <- 0
    Y <- matrix(NA_real_, n_steps + 1L, 3); Y[1, ] <- p

    for (ii in 1:n_steps) {
      dt <- dt_def
      if ((signT > 0 && (t[ii] + dt) > t_final) ||
          (signT < 0 && (t[ii] + dt) < t_final)) {
        dt <- t_final - t[ii]
      }
      x <- Y[ii, 1]; y <- Y[ii, 2]; z <- Y[ii, 3]; ti <- t[ii]
      k1 <- field_eval(x, y, z, ti)
      k2 <- field_eval(x + 0.5 * dt * k1[1],
                       y + 0.5 * dt * k1[2],
                       z + 0.5 * dt * k1[3],
                       ti + 0.5 * dt)
      k3 <- field_eval(x + 0.5 * dt * k2[1],
                       y + 0.5 * dt * k2[2],
                       z + 0.5 * dt * k2[3],
                       ti + 0.5 * dt)
      k4 <- field_eval(x + dt * k3[1],
                       y + dt * k3[2],
                       z + dt * k3[3],
                       ti + dt)
      incr <- (k1 + 2 * k2 + 2 * k3 + k4) / 6
      Y[ii + 1, ] <- c(x, y, z) + dt * incr
      t[ii + 1]   <- t[ii] + dt
    }

    speed <- vapply(seq_along(t), function(jj) {
      fj <- field_eval(Y[jj, 1], Y[jj, 2], Y[jj, 3], t[jj])
      sqrt(sum(fj * fj))
    }, numeric(1))

    traj <- data.frame(
      t = t,
      x = Y[, 1], y = Y[, 2], z = Y[, 3],
      speed = speed
    )
  }

  # --- figure
  plt <- plotly::plot_ly()

  # field shafts
  if (arrows %in% c("both", "line", "cone")) {
    xs <- as.numeric(rbind(field_segments$x0, field_segments$x1, NA))
    ys <- as.numeric(rbind(field_segments$y0, field_segments$y1, NA))
    zs <- as.numeric(rbind(field_segments$z0, field_segments$z1, NA))
    plt <- plotly::add_trace(
      plt,
      x = xs, y = ys, z = zs,
      type = "scatter3d", mode = "lines",
      line = list(color = normal_color, width = normal_width),
      hoverinfo = "none", showlegend = FALSE
    )
  }

  # field heads
  if (arrows %in% c("both", "cone")) {
    if (isTRUE(has_cone)) {
      cone_len <- max(1e-8, arrow_scale * span * arrow_size)
      cs_cone  <- list(list(0, arrow_color), list(1, arrow_color))
      plt <- plotly::add_trace(
        plt,
        type = "cone",
        x = field_segments$x1,
        y = field_segments$y1,
        z = field_segments$z1,
        u = diru[, 1], v = diru[, 2], w = diru[, 3],
        anchor     = "tip",
        sizemode   = "absolute",
        sizeref    = cone_len,
        colorscale = cs_cone,
        showscale  = FALSE,
        opacity    = arrow_opacity
      )
    } else {
      head_len <- max(1e-8, arrow_scale * span * arrow_size)
      head_w   <- head_len * 0.5
      hx <- hy <- hz <- numeric(0)
      for (i in seq_len(nrow(field_segments))) {
        tip  <- c(field_segments$x1[i],
                  field_segments$y1[i],
                  field_segments$z1[i])
        n    <- diru[i, ]
        if (all(n == 0)) next
        base <- tip - head_len * n
        pperp <- perp_unit(n)
        left  <- base + head_w * pperp
        right <- base - head_w * pperp
        hx <- c(hx, tip[1], left[1], NA, tip[1], right[1], NA)
        hy <- c(hy, tip[2], left[2], NA, tip[2], right[2], NA)
        hz <- c(hz, tip[3], left[3], NA, tip[3], right[3], NA)
      }
      plt <- plotly::add_trace(
        plt,
        x = hx, y = hy, z = hz,
        type = "scatter3d", mode = "lines",
        line = list(color = arrow_color, width = max(1, normal_width)),
        hoverinfo = "none", showlegend = FALSE
      )
    }
  }

  # streamline
  plt <- plotly::add_trace(
    plt,
    x = traj$x, y = traj$y, z = traj$z,
    type = "scatter3d",
    mode = if (isTRUE(traj_markers)) "lines+markers" else "lines",
    line = list(color = traj_color, width = traj_width),
    marker = list(size = traj_marker_size, color = traj_color),
    showlegend = FALSE, name = "streamline"
  )

  # start/end markers
  n_last <- nrow(traj)
  plt <- plotly::add_markers(
    plt,
    x = traj$x[1], y = traj$y[1], z = traj$z[1],
    marker = list(size = traj_marker_size + 2, color = "forestgreen"),
    showlegend = FALSE, hoverinfo = "none"
  )
  plt <- plotly::add_markers(
    plt,
    x = traj$x[n_last], y = traj$y[n_last], z = traj$z[n_last],
    marker = list(size = traj_marker_size + 2, color = "firebrick"),
    showlegend = FALSE, hoverinfo = "none"
  )

  plt <- plotly::layout(
    plt,
    title = "Vector field + streamline (combined)",
    scene = scene,
    paper_bgcolor = bg$paper,
    plot_bgcolor  = bg$plot
  )

  print(plt)

  list(
    field_points   = field_points,
    field_segments = field_segments,
    traj = traj,
    fig = plt
  )
}
