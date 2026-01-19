#' Riemann rectangular prisms over a planar region
#'
#' @description
#' Approximates the double integral of a scalar function over a planar region
#' using an N-by-M rectangular partition and rectangular prisms with constant
#' height on each cell. The region is defined by an x-interval and two functions
#' giving the lower and upper y-limits. Each valid cell produces a prism whose
#' height corresponds to a chosen estimate of the function on that cell:
#' lower value, upper value, or mean value.
#'
#' When \code{plot = TRUE}, a 3D visualization of the prisms is produced using
#' \pkg{plotly}. Optionally, the actual surface \code{z = F(x, y)} can also be
#' drawn over the rectangular bounding box that contains the region.
#'
#' @param f A function \code{f(x, y)} returning a numeric scalar.
#' @param a,b Numeric endpoints of the x-interval. Must satisfy \code{a < b}.
#' @param f1,f2 Functions returning the lower and upper y-boundaries for each x.
#'   The valid y-range at each x is the interval between the minimum and maximum
#'   of these two functions.
#' @param N,M Integer numbers of subdivisions in x and y for the rectangular
#'   partition.
#' @param plot Logical. If \code{TRUE}, the 3D visualization is generated.
#' @param estimate Character. One of \code{"lower"}, \code{"upper"},
#'   \code{"mean"}, or \code{"all"}, indicating which estimate to highlight.
#' @param sample_n Integer. Number of evaluation points per direction inside
#'   each cell when computing lower, upper, and mean values.
#'
#' @param show_surface Logical. If \code{TRUE}, draws the surface
#'   \code{z = F(x, y)} over the entire rectangular bounding box.
#' @param surface_colorscale Colorscale used for the surface.
#' @param surface_opacity Opacity for the surface.
#' @param show_surface_grid Logical. If \code{TRUE}, draws a grid on the surface.
#' @param surface_grid_color Color of the grid lines.
#' @param surface_grid_width Width of the grid lines.
#'
#' @param color_by Character. Determines the value used to color the top of
#'   each prism: \code{"mean"}, \code{"lower"}, or \code{"upper"}.
#' @param top_colorscale Colorscale for the prism tops.
#' @param top_opacity Opacity of the prism tops.
#' @param side_color Color for the vertical faces of the prisms.
#' @param side_opacity Opacity of the prism sides.
#' @param frame_color Color for the prism frame lines.
#' @param frame_width Width of the frame lines.
#'
#' @param scene A list with plotly 3D scene settings.
#' @param bg A list with background colors for the figure.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{sum_lower}}{Lower Riemann sum for the chosen partition.}
#'   \item{\code{sum_upper}}{Upper Riemann sum.}
#'   \item{\code{sum_mean}}{Mean-value Riemann sum.}
#'   \item{\code{cells}}{A tibble describing all valid cells.}
#'   \item{\code{fig}}{A \pkg{plotly} object if \code{plot = TRUE}, otherwise \code{NULL}.}
#'   \item{\code{estimate}}{When \code{estimate != "all"}, the selected
#'     Riemann sum (lower/upper/mean) is repeated here for convenience.}
#' }
#'
#' @examples
#' f  <- function(x, y) x * y
#' f1 <- function(x) 0
#' f2 <- function(x) 1 - x
#' riemann_prisms3d(f, f1, f2, 0, 1, N = 10, M = 10, plot = FALSE)
#'
#' @export
riemann_prisms3d <- function(
    f, f1, f2,
    a, b, N, M,
    plot = TRUE,
    estimate = c("lower","upper","mean","all"),
    sample_n = 6,
    # visual options
    show_surface = FALSE,
    surface_colorscale = "Viridis",
    surface_opacity    = 0.35,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    color_by = c("mean","lower","upper"),
    top_colorscale = "YlOrBr",
    top_opacity    = 0.85,
    side_color     = "rgba(60,60,60,0.25)",
    side_opacity   = 0.35,
    frame_color    = "rgba(0,0,0,0.55)",
    frame_width    = 1.5,
    scene = list(aspectmode="data",
                 xaxis=list(title="x"),
                 yaxis=list(title="y"),
                 zaxis=list(title="z")),
    bg = list(paper="white", plot="white")
){
  estimate <- match.arg(estimate)
  color_by <- match.arg(color_by)

  # --- tiny helpers (colorscale handling)
  is_color <- function(x) {
    if (!is.character(x) || length(x) != 1L) return(FALSE)
    if (grepl("^rgba?\\(", x, TRUE)) return(TRUE)
    if (grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)) return(TRUE)
    ok <- TRUE
    tryCatch(grDevices::col2rgb(x), error = function(...) ok <<- FALSE)
    ok
  }
  as_rgba <- function(col, alpha = NULL) {
    if (grepl("^rgba?\\(", col, TRUE)) return(col)
    rgb <- grDevices::col2rgb(col) / 255
    a   <- if (is.null(alpha)) 1 else max(0, min(1, alpha))
    sprintf("rgba(%g,%g,%g,%g)", 255*rgb[1], 255*rgb[2], 255*rgb[3], a)
  }
  as_colorscale <- function(x, alpha = NULL) {
    if (is.null(x)) return(NULL)
    if (is.list(x) && length(x) >= 2L && is.numeric(x[[1]][[1]])) return(x)
    if (is.character(x) && length(x) == 1L) {
      if (is_color(x)) {
        ccol <- as_rgba(x, alpha)
        return(list(list(0, ccol), list(1, ccol)))
      }
      return(x) # named Plotly scale
    }
    if (is.character(x) && length(x) > 1L) {
      cols <- vapply(x, as_rgba, character(1), alpha = alpha)
      pos  <- seq(0, 1, length.out = length(cols))
      return(lapply(seq_along(cols), function(i) list(pos[i], cols[i])))
    }
    stop("Unrecognized colorscale format.", call. = FALSE)
  }

  # --- validation
  if (!is.function(f) || !is.function(f1) || !is.function(f2)) {
    stop("f, f1 and f2 must be functions.", call. = FALSE)
  }
  if (!is.numeric(a) || !is.numeric(b) ||
      length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || b <= a) {
    stop("'a' and 'b' must be finite numeric scalars with a < b.", call. = FALSE)
  }

  N <- as.integer(N)
  M <- as.integer(M)
  sample_n <- as.integer(sample_n)

  for (nm in c("N","M","sample_n")) {
    val <- get(nm)
    if (!is.finite(val) || val < 1L) {
      stop(nm, " must be a positive integer.", call. = FALSE)
    }
  }

  # --- base partition in x, uniform
  hx <- (b - a) / N
  xs_nodes <- seq(a, b, length.out = N + 1L)

  # global y-span used for the rectangular grid track;
  # we still clip each cell against the y-intervals at x0 and x1 to stay inside Omega.
  ylo_nodes <- pmin(vapply(xs_nodes, f1, numeric(1)), vapply(xs_nodes, f2, numeric(1)))
  yhi_nodes <- pmax(vapply(xs_nodes, f1, numeric(1)), vapply(xs_nodes, f2, numeric(1)))
  cmin <- min(ylo_nodes)
  dmax <- max(yhi_nodes)
  hy <- (dmax - cmin) / M
  y_nodes <- seq(cmin, dmax, length.out = M + 1L)

  # iterate rectangular base cells [x0,x1] x [y0,y1], then clip to intersection
  cells <- list()
  eval_cell_stats <- function(x0, x1, y0, y1) {
    xs <- seq(x0, x1, length.out = sample_n + 1L)
    ys <- seq(y0, y1, length.out = sample_n + 1L)
    Z  <- outer(ys, xs, function(yy, xx) f(xx, yy))
    c(
      z_lower = min(Z, na.rm = TRUE),
      z_upper = max(Z, na.rm = TRUE),
      z_mean  = mean(Z, na.rm = TRUE)
    )
  }

  for (i in seq_len(N)) {
    x0 <- xs_nodes[i]
    x1 <- xs_nodes[i + 1L]
    # y interval allowed at x0 and x1
    yl0 <- min(f1(x0), f2(x0)); yh0 <- max(f1(x0), f2(x0))
    yl1 <- min(f1(x1), f2(x1)); yh1 <- max(f1(x1), f2(x1))
    # conservative rectangle: intersect the two y-intervals
    yL <- max(yl0, yl1)
    yH <- min(yh0, yh1)
    if (yL >= yH) next  # strip has no feasible y at both ends

    for (j in seq_len(M)) {
      yy0 <- y_nodes[j]
      yy1 <- y_nodes[j + 1L]
      # intersect this cell with [yL,yH]
      y0 <- max(yy0, yL)
      y1 <- min(yy1, yH)
      if (y0 >= y1) next

      zstat <- eval_cell_stats(x0, x1, y0, y1)
      cells[[length(cells) + 1L]] <- c(
        x0 = x0, y0 = y0,
        hx = x1 - x0, hy = y1 - y0,
        zstat
      )
    }
  }

  if (!length(cells)) {
    warning("No valid cells for these parameters.", call. = FALSE)
    return(list(
      cells     = tibble::tibble(),
      sum_lower = 0,
      sum_upper = 0,
      sum_mean  = 0,
      fig       = NULL
    ))
  }

  cells <- tibble::as_tibble(do.call(rbind, cells))

  # --- sums
  area <- cells$hx * cells$hy
  sum_lower <- sum(cells$z_lower * area)
  sum_upper <- sum(cells$z_upper * area)
  sum_mean  <- sum(cells$z_mean  * area)

  # --- plotting
  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Plotting requires 'plotly' installed.", call. = FALSE)
    } else {
      cs_top <- as_colorscale(top_colorscale, alpha = top_opacity)

      add_rect_surface <- function(fig, x0, x1, y0, y1, z, colorscale, opacity) {
        X <- rbind(c(x0, x1), c(x0, x1))
        Y <- rbind(c(y0, y0), c(y1, y1))
        Z <- matrix(z, nrow = 2L, ncol = 2L)
        fig |>
          plotly::add_surface(
            x = X, y = Y, z = Z,
            colorscale = colorscale,
            showscale  = FALSE,
            opacity    = opacity
          )
      }
      add_side <- function(fig, x, y, z0, z1, horiz = c("x","y")) {
        horiz <- match.arg(horiz)
        if (horiz == "x") {
          X <- rbind(c(x[1], x[2]), c(x[1], x[2]))
          Y <- rbind(c(y, y),       c(y, y))
        } else {
          X <- rbind(c(x, x),       c(x, x))
          Y <- rbind(c(y[1], y[2]), c(y[1], y[2]))
        }
        Z <- rbind(c(z0, z0), c(z1, z1))
        plotly::add_surface(
          fig, x = X, y = Y, z = Z,
          colorscale = list(list(0, side_color), list(1, side_color)),
          showscale  = FALSE,
          opacity    = side_opacity
        )
      }
      add_edge <- function(fig, xs, ys, zs) {
        if (is.na(frame_color)) return(fig)
        plotly::add_trace(
          fig, x = xs, y = ys, z = zs,
          type = "scatter3d", mode = "lines",
          line = list(color = frame_color, width = frame_width),
          hoverinfo = "none", showlegend = FALSE
        )
      }

      fig <- plotly::plot_ly()

      # optional true surface over rectangular bbox [a,b] x [cmin,dmax]
      if (isTRUE(show_surface)) {
        nx <- max(40L, N * 6L)
        ny <- max(40L, M * 6L)
        xs <- seq(a, b, length.out = nx)
        ys <- seq(cmin, dmax, length.out = ny)
        Xs <- matrix(rep(xs, each = ny), nrow = ny)
        Ys <- matrix(rep(ys, times = nx), nrow = ny)
        Zs <- matrix(NA_real_, nrow = ny, ncol = nx)
        for (j in seq_len(nx)) {
          Zs[, j] <- vapply(ys, function(yy) f(xs[j], yy), numeric(1))
        }

        contours_arg <- if (isTRUE(show_surface_grid)) list(
          x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
          y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
          z = list(show = FALSE)
        ) else NULL

        fig <- fig |>
          plotly::add_surface(
            x = Xs, y = Ys, z = Zs,
            colorscale = as_colorscale(surface_colorscale),
            showscale  = FALSE,
            opacity    = surface_opacity,
            contours   = contours_arg
          )
      }

      # draw prisms (top + 4 sides + edges)
      z_col <- switch(
        color_by,
        lower = cells$z_lower,
        upper = cells$z_upper,
        mean  = cells$z_mean
      )

      for (i in seq_len(nrow(cells))) {
        x0 <- cells$x0[i]; x1 <- x0 + cells$hx[i]
        y0 <- cells$y0[i]; y1 <- y0 + cells$hy[i]
        zt <- z_col[i]

        # top (flat)
        fig <- add_rect_surface(fig, x0, x1, y0, y1, zt, cs_top, top_opacity)
        # sides (vertical rectangles)
        fig <- add_side(fig, x = c(x0, x1), y = y0, z0 = 0, z1 = zt, horiz = "x") # front
        fig <- add_side(fig, x = c(x0, x1), y = y1, z0 = 0, z1 = zt, horiz = "x") # back
        fig <- add_side(fig, x = x0,       y = c(y0, y1), z0 = 0, z1 = zt, horiz = "y") # left
        fig <- add_side(fig, x = x1,       y = c(y0, y1), z0 = 0, z1 = zt, horiz = "y") # right

        # edges (bottom rectangle @ z=0 + top rectangle @ z=zt + 4 verticals)
        fig <- add_edge(fig, c(x0, x1), c(y0, y0), c(0, 0))
        fig <- add_edge(fig, c(x0, x1), c(y1, y1), c(0, 0))
        fig <- add_edge(fig, c(x0, x0), c(y0, y1), c(0, 0))
        fig <- add_edge(fig, c(x1, x1), c(y0, y1), c(0, 0))

        fig <- add_edge(fig, c(x0, x1), c(y0, y0), c(zt, zt))
        fig <- add_edge(fig, c(x0, x1), c(y1, y1), c(zt, zt))
        fig <- add_edge(fig, c(x0, x0), c(y0, y1), c(zt, zt))
        fig <- add_edge(fig, c(x1, x1), c(y0, y1), c(zt, zt))

        fig <- add_edge(fig, c(x0, x0), c(y0, y0), c(0, zt))
        fig <- add_edge(fig, c(x1, x1), c(y0, y0), c(0, zt))
        fig <- add_edge(fig, c(x0, x0), c(y1, y1), c(0, zt))
        fig <- add_edge(fig, c(x1, x1), c(y1, y1), c(0, zt))
      }

      fig <- fig |>
        plotly::layout(
          title = "Riemann sum with rectangular prisms",
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )
      print(fig)
    }
  }

  out <- list(
    cells     = cells,
    sum_lower = sum_lower,
    sum_upper = sum_upper,
    sum_mean  = sum_mean,
    fig       = fig
  )
  if (estimate != "all") {
    key <- switch(estimate,
                  lower = "sum_lower",
                  upper = "sum_upper",
                  mean  = "sum_mean")
    out$estimate <- out[[key]]
  }
  out
}
