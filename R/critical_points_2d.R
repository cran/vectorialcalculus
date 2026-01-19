#' Critical points of a two-variable function using gradient and Hessian
#'
#' Identifies stationary points of a function of two variables over a given
#' rectangular domain. A set of initial points is generated on a regular grid
#' together with additional random starting points. Each start is refined using
#' numerical optimization applied to the squared gradient norm. Points with a
#' sufficiently small gradient norm are kept and then merged if they are too
#' close to each other.
#'
#' Each surviving point is classified based on the eigenvalues of the numerical
#' Hessian matrix. The Hessian classification uses four categories:
#' \itemize{
#'   \item \code{"minimum"}  - both eigenvalues positive,
#'   \item \code{"maximum"}  - both eigenvalues negative,
#'   \item \code{"saddle"}   - mixed signs,
#'   \item \code{"flat"}     - small eigenvalues, inconclusive classification.
#' }
#'
#' Optionally, a 3D surface with the detected critical points can be displayed
#' using \pkg{plotly}.
#'
#' @param f Function of two variables \code{f(x, y)} returning a numeric scalar.
#' @param xlim Numeric vector \code{c(xmin, xmax)} defining the domain in the
#'   x-direction.
#' @param ylim Numeric vector \code{c(ymin, ymax)} defining the domain in the
#'   y-direction.
#' @param start_n Integer vector of length two. Number of regular starting
#'   points per axis for the grid.
#' @param n_rand Integer. Number of additional random starting points inside
#'   the domain.
#' @param h Numeric step size for finite-difference gradients and Hessians.
#'   If \code{NULL}, an automatic step size is used.
#' @param tol_grad Numeric tolerance. Points with gradient norm below this
#'   threshold are treated as stationary.
#' @param tol_merge Numeric tolerance for merging nearby solutions.
#' @param tol_eig Numeric tolerance used when classifying eigenvalues of the
#'   Hessian.
#' @param maxit Maximum number of iterations permitted for numerical
#'   optimization.
#' @param optim_method Character string naming an optimization method supported
#'   by \code{stats::optim}, such as \code{"BFGS"} or \code{"Nelder-Mead"}.
#' @param plot Logical. If \code{TRUE}, a \pkg{plotly} surface with critical
#'   points is drawn.
#' @param grid_plot Integer vector of length two. Resolution of the grid used
#'   when drawing the surface.
#' @param surface_colorscale Character. Name of the Plotly colorscale for the
#'   surface.
#' @param surface_opacity Numeric between 0 and 1 giving the opacity of the
#'   surface.
#' @param cp_colors Named list mapping each critical point type to a color.
#'   Expected names: \code{"minimum"}, \code{"maximum"}, \code{"saddle"},
#'   \code{"flat"}.
#' @param cp_size Numeric. Size of the point markers.
#' @param scene List of Plotly scene options (axis titles, aspect mode, and so
#'   on).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{critical_points}}{
#'       A data frame with columns \code{x}, \code{y}, \code{z}, the gradient
#'       norm, and the classification label.
#'   }
#'   \item{\code{fig}}{
#'       A \pkg{plotly} object if \code{plot = TRUE}, or \code{NULL} otherwise.
#'   }
#' }
#'
#' @examples
#' \donttest{
#' f <- function(x, y) x^2 + y^2
#' critical_points_2d(f, xlim = c(-2, 2), ylim = c(-2, 2), plot = FALSE)
#' }
#'
#' @export
critical_points_2d <- function(
    f,
    xlim, ylim,
    start_n = c(7L, 7L),
    n_rand = 40L,
    h = NULL,
    tol_grad = 1e-6,
    tol_merge = 1e-3,
    tol_eig = 1e-6,
    maxit = 200,
    optim_method = c("BFGS", "Nelder-Mead"),
    plot = TRUE,
    grid_plot = c(60L, 60L),
    surface_colorscale = "YlGnBu",
    surface_opacity = 0.85,
    cp_colors = list(
      minimum = "#2ca02c",
      maximum = "#d62728",
      saddle  = "#1f77b4",
      flat    = "#ff7f0e"
    ),
    cp_size = 6,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "f(x,y)")
    )
) {
  x <- y <- grad_norm <- NULL
  optim_method <- match.arg(optim_method)
  stopifnot(is.function(f))
  if (length(xlim) != 2 || length(ylim) != 2 ||
      xlim[2] <= xlim[1] || ylim[2] <= ylim[1]) {
    stop("'xlim' and 'ylim' must be c(min,max) with min < max.", call. = FALSE)
  }

  # ---- steps per coordinate
  .steps <- function(x, y) {
    if (is.null(h)) {
      c(1e-4 * (1 + abs(x)), 1e-4 * (1 + abs(y)))
    } else if (length(h) == 1L) {
      rep(as.numeric(h), 2L)
    } else {
      as.numeric(h[1:2])
    }
  }

  # central gradient (order 2)
  .grad <- function(x, y) {
    st <- .steps(x, y); hx <- st[1]; hy <- st[2]
    fxp <- f(x + hx, y); fxm <- f(x - hx, y)
    fyp <- f(x, y + hy); fym <- f(x, y - hy)
    c((fxp - fxm) / (2 * hx), (fyp - fym) / (2 * hy))
  }

  # central Hessian (order 2)
  .hess <- function(x, y) {
    st <- .steps(x, y); hx <- st[1]; hy <- st[2]
    fv  <- f(x, y)
    fxx <- (f(x + hx, y) - 2 * fv + f(x - hx, y)) / (hx^2)
    fyy <- (f(x, y + hy) - 2 * fv + f(x, y - hy)) / (hy^2)
    fxy <- (f(x + hx, y + hy) - f(x + hx, y - hy) -
              f(x - hx, y + hy) + f(x - hx, y - hy)) / (4 * hx * hy)
    matrix(c(fxx, fxy, fxy, fyy), 2, 2)
  }

  # ||grad||^2 objective
  .gnorm2 <- function(p) {
    g <- .grad(p[1], p[2])
    sum(g * g)
  }

  # classification from Hessian eigenvalues
  .classify <- function(H) {
    ev  <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    pos <- ev >  tol_eig
    neg <- ev < -tol_eig
    if (all(pos)) {
      "minimum"
    } else if (all(neg)) {
      "maximum"
    } else if (any(pos) && any(neg)) {
      "saddle"
    } else {
      "flat"
    }
  }

  # merge nearby points
  .merge_pts <- function(P, r) {
    if (nrow(P) <= 1) return(P)
    keep <- rep(TRUE, nrow(P))
    for (i in seq_len(nrow(P))) {
      if (!keep[i]) next
      di <- sqrt((P$x[i] - P$x)^2 + (P$y[i] - P$y)^2)
      dup <- which(di < r)
      dup <- dup[dup > i]
      keep[dup] <- FALSE
    }
    P[keep, , drop = FALSE]
  }

  # clamp to domain
  .clip <- function(p) {
    c(
      min(max(p[1], xlim[1]), xlim[2]),
      min(max(p[2], ylim[1]), ylim[2])
    )
  }

  # safe optim with fallback
  .safe_optim <- function(par0) {
    par0 <- .clip(par0)
    fn   <- function(p) .gnorm2(p)
    ctrl <- list(maxit = maxit)
    res  <- try(
      stats::optim(par0, fn, method = optim_method, control = ctrl),
      silent = TRUE
    )
    if (inherits(res, "try-error")) {
      res <- stats::optim(par0, fn, method = "Nelder-Mead", control = ctrl)
    }
    res
  }

  # ---- multi-start (grid + random)
  xs <- seq(xlim[1], xlim[2], length.out = start_n[1])
  ys <- seq(ylim[1], ylim[2], length.out = start_n[2])
  grid_starts <- as.data.frame(expand.grid(x = xs, y = ys))
  rand_starts <- data.frame(
    x = stats::runif(n_rand, min = xlim[1], max = xlim[2]),
    y = stats::runif(n_rand, min = ylim[1], max = ylim[2])
  )
  starts <- rbind(grid_starts, rand_starts)

  # ---- minimize g from each start
  sols <- lapply(
    seq_len(nrow(starts)),
    function(i) .safe_optim(as.numeric(starts[i, ]))
  )
  pars <- do.call(rbind, lapply(sols, function(s) s$par))
  vals <- vapply(sols, function(s) s$value, numeric(1))
  df   <- data.frame(x = pars[, 1], y = pars[, 2], g = vals)

  # filter by grad norm and domain
  df$grad_norm <- sqrt(df$g)
  df <- subset(
    df,
    is.finite(grad_norm) & grad_norm <= tol_grad &
      x >= xlim[1] & x <= xlim[2] &
      y >= ylim[1] & y <= ylim[2]
  )

  # order and merge duplicates
  if (nrow(df) > 0) {
    df <- df[order(df$grad_norm), , drop = FALSE]
    df <- .merge_pts(df[, c("x", "y", "grad_norm")], r = tol_merge)
  }

  # evaluate f and classify
  if (nrow(df) > 0) {
    df$z <- mapply(f, df$x, df$y)
    cls  <- character(nrow(df))
    for (i in seq_len(nrow(df))) {
      H       <- .hess(df$x[i], df$y[i])
      cls[i]  <- .classify(H)
    }
    df$class <- cls
    df <- df[, c("x", "y", "z", "grad_norm", "class")]
  } else {
    df <- data.frame(
      x         = numeric(0),
      y         = numeric(0),
      z         = numeric(0),
      grad_norm = numeric(0),
      class     = character(0)
    )
  }

  # ---- plot
  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Package 'plotly' is required for plotting.", call. = FALSE)
    } else {
      gx <- seq(xlim[1], xlim[2], length.out = grid_plot[1])
      gy <- seq(ylim[1], ylim[2], length.out = grid_plot[2])
      # plotly::add_surface expects z[y, x]
      Z <- outer(gy, gx, function(yy, xx) f(xx, yy))

      fig <- plotly::plot_ly() |>
        plotly::add_surface(
          x = gx, y = gy, z = Z,
          colorscale = surface_colorscale,
          opacity = surface_opacity,
          showscale = FALSE
        ) |>
        plotly::layout(
          title = "Surface f(x,y) with critical points",
          scene = scene
        )

      add_cls <- function(p, cls_name, sym) {
        sub <- df[df$class == cls_name, , drop = FALSE]
        if (nrow(sub) == 0) return(p)
        plotly::add_markers(
          p,
          x = sub$x, y = sub$y, z = sub$z,
          type = "scatter3d", mode = "markers",
          marker = list(
            size  = cp_size,
            color = if (!is.null(cp_colors[[cls_name]])) {
              cp_colors[[cls_name]]
            } else {
              "black"
            },
            line  = list(width = 1, color = "black")
          ),
          symbol = sym,
          name = switch(
            cls_name,
            minimum = "minimum",
            maximum = "maximum",
            saddle  = "saddle",
            flat    = "flat / indeterminate"
          ),
          hovertemplate = paste0(
            "<b>", toupper(cls_name), "</b><br>",
            "x=%{x:.4f}<br>y=%{y:.4f}<br>f=%{z:.6f}<br>",
            "||grad||=%{text:.2e}<extra></extra>"
          ),
          text = sub$grad_norm,
          showlegend = TRUE
        )
      }

      fig <- add_cls(fig, "minimum", "circle")
      fig <- add_cls(fig, "maximum", "diamond")
      fig <- add_cls(fig, "saddle",  "cross")
      fig <- add_cls(fig, "flat",    "square")

      print(fig)
    }
  }

  list(critical_points = df, fig = fig)
}
