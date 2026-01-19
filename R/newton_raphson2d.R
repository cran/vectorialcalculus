#' Newton-Raphson method for systems in R^2 with animation (Plotly)
#'
#' Applies the Newton-Raphson method to solve a nonlinear system in two
#' variables:
#' \deqn{
#'   \mathbf{F}(x,y)=\mathbf{0},
#' }
#' where
#' \deqn{
#'   \mathbf{F}(x,y)=\begin{pmatrix}f_1(x,y)\\ f_2(x,y)\end{pmatrix}.
#' }
#' At an iterate
#' \deqn{
#'   \mathbf{x}_n=(x_n,y_n),
#' }
#' the Newton update is
#' \deqn{
#'   \mathbf{x}_{n+1}=\mathbf{x}_n - J(\mathbf{x}_n)^{-1}\mathbf{F}(\mathbf{x}_n),
#' }
#' where the Jacobian matrix is
#' \deqn{
#'   J(x,y)=\begin{pmatrix}
#'   \frac{\partial f_1}{\partial x}(x,y) & \frac{\partial f_1}{\partial y}(x,y)\\
#'   \frac{\partial f_2}{\partial x}(x,y) & \frac{\partial f_2}{\partial y}(x,y)
#'   \end{pmatrix}.
#' }
#' If a Jacobian function is not provided, partial derivatives are approximated
#' numerically by central differences.
#'
#' The animation shows the two zero level curves
#' \deqn{
#'   f_1(x,y)=0 \quad \text{and} \quad f_2(x,y)=0
#' }
#' together with the Newton iterates and the step segments
#' \deqn{
#'   \mathbf{x}_n \to \mathbf{x}_{n+1}.
#' }
#'
#' @param f Function. A function f(x,y) that returns a numeric vector of length 2:
#'   c(f1(x,y), f2(x,y)).
#' @param x0 Numeric vector of length 2. Initial guess c(x0, y0).
#' @param J Optional function. A function J(x,y) returning a 2x2 numeric matrix
#'   (the Jacobian). If NULL, a numerical Jacobian is used.
#' @param h Numeric scalar. Step size for numerical partial derivatives (when J is NULL).
#' @param max_iter Integer. Maximum number of Newton iterations.
#' @param tol Numeric scalar. Stopping tolerance based on the infinity norm:
#'   max(|f1|,|f2|) <= tol.
#' @param xlim Numeric vector of length 2. Plot range for x. If NULL, it is chosen
#'   from the iterates.
#' @param ylim Numeric vector of length 2. Plot range for y. If NULL, it is chosen
#'   from the iterates.
#' @param n_grid Integer. Grid size per axis used to draw the zero level curves.
#' @param frame_ms Integer. Frame duration in milliseconds.
#' @param transition_ms Integer. Transition duration in milliseconds.
#' @param title Character. Plot title. If NULL, a default title is used.
#' @param safe_mode Logical. If TRUE, uses calmer animation defaults intended to
#'   reduce flicker and visual stress.
#'
#' @return A list with components:
#' \describe{
#' \item{plot}{A plotly object (htmlwidget) with animation frames.}
#' \item{iterates}{Data frame of iterates (n, x, y, f1, f2, norm_inf).}
#' \item{root}{Numeric vector c(x, y) with the last iterate.}
#' \item{converged}{Logical. TRUE if convergence was detected within max_iter.}
#' }
#'
#' @examples
#' \donttest{
#' library(plotly)
#'
#' # Example system:
#' # f1(x,y)=x^2+y^2-1 (unit circle)
#' # f2(x,y)=x-y (line)
#' f <- function(x, y) c(x^2 + y^2 - 1, x - y)
#'
#' out <- newton_raphson2d(f, x0 = c(0.8, 0.2), max_iter = 8)
#' out$plot
#' out$iterates
#' out$converged
#' out$root
#' }
#'
#' @export
newton_raphson2d <- function(
    f,
    x0,
    J = NULL,
    h = 1e-4,
    max_iter = 10L,
    tol = 1e-8,
    xlim = NULL,
    ylim = NULL,
    n_grid = 120L,
    frame_ms = 700L,
    transition_ms = 450L,
    title = NULL,
    safe_mode = TRUE
) {

  if (!is.function(f)) stop("'f' must be a function.")
  if (!is.numeric(x0) || length(x0) != 2L) stop("'x0' must be a numeric vector of length 2.")
  if (!is.null(J) && !is.function(J)) stop("'J' must be NULL or a function.")
  if (!is.numeric(h) || length(h) != 1L || h <= 0) stop("'h' must be a positive numeric scalar.")
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter < 1) stop("'max_iter' must be an integer >= 1.")
  if (!is.numeric(tol) || length(tol) != 1L || tol <= 0) stop("'tol' must be a positive numeric scalar.")
  if (!is.numeric(n_grid) || length(n_grid) != 1L || n_grid < 50) stop("'n_grid' must be an integer >= 50.")

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Please install it.")
  }

  max_iter <- as.integer(max_iter)
  n_grid <- as.integer(n_grid)

  if (isTRUE(safe_mode)) {
    frame_ms <- max(as.integer(frame_ms), 600L)
    transition_ms <- max(as.integer(transition_ms), 350L)
  }

  # Helper: check f output
  f_eval <- function(x, y) {
    v <- f(x, y)
    if (!is.numeric(v) || length(v) != 2L || any(!is.finite(v))) {
      stop("'f' must return a finite numeric vector of length 2.")
    }
    v
  }

  # Numerical Jacobian if J is NULL: central differences
  J_num <- function(x, y) {
    f_xph <- f_eval(x + h, y)
    f_xmh <- f_eval(x - h, y)
    f_yph <- f_eval(x, y + h)
    f_ymh <- f_eval(x, y - h)
    dfdx <- (f_xph - f_xmh) / (2*h)
    dfdy <- (f_yph - f_ymh) / (2*h)
    matrix(c(dfdx[1], dfdy[1], dfdx[2], dfdy[2]), nrow = 2, byrow = TRUE)
  }

  J_eval <- function(x, y) {
    if (is.null(J)) {
      J_num(x, y)
    } else {
      A <- J(x, y)
      if (!is.matrix(A) || any(dim(A) != c(2L, 2L)) || any(!is.finite(A))) {
        stop("'J' must return a finite 2x2 numeric matrix.")
      }
      A
    }
  }

  # Iteration storage
  xs <- numeric(max_iter + 1L)
  ys <- numeric(max_iter + 1L)
  f1s <- numeric(max_iter + 1L)
  f2s <- numeric(max_iter + 1L)
  norms <- numeric(max_iter + 1L)

  xs[1] <- x0[1]
  ys[1] <- x0[2]
  F0 <- f_eval(xs[1], ys[1])
  f1s[1] <- F0[1]
  f2s[1] <- F0[2]
  norms[1] <- max(abs(F0))

  converged <- is.finite(norms[1]) && norms[1] <= tol
  n_stop <- 0L

  if (!converged) {
    for (k in 1:max_iter) {
      A <- J_eval(xs[k], ys[k])
      Fk <- c(f1s[k], f2s[k])

      # Solve A * s = Fk, then x_{k+1} = x_k - s
      s <- tryCatch(
        solve(A, Fk),
        error = function(e) NA_real_
      )
      if (length(s) != 2L || any(!is.finite(s))) {
        n_stop <- k - 1L
        break
      }

      xs[k + 1] <- xs[k] - s[1]
      ys[k + 1] <- ys[k] - s[2]

      F_next <- f_eval(xs[k + 1], ys[k + 1])
      f1s[k + 1] <- F_next[1]
      f2s[k + 1] <- F_next[2]
      norms[k + 1] <- max(abs(F_next))

      n_stop <- k
      if (is.finite(norms[k + 1]) && norms[k + 1] <= tol) {
        converged <- TRUE
        break
      }
    }
  }

  iter_df <- data.frame(
    n = 0:n_stop,
    x = xs[1:(n_stop + 1L)],
    y = ys[1:(n_stop + 1L)],
    f1 = f1s[1:(n_stop + 1L)],
    f2 = f2s[1:(n_stop + 1L)],
    norm_inf = norms[1:(n_stop + 1L)]
  )

  # Choose plot ranges if NULL
  if (is.null(xlim) || is.null(ylim)) {
    xr <- range(iter_df$x, finite = TRUE)
    yr <- range(iter_df$y, finite = TRUE)
    sx <- xr[2] - xr[1]; if (!is.finite(sx) || sx == 0) sx <- 1
    sy <- yr[2] - yr[1]; if (!is.finite(sy) || sy == 0) sy <- 1
    if (is.null(xlim)) xlim <- c(xr[1] - 0.9*sx, xr[2] + 0.9*sx)
    if (is.null(ylim)) ylim <- c(yr[1] - 0.9*sy, yr[2] + 0.9*sy)
  }

  if (!is.numeric(xlim) || length(xlim) != 2L || !(xlim[2] > xlim[1])) stop("'xlim' must be length 2 with xlim[2] > xlim[1].")
  if (!is.numeric(ylim) || length(ylim) != 2L || !(ylim[2] > ylim[1])) stop("'ylim' must be length 2 with ylim[2] > ylim[1].")

  # Grid for zero level curves
  xg <- seq(xlim[1], xlim[2], length.out = n_grid)
  yg <- seq(ylim[1], ylim[2], length.out = n_grid)

  z1 <- outer(xg, yg, Vectorize(function(a, b) f_eval(a, b)[1]))
  z2 <- outer(xg, yg, Vectorize(function(a, b) f_eval(a, b)[2]))

  # Frames
  frame_id <- sprintf("it%02d", iter_df$n + 1L)

  # Animated points (iterate)
  df_pt <- data.frame(
    x = iter_df$x,
    y = iter_df$y,
    frame = frame_id,
    txt = paste0(
      "n = ", iter_df$n,
      "<br>(x_n, y_n) = (", formatC(iter_df$x, digits = 8, format = "f"),
      ", ", formatC(iter_df$y, digits = 8, format = "f"), ")",
      "<br>norm = ", formatC(iter_df$norm_inf, digits = 3, format = "e")
    )
  )

  # Step segment x_n -> x_{n+1} (last frame repeats the point)
  df_step <- lapply(seq_len(nrow(iter_df)), function(i) {
    x1 <- iter_df$x[i]
    y1 <- iter_df$y[i]
    if (i < nrow(iter_df)) {
      x2 <- iter_df$x[i + 1L]
      y2 <- iter_df$y[i + 1L]
    } else {
      x2 <- x1; y2 <- y1
    }
    data.frame(
      x = c(x1, x2, NA_real_),
      y = c(y1, y2, NA_real_),
      frame = frame_id[i]
    )
  })
  df_step <- do.call(rbind, df_step)

  # On-screen label
  x_label <- xlim[1] + 0.03 * (xlim[2] - xlim[1])
  y_label <- ylim[2] - 0.05 * (ylim[2] - ylim[1])

  df_label <- data.frame(
    x = rep(x_label, nrow(iter_df)),
    y = rep(y_label, nrow(iter_df)),
    frame = frame_id,
    txt = paste0(
      "n = ", iter_df$n,
      "    norm = ", formatC(iter_df$norm_inf, digits = 3, format = "e"),
      if (converged) "" else ""
    )
  )

  if (is.null(title)) title <- "Newton-Raphson for systems in R^2"

  plot <- plotly::plot_ly() |>
    # f1(x,y)=0
    plotly::add_trace(
      x = xg, y = yg, z = z1,
      type = "contour",
      showscale = FALSE,
      contours = list(start = 0, end = 0, size = 1, coloring = "lines"),
      name = "f1(x,y)=0",
      hoverinfo = "skip"
    ) |>
    # f2(x,y)=0
    plotly::add_trace(
      x = xg, y = yg, z = z2,
      type = "contour",
      showscale = FALSE,
      contours = list(start = 0, end = 0, size = 1, coloring = "lines"),
      name = "f2(x,y)=0",
      hoverinfo = "skip"
    ) |>
    # step segment
    plotly::add_trace(
      data = df_step,
      x = ~x, y = ~y,
      frame = ~frame,
      type = "scatter", mode = "lines",
      name = "step",
      hoverinfo = "skip",
      line = list(width = 3)
    ) |>
    # iterate point
    plotly::add_trace(
      data = df_pt,
      x = ~x, y = ~y,
      frame = ~frame,
      type = "scatter", mode = "markers",
      name = "iterate",
      text = ~txt,
      hoverinfo = "text",
      marker = list(size = 10)
    ) |>
    # label
    plotly::add_trace(
      data = df_label,
      x = ~x, y = ~y,
      frame = ~frame,
      type = "scatter", mode = "text",
      text = ~txt,
      textposition = "top left",
      textfont = list(size = 14, color = "black"),
      hoverinfo = "skip",
      showlegend = FALSE
    ) |>
    plotly::layout(
      title = title,
      xaxis = list(title = "x", range = xlim),
      yaxis = list(title = "y", range = ylim)
    ) |>
    plotly::animation_opts(
      frame = as.integer(frame_ms),
      transition = as.integer(transition_ms),
      easing = "linear",
      redraw = FALSE
    ) |>
    plotly::animation_button(
      x = 1, xanchor = "right",
      y = 1, yanchor = "top"
    ) |>
    plotly::animation_slider(currentvalue = list(prefix = "iter: "))

  list(
    plot = plot,
    iterates = iter_df,
    root = c(iter_df$x[nrow(iter_df)], iter_df$y[nrow(iter_df)]),
    converged = converged
  )
}
