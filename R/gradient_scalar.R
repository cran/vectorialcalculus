#' Gradient of a scalar field in R^n
#'
#' Computes a numerical approximation of the gradient of a scalar function
#' at a given point using central finite differences. The function \code{f}
#' is assumed to take a numeric vector as input and return a scalar.
#'
#' Optionally, if the input point has length 2 or 3 and \code{plot = TRUE},
#' a simple visualization of the gradient vector is produced using
#' \pkg{plotly}.
#'
#' @param f Function of a numeric vector \code{f(x)} returning a numeric scalar.
#' @param x0 Numeric vector giving the evaluation point.
#' @param h Numeric step size for finite differences. Can be:
#'   \itemize{
#'     \item \code{NULL} (default): a step is chosen as \code{1e-4 * (1 + abs(x0))} for each component;
#'     \item A scalar, used for all components;
#'     \item A numeric vector of the same length as \code{x0}.
#'   }
#' @param plot Logical; if \code{TRUE} and \code{length(x0)} is 2 or 3,
#'   draws the gradient vector with \pkg{plotly}.
#'
#' @return A numeric vector of the same length as \code{x0} with the
#'   components of the gradient.
#'
#' @examples
#' f <- function(v) exp(-(v[1]^2 + v[2]^2)) + 0.3 * sin(2 * v[1] * v[2])
#' gradient_scalar(f, c(0.6, -0.4))
#'
#' @export
gradient_scalar <- function(f, x0, h = NULL, plot = FALSE) {
  if (!is.function(f)) {
    stop("'f' must be a function(x).", call. = FALSE)
  }
  if (!is.numeric(x0) || any(!is.finite(x0))) {
    stop("'x0' must be a finite numeric vector.", call. = FALSE)
  }

  x0 <- as.numeric(x0)
  n  <- length(x0)
  if (n < 1L) {
    stop("'x0' must have length at least 1.", call. = FALSE)
  }

  # step sizes
  if (is.null(h)) {
    h <- 1e-4 * (1 + abs(x0))
  } else {
    if (!is.numeric(h) || any(!is.finite(h))) {
      stop("'h' must be numeric and finite.", call. = FALSE)
    }
    if (length(h) == 1L) {
      h <- rep(as.numeric(h), n)
    } else if (length(h) != n) {
      stop("'h' must be scalar, vector of length n, or NULL.", call. = FALSE)
    }
  }
  if (any(h <= 0)) {
    stop("Step sizes 'h' must be positive.", call. = FALSE)
  }

  # robust evaluator: f(x) must return finite scalar
  eval_f <- function(xx) {
    xx <- as.numeric(xx)
    val <- f(xx)
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
      stop("f(x) must return a finite numeric scalar.", call. = FALSE)
    }
    val
  }

  grad <- numeric(n)
  for (i in seq_len(n)) {
    ei <- rep(0, n)
    ei[i] <- h[i]
    fp <- eval_f(x0 + ei)
    fm <- eval_f(x0 - ei)
    grad[i] <- (fp - fm) / (2 * h[i])
  }

  # Optional plot
  if (isTRUE(plot) && (n %in% c(2L, 3L))) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Package 'plotly' is required for plotting.", call. = FALSE)
    } else {
      if (n == 2L) {
        p <- plotly::plot_ly() |>
          plotly::add_trace(
            x = x0[1], y = x0[2],
            type = "scatter", mode = "markers",
            marker = list(color = "blue", size = 6),
            name = "Point"
          ) |>
          plotly::add_trace(
            x = c(x0[1], x0[1] + grad[1]),
            y = c(x0[2], x0[2] + grad[2]),
            type = "scatter", mode = "lines+markers",
            line = list(color = "red", width = 3),
            marker = list(size = 4, color = "red"),
            name = "Gradient"
          ) |>
          plotly::layout(
            title = "Gradient in 2D",
            xaxis = list(title = "x"),
            yaxis = list(title = "y"),
            showlegend = TRUE
          )
        print(p)
      }
      if (n == 3L) {
        p <- plotly::plot_ly() |>
          plotly::add_markers(
            x = x0[1], y = x0[2], z = x0[3],
            marker = list(color = "blue", size = 4),
            name = "Point"
          ) |>
          plotly::add_trace(
            x = c(x0[1], x0[1] + grad[1]),
            y = c(x0[2], x0[2] + grad[2]),
            z = c(x0[3], x0[3] + grad[3]),
            type = "scatter3d", mode = "lines+markers",
            line = list(color = "red", width = 6),
            marker = list(size = 3, color = "red"),
            name = "Gradient"
          ) |>
          plotly::layout(
            title = "Gradient in 3D",
            scene = list(
              xaxis = list(title = "x"),
              yaxis = list(title = "y"),
              zaxis = list(title = "z")
            )
          )
        print(p)
      }
    }
  }

  grad
}
