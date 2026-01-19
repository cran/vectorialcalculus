#' @title Numerical Triple Integration over a General Region
#'
#' @description
#' Calculates the definite triple integral of a function f(x, y, z) over a
#' general region W defined by the order dz dy dx.
#' The region W is defined by constant limits for the outer integral (x),
#' functional limits depending on x for the middle integral (y), and
#' functional limits depending on both x and y for the inner integral (z).
#' Uses the Composite Simpson's Rule for numerical approximation.
#'
#' @param f A function in R of three variables, \code{f(x, y, z)}, returning a numeric value.
#' @param x_min The constant lower limit for the outermost integral (x = a).
#' @param x_max The constant upper limit for the outermost integral (x = b).
#' @param y_limit1 A function in R of one variable defining the middle integral's lower limit (y = h1(x)).
#' @param y_limit2 A function in R of one variable defining the middle integral's upper limit (y = h2(x)).
#' @param z_limit1 A function in R of two variables defining the inner integral's lower limit (z = g1(x, y)).
#' @param z_limit2 A function in R of two variables defining the inner integral's upper limit (z = g2(x, y)).
#' @param n_outer Number of subintervals for the outermost integral (x). Must be even. Default is 50.
#' @param n_middle Number of subintervals for the middle integral (y). Must be even. Default is 50.
#' @param n_inner Number of subintervals for the innermost integral (z). Must be even. Default is 50.
#' @param plot_xy_domain Logical. If TRUE, generates a ggplot2 plot of the projection
#'   of the domain W onto the xy-plane. Default is TRUE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{integral_value}: The calculated approximation of the integral.
#'     \item \code{domain_plot}: The ggplot2 object representing the xy-projection domain
#'           (if \code{plot_xy_domain = TRUE}).
#'   }
#'
#' @export
integrate_triple_general <- function(
    f,
    x_min, x_max,
    y_limit1, y_limit2,
    z_limit1, z_limit2,
    n_outer = 50,
    n_middle = 50,
    n_inner = 50,
    plot_xy_domain = TRUE
) {
  # --- 1. Validation and setup ----------------------------------------------
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("The 'purrr' package is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning(
      "The 'ggplot2' package is required for plotting but is not installed. 'plot_xy_domain' set to FALSE.",
      call. = FALSE
    )
    plot_xy_domain <- FALSE
  }

  if (!is.function(f)) {
    stop("'f' must be a function of three variables.", call. = FALSE)
  }
  if (!is.function(y_limit1) || !is.function(y_limit2)) {
    stop("'y_limit1' and 'y_limit2' must be functions of one variable.", call. = FALSE)
  }
  if (!is.function(z_limit1) || !is.function(z_limit2)) {
    stop("'z_limit1' and 'z_limit2' must be functions of two variables.", call. = FALSE)
  }

  if (!is.numeric(x_min) || !is.numeric(x_max) ||
      length(x_min) != 1L || length(x_max) != 1L ||
      !is.finite(x_min) || !is.finite(x_max)) {
    stop("'x_min' and 'x_max' must be finite numeric scalars.", call. = FALSE)
  }
  if (x_max <= x_min) {
    stop("'x_max' must be greater than 'x_min'.", call. = FALSE)
  }

  if (!is.numeric(n_outer) || !is.numeric(n_middle) || !is.numeric(n_inner) ||
      length(n_outer) != 1L || length(n_middle) != 1L || length(n_inner) != 1L) {
    stop("'n_outer', 'n_middle' and 'n_inner' must be numeric scalars.", call. = FALSE)
  }

  n_outer  <- as.integer(n_outer)
  n_middle <- as.integer(n_middle)
  n_inner  <- as.integer(n_inner)

  if (n_outer <= 0L || n_middle <= 0L || n_inner <= 0L ||
      n_outer %% 2L != 0L || n_middle %% 2L != 0L || n_inner %% 2L != 0L) {
    stop("All 'n' parameters must be positive even integers for Simpson's Rule.", call. = FALSE)
  }

  # Simpson coefficients for each level
  make_coef <- function(n) c(1, rep(c(4, 2), n / 2L - 1L), 4, 1)
  coef_inner  <- make_coef(n_inner)
  coef_middle <- make_coef(n_middle)
  coef_outer  <- make_coef(n_outer)

  # --- 2. Innermost integral in z -------------------------------------------
  integral_inner_z <- function(val_x, val_y, n_inner) {
    z_min <- z_limit1(val_x, val_y)
    z_max <- z_limit2(val_x, val_y)

    if (!is.numeric(z_min) || !is.numeric(z_max) ||
        length(z_min) != 1L || length(z_max) != 1L ||
        !is.finite(z_min) || !is.finite(z_max)) {
      stop("Inner limits 'z_limit1' and 'z_limit2' must return finite numeric scalars.", call. = FALSE)
    }

    if (z_max <= z_min) {
      return(0)
    }

    h_inner <- (z_max - z_min) / n_inner
    z_vals  <- seq(z_min, z_max, length.out = n_inner + 1L)

    f_vals <- purrr::map_dbl(z_vals, function(z) f(val_x, val_y, z))
    if (any(!is.finite(f_vals))) {
      stop("Function 'f' returned non-finite values in the z-integration.", call. = FALSE)
    }

    sum(f_vals * coef_inner) * (h_inner / 3)
  }

  # --- 3. Middle integral in y ----------------------------------------------
  integral_middle_y <- function(val_x, n_middle, n_inner) {
    y_min <- y_limit1(val_x)
    y_max <- y_limit2(val_x)

    if (!is.numeric(y_min) || !is.numeric(y_max) ||
        length(y_min) != 1L || length(y_max) != 1L ||
        !is.finite(y_min) || !is.finite(y_max)) {
      stop("Middle limits 'y_limit1' and 'y_limit2' must return finite numeric scalars.", call. = FALSE)
    }

    if (y_max <= y_min) {
      return(0)
    }

    h_middle <- (y_max - y_min) / n_middle
    y_vals   <- seq(y_min, y_max, length.out = n_middle + 1L)

    z_integral_results <- purrr::map_dbl(
      y_vals,
      function(y) integral_inner_z(val_x, y, n_inner = n_inner)
    )
    if (any(!is.finite(z_integral_results))) {
      stop("Non-finite values encountered in y-integration.", call. = FALSE)
    }

    sum(z_integral_results * coef_middle) * (h_middle / 3)
  }

  # --- 4. Outermost integral in x -------------------------------------------
  h_outer <- (x_max - x_min) / n_outer
  x_vals  <- seq(x_min, x_max, length.out = n_outer + 1L)

  integral_outer_results <- purrr::map_dbl(
    x_vals,
    integral_middle_y,
    n_middle = n_middle,
    n_inner  = n_inner
  )
  if (any(!is.finite(integral_outer_results))) {
    stop("Non-finite values encountered in x-integration.", call. = FALSE)
  }

  integral_final <- sum(integral_outer_results * coef_outer) * (h_outer / 3)

  # --- 5. Domain plotting (xy projection) -----------------------------------
  domain_plot <- NULL
  if (isTRUE(plot_xy_domain)) {
    x_seq <- seq(x_min, x_max, length.out = 100L)
    y_min_vec <- y_limit1(x_seq)
    y_max_vec <- y_limit2(x_seq)

    if (!is.numeric(y_min_vec) || !is.numeric(y_max_vec) ||
        length(y_min_vec) != length(x_seq) ||
        length(y_max_vec) != length(x_seq)) {
      stop("Functions 'y_limit1' and 'y_limit2' must be vectorizable over x.", call. = FALSE)
    }

    data_plot <- data.frame(
      x    = x_seq,
      y_min = y_min_vec,
      y_max = y_max_vec
    )

    domain_plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = x)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = y_min, ymax = y_max),
        fill = "darkgreen", alpha = 0.4
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = y_min),
        color = "green4", linetype = "dashed"
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = y_max),
        color = "green4", linetype = "dashed"
      ) +
      ggplot2::labs(
        title = "Domain projection (xy-plane) for triple integral",
        x = "x (outer variable)",
        y = "y (middle variable)"
      ) +
      ggplot2::theme_minimal()
  }

  # --- 6. Return -------------------------------------------------------------
  list(
    integral_value = integral_final,
    domain_plot    = domain_plot
  )
}
