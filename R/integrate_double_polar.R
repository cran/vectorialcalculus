#' @title Numerical Double Integration in Polar Coordinates
#'
#' @description
#' Calculates the definite double integral of a function f(x, y) over a
#' polar region R, defined by constant limits for theta and functional limits
#' for the radius r, R: \eqn{g1(theta) <= r <= g2(theta)}, \eqn{alpha <= theta <= beta}.
#' The integration order used is r dr dtheta.
#' Uses the Composite Simpson's Rule for numerical approximation.
#'
#' @param f A function R of two variables, \code{f(x, y)}, returning a numeric value (the original integrand).
#' @param theta_min The constant lower limit for the outer integral (alpha).
#' @param theta_max The constant upper limit for the outer integral (beta).
#' @param r_limit1 A function R of one variable defining the inner integral's lower limit (r = g1(theta)).
#' @param r_limit2 A function R of one variable defining the inner integral's upper limit (r = g2(theta)).
#' @param n_theta Number of subintervals for the outer integration (theta). Must be even. Default is 100.
#' @param n_r Number of subintervals for the inner integration (r). Must be even. Default is 100.
#' @param plot_domain Logical. If TRUE, generates a ggplot2 plot of the integration domain in the Cartesian plane. Default is TRUE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{integral_value}: The calculated approximation of the integral.
#'     \item \code{domain_plot}: The ggplot2 object representing the domain (if plot_domain = TRUE).
#'   }
#'
#' @export
integrate_double_polar <- function(f, theta_min, theta_max, r_limit1, r_limit2,
                                   n_theta = 100, n_r = 100, plot_domain = TRUE) {

  # --- 1. VALIDATION AND SETUP ---
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("The 'purrr' package is required. Please install it.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("The 'ggplot2' package is required for plotting, but not installed. Plotting disabled.")
    plot_domain <- FALSE
  }
  if (theta_max <= theta_min) stop("'theta_max' must be greater than 'theta_min'.")
  if (n_theta %% 2 != 0 || n_r %% 2 != 0) stop("'n_theta' and 'n_r' must be even for Simpson's Rule.")

  # Pre-calculate Simpson's coefficients
  coef_simp_inner <- c(1, rep(c(4, 2), n_r / 2), 1)
  coef_simp_outer <- c(1, rep(c(4, 2), n_theta / 2), 1)

  # --- 2. INTEGRATION CORE LOGIC (Simpson's Rule) ---

  # Function for the INNER INTEGRAL (w.r.t r)
  integral_inner_r <- function(val_theta, n_r) {
    r_min <- r_limit1(val_theta)
    r_max <- r_limit2(val_theta)

    if (r_max <= r_min) {
      return(0)
    }

    # Discretization for the inner integral (r)
    h_r <- (r_max - r_min) / n_r
    r_vals <- seq(r_min, r_max, length.out = n_r + 1)

    # Define the new integrand g(r, theta) = f(r cos(theta), r sin(theta)) * r
    g_vals <- purrr::map_dbl(r_vals, function(r_val) {
      # Convert polar (r, theta) to Cartesian (x, y)
      x <- r_val * cos(val_theta)
      y <- r_val * sin(val_theta)
      # Calculate the value of the new integrand: f(x, y) * r
      return(f(x, y) * r_val)
    })

    # Apply Simpson's Rule for the inner integral (r)
    result_inner <- sum(g_vals * coef_simp_inner) * (h_r / 3)
    return(result_inner)
  }

  # --- 3. OUTER INTEGRATION (Aggregation, w.r.t theta) ---
  h_theta <- (theta_max - theta_min) / n_theta
  theta_vals <- seq(theta_min, theta_max, length.out = n_theta + 1)

  # Calculate the inner integral result for every point in the outer integral
  integral_outer_results <- purrr::map_dbl(theta_vals, integral_inner_r, n_r = n_r)

  # Apply Simpson's Rule for the outer integral (theta)
  integral_final <- sum(integral_outer_results * coef_simp_outer) * (h_theta / 3)

  # --- 4. DOMAIN PLOTTING (ggplot2) ---
  domain_plot <- NULL
  if (plot_domain) {
    # Generate points for plotting
    theta_seq <- seq(theta_min, theta_max, length.out = 100)

    # Generate the boundary points for the outer limit (r_max)
    r2 <- r_limit2(theta_seq)
    x2 <- r2 * cos(theta_seq)
    y2 <- r2 * sin(theta_seq)

    # Generate the boundary points for the inner limit (r_min)
    r1 <- r_limit1(rev(theta_seq)) # Reverse order for plotting polygon boundary
    x1 <- r1 * cos(rev(theta_seq))
    y1 <- r1 * sin(rev(theta_seq))

    # Combine points to form a closed polygon for the region D
    data_plot <- data.frame(
      x = c(x2, x1),
      y = c(y2, y1)
    )

    domain_plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_polygon(fill = "gold", alpha = 0.6, color = "darkgoldenrod4") +
      # Ensure aspect ratio is 1:1 for correct geometric representation
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::labs(
        title = "Domain Plot in Cartesian Coordinates (Polar Region)",
        x = "x = r cos(theta)",
        y = "y = r sin(theta)"
      ) + ggplot2::theme_minimal()
  }

  # --- 5. RETURN RESULTS ---
  return(list(
    integral_value = integral_final,
    domain_plot = domain_plot
  ))
}
# -------------------------------------------------------------------
# (3) BLOQUE DE VARIABLES GLOBALES (¡Agregado aquí!)
# Debe estar fuera de cualquier función.

# Declaración de variables globales para satisfacer el requisito de CRAN
if (getRversion() >= "2.15.1") utils::globalVariables(c(
  "x", "y", "outer_var", "inner_min", "inner_max", "x_min", "x_max", "y_min", "y_max"
))
