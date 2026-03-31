#' Fit a single site's tilted distribution and compute EP site update
#'
#' Internal function called by [fit_ep()] for each site at each EP iteration.
#' Computes the cavity distribution by removing the current site's contribution
#' from the global approximation, fits the tilted distribution using the
#' user-supplied `fit_model` function, and returns the change in site natural
#' parameters needed to update the global approximation.
#'
#' @param site_data Data frame containing the observations for this site.
#' @param d Integer. Dimension of the shared parameter vector \eqn{\phi}.
#' @param Q_global Numeric \eqn{d \times d} matrix. Current global precision
#'   matrix.
#' @param r_global Numeric vector of length \eqn{d}. Current global precision
#'   mean.
#' @param Q_site Numeric \eqn{d \times d} matrix. Current site precision matrix.
#' @param r_site Numeric vector of length \eqn{d}. Current site precision mean.
#' @param eta Numeric in \eqn{(0, 1]}. Power EP parameter. \eqn{\eta = 1} (default)
#'   corresponds to standard EP with KL divergence minimization. Values less
#'   than one correspond to power EP with \eqn{\alpha}-divergence.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{delta_Q}{Numeric \eqn{d \times d} matrix. Change in site precision
#'     matrix. Zero matrix if site failed or was skipped.}
#'   \item{delta_r}{Numeric vector of length \eqn{d}. Change in site precision
#'     mean. Zero vector if the site failed or was skipped.}
#'   \item{status}{Character string: `"success"`, `"skipped"` (improper cavity
#'     distribution), or `"failed"` (`fit_model` returned `NULL`).}
#'   \item{fit}{The full return value from `fit_model`, or `NULL` if skipped/failed.}
#'   \item{phi_draws}{Numeric matrix of shared parameter draws from the tilted
#'     distribution (\eqn{n} rows \eqn{\times d} columns), or `NULL` if
#'     skipped/failed.}
#' }
#'
#' @seealso [fit_ep()] which calls this function iteratively for each site.
#'
#' @importFrom stats cov
#'
#' @keywords internal
fit_site <- function(site_data, d, Q_global, r_global, Q_site, r_site, eta = 1, fit_model) {
  ##### Check arguments #####
  checkmate::assert_data_frame(site_data)
  checkmate::assert_count(d, positive = TRUE)
  checkmate::assert_matrix(Q_global, mode = "numeric", nrows = d, ncols = d)
  checkmate::assert_numeric(r_global, len = d)
  checkmate::assert_matrix(Q_site, mode = "numeric", nrows = d, ncols = d)
  checkmate::assert_numeric(r_site, len = d)
  checkmate::assert_number(eta, lower = 0, upper = 1)
  checkmate::assert_function(fit_model)

  ###### Form cavity distribution #####
  Q_cavity <- Q_global - eta * Q_site
  r_cavity <- r_global - eta * r_site

  # Check that cavity distribution is proper
  eig_min <- min(eigen(Q_cavity, symmetric = TRUE, only.values = TRUE)$values)
  if (eig_min <= 1e-8) {
    message(sprintf("Improper cavity distribution (min eigenvalue = %0.4f), skipping", eig_min))
    return(list(
      delta_Q = matrix(0, d, d),
      delta_r = rep(0, d),
      status = "skipped",
      group_index = NULL,
      fit = NULL,
      phi_draws = NULL
    ))
  }

  cavity_Sigma <- solve(Q_cavity)
  cavity_mu <- as.vector(cavity_Sigma %*% r_cavity)

  result <- fit_model(
    data = site_data,
    use_cavity = TRUE,
    cavity_mu = cavity_mu,
    cavity_Sigma = cavity_Sigma
  )

  if (is.null(result)) {
    return(list(
      delta_Q = matrix(0, d, d),
      delta_r = rep(0, d),
      status = "failed",
      fit = NULL,
      group_index = NULL,
      phi_draws = NULL
    ))
  }

  ##### Extract shared parameters from fit #####

  phi_draws <- as.matrix(result$phi)
  n_draws <- nrow(phi_draws)

  if (ncol(phi_draws) != d) {
    warning(sprintf("Expected %d shared parameters but got %d from phi draws", d, ncol(phi_draws)))
  }

  phi_mean <- colMeans(phi_draws)
  phi_cov <- stats::cov(phi_draws)

  Q_tilted <- (n_draws - d - 2) * solve(phi_cov * (n_draws - 1))
  r_tilted <- as.vector(Q_tilted %*% phi_mean)

  delta_Q <- (1 / eta) * (Q_tilted - Q_cavity) - Q_site
  delta_r <- (1 / eta) * (r_tilted - r_cavity) - r_site

  list(
    delta_Q = delta_Q,
    delta_r = delta_r,
    status = "success",
    fit = result,
    phi_draws = phi_draws
  )
}
