#' Fit a hierarchical Bayesian model using distributed Expectation Propagation
#'
#' Approximates the posterior distribution of shared (hierarchical) parameters
#' in a Bayesian hierarchical model by iteratively passing messages between
#' independently fitted subsets of the data. This implements the EP framework
#' for described by Vehtari et al. (2020).
#'
#' @section Overview:
#' In a hierarchical model with shared parameters \eqn{\phi} and group specific
#' parameters \eqn{\alpha_1, \ldots, \alpha_C}, the joint posterior is
#' \deqn{p(\phi, \alpha | y) \propto p(\phi) \prod_{c=1}^C p(y_c | \alpha_c, \phi) p(\alpha_c | \phi)}
#'
#' Rather than fitting this model jointly (which may be computationally
#' prohibitive for large \eqn{C}), `fit_ep` partitions the groups into \eqn{K}
#' sites and iteratively:
#' 1. Computes a **cavity distribution** for each site by removing that site's
#'    contribution from the current global approximation
#' 2. Fits the model at each site using only that site's data and the cavity
#'    distribution as a prior on the shared parameters
#' 3. Updates the global approximation based on the site's resutls
#'
#' The use of the cavity distribution ensures that each site's contribution is
#' counted exactly once, with information from all other sites regularizing
#' local inference.
#'
#' @section The fit_model contract:
#' The `fit_model` argument must be a user-supplied function that fits the
#' hierarchical model to a subset of the data, using either the original priors
#' or a cavity distribution as the prior on the shared parameters.
#'
#' **Signature:**
#' ```
#' fit_model(data, use_cavity, cavity_mu, cavity_Sigma)
#' ````
#'
#' **Arguments:**
#' \describe{
#'   \item{data}{A data frame containing the observations for a single site.
#'     This will be a subset of the `data` argument passed to `fit_ep`,
#'     filtered to include only the groups assigned to this site.}
#'   \item{use_cavity}{Logical. If `TRUE`, the model should use the cavity
#'     distribution as the prior on the shared parameters. If `FALSE`, the
#'     model should use the original priors (e.g., for a full joint fit.)}
#'   \item{cavity_mu}{Numeric vector of length `d`. The mean of the multivariate
#'     normal cavity distribution over the shared parameters. When
#'     `use_cavity = FALSE`, this will be `NULL`.}
#'   \item{cavity_Sigma}{A `d x d` positive definite covariance matrix of the
#'     multivariate normal cavity distribution over the shared parameters. When
#'     `use_cavity = FALSE`, this will be `NULL`.}
#' }
#'
#' **Return value:**
#' On success, `fit_model` must return a named list containing (at minimum):
#' \describe{
#'   \item{phi}{An object coercible to a matrix via `as.matrix`, where
#'     rows are posterior draws and columns are the `d` shared parameters.
#'     This is the only element inspected by the EP algorithm. A
#'     `draws_matrix` from the `posterior` package, or a plain numeric matrix,
#'     will work. The column ordering must be consistent across all calls.}
#' }
#'
#' Any additional elements in the returned list (e.g., the full fit object,
#' group indices, diagnostics) are passed through and stored in the `tilted_fits`
#' element of the `fit_ep` return value for the user's convenience.
#'
#' On failure, `fit_model` should return NULL. Alternatively, it may throw an
#' error, which will be caught internally and treated as a failure. Failed
#' sites are skipped for that iteration without affecting other sites.
#'
#' @section Shared parameter vector (phi):
#' The EP algorithm approximates the marginal posterior of the shared parameters
#' as a multivariate normal distribution parameters by its natural parameters:
#' a precision matrix \eqn{Q} and a precision mean vector \eqn{r}.
#'
#' The shared parameter vector `phi` return by `fit_model` must have a consistent
#' ordering across all calls.
#'
#' Shared parameters should be defined on an unconstrained scale (e.g., use
#' log-transformed standard deviations rather than raw standard deviations)
#' so that the multivariate normal approximation is reasonable.
#'
#' @section Damping:
#' To stabilize convergence when updating all sites in parallel within each
#' iteration, a damping factor \eqn{\delta \in (0, 1]} is applied to the site
#' parameter updates. The default schedule starts at \eqn{\delta = 0.5}
#' for the first iteration and decays exponentially toward \eqn{\min(1/K, 0.2)}
#' over subsequent iterations (Vehtari et al. 2020).
#'
#' @section Convergence:
#' The algorithm terminates when the maximum absolute change in the site natural
#' parameters falls below `conv_tol` for all sites, or when `max_iter` is reached
#' (Vehtari et al. 2020).
#'
#' @param data A data frame containing all observations across all groups.
#' @param group_column A string giving the name of the column in `data` that
#'   identifies group membership.
#' @param K Integer. The number of EP sites.
#' @param d Integer. The dimension of the shared parameter vector `phi`.
#' @param fit_model A function implementing the model fitting. See the
#'   **fit_model contract** section for details.
#' @param max_iter Integer. Maximum number of EP iterations (default: 20).
#' @param conv_tol Numeric. Convergence tolerance (default: 0.1).
#' @param verbose Logical. If `TRUE` (default), prints progress messages.
#' @param save_all_tilted_fits Logical. If `TRUE`, the tilted distribution fit
#'   results from every iteration is stored. If `FALSE` (default), only the
#'   final iteration fits are stored.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{mu}{Numeric vector of length `d`. The mean of the EP approximation to
#'     the marginal posterior of the shared parameters.}
#'   \item{Sigma}{A `d x d` covariance matrix of the EP approximation to the
#'     marginal posterior of the shared parameters.}
#'   \item{ep_log}{A list of length equal to the number of completed iterations.
#'     Each element is a list containing:
#'     \describe{
#'       \item{iteration}{Iteration number.}
#'       \item{mu}{Posterior mean at this iteration.}
#'       \item{Sigma}{Posterior covariance at this iteration.}
#'       \item{max_delta_Q}{Maximum absolute change in site precision matrices.}
#'       \item{max_delta_r}{Maximum absolute change in site precision means.}
#'       \item{n_success}{Number of sites successfully updated.}
#'       \item{n_failed}{Number of sites where `fit_model` returned `NULL`.}
#'       \item{n_skipped}{Number of sites skipped due to improper cavity
#'         distributions.}
#'     }}
#'   \item{tilted_fits}{The fit results from `fit_site` for each site. If
#'     `save_all_tilted_fits = TRUE`, this is a list with length equal to the
#'     number of completed iterations, where each list element is a length `K`
#'     list containing the results from that iteration. If
#'     `save_all_tilted_fits = FALSE`, this is a list with only one element,
#'     a length `K` list containing the results from the final iteration.
#'   }
#'   \item{site_assignments}{A length `K` list containing the groups assigned
#'     to each site.}
#' }
#'
#' @examples
#' \dontrun{
#' # Define a fit_model function for use with cmdstanr
#' my_fitter <- function(data, use_cavity = FALSE, cavity_mu = NULL, cavity_Sigma = NULL) {
#'   if(is.null(cavity_mu)) cavity_mu <- rep(0, 3)
#'   if(is.null(cavity_Sigma)) cavity_Sigma <- diag(3)
#'
#'   stan_data <- list(
#'     N            = nrow(data),
#'     y            = data$y,
#'     group        = data$group_index,
#'     C            = length(unique(data$group_index)),
#'     use_cavity   = as.integer(use_cavity),
#'     cavity_mu    = cavity_mu,
#'     cavity_Sigma = cavity_Sigma
#'   )
#'
#'   samples <- model$sample(data = stan_data, chains = 4, refresh = 0)
#'
#'   list(
#'     phi = posterior::as_draws_matrix(samples$draws("phi")),
#'     fit = fit
#'   )
#' }
#'
#' result <- fit_ep(
#'   data = my_data,
#'   group_column = "group",
#'   K = 10,
#'   d = 3,
#'   fit_model = my_fitter,
#'   max_iter = 20,
#'   conv_tol = 1
#' )
#' }
#'
#' @references
#' Vehtari, A., Gelman, A., Sivula, T., Jylänki, P., Tran, D., Sahai, S.,
#' Blomstedt, P., Cunningham, J. P., Schiminovich, D., and Robert, C. P. (2020).
#' Expectation Propagation as a Way of Life: A Framework for Bayesian Inference
#' on Partitioned Data. Journal of Machine Learning Research, 21(17), 1–53.
#'
#' @export
fit_ep <- function(
  data,
  group_column,
  K,
  d,
  fit_model,
  max_iter = 20,
  conv_tol = 0.1,
  verbose = TRUE,
  save_all_tilted_fits = FALSE,
  workers = NULL
) {
  ##### Argument Checks #####
  checkmate::assert_data_frame(data, min.rows = 1)
  checkmate::assert_string(group_column)
  checkmate::assert_choice(group_column, names(data))
  checkmate::assert_count(K, positive = TRUE)
  checkmate::assert_count(d, positive = TRUE)
  checkmate::assert_function(
    fit_model,
    args = c("data", "use_cavity", "cavity_mu", "cavity_Sigma")
  )
  checkmate::assert_count(max_iter, positive = TRUE)
  checkmate::assert_number(conv_tol, lower = 0, finite = TRUE)
  checkmate::assert_flag(verbose)
  checkmate::assert_flag(save_all_tilted_fits)
  checkmate::assert_class(workers, c("Pool", "R6"), null.ok = TRUE)
  if (!is.null(workers)) {
    if (!requireNamespace("clustermq", quiety = TRUE)) {
      stop("The 'clustermq' package is required when providing a worker pool.")
    }
  }

  # For now eta is hardcoded
  eta <- 1

  ##### Assign groups to sites #####
  group_names <- sort(unique(data[[group_column]]))
  C <- length(group_names)
  checkmate::assert_true(C >= 2, .var.name = "number of groups (must be >= 2)")

  if (K > C) {
    warning("Number of sites greater than number of groups. Setting K = C.")
    K <- C
  }

  site_assignments <- assign_groups_to_sites(C, K)
  site_data_list <- lapply(1:K, function(k) {
    data[data[[group_column]] %in% group_names[site_assignments[[k]]], ]
  })

  ##### Initialize EP #####
  Q_prior <- diag(d)
  r_prior <- rep(0, d)

  Q_global <- Q_prior
  r_global <- r_prior

  Q_sites <- lapply(1:K, function(k) matrix(0, d, d))
  r_sites <- lapply(1:K, function(k) rep(0, d))

  ep_log <- list()
  tilted_fits <- list(list())

  ##### Send data to workers if parallel #####
  if (!is.null(workers)) {
    if (verbose) {
      cat("Sending data and functions to workers...\n")
    }

    workers$env(list(
      fit_site = fit_site,
      fit_model = fit_model,
      d = d,
      eta = eta
    ))
  }

  ##### EP Loop #####
  iter <- 1
  converged <- FALSE
  while (converged == FALSE && iter <= max_iter) {
    if (save_all_tilted_fits) {
      tilted_fits[[iter]] <- list()
    }

    if (iter == 1) {
      delta <- 0.5
    } else {
      delta_target <- min(1 / K, 0.2)
      delta <- delta_target + (0.5 - delta_target) * exp(-3 * (iter - 1) / K)
    }

    if (verbose == TRUE) {
      cat(sprintf(
        "--- Starting iteration %d (damping = %.4f) ---\n",
        iter,
        delta
      ))
    }

    # Fit all the sites
    site_results <- ep_iteration(
      site_data_list = site_data_list,
      d = d,
      Q_global = Q_global,
      r_global = r_global,
      Q_sites = Q_sites,
      r_sites = r_sites,
      eta = 1,
      fit_model = fit_model,
      workers = workers,
      verbose = verbose
    )

    # Aggregate site results
    total_delta_Q <- matrix(0, d, d)
    total_delta_r <- rep(0, d)
    n_success <- 0
    n_skipped <- 0
    n_failed <- 0

    for (k in 1:K) {
      result <- site_results[[k]]

      if (result$status == "success") {
        n_success <- n_success + 1
        total_delta_Q <- total_delta_Q + result$delta_Q
        total_delta_r <- total_delta_r + result$delta_r
        Q_sites[[k]] <- Q_sites[[k]] + delta * result$delta_Q
        r_sites[[k]] <- r_sites[[k]] + delta * result$delta_r
      } else if (result$status == "skipped") {
        n_skipped <- n_skipped + 1
      } else if (result$status == "failed") {
        n_failed <- n_failed + 1
      }
      if (save_all_tilted_fits) {
        tilted_fits[[iter]][[k]] <- result
      } else {
        tilted_fits[[1]][[k]] <- result
      }
    }

    Q_new <- Q_global + delta * total_delta_Q
    r_new <- r_global + delta * total_delta_r

    eig_min <- min(eigen(Q_new, symmetric = TRUE, only.values = TRUE)$values)

    if (eig_min <= 1e-8) {
      cat(sprintf("Global Q not positive definite (min eig = %0.4f", eig_min))
    }

    Q_global <- Q_new
    r_global <- r_new

    Sigma_global <- solve(Q_global)
    mu_global <- as.vector(Sigma_global %*% r_global)

    max_delta_Q <- max(abs(total_delta_Q))
    max_delta_r <- max(abs(total_delta_r))

    ep_log[[iter]] <- list(
      iteration = iter,
      mu = mu_global,
      Sigma = Sigma_global,
      max_delta_Q = max_delta_Q,
      max_delta_r = max_delta_r,
      n_success = n_success,
      n_failed = n_failed,
      n_skipped = n_skipped
    )

    if (verbose == TRUE) {
      cat(sprintf(
        "--- Finished iteration %d (%d success, %d skipped, %d failed, max_delta_Q = %.1f, max_delta_r = %0.1f) ---\n",
        iter,
        n_success,
        n_skipped,
        n_failed,
        max_delta_Q,
        max_delta_r
      ))
    }

    if (max_delta_Q < conv_tol && max_delta_r < conv_tol) {
      if (verbose == TRUE) {
        cat(sprintf("=== Converged (conv_tol = %0.1f) ===\n", conv_tol))
      }
      converged <- TRUE
    }
    iter <- iter + 1
  }

  return(list(
    site_assignments = site_assignments,
    ep_log = ep_log,
    mu = mu_global,
    Sigma = Sigma_global,
    tilted_fits = tilted_fits
  ))
}


# Fit all sites for a single EP iteration
ep_iteration <- function(site_data_list, d, Q_global, r_global, Q_sites, r_sites, eta, fit_model, workers, verbose) {
  K <- length(site_data_list)

  if (is.null(workers)) {
    # Sequential mode
    if (verbose) {
      cat("    > Site: ")
    }
    results <- lapply(1:K, function(k) {
      if (verbose == TRUE) {
        cat(sprintf("%d ", k))
      }
      result <- fit_site(
        site_data = site_data_list[[k]],
        d = d,
        Q_global = Q_global,
        r_global = r_global,
        Q_site = Q_sites[[k]],
        r_site = r_sites[[k]],
        eta = eta,
        fit_model = fit_model
      )
    })

    if (verbose) {
      cat("\n")
    }

    results
  } else {
    # Parallel mode using clustermq workers
    results <- tryCatch(
      {
        clustermq::Q(
          fun = fit_site,
          site_data = site_data_list,
          Q_site = Q_sites,
          r_site = r_sites,
          const = list(
            d = d,
            Q_global = Q_global,
            r_global = r_global,
            eta = eta,
            fit_model = fit_model
          ),
          workers = workers,
          fail_on_error = FALSE
        )
      },
      error = function(e) {
        stop(sprintf("clustermq error: %s.", e$message))
      }
    )

    # Handle any individual failed jobs and convert to our format for failed jobs
    results <- lapply(results, function(r) {
      if (inherits(r, "error") || inherits(r, "condition")) {
        warning(sprintf("Worker error: %s", r))
        list(
          delta_Q = matrix(0, d, d),
          delta_r = rep(0, d),
          status = "failed",
          fit = NULL,
          phi_draws = NULL
        )
      } else {
        r
      }
    })

    results
  }
}
