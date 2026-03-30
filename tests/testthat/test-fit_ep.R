##### Helpers ####

simulate_test_data <- function(C = 6, N = 20, seed = 1) {
  set.seed(seed)
  true_beta <- rnorm(C, mean = 2, sd = 0.5)
  data.frame(
    group = rep(1:C, each = N),
    y = unlist(lapply(1:C, function(c) rnorm(N, true_beta[c], sd = 1)))
  )
}

make_fit_model <- function(d = 3, n_draws = 500, noise_sd = 0.1) {
  function(data, use_cavity = FALSE, cavity_mu = NULL, cavity_Sigma = NULL) {
    if (use_cavity && !is.null(cavity_mu)) {
      phi <- MASS::mvrnorm(
        n_draws,
        mu = cavity_mu,
        Sigma = cavity_Sigma * noise_sd
      )
    } else {
      phi <- MASS::mvrnorm(n_draws, mu = rep(0, d), Sigma = diag(d))
    }

    list(phi = phi)
  }
}

make_failing_fit_model <- function(d = 3) {
  function(data, use_cavity = FALSE, cavity_mu = NULL, cavity_Sigma = NULL) {
    NULL
  }
}

skip_if_no_stan <- function() {
  if (!requireNamespace("cmdstanr", quiety = TRUE)) {
    skip("cmdstanr not available")
  }
  tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) skip("CmdStan not installed")
  )
}


##### Check argument validation #####

test_that("fit_ep rejects invalid data argument", {
  fm <- make_fit_model()
  expect_error(fit_ep(data = "invalid", group_column = "group", K = 2, d = 3, fit_model = fm))
})

test_that("fit_ep rejects invalid group_column", {
  dat <- simulate_test_data()
  fm <- make_fit_model()

  expect_error(fit_ep(data = dat, group_column = "nonexistent", K = 2, d = 3, fit_model = fm))
  expect_error(fit_ep(data = dat, group_column = 123, K = 2, d = 3, fit_model = fm))
})

test_that("fit_ep rejects invaliad K and d", {
  dat <- simulate_test_data()
  fm <- make_fit_model()

  expect_error(fit_ep(data = dat, group_column = "group", K = 0, d = 3, fit_model = fm))
  expect_error(fit_ep(data = dat, group_column = "group", K = -1, d = 3, fit_model = fm))
  expect_error(fit_ep(data = dat, group_column = "group", K = 2.5, d = 3, fit_model = fm))
  expect_error(fit_ep(data = dat, group_column = "group", K = 2.5, d = 0, fit_model = fm))
})


test_that("fit_ep rejects invalid fit_model", {
  dat <- simulate_test_data()
  expect_error(fit_ep(data = dat, group_column = "group", K = 2, d = 0, fit_model = "not a function"))

  wrong_argument_names <- function(x, y, z) NULL
  expect_error(fit_ep(data = dat, group_column = "group", K = 2, d = 3, fit_model = wrong_argument_names))

  additional_argument_names <- function(data, use_cavity, cavity_mu, cavity_Sigma, additional_arguments) NULL
  expect_no_error(fit_ep(data = dat, group_column = "group", K = 2, d = 3, fit_model = additional_argument_names, max_iter = 1))
})

test_that("fit_ep rejects data with fewer than 2 groups", {
  one_group <- data.frame(group = rep(1, 10), y = rnorm(10))
  fm <- make_fit_model()
  expect_error(fit_ep(
    data = one_group,
    group_column = "group",
    K = 1,
    d = 3,
    fit_model = fm
  ))
})

##### Test return structure #####

test_that("fit_ep returns correct structure", {
  dat <- simulate_test_data(C = 4)
  fm <- make_fit_model()
  result <- fit_ep(data = dat, group_column = "group", K = 2, d = 3, fit_model = fm, max_iter = 2, verbose = FALSE)

  expect_type(result, "list")
  expect_true("mu" %in% names(result))
  expect_true("Sigma" %in% names(result))
  expect_true("ep_log" %in% names(result))
  expect_true("tilted_fits" %in% names(result))

  # Check dimensions equal input dimension
  expect_length(result$mu, 3)
  expect_equal(dim(result$Sigma), c(3, 3))
  expect_true(length(result$ep_log) >= 1)
  expect_true(length(result$ep_log) <= 2)
})

test_that("fit_ep Sigma is symmetric positive definite", {
  dat <- simulate_test_data(C = 4)
  fm <- make_fit_model()
  result <- fit_ep(data = dat, group_column = "group", K = 2, d = 3, fit_model = fm, max_iter = 2, verbose = FALSE)

  expect_equal(result$Sigma, t(result$Sigma))
  eigs <- eigen(result$Sigma, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

##### Test convergence behavior #####

test_that("fit_ep stops early when convergence is reached", {
  dat <- simulate_test_data(C = 4)
  fm <- make_fit_model(noise_sd = 1e-3)

  result <- fit_ep(data = dat, group_column = "group", K = 2, d = 3, fit_model = fm, max_iter = 50, conv_tol = 2e3, verbose = FALSE)
  expect_true(length(result$ep_log) < 50)
})

test_that("fit_ep runs all iterations if conv_tol is small", {
  dat <- simulate_test_data(C = 4)
  fm <- make_fit_model(noise_sd = 1e-3)
  max_iter <- 3
  result <- fit_ep(data = dat, group_column = "group", K = 2, d = 3, fit_model = fm, max_iter = max_iter, conv_tol = 1e-3, verbose = FALSE)
  expect_length(result$ep_log, max_iter)
})

##### Test failing fit_model #####
test_that("fit_ep handles failing fit_model", {
  dat <- simulate_test_data(C = 4)
  fm <- make_failing_fit_model()
  result <- fit_ep(data = dat, group_column = "group", K = 2, d = 3, fit_model = fm, max_iter = 2, verbose = FALSE)

  expect_type(result, "list")
  expect_equal(result$ep_log[[1]]$n_failed, 2)
  expect_equal(result$ep_log[[1]]$n_success, 0)
  expect_equal(result$ep_log[[1]]$n_skipped, 0)
  expect_equal(result$mu, rep(0, 3))
})

##### Test group handling #####
test_that("fit_ep works with character group labels", {
  dat <- data.frame(
    group = rep(c("A", "B", "C", "D"), each = 15),
    y = rnorm(60)
  )
  fm <- make_fit_model()
  result <- fit_ep(data = dat, group_column = "group", K = 2, d = 3, fit_model = fm, max_iter = 1, verbose = FALSE)

  expect_type(result, "list")
  expect_length(result$mu, 3)
})

##### Test integration with Stan #####
test_that("fip ep recovers approximate posterior with Stan", {
  skip_on_ci()
  skip_if_no_stan()

  set.seed(10016)
  C <- 10
  N <- 30
  true_beta_mu <- 2.0
  true_beta_sigma <- 0.8
  true_sigma <- 1.0
  true_beta <- rnorm(C, true_beta_mu, true_beta_sigma)

  dat <- data.frame(
    group = rep(1:C, each = N),
    y = unlist(lapply(1:C, function(c) rnorm(N, true_beta[c], true_sigma)))
  )

  stan_model_code <- "
  data {
    int<lower=0> C;
    int<lower=0> N;
    vector[N] y;
    array[N] int<lower=1, upper=C> group_index;

    // Flag: 0 = use original priors, 1 = use cavity distribution
    int<lower=0, upper=1> use_cavity;

    vector[3] cavity_mu;
    matrix[3, 3] cavity_Sigma;
  }
  parameters {
    vector[C] raw_beta;
    real beta_mu;
    real log_beta_sigma;
    real log_sigma;
  }
  transformed parameters {
    real sigma = exp(log_sigma);
    real beta_sigma = exp(log_beta_sigma);

    vector[C] beta = beta_mu + beta_sigma * raw_beta;

    vector[3] phi;
    phi[1] = beta_mu;
    phi[2] = log_beta_sigma;
    phi[3] = log_sigma;
  }
  model {
    raw_beta ~ std_normal();

    if(use_cavity == 1) {
      phi ~ multi_normal(cavity_mu, cavity_Sigma);
    }
    else {
      beta_mu ~ std_normal();
      log_beta_sigma ~ std_normal();
      log_sigma ~ std_normal();
    }

    for(n in 1:N) {
      y[n] ~ normal(beta[group_index[n]], sigma);
    }
  }
  "

  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stan_model_code))

  fitter <- function(data, use_cavity = FALSE, cavity_mu = NULL, cavity_Sigma = NULL) {
    if (is.null(cavity_mu)) {
      cavity_mu <- rep(0, 3)
      cavity_Sigma <- diag(3)
    }

    group_index <- data.frame(group = unique(data$group))
    group_index$ci <- 1:nrow(group_index)

    data$ci <- group_index$ci[base::match(data$group, group_index$group)]

    stan_data <- list(
      N = nrow(data),
      y = data$y,
      group_index = data$ci,
      C = nrow(group_index),
      use_cavity = as.integer(use_cavity),
      cavity_mu = cavity_mu,
      cavity_Sigma = cavity_Sigma
    )

    fit <- model$sample(
      data = stan_data,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 500,
      iter_sampling = 500,
      refresh = 0,
      show_messages = FALSE,
      show_exceptions = FALSE,
      seed = 10016
    )

    list(
      phi = posterior::as_draws_matrix(fit$draws("phi")),
      fit = fit,
      group_index = group_index
    )
  }

  full_result <- fitter(dat, use_cavity = FALSE)
  full_phi <- as.matrix(full_result$phi)
  full_mean <- colMeans(full_phi)
  full_sd <- apply(full_phi, 2, sd)

  ep_result <- fit_ep(
    data = dat,
    group_column = "group",
    K = 2,
    d = 3,
    fit_model = fitter,
    max_iter = 15,
    conv_tol = 1,
    verbose = TRUE
  )

  expect_equal(ep_result$mu, unname(full_mean), tolerance = 0.1)

  ep_sd <- sqrt(diag(ep_result$Sigma))
  for (i in 1:3) {
    expect_true(ep_sd[i] > full_sd[i] * 0.3, label = sprintf("EP SD[%d] not too small", i))
    expect_true(ep_sd[i] < full_sd[i] * 3.0, label = sprintf("EP SD[%d] not too large", i))
  }
})
