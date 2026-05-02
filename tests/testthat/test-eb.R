test_that("eb_posterior_snr moments are reasonable on synthetic z", {
  set.seed(1)
  # For z=0, by symmetry the posterior mean on SNR should be ~0, and the
  # overall shrinkage B_bar should lie strictly in (0, 1).
  fit0 <- eb_posterior_snr(0)
  expect_true(is.finite(fit0$mean))
  expect_equal(fit0$mean, 0, tolerance = 1e-10)
  expect_gt(fit0$B_bar, 0)
  expect_lt(fit0$B_bar, 1)

  # For a large positive z the posterior mean should be positive but
  # strictly less than z (shrinkage toward 0).
  z <- 3
  fit <- eb_posterior_snr(z)
  expect_gt(fit$mean, 0)
  expect_lt(fit$mean, z)

  # Sampler should produce draws whose empirical mean / SD are close
  # to the analytic ones.
  dr <- fit$sampler(50000)
  expect_equal(mean(dr), fit$mean, tolerance = 0.05)
  expect_equal(sd(dr),   fit$sd,   tolerance = 0.05)
})

test_that("eb_marginal_logRR_prior scales with SE correctly", {
  pr1 <- eb_marginal_logRR_prior(1)
  pr2 <- eb_marginal_logRR_prior(2)
  expect_equal(pr1["mu0"], c(mu0 = 0))
  expect_equal(pr2["tau"] / pr1["tau"], c(tau = 2), tolerance = 1e-10)
})

test_that("eb_posterior_logRR_draws shrinks toward 0", {
  set.seed(7)
  y <- log(0.61); sigma <- 0.2
  draws <- eb_posterior_logRR_draws(y, sigma, ndraws = 20000)
  # Median log-RR should lie between y and 0 (i.e. shrunk toward null).
  m <- median(draws)
  expect_lt(m, 0)
  expect_gt(m, y)
})
