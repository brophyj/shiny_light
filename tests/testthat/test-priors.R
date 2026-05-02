test_that("sigma_from_ci matches the normal CI formula", {
  # For a 95% CI of [0.5, 2] on RR:
  # sigma = (log(2) - log(0.5)) / (2 * 1.96)
  expect_equal(sigma_from_ci(0.5, 2.0),
               (log(2) - log(0.5)) / (2 * qnorm(0.975)),
               tolerance = 1e-12)
  # symmetric CI on log-scale: sigma should be strictly positive
  expect_gt(sigma_from_ci(0.8, 1.25), 0)
})

test_that("post_from_prior_normal returns conjugate Normal update", {
  # Prior N(0, 2) combined with observation y=0.2, sigma=0.1
  y <- 0.2; sigma <- 0.1; mu0 <- 0; tau <- 2
  w_data <- 1 / sigma^2; w_prio <- 1 / tau^2
  expected_mu <- (w_prio * mu0 + w_data * y) / (w_prio + w_data)
  expected_sd <- sqrt(1 / (w_prio + w_data))
  out <- post_from_prior_normal(y, sigma, mu0, tau)
  expect_equal(out$mu, expected_mu, tolerance = 1e-10)
  expect_equal(out$sd, expected_sd, tolerance = 1e-10)
  # Tighter prior pulls the posterior mean toward mu0
  tight <- post_from_prior_normal(y, sigma, 0, 0.01)
  expect_lt(abs(tight$mu), abs(out$mu))
})

test_that("map_logRR_prior_to_alpha returns finite and sensible values", {
  res <- map_logRR_prior_to_alpha(mu_eta = log(0.8), sd_eta = 0.2, p2 = 0.1)
  expect_true(all(is.finite(res)))
  expect_true(res["s_alpha"] > 0)
})
