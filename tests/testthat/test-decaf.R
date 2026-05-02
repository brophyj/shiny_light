test_that("EB posterior HR median for DECAF lies in [0.65, 0.90]", {
  # DECAF trial summary from eb2.qmd: HR 0.61 (95% CI 0.42-0.89).
  # SE[log HR] is derived from the CI; the EB mixture (van Zwet et al. 2021/2025)
  # shrinks the raw log-HR toward 0, yielding a posterior HR closer to 1.
  # We expect the posterior median HR to land well inside [0.65, 0.90].
  set.seed(2026)
  y  <- log(0.61)
  se <- (log(0.89) - log(0.42)) / (2 * qnorm(0.975))
  dr <- eb_posterior_logRR_draws(y, se, ndraws = 200000)
  med_hr <- median(exp(dr))
  expect_gte(med_hr, 0.65)
  expect_lte(med_hr, 0.90)
})
