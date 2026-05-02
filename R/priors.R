# R/priors.R
#
# Helpers for building / summarising prior distributions used by the Shiny app.
# All priors here are expressed on the log(RR) scale as Normal(mu0, tau^2),
# so densities can be put on the RR scale with a simple change of variable.
#
# Functions:
#   sigma_from_ci()         : derive sigma (log-RR SD) from a 95% CI on RR
#   post_from_prior_normal(): Normal-Normal conjugate update on log-RR
#   logit()/invlogit()      : small numeric helpers shared with other modules

#' Derive SD on log(RR) from a symmetric 95% CI on RR
#'
#' @param lower Lower bound of a 95% CI on the risk-ratio scale (>0).
#' @param upper Upper bound of a 95% CI on the risk-ratio scale (> lower).
#' @return A scalar SD on the log(RR) scale implied by a Normal CI.
#' @details Uses the classical back-calculation
#'   sigma = (log(upper) - log(lower)) / (2 * z_{0.975}).
sigma_from_ci <- function(lower, upper) {
  (log(upper) - log(lower)) / (2 * qnorm(0.975))
}

#' Conjugate Normal update for a Normal prior on log(RR)
#'
#' @param y Observed log(RR) (scalar).
#' @param sigma SE of the observed log(RR) (scalar, > 0).
#' @param mu0 Prior mean on log(RR).
#' @param tau Prior SD on log(RR) (> 0).
#' @return A list with posterior mean (mu) and posterior SD (sd) for log(RR).
post_from_prior_normal <- function(y, sigma, mu0, tau) {
  w_data <- 1 / sigma^2
  w_prio <- 1 / tau^2
  mu_post <- (w_prio * mu0 + w_data * y) / (w_prio + w_data)
  sd_post <- sqrt(1 / (w_prio + w_data))
  list(mu = mu_post, sd = sd_post)
}

#' Numerically safe logit / inverse-logit
#'
#' @param p Probability in (0, 1) for logit(); real number for invlogit().
#' @return Transformed value.
logit <- function(p) log(p / (1 - p))
invlogit <- function(x) 1 / (1 + exp(-x))

#' Translate a Normal prior on log(RR) to a Normal prior on the arm-logit
#' contrast (alpha).
#'
#' Used by Laplace / MCMC paths which parameterise in logit-space.  Works by
#' first mapping mu via p1 = p2 * exp(mu), then applying a first-order delta
#' transform for the variance:  Var(logit(p1)) approximately equal to
#' (1/(1-p1))^2 * Var(log(p1)), and similarly for the contrast.
#'
#' @param mu_eta Prior mean on log(RR).
#' @param sd_eta Prior SD on log(RR).
#' @param p2     Control-arm event rate used to anchor the delta approximation.
#' @return Named numeric vector c(m_alpha, s_alpha): mean and SD of the
#'   induced Normal prior on the alpha (arm-effect) logit coefficient.
map_logRR_prior_to_alpha <- function(mu_eta, sd_eta, p2) {
  eps <- 1e-8
  p2  <- min(max(p2, eps), 1 - eps)
  p1  <- min(max(p2 * exp(mu_eta), eps), 1 - eps)
  mu_alpha <- logit(p1) - logit(p2)
  sd_alpha <- sd_eta / (1 - p1)   # first-order delta method
  sd_alpha <- max(sd_alpha, 1e-6)
  c(m_alpha = mu_alpha, s_alpha = sd_alpha)
}

#' Derive SD on the risk-difference (proportion) scale from a symmetric
#' 95% CI.
#'
#' @param lower Lower bound of a 95% CI on the RD scale.
#' @param upper Upper bound (> lower).
#' @return Scalar SD on the RD (proportion) scale.
sigma_from_ci_rd <- function(lower, upper) {
  (upper - lower) / (2 * qnorm(0.975))
}

#' Translate a Normal prior on RD (proportion scale) to a Normal prior
#' on the arm-logit contrast (alpha).
#'
#' Analogous to map_logRR_prior_to_alpha() but for the additive
#' risk-difference parameterisation.  Uses a first-order delta-method
#' approximation:  dRD/dalpha = p1*(1-p1) where p1 = p2 + mu_rd.
#'
#' @param mu_rd Prior mean on RD (proportion scale).
#' @param sd_rd Prior SD on RD (proportion scale).
#' @param p2    Control-arm event rate.
#' @return Named numeric vector c(m_alpha, s_alpha).
map_RD_prior_to_alpha <- function(mu_rd, sd_rd, p2) {
  eps <- 1e-8
  p2 <- min(max(p2, eps), 1 - eps)
  p1 <- min(max(p2 + mu_rd, eps), 1 - eps)
  mu_alpha <- logit(p1) - logit(p2)
  # delta method: dRD/dalpha = p1*(1-p1)
  sd_alpha <- sd_rd / (p1 * (1 - p1))
  sd_alpha <- max(sd_alpha, 1e-6)
  c(m_alpha = mu_alpha, s_alpha = sd_alpha)
}
