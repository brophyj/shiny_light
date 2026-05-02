# R/eb_mixture.R
#
# Empirical-Bayes prior on the signal-to-noise ratio (SNR = beta / SE),
# following van Zwet & Gelman and the "Shrinkage Trilogy" series:
#
#   - van Zwet, Schwab & Senn (2021) Stat in Med 40:6107 (4-component
#     zero-mean Normal mixture on the z-scale; this is the principal
#     reference for the numeric parameters).
#   - van Zwet, Wiecek & Gelman (2025) Stat Methods Med Res 34(12):2302
#     "Meta-analysis with a single study" (qualitative framework).
#
# Note on parameters:
#   The `sd_z`/`p_mix` values below match those used in eb2.qmd:
#     sd_z  = c(1.17, 1.74, 2.38, 5.73)
#     p_mix = c(0.32, 0.31, 0.30, 0.07)
#   These are very close to the parameters published in van Zwet et al. 2021
#   Table 1 (sd_z ~ c(1.19, 1.71, 2.40, 5.65), same mixing weights).
#   The full 2025 paper (doi:10.1177/09622802251380628) is paywalled at
#   publication time and we were unable to verify whether it reports a
#   materially different numerical fit.  If it does, update `EB_DEFAULTS`
#   below; until then we retain the 2021 / eb2.qmd parameterisation used
#   throughout the deck-level analysis.

#' Default Cochrane EB mixture parameters (z-scale).
#'
#' @format List with elements sd_z (component SDs) and p_mix (weights).
EB_DEFAULTS <- list(
  sd_z  = c(1.17, 1.74, 2.38, 5.73),
  p_mix = c(0.32, 0.31, 0.30, 0.07)
)

#' Posterior distribution of SNR given observed z, under the 4-component
#' EB mixture prior
#'
#' SNR ~ sum_k p_k N(0, s_k^2) with s_k^2 = sd_z_k^2 - 1 on the z scale
#' (deconvolution subtracts unit noise).  Conditional on component k,
#' SNR | z ~ N(B_k z, B_k) with B_k = s_k^2 / (s_k^2 + 1).  The posterior
#' weights update via the marginal z-likelihood in each component.
#'
#' @param z Observed z-value (scalar or length-1 numeric).
#' @param sd_z Component SDs on the z-scale (length-4 default).
#' @param p_mix Component weights (length-4 default, summing to ~ 1).
#' @param n Optional number of posterior SNR draws to pre-compute (for
#'   callers that want the draws directly rather than via `sampler`).  If
#'   left NULL, no draws are generated.
#' @return list with:
#'   * components: tibble of component-level posterior weights, means, SDs
#'   * mean, sd:   pooled posterior mean and SD on the SNR scale (PR1 API)
#'   * post_mean, post_sd: aliases of mean/sd, matching the pipeline API
#'   * B_bar:      pooled shrinkage factor
#'   * sampler:    function(n) that returns n draws from the posterior mixture
#'   * snr:        pre-computed draws if `n` was supplied, else NULL
eb_posterior_snr <- function(z,
                             sd_z  = EB_DEFAULTS$sd_z,
                             p_mix = EB_DEFAULTS$p_mix,
                             n     = NULL) {
  stopifnot(length(sd_z) == length(p_mix),
            all(sd_z > 1), all(p_mix > 0))
  p_mix  <- p_mix / sum(p_mix)
  sd_snr <- sqrt(sd_z^2 - 1)
  like   <- dnorm(z, 0, sqrt(sd_snr^2 + 1))
  post_w <- p_mix * like
  post_w <- post_w / sum(post_w)
  Bk      <- sd_snr^2 / (sd_snr^2 + 1)
  post_mu <- Bk * z
  post_sd <- sqrt(Bk)

  mean_mix <- sum(post_w * post_mu)
  var_mix  <- sum(post_w * (post_sd^2 + (post_mu - mean_mix)^2))
  sd_mix   <- sqrt(var_mix)
  B_bar    <- sum(post_w * Bk)

  sampler <- function(n) {
    comp <- sample.int(length(post_w), n, replace = TRUE, prob = post_w)
    rnorm(n, post_mu[comp], post_sd[comp])
  }
  snr_draws <- if (!is.null(n) && is.finite(n) && n > 0) sampler(n) else NULL
  list(
    components = data.frame(weight = post_w, mean = post_mu, sd = post_sd,
                            B = Bk),
    mean      = mean_mix,
    sd        = sd_mix,
    post_mean = mean_mix,   # alias used by rct_pipeline.qmd
    post_sd   = sd_mix,     # alias used by rct_pipeline.qmd
    B_bar     = B_bar,
    sampler   = sampler,
    snr       = snr_draws
  )
}

#' Derive a Normal prior on log(RR) that approximates the EB prior
#' (marginal, before conditioning on the data).
#'
#' Used for Normal-likelihood calibration and as the starting point for
#' mapping to the logit-contrast scale.  Marginal SNR variance is
#' V = sum_k p_k (sd_z_k^2 - 1); scaling by SE_logRR yields a Normal prior
#' on log(RR) with mean 0 and SD sqrt(V) * SE_logRR.
#'
#' @param se_logRR Observed SE of log(RR) (scalar, > 0).
#' @param sd_z,p_mix Mixture parameters (defaults match EB_DEFAULTS).
#' @return Named numeric c(mu0, tau) on the log(RR) scale.
eb_marginal_logRR_prior <- function(se_logRR,
                                    sd_z  = EB_DEFAULTS$sd_z,
                                    p_mix = EB_DEFAULTS$p_mix) {
  p_mix  <- p_mix / sum(p_mix)
  var_snr <- sum(p_mix * (sd_z^2 - 1))
  c(mu0 = 0, tau = sqrt(var_snr) * se_logRR)
}

#' Draw posterior log(RR) samples under the EB prior using the Normal
#' likelihood approximation (z-path).
#'
#' @param y       Observed log(RR).
#' @param sigma   SE of log(RR).
#' @param ndraws  Number of draws.
#' @param sd_z,p_mix  Mixture parameters.
#' @return Numeric vector of posterior log(RR) draws (i.e. SNR * sigma).
eb_posterior_logRR_draws <- function(y, sigma, ndraws = 40000,
                                     sd_z  = EB_DEFAULTS$sd_z,
                                     p_mix = EB_DEFAULTS$p_mix) {
  z   <- y / sigma
  fit <- eb_posterior_snr(z, sd_z = sd_z, p_mix = p_mix)
  snr <- fit$sampler(ndraws)
  snr * sigma
}

#' Derive a Normal prior on RD (proportion scale) that approximates the
#' EB prior (marginal, before conditioning on the data).
#'
#' Analogous to eb_marginal_logRR_prior() but for the additive RD scale.
#' The z-score SNR structure is scale-invariant, so the same mixture
#' parameters apply.
#'
#' @param se_rd Observed SE of the risk difference (proportion scale).
#' @param sd_z,p_mix Mixture parameters (defaults match EB_DEFAULTS).
#' @return Named numeric c(mu0, tau) on the RD (proportion) scale.
eb_marginal_RD_prior <- function(se_rd,
                                 sd_z  = EB_DEFAULTS$sd_z,
                                 p_mix = EB_DEFAULTS$p_mix) {
  p_mix  <- p_mix / sum(p_mix)
  var_snr <- sum(p_mix * (sd_z^2 - 1))
  c(mu0 = 0, tau = sqrt(var_snr) * se_rd)
}

#' Draw posterior RD samples under the EB prior using the Normal
#' likelihood approximation (z-path).
#'
#' @param rd     Observed risk difference (proportion scale).
#' @param se_rd  SE of the risk difference.
#' @param ndraws Number of draws.
#' @param sd_z,p_mix Mixture parameters.
#' @return Numeric vector of posterior RD draws (proportion scale).
eb_posterior_RD_draws <- function(rd, se_rd, ndraws = 40000,
                                  sd_z  = EB_DEFAULTS$sd_z,
                                  p_mix = EB_DEFAULTS$p_mix) {
  z   <- rd / se_rd
  fit <- eb_posterior_snr(z, sd_z = sd_z, p_mix = p_mix)
  snr <- fit$sampler(ndraws)
  snr * se_rd
}
