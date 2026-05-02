# R/summaries.R
#
# Posterior summary helpers.  These operate on long-format tibbles with
# columns (prior, mu, RR) as produced by the posterior_df() reactive.

#' Summary of posterior RR draws per prior, including MCID probability
#'
#' @param df  tibble with columns (prior, RR).
#' @param mcid numeric, the MCID threshold on the RR scale (e.g. 0.80).
#' @return tibble with mean/CrI and probability columns.  Includes
#'   P(RR <= MCID) and P(RR < 1) as requested.
summarise_posterior_rr <- function(df, mcid = 0.80) {
  stopifnot(all(c("prior", "RR") %in% names(df)))
  df %>%
    dplyr::group_by(.data$prior) %>%
    dplyr::summarise(
      `Mean [95% CrI]` = sprintf(
        "%.2f [%.2f, %.2f]",
        mean(.data$RR),
        stats::quantile(.data$RR, 0.025),
        stats::quantile(.data$RR, 0.975)
      ),
      `P(RR <= MCID)` = mean(.data$RR <= mcid),
      `P(RR < 1)`     = mean(.data$RR < 1),
      `P(RR < 0.80)`  = mean(.data$RR < 0.80),
      `P(RR < 0.90)`  = mean(.data$RR < 0.90),
      `P(RR < 1.15)`  = mean(.data$RR < 1.15),
      .groups = "drop"
    )
}

#' Posterior predictive summary
#'
#' @param df tibble with columns (prior, RR_pred)
#' @return tibble with 2.5 / 50 / 97.5% quantiles
summarise_predictive_rr <- function(df) {
  df %>%
    dplyr::group_by(.data$prior) %>%
    dplyr::summarise(
      `Predictive 2.5%`  = stats::quantile(.data$RR_pred, 0.025),
      `Predictive 50%`   = stats::quantile(.data$RR_pred, 0.50),
      `Predictive 97.5%` = stats::quantile(.data$RR_pred, 0.975),
      .groups = "drop"
    )
}

# ---- Risk Difference (RD) summaries ----
# Draws are on the proportion scale; display values are in percentage
# points (pp) so readers see e.g. -1.3 pp rather than -0.013.

#' Summary of posterior RD draws per prior, including MCID probability
#'
#' @param df   tibble with columns (prior, RD) on the proportion scale.
#' @param mcid_rd numeric MCID threshold on the RD proportion scale
#'   (e.g. -0.02 for a 2 pp absolute reduction).
#' @return tibble with mean/CrI (in pp) and probability columns.
summarise_posterior_rd <- function(df, mcid_rd = -0.02) {
  stopifnot(all(c("prior", "RD") %in% names(df)))
  df %>%
    dplyr::group_by(.data$prior) %>%
    dplyr::summarise(
      `Mean [95% CrI] (pp)` = sprintf(
        "%.2f [%.2f, %.2f]",
        mean(.data$RD) * 100,
        stats::quantile(.data$RD, 0.025) * 100,
        stats::quantile(.data$RD, 0.975) * 100
      ),
      `P(RD <= MCID)` = mean(.data$RD <= mcid_rd),
      `P(RD < 0)`     = mean(.data$RD < 0),
      .groups = "drop"
    )
}

#' Posterior predictive summary on the RD scale
#'
#' @param df tibble with columns (prior, RD_pred) on proportion scale.
#' @return tibble with 2.5 / 50 / 97.5% quantiles in percentage points.
summarise_predictive_rd <- function(df) {
  df %>%
    dplyr::group_by(.data$prior) %>%
    dplyr::summarise(
      `Predictive 2.5% (pp)`  = stats::quantile(.data$RD_pred, 0.025) * 100,
      `Predictive 50% (pp)`   = stats::quantile(.data$RD_pred, 0.50)  * 100,
      `Predictive 97.5% (pp)` = stats::quantile(.data$RD_pred, 0.975) * 100,
      .groups = "drop"
    )
}

#' Data-driven x-axis limits from posterior draws
#'
#' Uses quantiles of the pooled draws with a small padding.  Floors the
#' lower bound at 0 since RR must be positive.
#'
#' @param x numeric vector of RR draws.
#' @param q two-element vector of lower/upper probabilities for the range.
#' @param pad multiplicative padding on each side (on the log scale).
#' @return numeric c(lo, hi) for use with coord_cartesian.
rr_plot_limits <- function(x, q = c(0.005, 0.995), pad = 0.10) {
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(c(0.4, 1.6))
  lim <- stats::quantile(x, q)
  llim <- log(lim)
  span <- diff(llim)
  c(exp(llim[1] - pad * span), exp(llim[2] + pad * span))
}

#' Data-driven x-axis limits for RD draws (percentage-point scale).
#'
#' @param x numeric vector of RD draws (proportion scale).
#' @param q two-element vector of lower/upper probabilities.
#' @param pad additive padding in proportion units on each side.
#' @return numeric c(lo, hi) in PROPORTION scale for coord_cartesian.
rd_plot_limits <- function(x, q = c(0.005, 0.995), pad = 0.005) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(-0.05, 0.05))
  lim  <- stats::quantile(x, q)
  span <- diff(lim)
  c(lim[1] - pad - 0.10 * span, lim[2] + pad + 0.10 * span)
}

#' Translate an MCID specified implicitly via a power calculation into
#' an RR.  If the user specifies an explicit target_rr, it is returned
#' directly.  Otherwise a rough effect size implied by the power target
#' is derived from Normal-approx power.
#'
#' @param mcid_rr user-supplied MCID (numeric or NA).
#' @param n_per_arm per-arm sample size for the implied power calculation.
#' @param alpha significance level (two-sided).
#' @param target_rr hypothesised RR used as a fallback.
#' @param power target power (e.g. 0.80).
#' @param p_ctrl control-arm event rate (defaults to 0.10).
#' @return numeric MCID on the RR scale.
resolve_mcid <- function(mcid_rr, n_per_arm = NA, alpha = 0.05,
                         target_rr = NA, power = 0.80, p_ctrl = 0.10) {
  if (!is.na(mcid_rr) && is.finite(mcid_rr) && mcid_rr > 0) return(mcid_rr)
  if (!is.na(target_rr) && is.finite(target_rr) && target_rr > 0) return(target_rr)
  if (is.na(n_per_arm) || !is.finite(n_per_arm) || n_per_arm <= 0) return(0.80)
  # rough derivation: find RR such that normal-approx two-sample test has
  # the requested power at the control-arm event rate p_ctrl.
  z_a <- qnorm(1 - alpha / 2)
  z_b <- qnorm(power)
  se_per_arm <- sqrt((1 - p_ctrl) / (p_ctrl * n_per_arm) * 2)
  log_rr <- -(z_a + z_b) * se_per_arm
  exp(log_rr)
}
