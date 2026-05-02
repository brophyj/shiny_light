# R/pipeline.R
#
# Pipeline-specific helpers for `rct_pipeline.qmd` / `run_pipeline.R`.
# These complement R/eb_mixture.R, R/summaries.R, R/plots.R (which are
# shared with the Shiny app) with functions that are only useful in the
# systematic CV-RCT re-analysis:
#
#   * Frequentist-calibration layer: type_m_s(), s_value()
#   * Bayesian fit layer: choose_brms_backend(), fit_brms_study(),
#     build_priors()
#   * Pipeline-specific summary: summarise_posterior()  (per-study reduction
#     over draws with columns prior, RR)
#   * Pipeline-specific plots: plot_study_densities(), plot_forest(),
#     plot_heatmap_mcid(), plot_shrinkage_scatter()
#
# The Shiny app does NOT source this file.

# --- Frequentist-calibration layer -------------------------------------------

#' Gelman / Carlin retrodesign-style Type M and Type S errors.
#'
#' Simulates a hypothetical design with the given `power` to detect
#' `target_effect` at alpha = 0.05 (two-sided).  Returns the probability of
#' a sign error conditional on statistical significance (Type S) and the
#' median exaggeration ratio among significant estimates (Type M).
#'
#' Implementation follows Gelman & Carlin 2014 (Perspect. Psych. Sci.).
#'
#' @param power        nominal power (0 < power < 1).
#' @param target_effect assumed true effect size on the SNR (unit-SE) scale.
#' @param n_sim        Monte Carlo replications.
#' @param alpha        two-sided significance threshold.
#' @return list with `type_s`, `type_m` (median |ratio|), `type_m_mean`, `se`.
type_m_s <- function(power, target_effect = 1, n_sim = 1e5, alpha = 0.05) {
  stopifnot(power > 0, power < 1)
  z_crit <- qnorm(1 - alpha / 2)
  solve_se <- function(target_effect, power) {
    f <- function(se) {
      mu <- target_effect / se
      pnorm(mu - z_crit) + pnorm(-mu - z_crit) - power
    }
    tryCatch(stats::uniroot(f, interval = c(1e-4, 1e4))$root,
             error = function(e) NA_real_)
  }
  se <- solve_se(target_effect, power)
  if (!is.finite(se)) {
    return(list(type_s = NA_real_, type_m = NA_real_,
                type_m_mean = NA_real_, se = NA_real_))
  }
  est <- rnorm(n_sim, mean = target_effect, sd = se)
  sig <- abs(est) > z_crit * se
  if (!any(sig)) {
    return(list(type_s = 0, type_m = NA_real_, type_m_mean = NA_real_, se = se))
  }
  sig_est <- est[sig]
  list(
    type_s      = mean(sign(sig_est) != sign(target_effect)),
    type_m      = median(abs(sig_est / target_effect)),
    type_m_mean = mean(abs(sig_est / target_effect)),
    se          = se
  )
}

#' Shannon surprisal / S-value: `-log2(p)`.
#' @param p numeric p-value in (0, 1].
s_value <- function(p) -log2(p)

# --- Bayesian fit layer ------------------------------------------------------

# Module-level cache for the compiled brms model so per-study loops do not
# recompile Stan each call.
.brms_cache <- new.env(parent = emptyenv())
.brms_cache$fit     <- NULL
.brms_cache$backend <- NULL

#' Pick an available brms backend (cmdstanr preferred, else rstan).
choose_brms_backend <- function() {
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    ok <- tryCatch({
      v <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
      !is.null(v) && !is.na(v)
    }, error = function(e) FALSE)
    if (ok) return("cmdstanr")
  }
  if (requireNamespace("rstan", quietly = TRUE)) return("rstan")
  stop("Neither cmdstanr (with cmdstan installed) nor rstan is available.")
}

#' Fit a brms normal-pseudo-likelihood model for a single study.
#'
#' Data is expected to be a single-row tibble with columns `effect` (log-HR)
#' and `se` (SE of log-HR).  Model: `effect | se(se) ~ 1` with a user-supplied
#' Normal prior on the Intercept.
#'
#' The compiled Stan model is cached across calls, so looping through
#' (study, prior) combinations reuses the underlying binary.
#'
#' @return tibble with columns (prior, mu, RR).
fit_brms_study <- function(df, priors, mcid = NULL, chains = 2, iter = 1000,
                           seed = 2026, refresh = 0, backend = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is required.")
  stopifnot(nrow(df) == 1, all(c("effect", "se") %in% names(df)))
  if (is.null(backend)) backend <- choose_brms_backend()

  base_fit     <- .brms_cache$fit
  same_backend <- identical(.brms_cache$backend, backend)

  out <- list()
  for (pname in names(priors)) {
    pr <- priors[[pname]]
    if (is.null(base_fit) || !same_backend) {
      fit <- brms::brm(
        formula = brms::bf(effect | se(se) ~ 1),
        data    = df, family = gaussian(), prior = pr,
        chains  = chains, iter = iter, seed = seed,
        refresh = refresh, backend = backend, silent = 2
      )
      .brms_cache$fit     <- fit
      .brms_cache$backend <- backend
      base_fit <- fit; same_backend <- TRUE
    } else {
      fit <- tryCatch(
        brms::update.brmsfit(base_fit, newdata = df, prior = pr,
                             chains = chains, iter = iter, seed = seed,
                             refresh = refresh, silent = 2, recompile = FALSE),
        error = function(e) brms::brm(
          formula = brms::bf(effect | se(se) ~ 1),
          data = df, family = gaussian(), prior = pr,
          chains = chains, iter = iter, seed = seed,
          refresh = refresh, backend = backend, silent = 2
        )
      )
    }
    draws <- as.data.frame(brms::as_draws_df(fit))
    mu_draws <- draws[["b_Intercept"]]
    out[[pname]] <- tibble::tibble(
      prior = pname, mu = mu_draws, RR = exp(mu_draws)
    )
  }
  dplyr::bind_rows(out)
}

#' Pipeline per-study summary reducer.
#'
#' @param draws tibble with columns (prior, RR).
#' @param mcid  MCID on the RR scale (for `P_RR_le_mcid`).
#' @return tibble with one row per prior: median_RR, lo95, hi95,
#'   P_RR_lt_1, P_RR_le_mcid.
summarise_posterior <- function(draws, mcid = NULL) {
  stopifnot(all(c("prior", "RR") %in% names(draws)))
  grp <- split(draws, draws$prior)
  rows <- lapply(names(grp), function(p) {
    rr <- grp[[p]]$RR
    tibble::tibble(
      prior        = p,
      median_RR    = stats::median(rr),
      lo95         = unname(stats::quantile(rr, 0.025)),
      hi95         = unname(stats::quantile(rr, 0.975)),
      P_RR_lt_1    = mean(rr < 1),
      P_RR_le_mcid = if (!is.null(mcid)) mean(rr <= mcid) else NA_real_
    )
  })
  dplyr::bind_rows(rows)
}

#' Build the five standard priors on log-RR (class = Intercept for brms).
#'
#' - Weak:            N(0, 10^2)
#' - Skeptical:       N(0, 0.354^2)   (~ 95% HR in [0.5, 2])
#' - Enthusiastic:    N(log(mcid_rr), 0.354^2)    -- observed-direction MCID
#' - EB:              N(eb_mu_log, eb_sd_log^2)
#' - DesignInformed:  N(log(design_hr), SE_logHR^2) -- trial-design-direction
#'
#' `mcid_rr` and `design_hr` are usually identical, but for "twist" trials
#' (e.g. DECAF: observed benefit in the opposite direction to the design)
#' they can differ. `design_hr` defaults to `mcid_rr` when not supplied.
build_priors <- function(se_logHR, mcid_rr, eb_mu_log, eb_sd_log,
                         design_hr = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms required.")
  if (is.null(design_hr) || !is.finite(design_hr) || design_hr <= 0) {
    design_hr <- mcid_rr
  }
  list(
    Weak           = brms::prior_string("normal(0, 10)",    class = "Intercept"),
    Skeptical      = brms::prior_string("normal(0, 0.354)", class = "Intercept"),
    Enthusiastic   = brms::prior_string(
      sprintf("normal(%f, 0.354)", log(mcid_rr)),           class = "Intercept"),
    EB             = brms::prior_string(
      sprintf("normal(%f, %f)", eb_mu_log, max(eb_sd_log, 1e-3)),
                                                            class = "Intercept"),
    DesignInformed = brms::prior_string(
      sprintf("normal(%f, %f)", log(design_hr), se_logHR),  class = "Intercept")
  )
}

# --- Pipeline plots ----------------------------------------------------------

#' Faceted density plot of posterior HR draws by prior for a single study.
#'
#' The x-axis defaults to (0.3, 1.2) so benefit-side distributions are
#' readable.  If any prior's 95 % HR quantiles extend past that window,
#' the limits are widened to accommodate (plus a small padding).
plot_study_densities <- function(draws, mcid_rr, study = "",
                                 xlim_default = c(0.3, 1.2)) {
  draws$prior <- factor(draws$prior,
                        levels = c("Weak", "Skeptical", "Enthusiastic",
                                   "EB", "DesignInformed"))
  qs <- stats::quantile(draws$RR, c(0.01, 0.99), na.rm = TRUE)
  xlim <- c(min(xlim_default[1], floor(qs[1] * 10) / 10),
            max(xlim_default[2], ceiling(qs[2] * 10) / 10))
  ggplot2::ggplot(draws, ggplot2::aes(RR)) +
    ggplot2::geom_density(fill = "#3182bd", color = "#3182bd", alpha = 0.3) +
    ggplot2::geom_vline(xintercept = 1,       linetype = "dashed", color = "gray40") +
    ggplot2::geom_vline(xintercept = mcid_rr, linetype = "dotdash", color = "#de2d26") +
    ggplot2::facet_wrap(~ prior, ncol = 3, scales = "free_y") +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(
      title    = study,
      subtitle = sprintf("Dashed = 1 (null); dot-dash = MCID (RR=%.2f)", mcid_rr),
      x = "HR (posterior)", y = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Cross-study forest / dot-whisker of posterior medians and 95 % CrIs.
#'
#' Uses a log-scale x-axis.  Limits default to (0.3, 1.2) and widen only
#' if any study's 95 % CrI extends past them.
plot_forest <- function(summary_df, xlim_default = c(0.3, 1.2)) {
  summary_df$prior <- factor(summary_df$prior,
                             levels = c("Weak", "Skeptical", "Enthusiastic",
                                        "EB", "DesignInformed"))
  xlim <- c(min(xlim_default[1], floor(min(summary_df$lo95, na.rm = TRUE) * 10) / 10),
            max(xlim_default[2], ceiling(max(summary_df$hi95, na.rm = TRUE) * 10) / 10))
  ggplot2::ggplot(summary_df, ggplot2::aes(x = median_RR, y = study,
                                           color = prior)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6),
                        size = 2) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = lo95, xmax = hi95),
      position = ggplot2::position_dodge(width = 0.6), height = 0.2
    ) +
    ggplot2::scale_x_log10(limits = xlim) +
    ggplot2::labs(title = "Posterior HR by prior (forest)",
                  x = "HR (log scale)", y = NULL, color = "Prior") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

#' Heatmap of P(HR <= MCID) by study x prior.
plot_heatmap_mcid <- function(summary_df) {
  summary_df$prior <- factor(summary_df$prior,
                             levels = c("Weak", "Skeptical", "Enthusiastic",
                                        "EB", "DesignInformed"))
  ggplot2::ggplot(summary_df, ggplot2::aes(x = prior, y = study,
                                           fill = P_RR_le_mcid)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", P_RR_le_mcid)),
                       color = "black", size = 3) +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08306b",
                                 limits = c(0, 1), name = "P(HR \u2264 MCID)") +
    ggplot2::labs(title = "Posterior P(HR \u2264 MCID)", x = "Prior", y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}

#' Scatter of observed |z| vs EB-shrunk |SNR posterior mean|.
plot_shrinkage_scatter <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = abs_z, y = abs_snr_post)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "gray60") +
    ggplot2::geom_point(size = 3, color = "#08519c") +
    .ggrepel_or_text(df) +
    ggplot2::labs(title = "EB shrinkage of |z|",
                  subtitle = "y = x is no-shrinkage; below y=x is shrinkage toward 0.",
                  x = "Observed |z|", y = "EB posterior |SNR| (mean)") +
    ggplot2::theme_minimal(base_size = 12)
}

.ggrepel_or_text <- function(df) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(ggplot2::aes(label = study),
                             size = 3, max.overlaps = Inf)
  } else {
    ggplot2::geom_text(ggplot2::aes(label = study), size = 3, vjust = -0.8)
  }
}
