# R/plots.R
#
# ggplot2 helpers for prior/likelihood/posterior overlays and forest plots.
# All helpers take tidy inputs.  Axis limits are data-driven via
# rr_plot_limits() (see R/summaries.R).

#' Prior-density tibble on the RR scale from a Normal(mu0, tau) on log(RR).
#'
#' @param priors named list of c(mu0, tau) pairs.
#' @param rr_grid numeric vector of RR points.
#' @return long tibble with columns prior / RR / density.
prior_density_rr <- function(priors, rr_grid) {
  dplyr::bind_rows(lapply(names(priors), function(nm) {
    pr <- priors[[nm]]
    tibble::tibble(
      prior   = nm,
      RR      = rr_grid,
      density = stats::dnorm(log(rr_grid), pr["mu0"], pr["tau"]) / rr_grid
    )
  }))
}

#' Likelihood density tibble on the RR scale from observed y, sigma.
#'
#' @param y log(RR) point estimate.
#' @param sigma SE[log(RR)].
#' @param rr_grid numeric vector of RR points.
#' @return tibble with columns RR / density.
likelihood_density_rr <- function(y, sigma, rr_grid) {
  tibble::tibble(
    RR      = rr_grid,
    density = stats::dnorm(log(rr_grid), y, sigma) / rr_grid
  )
}

#' Build the "Prior vs Likelihood vs Posterior" facet plot (per prior).
#'
#' @param posterior_df tibble (prior, RR) of posterior draws.
#' @param priors named list of c(mu0, tau) as in prior_density_rr().
#' @param obs_y,obs_sigma observed log(RR) and SE (for likelihood overlay).
#' @param mcid numeric MCID threshold for shading and reference line.
#' @return ggplot2 object faceted by prior.
plv_plot <- function(posterior_df, priors, obs_y, obs_sigma, mcid = 0.80) {
  lims <- rr_plot_limits(posterior_df$RR)
  rr_grid <- exp(seq(log(lims[1]), log(lims[2]), length.out = 400))

  pr_df <- prior_density_rr(priors, rr_grid)
  lk_df <- if (is.finite(obs_y) && is.finite(obs_sigma))
    likelihood_density_rr(obs_y, obs_sigma, rr_grid) else NULL

  # posterior density per prior via kernel density on RR draws
  post_df <- posterior_df %>%
    dplyr::group_by(.data$prior) %>%
    dplyr::group_modify(~{
      d <- stats::density(.x$RR, from = lims[1], to = lims[2], n = 400)
      tibble::tibble(RR = d$x, density = d$y)
    }) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = 0, xmax = mcid, ymin = 0, ymax = Inf,
                      alpha = 0.08, fill = "steelblue") +
    ggplot2::geom_line(data = pr_df,
                       ggplot2::aes(.data$RR, .data$density, colour = "Prior"),
                       linewidth = 0.8) +
    ggplot2::geom_line(data = post_df,
                       ggplot2::aes(.data$RR, .data$density, colour = "Posterior"),
                       linewidth = 0.9)

  if (!is.null(lk_df)) {
    p <- p + ggplot2::geom_line(
      data = lk_df,
      ggplot2::aes(.data$RR, .data$density, colour = "Likelihood"),
      linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE
    )
  }

  p +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = mcid, linetype = "dotted",
                        colour = "steelblue") +
    ggplot2::facet_wrap(~ .data$prior) +
    ggplot2::scale_colour_manual(
      values = c(Prior = "#756bb1", Likelihood = "#e6550d",
                 Posterior = "#31a354"),
      name = NULL
    ) +
    ggplot2::coord_cartesian(xlim = lims) +
    ggplot2::labs(title = "Prior vs Likelihood vs Posterior",
                  x = "RR", y = "Density") +
    ggplot2::theme_minimal()
}

#' Forest-style plot of posterior median + 95% CrI per prior.
#'
#' @param posterior_df tibble (prior, RR).
#' @param mcid MCID threshold.
#' @return ggplot2 object.
forest_plot <- function(posterior_df, mcid = 0.80) {
  summ <- posterior_df %>%
    dplyr::group_by(.data$prior) %>%
    dplyr::summarise(
      med = stats::median(.data$RR),
      lo  = stats::quantile(.data$RR, 0.025),
      hi  = stats::quantile(.data$RR, 0.975),
      .groups = "drop"
    )
  lims <- rr_plot_limits(posterior_df$RR)
  ggplot2::ggplot(summ,
                  ggplot2::aes(x = .data$med, y = .data$prior)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = mcid, linetype = "dotted",
                        colour = "steelblue") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
      height = 0.15, linewidth = 0.7
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::coord_cartesian(xlim = lims) +
    ggplot2::labs(x = "RR (posterior median, 95% CrI)", y = NULL,
                  title = "Forest: posterior RR by prior") +
    ggplot2::theme_minimal()
}

# ---------- Risk Difference plot helpers ----------

#' Prior-density tibble on the RD (proportion) scale.
#'
#' Priors are Normal(mu0, tau) directly on the RD proportion scale, so
#' no Jacobian transform is needed.
#'
#' @param priors named list of c(mu0, tau) pairs on the RD proportion scale.
#' @param rd_grid numeric vector of RD points (proportion scale).
#' @return long tibble with columns prior / RD / density.
prior_density_rd <- function(priors, rd_grid) {
  dplyr::bind_rows(lapply(names(priors), function(nm) {
    pr <- priors[[nm]]
    tibble::tibble(
      prior   = nm,
      RD      = rd_grid,
      density = stats::dnorm(rd_grid, pr["mu0"], pr["tau"])
    )
  }))
}

#' Likelihood density on the RD (proportion) scale.
#'
#' @param rd  Observed risk difference (proportion scale).
#' @param se_rd SE of the risk difference.
#' @param rd_grid numeric vector of RD points.
#' @return tibble with columns RD / density.
likelihood_density_rd <- function(rd, se_rd, rd_grid) {
  tibble::tibble(
    RD      = rd_grid,
    density = stats::dnorm(rd_grid, rd, se_rd)
  )
}

#' Prior vs Likelihood vs Posterior facet plot on the RD scale.
#'
#' @param posterior_df tibble (prior, RD) of posterior draws (proportion).
#' @param priors named list of c(mu0, tau) on the proportion scale.
#' @param obs_rd,obs_se_rd observed RD and SE (proportion scale).
#' @param mcid_rd MCID threshold (proportion scale, e.g. -0.02).
#' @return ggplot2 object faceted by prior.
plv_plot_rd <- function(posterior_df, priors, obs_rd, obs_se_rd,
                        mcid_rd = -0.02) {
  lims <- rd_plot_limits(posterior_df$RD)
  rd_grid <- seq(lims[1], lims[2], length.out = 400)

  pr_df <- prior_density_rd(priors, rd_grid)
  lk_df <- if (is.finite(obs_rd) && is.finite(obs_se_rd))
    likelihood_density_rd(obs_rd, obs_se_rd, rd_grid) else NULL

  post_df <- posterior_df %>%
    dplyr::group_by(.data$prior) %>%
    dplyr::group_modify(~{
      d <- stats::density(.x$RD, from = lims[1], to = lims[2], n = 400)
      tibble::tibble(RD = d$x, density = d$y)
    }) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = -Inf, xmax = mcid_rd,
                      ymin = 0, ymax = Inf,
                      alpha = 0.08, fill = "steelblue") +
    ggplot2::geom_line(data = pr_df,
                       ggplot2::aes(.data$RD, .data$density, colour = "Prior"),
                       linewidth = 0.8) +
    ggplot2::geom_line(data = post_df,
                       ggplot2::aes(.data$RD, .data$density, colour = "Posterior"),
                       linewidth = 0.9)

  if (!is.null(lk_df)) {
    p <- p + ggplot2::geom_line(
      data = lk_df,
      ggplot2::aes(.data$RD, .data$density, colour = "Likelihood"),
      linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE
    )
  }

  p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = mcid_rd, linetype = "dotted",
                        colour = "steelblue") +
    ggplot2::facet_wrap(~ .data$prior) +
    ggplot2::scale_colour_manual(
      values = c(Prior = "#756bb1", Likelihood = "#e6550d",
                 Posterior = "#31a354"),
      name = NULL
    ) +
    ggplot2::coord_cartesian(xlim = lims) +
    ggplot2::labs(title = "Prior vs Likelihood vs Posterior (RD)",
                  x = "Risk Difference", y = "Density") +
    ggplot2::theme_minimal()
}

#' Forest-style plot of posterior median + 95% CrI per prior (RD scale).
#'
#' @param posterior_df tibble (prior, RD) on the proportion scale.
#' @param mcid_rd MCID threshold (proportion scale).
#' @return ggplot2 object.
forest_plot_rd <- function(posterior_df, mcid_rd = -0.02) {
  summ <- posterior_df %>%
    dplyr::group_by(.data$prior) %>%
    dplyr::summarise(
      med = stats::median(.data$RD),
      lo  = stats::quantile(.data$RD, 0.025),
      hi  = stats::quantile(.data$RD, 0.975),
      .groups = "drop"
    )
  lims <- rd_plot_limits(posterior_df$RD)
  ggplot2::ggplot(summ,
                  ggplot2::aes(x = .data$med, y = .data$prior)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = mcid_rd, linetype = "dotted",
                        colour = "steelblue") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
      height = 0.15, linewidth = 0.7
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::coord_cartesian(xlim = lims) +
    ggplot2::labs(x = "RD (posterior median, 95% CrI)", y = NULL,
                  title = "Forest: posterior RD by prior") +
    ggplot2::theme_minimal()
}
