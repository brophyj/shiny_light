# app_bayes_alone.R -- Bayesian-only re-analysis of a single RCT
#
# Slimmer counterpart to `app_full.R`.  Use this when you want the
# Bayesian re-analysis without the frequentist calibration extras
# (S-values, Type M/S retrodesign, per-prior posterior calibration).
# The empirical-Bayes prior is *retained* here -- it is itself a
# Bayesian construction.
#
# Features:
#   * Four priors: Weak / Enthusiastic / Skeptical / Empirical-Bayes
#     (van Zwet et al. 2021/2025; see R/eb_mixture.R).
#   * Three likelihoods: Normal (analytic), Binomial (Laplace), MCMC
#     (brms + cmdstanr) on demand.
#   * MCID-aware posterior tables and plots, including a
#     "Prior vs Likelihood vs Posterior" facet tab and a forest tab.
#   * Sidebar shows the basic frequentist observed RR / 95% CI / z / p
#     for orientation, but no S-values or Type M/S calibration.
#
# For the full-feature build with calibration metrics and the
# Design-informed prior, see `app_full.R`.
#
# All posterior draws are returned in a long tibble with columns
# (prior, mu, RR) so downstream summary/plot code is likelihood-agnostic.

library(shiny)
library(ggplot2)
library(dplyr)
library(ggdist)
library(tibble)
library(bslib)
library(rlang)

# ---------- source helpers ----------
local({
  here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
  if (is.null(here) || !nzchar(here)) here <- getwd()
  for (f in c("priors.R", "likelihoods.R", "eb_mixture.R",
              "summaries.R", "plots.R")) {
    src <- file.path(here, "R", f)
    if (!file.exists(src)) src <- file.path("R", f)
    source(src, local = FALSE)
  }
})

# ---------- optional MCMC backend ----------
# We pre-compile a generic binomial regression once at app startup so
# subsequent fits reuse the Stan binary via brms::update().
has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
has_brms     <- requireNamespace("brms",     quietly = TRUE)
has_memoise  <- requireNamespace("memoise",  quietly = TRUE)
has_rstan    <- requireNamespace("rstan",    quietly = TRUE)

#' Ensure a working cmdstan installation exists; safe to call repeatedly.
#' Returns TRUE if cmdstan can be located, FALSE otherwise.
ensure_cmdstan <- local({
  checked <- FALSE
  ok <- FALSE
  function() {
    if (checked) return(ok)
    checked <<- TRUE
    if (!has_cmdstanr) { ok <<- FALSE; return(FALSE) }
    ok <<- tryCatch({
      path <- cmdstanr::cmdstan_path()
      !is.null(path) && nzchar(path) && dir.exists(path)
    }, error = function(e) FALSE)
    if (!ok) {
      ok <<- tryCatch({
        cmdstanr::install_cmdstan(cores = 2, overwrite = FALSE)
        path <- cmdstanr::cmdstan_path()
        !is.null(path) && nzchar(path) && dir.exists(path)
      }, error = function(e) FALSE)
    }
    ok
  }
})

mcmc_available <- function() {
  has_brms && ((has_cmdstanr && isTRUE(ensure_cmdstan())) || has_rstan)
}

# Pre-compiled brms model (fitted once against a tiny placeholder dataset).
.brms_base <- NULL
init_brms_base <- function() {
  if (!has_brms) return(NULL)
  if (!is.null(.brms_base)) return(.brms_base)
  backend <- if (has_cmdstanr && ensure_cmdstan()) "cmdstanr"
             else if (has_rstan) "rstan" else return(NULL)
  placeholder <- data.frame(
    events = c(1, 1),
    n      = c(10, 10),
    arm    = factor(c("control", "treatment"),
                    levels = c("control", "treatment"))
  )
  fit <- tryCatch(
    brms::brm(
      events | trials(n) ~ arm,
      data = placeholder,
      family = brms::brmsfamily("binomial", link = "logit"),
      backend = backend,
      chains = 1, iter = 500, warmup = 200, refresh = 0, silent = 2,
      seed = 1
    ),
    error = function(e) NULL
  )
  .brms_base <<- fit
  fit
}

#' Run brms on the provided RCT given a named prior on the arm-effect.
#' @return tibble with columns (prior, mu, RR).
run_mcmc_fit <- function(e_t, n_t, e_c, n_c, prior_name,
                         m_alpha, s_alpha, chains = 2, iter = 2000,
                         base = NULL) {
  if (is.null(base)) base <- .brms_base
  if (is.null(base)) stop("brms base model not compiled")
  dat <- data.frame(
    events = c(e_c, e_t),
    n      = c(n_c, n_t),
    arm    = factor(c("control", "treatment"),
                    levels = c("control", "treatment"))
  )
  pr <- brms::prior_string(
    sprintf("normal(%.6f, %.6f)", m_alpha, max(s_alpha, 1e-3)),
    class = "b"
  )
  fit <- brms::update(base,
                      newdata = dat,
                      prior = pr,
                      chains = chains, iter = iter, warmup = floor(iter / 2),
                      refresh = 0, silent = 2)
  dr <- posterior::as_draws_df(fit)
  cand <- grep("^b_arm", colnames(dr), value = TRUE)
  if (!length(cand)) stop("MCMC: could not find arm-effect coefficient in draws")
  alpha_draws <- dr[[cand[1]]]
  b_draws     <- dr[["b_Intercept"]]
  p_c <- plogis(b_draws)
  p_t <- plogis(b_draws + alpha_draws)
  RR  <- p_t / p_c
  RD  <- p_t - p_c
  tibble(prior = prior_name, mu = log(RR), RR = RR, RD = RD)
}

mcmc_fit <- if (has_memoise) memoise::memoise(run_mcmc_fit) else run_mcmc_fit

# ---------- UI ----------
ui <- fluidPage(
  theme = bs_theme(version = 5),
  titlePanel("Bayesian Interpretation of a Randomized Clinical Trial"),
  sidebarLayout(
    sidebarPanel(
      h4("Likelihood"),
      radioButtons(
        "likelihood", NULL,
        choices = c(
          "Normal (analytic)" = "normal",
          "Binomial (Laplace, no normality)" = "binom",
          "MCMC (brms + cmdstanr)" = "mcmc"
        ),
        selected = "normal"
      ),
      uiOutput("mcmc_note"),
      conditionalPanel(
        condition = "input.likelihood == 'mcmc'",
        numericInput("mcmc_chains", "MCMC chains", value = 2, min = 1, max = 8, step = 1),
        numericInput("mcmc_iter",   "MCMC iterations (per chain)",
                     value = 2000, min = 500, step = 500),
        actionButton("run_mcmc", "Run MCMC", class = "btn-primary")
      ),
      hr(),
      h4("Observed RCT counts"),
      numericInput("n_t", "Treatment n1", value = 3528, min = 2, step = 1),
      numericInput("e_t", "Treatment events e1", value = 322, min = 0, step = 1),
      numericInput("n_c", "Control n2",   value = 3534, min = 2, step = 1),
      numericInput("e_c", "Control events e2", value = 327, min = 0, step = 1),
      wellPanel(
        style = "padding: 8px 12px; margin-top: 6px;",
        strong("Frequentist summary"),
        htmlOutput("obs_freq")
      ),
      hr(),
      h4("Effect measure & priors"),
      tabsetPanel(
        id = "effect_measure", type = "pills",
        tabPanel("RR",
          radioButtons("prior_mode", NULL,
                       choices = c("RR & 95% CI" = "rrci",
                                   "Mean log(RR) & SD" = "logsd"),
                       selected = "rrci", inline = TRUE),
          h5("Weak prior"),         uiOutput("weak_ui"),
          h5("Enthusiastic prior"), uiOutput("enth_ui"),
          h5("Skeptical prior"),    uiOutput("skep_ui"),
          helpText("Empirical-Bayes prior (van Zwet et al. 2025) has ",
                   "no user-tunable parameters.")
        ),
        tabPanel("RD",
          radioButtons("rd_prior_mode", NULL,
                       choices = c("RD (pp) & 95% CI" = "rdci",
                                   "Mean RD (pp) & SD" = "rdsd"),
                       selected = "rdsd", inline = TRUE),
          helpText("All values in percentage points (pp)."),
          h5("Weak prior"),         uiOutput("weak_rd_ui"),
          h5("Enthusiastic prior"), uiOutput("enth_rd_ui"),
          h5("Skeptical prior"),    uiOutput("skep_rd_ui"),
          helpText("Empirical-Bayes prior (van Zwet et al. 2025) has ",
                   "no user-tunable parameters.")
        )
      ),
      hr(),
      h4("Clinically meaningful benefit (MCID)"),
      conditionalPanel(
        condition = "input.effect_measure == 'RR'",
        helpText(
          "The MCID is the smallest effect considered meaningful. ",
          "Posterior tables report ", tags$code("P(RR <= MCID)"), "."
        ),
        numericInput("mcid_rr", "MCID (RR)", value = 0.80, step = 0.01, min = 1e-3),
        helpText(
          tags$em("Optional:"), " leave blank/<=0 and fill in power-calc ",
          "inputs below to derive MCID from the original design."
        ),
        numericInput("pow_n_per_arm", "Power calc: n per arm", value = NA, min = 1, step = 1),
        numericInput("pow_alpha",     "Power calc: alpha",     value = 0.05, min = 1e-4, step = 0.01),
        numericInput("pow_target",    "Power calc: target power",
                     value = 0.80, min = 0.01, max = 0.999, step = 0.01)
      ),
      conditionalPanel(
        condition = "input.effect_measure == 'RD'",
        helpText(
          "The MCID is the smallest absolute risk difference considered ",
          "meaningful. Posterior tables report ", tags$code("P(RD <= MCID)"),
          ". Enter in percentage points (pp)."
        ),
        numericInput("mcid_rd", "MCID (RD, pp)", value = -2, step = 0.1)
      ),
      hr(),
      h4("Predictive settings"),
      conditionalPanel(
        condition = "input.effect_measure == 'RR'",
        helpText("SE for a future study's log(RR). Leave blank/<=0 to reuse observed SE."),
        numericInput("se_future", "SE_new (log RR)", value = NA, min = 0, step = 0.0001)
      ),
      conditionalPanel(
        condition = "input.effect_measure == 'RD'",
        helpText("SE for a future study's RD. Leave blank/<=0 to reuse observed SE(RD)."),
        numericInput("se_future_rd", "SE_new (RD, proportion)", value = NA, min = 0, step = 0.0001)
      ),
      hr(),
      downloadButton("download_csv", "Download summary CSV")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Observed / Summary",
                 h4("Observed (from counts)"),
                 verbatimTextOutput("obs_text"),
                 h5("MCID in use"), verbatimTextOutput("mcid_text")),
        tabPanel("Priors",
                 conditionalPanel("input.effect_measure == 'RR'",
                   plotOutput("priors_plot", height = 320)),
                 conditionalPanel("input.effect_measure == 'RD'",
                   plotOutput("priors_plot_rd", height = 320))),
        tabPanel("Posteriors",
                 conditionalPanel("input.effect_measure == 'RR'",
                   plotOutput("post_plot", height = 360),
                   br(), strong("Posterior numeric summaries"),
                   tableOutput("post_tbl_small"),
                   uiOutput("mcid_footnote_small")),
                 conditionalPanel("input.effect_measure == 'RD'",
                   plotOutput("post_plot_rd", height = 360),
                   br(), strong("Posterior numeric summaries"),
                   tableOutput("post_tbl_small_rd"),
                   uiOutput("mcid_footnote_small_rd"))),
        tabPanel("Predictive",
                 conditionalPanel("input.effect_measure == 'RR'",
                   plotOutput("pred_plot", height = 360),
                   br(), strong("Predictive numeric summaries"),
                   tableOutput("pred_tbl_small")),
                 conditionalPanel("input.effect_measure == 'RD'",
                   plotOutput("pred_plot_rd", height = 360),
                   br(), strong("Predictive numeric summaries"),
                   tableOutput("pred_tbl_small_rd"))),
        tabPanel("Prior vs Likelihood vs Posterior",
                 conditionalPanel("input.effect_measure == 'RR'",
                   plotOutput("plv_plot", height = 420)),
                 conditionalPanel("input.effect_measure == 'RD'",
                   plotOutput("plv_plot_rd", height = 420))),
        tabPanel("Forest",
                 conditionalPanel("input.effect_measure == 'RR'",
                   plotOutput("forest_plot", height = 320)),
                 conditionalPanel("input.effect_measure == 'RD'",
                   plotOutput("forest_plot_rd", height = 320))),
        tabPanel("Tables",
                 conditionalPanel("input.effect_measure == 'RR'",
                   h5("Posterior probabilities / intervals"),
                   tableOutput("post_tbl"),
                   uiOutput("mcid_footnote"),
                   h5("Predictive intervals"), tableOutput("pred_tbl")),
                 conditionalPanel("input.effect_measure == 'RD'",
                   h5("Posterior probabilities / intervals"),
                   tableOutput("post_tbl_rd"),
                   uiOutput("mcid_footnote_rd"),
                   h5("Predictive intervals"), tableOutput("pred_tbl_rd")))
      )
    )
  )
)

# ---------- server ----------
server <- function(input, output, session) {

  output$mcmc_note <- renderUI({
    if (mcmc_available()) {
      helpText("MCMC backend available.")
    } else {
      tagList(
        helpText(em("MCMC disabled: cmdstanr/rstan not available in this R session.")),
        tags$script("setTimeout(function(){var e=document.querySelector(\"input[value='mcmc']\");if(e)e.disabled=true;},50)")
      )
    }
  })

  # ----- dynamic prior UIs -----
  prior_inputs <- function(id_prefix, default_rr = 1.0, default_lo = 0.8, default_hi = 1.2,
                           default_log = 0, default_sd = 2) {
    if (identical(input$prior_mode, "rrci")) {
      tagList(
        numericInput(paste0(id_prefix, "_rr"), "RR (central)", value = default_rr, step = 0.01),
        numericInput(paste0(id_prefix, "_lo"), "Lower 95% RR", value = default_lo, step = 0.01),
        numericInput(paste0(id_prefix, "_hi"), "Upper 95% RR", value = default_hi, step = 0.01)
      )
    } else {
      tagList(
        numericInput(paste0(id_prefix, "_m"),  "Mean log(RR)", value = default_log, step = 0.01),
        numericInput(paste0(id_prefix, "_sd"), "SD of log(RR)", value = default_sd, min = 1e-4, step = 0.01)
      )
    }
  }
  output$weak_ui <- renderUI(prior_inputs("weak", 1.00, 0.50, 2.00, 0, 2))
  output$enth_ui <- renderUI(prior_inputs("enth", 0.80, 0.70, 0.92,
                                          log(0.80), sigma_from_ci(0.70, 0.92)))
  output$skep_ui <- renderUI(prior_inputs("skep", 1.10, 1.00, 1.21,
                                          log(1.10), sigma_from_ci(1.00, 1.21)))

  # ----- dynamic RD prior UIs (percentage points) -----
  rd_prior_inputs <- function(id_prefix, default_cen = 0, default_lo = -4, default_hi = 4,
                              default_mean = 0, default_sd = 2) {
    if (identical(input$rd_prior_mode, "rdci")) {
      tagList(
        numericInput(paste0(id_prefix, "_cen"), "RD (pp, central)", value = default_cen, step = 0.1),
        numericInput(paste0(id_prefix, "_lo"),  "Lower 95% (pp)",  value = default_lo,  step = 0.1),
        numericInput(paste0(id_prefix, "_hi"),  "Upper 95% (pp)",  value = default_hi,  step = 0.1)
      )
    } else {
      tagList(
        numericInput(paste0(id_prefix, "_m"),  "Mean RD (pp)", value = default_mean, step = 0.1),
        numericInput(paste0(id_prefix, "_sd"), "SD (pp)",      value = default_sd, min = 1e-4, step = 0.1)
      )
    }
  }
  output$weak_rd_ui <- renderUI(rd_prior_inputs("rd_weak", 0, -4, 4, 0, 2))
  output$enth_rd_ui <- renderUI(rd_prior_inputs("rd_enth", -1, -2, 0, -1, 0.5))
  output$skep_rd_ui <- renderUI(rd_prior_inputs("rd_skep",  1,  0, 2,  1, 0.5))

  # ----- observed -----
  obs <- reactive(obs_from_counts(input$e_t, input$n_t, input$e_c, input$n_c))

  # Derived frequentist summary: RR, 95% CI on RR, z, two-sided p.
  obs_freq_parts <- reactive({
    o <- obs()
    if (!is.finite(o$rr) || !is.finite(o$y) || !is.finite(o$sigma)) {
      return(list(rr = NA, lo = NA, hi = NA, z = NA, p = NA))
    }
    zcrit <- qnorm(0.975)
    lo <- exp(o$y - zcrit * o$sigma)
    hi <- exp(o$y + zcrit * o$sigma)
    z  <- o$y / o$sigma
    p  <- 2 * pnorm(-abs(z))
    list(rr = o$rr, lo = lo, hi = hi, z = z, p = p)
  })

  # Sidebar mini-summary (HTML, so we can use &ndash; and italic p).
  output$obs_freq <- renderUI({
    f <- obs_freq_parts(); o <- obs()
    if (!is.finite(f$rr)) return(HTML("<em>RR not defined (zero events in an arm).</em>"))
    rr_line <- sprintf(
      "RR = <strong>%.3f</strong> (95%%%% CI %.3f &ndash; %.3f)",
      f$rr, f$lo, f$hi
    )
    rd_line <- ""
    if (is.finite(o$rd) && is.finite(o$se_rd) && o$se_rd > 0) {
      zcrit <- qnorm(0.975)
      rd_line <- sprintf(
        "<br>RD = <strong>%.2f pp</strong> (95%%%% CI %.2f &ndash; %.2f)",
        o$rd * 100, (o$rd - zcrit * o$se_rd) * 100, (o$rd + zcrit * o$se_rd) * 100
      )
    }
    z_line <- sprintf(
      "<br><em>z</em> = %+.2f, <em>p</em> = %s",
      f$z, ifelse(f$p < 1e-4, sprintf("%.1e", f$p), sprintf("%.3f", f$p))
    )
    HTML(paste0(rr_line, rd_line, z_line))
  })

  # Main-panel observed block (verbatim) — always shows both RR and RD.
  output$obs_text <- renderText({
    o <- obs(); f <- obs_freq_parts()
    rr_str <- paste0(
      sprintf("Observed RR = %s",
              ifelse(is.finite(o$rr), sprintf("%.3f", o$rr), "NA")),
      ifelse(is.finite(f$lo) && is.finite(f$hi),
             sprintf("  (95%% CI %.3f - %.3f)", f$lo, f$hi), ""),
      "\n",
      sprintf("Observed log(RR) = %s",
              ifelse(is.finite(o$y), sprintf("%.3f", o$y), "NA")), "\n",
      sprintf("SE[log(RR)]      = %s",
              ifelse(is.finite(o$sigma), sprintf("%.5f", o$sigma), "NA"))
    )
    rd_str <- ""
    if (is.finite(o$rd) && is.finite(o$se_rd) && o$se_rd > 0) {
      zcrit <- qnorm(0.975)
      rd_str <- sprintf(
        "\nObserved RD = %.2f pp  (95%%%% CI %.2f - %.2f)\nSE[RD]      = %.5f",
        o$rd * 100, (o$rd - zcrit * o$se_rd) * 100,
        (o$rd + zcrit * o$se_rd) * 100, o$se_rd
      )
    }
    paste0(rr_str, rd_str, "\n",
      sprintf("z = %s,  p (2-sided) = %s",
              ifelse(is.finite(f$z), sprintf("%+.3f", f$z), "NA"),
              ifelse(is.finite(f$p),
                     ifelse(f$p < 1e-4, sprintf("%.2e", f$p), sprintf("%.4f", f$p)),
                     "NA"))
    )
  })

  # ----- MCID -----
  mcid <- reactive({
    resolve_mcid(
      mcid_rr   = input$mcid_rr,
      n_per_arm = input$pow_n_per_arm,
      alpha     = if (isTruthy(input$pow_alpha))  input$pow_alpha  else 0.05,
      power     = if (isTruthy(input$pow_target)) input$pow_target else 0.80,
      p_ctrl    = {
        o <- obs(); if (is.finite(o$p2)) o$p2 else 0.10
      }
    )
  })
  output$mcid_text <- renderText({
    if (identical(input$effect_measure, "RD")) {
      m <- mcid_rd_val()
      sprintf("%.2f pp (RD scale)", m * 100)
    } else {
      sprintf("%.3f (RR scale)", mcid())
    }
  })

  # ----- prior parameters on log(RR) -----
  prior_pars_logRR <- function(prefix) {
    if (identical(input$prior_mode, "rrci")) {
      rr <- input[[paste0(prefix, "_rr")]]
      lo <- input[[paste0(prefix, "_lo")]]
      hi <- input[[paste0(prefix, "_hi")]]
      validate(need(all(is.finite(c(rr, lo, hi))) && rr > 0 && lo > 0 && hi > lo,
                    "RR & CI must be positive and ordered"))
      mu0 <- log(rr); tau <- sigma_from_ci(lo, hi)
    } else {
      mu0 <- input[[paste0(prefix, "_m")]]
      tau <- input[[paste0(prefix, "_sd")]]
      validate(need(is.finite(mu0) && is.finite(tau) && tau > 0,
                    "Prior SD must be > 0"))
    }
    c(mu0 = mu0, tau = tau)
  }

  all_priors_logRR <- reactive({
    o <- obs()
    se <- if (is.finite(o$sigma)) o$sigma else 0.1
    list(
      Weak              = prior_pars_logRR("weak"),
      Enthusiastic      = prior_pars_logRR("enth"),
      Skeptical         = prior_pars_logRR("skep"),
      `Empirical Bayes` = eb_marginal_logRR_prior(se)
    )
  })

  # ----- priors plot (RR scale) -----
  output$priors_plot <- renderPlot({
    priors <- all_priors_logRR()
    taus  <- sapply(priors, function(p) unname(p["tau"]))
    rr_max <- exp(2 * max(taus, 0.4, na.rm = TRUE))
    lims <- c(max(1 / rr_max, 0.05), rr_max)
    rr_vals <- exp(seq(log(lims[1]), log(lims[2]), length.out = 600))
    pr_df <- prior_density_rr(priors, rr_vals)
    ggplot(pr_df, aes(RR, density, color = prior)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = 1, linetype = "dashed") +
      geom_vline(xintercept = mcid(), linetype = "dotted", colour = "steelblue") +
      labs(title = "Prior distributions (RR scale)", x = "RR", y = "Density") +
      coord_cartesian(xlim = lims) +
      theme_minimal()
  })

  # ----- posterior draws -----
  posterior_df <- reactive({
    o <- obs()
    like <- input$likelihood

    if (identical(like, "normal")) {
      validate(need(is.finite(o$y) && is.finite(o$sigma),
                    "Normal path needs nonzero events in both groups; otherwise use Binomial (Laplace)."))
      pars <- list(
        Weak         = prior_pars_logRR("weak"),
        Enthusiastic = prior_pars_logRR("enth"),
        Skeptical    = prior_pars_logRR("skep")
      )
      out <- lapply(names(pars), function(nm) {
        pr <- pars[[nm]]
        post <- post_from_prior_normal(o$y, o$sigma, pr["mu0"], pr["tau"])
        mu <- rnorm(40000, post$mu, post$sd)
        tibble(prior = nm, mu = mu, RR = exp(mu))
      })
      eb_mu <- eb_posterior_logRR_draws(o$y, o$sigma, ndraws = 40000)
      out$eb <- tibble(prior = "Empirical Bayes", mu = eb_mu, RR = exp(eb_mu))
      bind_rows(out)

    } else if (identical(like, "binom")) {
      p2 <- if (is.finite(o$p2)) o$p2 else (input$e_c + input$e_t) / (input$n_c + input$n_t)
      p2 <- min(max(p2, 1e-6), 1 - 1e-6)
      se_eff <- if (is.finite(o$sigma)) o$sigma else 0.1
      pars <- list(
        Weak         = prior_pars_logRR("weak"),
        Enthusiastic = prior_pars_logRR("enth"),
        Skeptical    = prior_pars_logRR("skep")
      )
      pars$`Empirical Bayes` <- eb_marginal_logRR_prior(se_eff)

      draw_one <- function(pr, lab) {
        map_pars <- map_logRR_prior_to_alpha(pr["mu0"], pr["tau"], p2)
        dr <- laplace_binomial(input$e_t, input$n_t, input$e_c, input$n_c,
                               m_alpha = map_pars["m_alpha"],
                               s_alpha = map_pars["s_alpha"],
                               m_beta  = qlogis(p2), s_beta = 2,
                               ndraws = 40000)
        out <- mutate(as_tibble(dr), prior = lab) %>% select(prior, mu, RR)
        # For EB, refine by forming a pseudo-z from the Laplace posterior mean/SD
        # on log-RR and re-shrinking with the EB mixture.
        if (identical(lab, "Empirical Bayes")) {
          post_mean <- mean(out$mu); post_sd <- stats::sd(out$mu)
          if (is.finite(post_mean) && is.finite(post_sd) && post_sd > 0) {
            fit <- eb_posterior_snr(post_mean / post_sd)
            snr <- fit$sampler(nrow(out))
            out$mu <- snr * post_sd
            out$RR <- exp(out$mu)
          }
        }
        out
      }
      bind_rows(lapply(names(pars), function(nm) draw_one(pars[[nm]], nm)))

    } else if (identical(like, "mcmc")) {
      validate(need(mcmc_available(),
                    "MCMC backend not available in this R session."))
      validate(need(input$run_mcmc > 0,
                    "Press 'Run MCMC' to generate draws."))
      isolate({
        if (is.null(.brms_base)) init_brms_base()
        validate(need(!is.null(.brms_base),
                      "Failed to compile brms base model."))

        p2 <- if (is.finite(o$p2)) o$p2 else (input$e_c + input$e_t) / (input$n_c + input$n_t)
        p2 <- min(max(p2, 1e-6), 1 - 1e-6)
        se_eff <- if (is.finite(o$sigma)) o$sigma else 0.1

        priors <- list(
          Weak         = prior_pars_logRR("weak"),
          Enthusiastic = prior_pars_logRR("enth"),
          Skeptical    = prior_pars_logRR("skep")
        )
        priors$`Empirical Bayes` <- eb_marginal_logRR_prior(se_eff)

        nm_list <- names(priors)
        shiny::withProgress(
          message = "Running MCMC", value = 0,
          {
            res <- list()
            for (i in seq_along(nm_list)) {
              nm <- nm_list[i]
              pr <- priors[[nm]]
              map_pars <- map_logRR_prior_to_alpha(pr["mu0"], pr["tau"], p2)
              res[[nm]] <- mcmc_fit(
                e_t = input$e_t, n_t = input$n_t,
                e_c = input$e_c, n_c = input$n_c,
                prior_name = nm,
                m_alpha = map_pars["m_alpha"],
                s_alpha = map_pars["s_alpha"],
                chains  = input$mcmc_chains,
                iter    = input$mcmc_iter
              )
              shiny::incProgress(1 / length(nm_list),
                                 detail = sprintf("%s (%d/%d)",
                                                  nm, i, length(nm_list)))
            }
            bind_rows(res)
          }
        )
      })
    }
  })

  # ----- predictive draws -----
  predictive_df <- reactive({
    df <- posterior_df()
    o  <- obs()
    se_new <- input$se_future
    if (is.na(se_new) || se_new <= 0) se_new <- if (is.finite(o$sigma)) o$sigma else 0.1
    mu_pred <- rnorm(nrow(df), mean = df$mu, sd = se_new)
    tibble(prior = df$prior, RR_pred = exp(mu_pred))
  })

  # ----- plots -----
  output$post_plot <- renderPlot({
    df <- posterior_df()
    lims <- rr_plot_limits(df$RR)
    ggplot(df, aes(x = RR, fill = prior, color = prior)) +
      annotate("rect", xmin = 0, xmax = mcid(), ymin = -Inf, ymax = Inf,
               alpha = 0.08, fill = "steelblue") +
      stat_slabinterval(aes(thickness = after_stat(pdf)), .width = 0.95,
                        alpha = 0.5, position = "identity") +
      geom_vline(xintercept = 1, linetype = "dashed") +
      geom_vline(xintercept = mcid(), linetype = "dotted", colour = "steelblue") +
      labs(title = "Posterior RR under each prior", x = "RR", y = "Density") +
      coord_cartesian(xlim = lims) +
      theme_minimal()
  })
  output$pred_plot <- renderPlot({
    df <- predictive_df()
    lims <- rr_plot_limits(df$RR_pred)
    ggplot(df, aes(x = RR_pred, fill = prior, color = prior)) +
      stat_slabinterval(aes(thickness = after_stat(pdf)), .width = 0.95,
                        alpha = 0.5) +
      geom_vline(xintercept = 1, linetype = "dashed") +
      geom_vline(xintercept = mcid(), linetype = "dotted", colour = "steelblue") +
      labs(title = "Posterior predictive RR for a future study",
           x = "RR", y = "Density") +
      coord_cartesian(xlim = lims) +
      theme_minimal()
  })
  output$plv_plot <- renderPlot({
    df <- posterior_df()
    priors <- all_priors_logRR()
    o <- obs()
    plv_plot(df, priors, obs_y = o$y, obs_sigma = o$sigma, mcid = mcid())
  })
  output$forest_plot <- renderPlot({
    df <- posterior_df()
    forest_plot(df, mcid = mcid())
  })

  # ----- summary tables -----
  post_summary_tbl <- reactive(summarise_posterior_rr(posterior_df(), mcid = mcid()))
  pred_summary_tbl <- reactive(summarise_predictive_rr(predictive_df()))

  output$post_tbl_small <- renderTable(post_summary_tbl(), digits = 3, striped = TRUE, hover = TRUE)
  output$pred_tbl_small <- renderTable(pred_summary_tbl(), digits = 2, striped = TRUE, hover = TRUE)
  output$post_tbl       <- renderTable(post_summary_tbl(), digits = 3, striped = TRUE, hover = TRUE)
  output$pred_tbl       <- renderTable(pred_summary_tbl(), digits = 2, striped = TRUE, hover = TRUE)

  # Active MCID provenance: "sidebar input", "derived from power calc",
  # or "default".  Used by the footnote and the downloaded-CSV header.
  mcid_source_text <- reactive({
    if (isTRUE(is.finite(input$mcid_rr) && input$mcid_rr > 0)) {
      "sidebar input MCID (RR)"
    } else if (isTruthy(input$pow_n_per_arm)) {
      sprintf("derived from power calc (n=%g/arm, alpha=%g, target power=%g)",
              input$pow_n_per_arm,
              ifelse(isTruthy(input$pow_alpha),  input$pow_alpha,  0.05),
              ifelse(isTruthy(input$pow_target), input$pow_target, 0.80))
    } else {
      "default 0.80"
    }
  })

  # Footnote explaining what MCID is and where the active value came from.
  mcid_footnote_text <- reactive({
    HTML(sprintf(
      paste0(
        "<small><em>Footnote.</em> MCID (Minimum Clinically Important ",
        "Difference) is the smallest effect considered worth acting on. ",
        "In this app the active MCID is <strong>%.3f</strong> (RR scale), ",
        "taken from %s. A common pragmatic choice when no formal MCID ",
        "exists is to use the effect size the original trial was powered ",
        "to detect &mdash; its design target.</small>"
      ),
      mcid(), mcid_source_text()
    ))
  })
  output$mcid_footnote_small <- renderUI(mcid_footnote_text())
  output$mcid_footnote       <- renderUI(mcid_footnote_text())

  # ====================================================================
  # ==================  RISK DIFFERENCE (RD) PATH  ====================
  # ====================================================================

  # ----- RD prior parameters (proportion scale) -----
  prior_pars_rd <- function(prefix) {
    if (identical(input$rd_prior_mode, "rdci")) {
      cen <- input[[paste0(prefix, "_cen")]]
      lo  <- input[[paste0(prefix, "_lo")]]
      hi  <- input[[paste0(prefix, "_hi")]]
      validate(need(all(is.finite(c(cen, lo, hi))) && hi > lo,
                    "RD & CI must be finite and ordered"))
      mu0 <- cen / 100; tau <- sigma_from_ci_rd(lo / 100, hi / 100)
    } else {
      mu0 <- input[[paste0(prefix, "_m")]]
      tau <- input[[paste0(prefix, "_sd")]]
      validate(need(is.finite(mu0) && is.finite(tau) && tau > 0,
                    "Prior SD must be > 0"))
      mu0 <- mu0 / 100; tau <- tau / 100   # pp -> proportion
    }
    c(mu0 = mu0, tau = tau)
  }

  all_priors_rd <- reactive({
    o <- obs()
    se <- if (is.finite(o$se_rd) && o$se_rd > 0) o$se_rd else 0.01
    list(
      Weak              = prior_pars_rd("rd_weak"),
      Enthusiastic      = prior_pars_rd("rd_enth"),
      Skeptical         = prior_pars_rd("rd_skep"),
      `Empirical Bayes` = eb_marginal_RD_prior(se)
    )
  })

  # ----- RD MCID (proportion scale) -----
  mcid_rd_val <- reactive({
    v <- input$mcid_rd
    if (is.na(v) || !is.finite(v)) return(-0.02)
    v / 100
  })

  # ----- RD posterior draws -----
  posterior_df_rd <- reactive({
    req(identical(input$effect_measure, "RD"))
    o <- obs()
    like <- input$likelihood

    if (identical(like, "normal")) {
      validate(need(is.finite(o$rd) && is.finite(o$se_rd) && o$se_rd > 0,
                    "Normal path needs finite RD and SE(RD)."))
      pars <- list(
        Weak         = prior_pars_rd("rd_weak"),
        Enthusiastic = prior_pars_rd("rd_enth"),
        Skeptical    = prior_pars_rd("rd_skep")
      )
      out <- lapply(names(pars), function(nm) {
        pr <- pars[[nm]]
        post <- post_from_prior_normal(o$rd, o$se_rd, pr["mu0"], pr["tau"])
        rd <- rnorm(40000, post$mu, post$sd)
        tibble(prior = nm, mu = rd, RD = rd)
      })
      eb_rd <- eb_posterior_RD_draws(o$rd, o$se_rd, ndraws = 40000)
      out$eb <- tibble(prior = "Empirical Bayes", mu = eb_rd, RD = eb_rd)
      bind_rows(out)

    } else if (identical(like, "binom")) {
      p2 <- if (is.finite(o$p2)) o$p2 else
        (input$e_c + input$e_t) / (input$n_c + input$n_t)
      p2 <- min(max(p2, 1e-6), 1 - 1e-6)
      se_eff <- if (is.finite(o$se_rd) && o$se_rd > 0) o$se_rd else 0.01
      pars <- list(
        Weak         = prior_pars_rd("rd_weak"),
        Enthusiastic = prior_pars_rd("rd_enth"),
        Skeptical    = prior_pars_rd("rd_skep")
      )
      pars$`Empirical Bayes` <- eb_marginal_RD_prior(se_eff)

      draw_one_rd <- function(pr, lab) {
        map_pars <- map_RD_prior_to_alpha(pr["mu0"], pr["tau"], p2)
        dr <- laplace_binomial(input$e_t, input$n_t, input$e_c, input$n_c,
                               m_alpha = map_pars["m_alpha"],
                               s_alpha = map_pars["s_alpha"],
                               m_beta  = qlogis(p2), s_beta = 2,
                               ndraws = 40000)
        rd_draws <- dr$p1 - dr$p2
        out <- tibble(prior = lab, mu = rd_draws, RD = rd_draws)
        if (identical(lab, "Empirical Bayes")) {
          post_mean <- mean(out$RD); post_sd <- stats::sd(out$RD)
          if (is.finite(post_mean) && is.finite(post_sd) && post_sd > 0) {
            fit <- eb_posterior_snr(post_mean / post_sd)
            snr <- fit$sampler(nrow(out))
            out$RD <- snr * post_sd
            out$mu <- out$RD
          }
        }
        out
      }
      bind_rows(lapply(names(pars), function(nm) draw_one_rd(pars[[nm]], nm)))

    } else if (identical(like, "mcmc")) {
      validate(need(mcmc_available(), "MCMC backend not available."))
      validate(need(input$run_mcmc > 0, "Press 'Run MCMC' to generate draws."))
      isolate({
        if (is.null(.brms_base)) init_brms_base()
        validate(need(!is.null(.brms_base), "Failed to compile brms base model."))
        p2 <- if (is.finite(o$p2)) o$p2 else
          (input$e_c + input$e_t) / (input$n_c + input$n_t)
        p2 <- min(max(p2, 1e-6), 1 - 1e-6)
        se_eff <- if (is.finite(o$se_rd) && o$se_rd > 0) o$se_rd else 0.01

        priors <- list(
          Weak         = prior_pars_rd("rd_weak"),
          Enthusiastic = prior_pars_rd("rd_enth"),
          Skeptical    = prior_pars_rd("rd_skep")
        )
        priors$`Empirical Bayes` <- eb_marginal_RD_prior(se_eff)

        nm_list <- names(priors)
        shiny::withProgress(
          message = "Running MCMC (RD)", value = 0,
          {
            res <- list()
            for (i in seq_along(nm_list)) {
              nm <- nm_list[i]
              pr <- priors[[nm]]
              map_pars <- map_RD_prior_to_alpha(pr["mu0"], pr["tau"], p2)
              fit_res <- mcmc_fit(
                e_t = input$e_t, n_t = input$n_t,
                e_c = input$e_c, n_c = input$n_c,
                prior_name = nm,
                m_alpha = map_pars["m_alpha"],
                s_alpha = map_pars["s_alpha"],
                chains  = input$mcmc_chains,
                iter    = input$mcmc_iter
              )
              res[[nm]] <- tibble(prior = nm, mu = fit_res$RD, RD = fit_res$RD)
              shiny::incProgress(1 / length(nm_list),
                                 detail = sprintf("%s (%d/%d)", nm, i, length(nm_list)))
            }
            bind_rows(res)
          }
        )
      })
    }
  })

  # ----- RD predictive draws -----
  predictive_df_rd <- reactive({
    df <- posterior_df_rd()
    o  <- obs()
    se_new <- input$se_future_rd
    if (is.na(se_new) || se_new <= 0) se_new <- if (is.finite(o$se_rd)) o$se_rd else 0.01
    rd_pred <- rnorm(nrow(df), mean = df$RD, sd = se_new)
    tibble(prior = df$prior, RD_pred = rd_pred)
  })

  # ----- RD plots -----
  output$priors_plot_rd <- renderPlot({
    req(identical(input$effect_measure, "RD"))
    priors <- all_priors_rd()
    taus <- sapply(priors, function(p) unname(p["tau"]))
    rd_max <- 3 * max(taus, 0.01, na.rm = TRUE)
    rd_grid <- seq(-rd_max, rd_max, length.out = 600)
    pr_df <- prior_density_rd(priors, rd_grid)
    m <- mcid_rd_val()
    ggplot(pr_df, aes(RD, density, color = prior)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = m, linetype = "dotted", colour = "steelblue") +
      labs(title = "Prior distributions (RD scale)",
           x = "Risk Difference", y = "Density") +
      theme_minimal()
  })

  output$post_plot_rd <- renderPlot({
    req(identical(input$effect_measure, "RD"))
    df <- posterior_df_rd()
    lims <- rd_plot_limits(df$RD)
    m <- mcid_rd_val()
    ggplot(df, aes(x = RD, fill = prior, color = prior)) +
      annotate("rect", xmin = -Inf, xmax = m, ymin = -Inf, ymax = Inf,
               alpha = 0.08, fill = "steelblue") +
      stat_slabinterval(aes(thickness = after_stat(pdf)), .width = 0.95,
                        alpha = 0.5, position = "identity") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = m, linetype = "dotted", colour = "steelblue") +
      labs(title = "Posterior RD under each prior",
           x = "Risk Difference", y = "Density") +
      coord_cartesian(xlim = lims) +
      theme_minimal()
  })

  output$pred_plot_rd <- renderPlot({
    req(identical(input$effect_measure, "RD"))
    df <- predictive_df_rd()
    lims <- rd_plot_limits(df$RD_pred)
    m <- mcid_rd_val()
    ggplot(df, aes(x = RD_pred, fill = prior, color = prior)) +
      stat_slabinterval(aes(thickness = after_stat(pdf)), .width = 0.95,
                        alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = m, linetype = "dotted", colour = "steelblue") +
      labs(title = "Posterior predictive RD for a future study",
           x = "Risk Difference", y = "Density") +
      coord_cartesian(xlim = lims) +
      theme_minimal()
  })

  output$plv_plot_rd <- renderPlot({
    req(identical(input$effect_measure, "RD"))
    df <- posterior_df_rd()
    priors <- all_priors_rd()
    o <- obs()
    plv_plot_rd(df, priors, obs_rd = o$rd, obs_se_rd = o$se_rd,
               mcid_rd = mcid_rd_val())
  })

  output$forest_plot_rd <- renderPlot({
    req(identical(input$effect_measure, "RD"))
    df <- posterior_df_rd()
    forest_plot_rd(df, mcid_rd = mcid_rd_val())
  })

  # ----- RD summary tables -----
  post_summary_tbl_rd <- reactive(
    summarise_posterior_rd(posterior_df_rd(), mcid_rd = mcid_rd_val())
  )
  pred_summary_tbl_rd <- reactive(summarise_predictive_rd(predictive_df_rd()))

  output$post_tbl_small_rd <- renderTable(post_summary_tbl_rd(),
                                          digits = 3, striped = TRUE, hover = TRUE)
  output$pred_tbl_small_rd <- renderTable(pred_summary_tbl_rd(),
                                          digits = 2, striped = TRUE, hover = TRUE)
  output$post_tbl_rd       <- renderTable(post_summary_tbl_rd(),
                                          digits = 3, striped = TRUE, hover = TRUE)
  output$pred_tbl_rd       <- renderTable(pred_summary_tbl_rd(),
                                          digits = 2, striped = TRUE, hover = TRUE)

  # ----- RD MCID footnote -----
  mcid_footnote_text_rd <- reactive({
    m <- mcid_rd_val()
    HTML(sprintf(
      paste0(
        "<small><em>Footnote.</em> MCID (Minimum Clinically Important ",
        "Difference) is the smallest effect considered worth acting on. ",
        "The active MCID is <strong>%.2f pp</strong> (RD scale, ",
        "= %.4f on proportion scale).</small>"
      ),
      m * 100, m
    ))
  })
  output$mcid_footnote_small_rd <- renderUI(mcid_footnote_text_rd())
  output$mcid_footnote_rd       <- renderUI(mcid_footnote_text_rd())

  # ----- download -----
  output$download_csv <- downloadHandler(
    filename = function() sprintf("bayes_rct_summary_%s.csv", Sys.Date()),
    content = function(file) {
      is_rd <- identical(input$effect_measure, "RD")
      if (is_rd) {
        post <- post_summary_tbl_rd()
        pred <- pred_summary_tbl_rd()
        header <- c(
          sprintf("# Bayesian RCT summary  (exported %s)", Sys.time()),
          sprintf("# Effect measure    : Risk Difference"),
          sprintf("# MCID (RD, pp)     : %.2f", mcid_rd_val() * 100),
          sprintf("# Likelihood        : %s", input$likelihood),
          sprintf("# Counts (t / c)    : %d / %d events out of %d / %d",
                  input$e_t, input$e_c, input$n_t, input$n_c),
          "#"
        )
      } else {
        post <- post_summary_tbl()
        pred <- pred_summary_tbl()
        header <- c(
          sprintf("# Bayesian RCT summary  (exported %s)", Sys.time()),
          sprintf("# Effect measure    : Relative Risk"),
          sprintf("# MCID (RR)         : %.4f", mcid()),
          sprintf("# MCID source       : %s", mcid_source_text()),
          sprintf("# Likelihood        : %s", input$likelihood),
          sprintf("# Counts (t / c)    : %d / %d events out of %d / %d",
                  input$e_t, input$e_c, input$n_t, input$n_c),
          "#"
        )
      }
      writeLines(header, con = file)
      suppressWarnings(
        utils::write.table(
          dplyr::left_join(post, pred, by = "prior"),
          file = file, sep = ",", row.names = FALSE,
          append = TRUE, col.names = TRUE, quote = TRUE
        )
      )
    }
  )
}

shinyApp(ui, server)
