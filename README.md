# Bayesian Re-analysis of Published Randomized Trials

This repository contains two interoperable tools for re-interpreting published
randomized controlled trials (RCTs) through a Bayesian / frequentist-calibration
lens:

1. **Shiny apps** for interactive single-trial analysis:
   - **`app_full.r`** — full-feature: Bayesian posteriors + frequentist
     calibration (S-values, Type M/S, Design-informed prior, Calibration tab).
   - **`app_bayes_alone.r`** — Bayesian-only (simpler; no calibration extras).
2. **`rct_pipeline.qmd`** (+ `run_pipeline.R`) — a **Quarto pipeline** that
   re-analyses a curated cohort of published cardiovascular RCTs systematically.

Both tools share the same analytical core: `R/eb_mixture.R`, `R/priors.R`,
`R/summaries.R`, `R/plots.R`, `R/likelihoods.R`. Pipeline-only helpers
(`type_m_s`, `s_value`, `fit_brms_study`, `build_priors`, pipeline plots) live
in `R/pipeline.R` and are **not** sourced by the Shiny apps.

## Why

Frequentist RCT reports — p-values, CIs, point estimates — can be misleading on
their own, especially for small trials. This project operationalises four ideas
for honest single-trial interpretation:

- **Bayesian posteriors under multiple priors** (Weak / Skeptical /
  Enthusiastic / Empirical-Bayes / Design-informed) so the effect of the prior
  is explicit, not hidden.
- **Empirical Bayes from ~23 551 Cochrane trials** (van Zwet, Więcek & Gelman,
  *Stat Methods Med Res* 2025; 34(12):2302–2312;
  [doi:10.1177/09622802251380628](https://doi.org/10.1177/09622802251380628))
  via a 4-component Normal mixture on the signal-to-noise (z) scale.
- **Type-M / Type-S retrodesign** (Gelman & Carlin 2014) on the published power
  calculation to quantify exaggeration and sign-error risk.
- **MCID-aware posterior readouts**: in addition to `P(HR < 1)` ("probability
  of any benefit"), report `P(HR ≤ MCID)` — the probability the true effect
  reaches the clinically-meaningful threshold the trial was powered for.

## Installation

### Prerequisites
- R ≥ 4.3
- [Quarto](https://quarto.org/) ≥ 1.4 (for the pipeline and the lecture decks
  `eb1.qmd`, `eb2.qmd`)
- A working C++ toolchain and either **cmdstanr** (with `cmdstan` installed)
  or **rstan** for the MCMC paths — **optional**; both Shiny app and pipeline
  gracefully fall back to analytic / Laplace approximations if Stan is missing.

### R packages
```r
install.packages(c(
  "shiny", "ggplot2", "dplyr", "tidyr", "tibble", "readr", "purrr",
  "ggdist", "bslib", "rlang", "memoise", "future", "patchwork",
  "rentrez", "xml2", "brms", "bayesplot", "posterior",
  "metafor", "meta", "bayesmeta", "bayestestR",
  "testthat", "kableExtra", "pdftools"
))
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/",
                                       getOption("repos")))
cmdstanr::install_cmdstan(cores = 2)   # one-time, ~5-15 min
```

## Running the Shiny apps

```sh
# Full-feature version (recommended):
Rscript -e 'shiny::runApp("app_full.r", launch.browser = TRUE)'

# Bayesian-only (no calibration metrics):
Rscript -e 'shiny::runApp("app_bayes_alone.r", launch.browser = TRUE)'
```

> Full parameter reference (plain-text, no rendered math): [`docs/shiny-app.md`](docs/shiny-app.md). Covers every sidebar input, the EB prior, an in-depth treatment of the `SE_new (log RR)` predictive parameter, every output tab, and common troubleshooting.

### App layout

| Area | Contents |
|------|----------|
| **Left sidebar** | All inputs: likelihood selector, RCT counts, frequentist summary, prior entry, MCID, Design-informed prior, `SE_new (log RR)`, download button |
| **Right main panel (tabs)** | *Observed / Summary*, *Priors (RR scale)*, *Posteriors (RR)*, *Predictive (RR)*, *Prior vs Likelihood vs Posterior*, *Forest*, *Tables*, *Calibration* (`app_full.r` only) |

### What the apps do (per single trial)
- Accept 2×2 counts (`e_t, n_t, e_c, n_c`) and three user-defined priors
  entered either as `RR ± 95% CI` or as `mean log(RR) ± SD`.
- Compute posterior draws under **four** priors (`app_full.r` adds a 5th):
  - `Weak`, `Enthusiastic`, `Skeptical` — user-specified.
  - `Empirical Bayes` — van Zwet 2025 mixture; no user parameters.
  - `Design-informed` (**`app_full.r` only**) — centred on the trial's
    pre-specified design HR.
- Three interchangeable **likelihoods** selectable at runtime:
  - **Normal (analytic)** — closed-form Normal–Normal conjugate update on
    log-RR (fastest).
  - **Binomial (Laplace)** — Laplace approximation to the joint posterior of
    `(log-OR, log-odds(control))`; the EB prior is first applied on log-RR and
    re-shrunk from the Laplace posterior via a pseudo-z.
  - **MCMC (brms + cmdstanr)** — full Bayesian binomial-logit fit, on demand
    via the **Run MCMC** button. The Stan model is pre-compiled at app startup
    and `brms::update()` is used for subsequent fits; results are memoised.
- Optional **MCID** input (or implied by an `n / α / power` triple); posterior
  tables include `P(RR ≤ MCID)` and `P(RR < 1)`, plots shade the
  `[0, MCID]` region.
- **`app_full.r` extras**: frequentist calibration block (S-values vs RR=1 and
  vs MCID; Gelman–Carlin Type-S/M retrodesign), and a dedicated Calibration
  tab with per-prior SNR metrics and EB mixture decomposition.

If `cmdstanr`/`rstan` is not available the MCMC radio is disabled and a note
is shown; the Normal and Laplace paths continue to work.

## Running the systematic pipeline

The pipeline re-analyses a curated cohort of cardiovascular RCTs (file
`data/rct_cohort.csv`) using the **same** EB prior and prior-comparison
framework as the Shiny app, but with proper MCMC.

```sh
# Render the full Quarto document (HTML, with per-study MCMC fits):
quarto render rct_pipeline.qmd

# Or run the non-Quarto driver to (re)build results/ only:
Rscript run_pipeline.R
```

Key parameters for the Quarto front-matter (`params:`):
- `n_mcmc_chains` (default 2), `n_mcmc_iter` (default 1 000) — bump for
  publication-quality runs.
- `n_eb_draws` (default 100 000) — Monte-Carlo size for the EB mixture posterior.

Output artifacts:
- `rct_pipeline.html` — rendered report (committed to the repo).
- `results/pipeline_summary.csv` — flat per-(study × prior) summary table.
- `results/forest.png` — dot-whisker of posterior HR medians + 95% CrIs.
- `results/heatmap_mcid.png` — study × prior heatmap of `P(HR ≤ MCID)`.
- `results/shrinkage_scatter.png` — observed |z| vs EB-shrunk |SNR mean|.

### The curated cohort (`data/rct_cohort.csv`)

Columns: `pmid, study, journal, year, endpoint, effect_type, hr, ci_lo, ci_hi,
n_t, e_t, n_c, e_c, alpha, power, mcid_rr, notes`.

v1 cohort (8 trials):
DECAF, PARADIGM-HF, EMPA-REG OUTCOME, DAPA-HF, ISCHEMIA, ORBITA,
EMPEROR-Reduced, SELECT. Two caveats:
- **EMPEROR-Reduced** and **SELECT** fill the JAMA Cardiology / Circulation
  slots but were both originally published in NEJM — noted in each row's
  `notes` column.
- **ORBITA** uses a continuous primary endpoint (exercise time) — stored with
  `effect_type = "MD"` and skipped by the HR-based Bayesian fit.

Editing the cohort is the supported way to extend the pipeline. Each new row
should come with published `(hr, ci_lo, ci_hi)`, arm sizes, and either a
stated target effect or enough power-calc info to derive `mcid_rr`.

### Automated PubMed helper

`R/pubmed.R` exposes `fetch_cv_rcts(journals, start, end)` (plus a thin
wrapper in `pubmed.R`) for pulling candidate studies from PubMed into
`data/pubmed_citations.csv` / `.ris`. The v1 pipeline uses a **hybrid**
strategy: PubMed retrieves metadata; the numeric fields needed for the
analyses (`hr`, `ci_lo`, `ci_hi`, arm counts, power) come from the
hand-curated `data/rct_cohort.csv`. A future iteration can close the loop
with a proper PDF/abstract parser.

## Methodology

### Observed-trial quantities
For a published HR / RR with 95% CI `(lo, hi)`:

- `log_effect = log(point)`
- `SE = (log(hi) − log(lo)) / (2 · 1.96)`
- `z = log_effect / SE`, `p = 2 · pnorm(−|z|)`
- `S-value = −log₂(p)` (bits of Shannon surprisal against H₀)

### Gelman–Carlin Type M / Type S (`type_m_s`)
Solves for the SE implied by a nominal `power` at `alpha = 0.05` (two-sided)
and a `target_effect` on the SNR scale, then simulates estimates under that
design:

- **Type S** = `P(sign(est) ≠ sign(true) | |est| > z_α/2 · SE)` — the
  conditional sign-error probability.
- **Type M** = median `|est / true|` among significant estimates — the
  exaggeration factor.

### Empirical-Bayes mixture (`eb_posterior_snr`)
Uses the 4-component Normal mixture on the z-scale (van Zwet 2025):

- `sd_z = c(1.17, 1.74, 2.38, 5.73)`
- `p_mix = c(0.32, 0.31, 0.30, 0.07)`

Deconvolves unit noise component-wise (`s_k² = sd_z_k² − 1`). Within
component `k`, conjugacy gives:

```
SNR | (z, k) ~ N(B_k · z, B_k),   B_k = s_k² / (s_k² + 1)
```

Posterior component weights are updated by the marginal z-likelihood
`N(0, s_k² + 1)`. The returned object exposes pooled `post_mean`, `post_sd`,
`B_bar` (overall shrinkage), a `sampler(n)` closure, and optional
pre-computed `snr` draws.

> **Parameter provenance.** The `sd_z` / `p_mix` values match those used in
> `eb2.qmd`, which come from the 2021 van Zwet paper. These are numerically
> very close to any values published in the 2025 follow-up; if you obtain the
> 2025 paper's supplementary materials and the numeric fit differs, update
> `EB_DEFAULTS` in `R/eb_mixture.R`.

### Bayesian fits
The Shiny app offers three interchangeable likelihoods (see above). The
pipeline uses **`fit_brms_study()`**, a `brms` normal-pseudo-likelihood model
on `(effect = log-HR, se = SE[log-HR])`:

```
effect | se(se) ~ 1,   Intercept ~ prior
```

Five priors built by **`build_priors()`**:

| Prior            | Mean (log-RR)    | SD (log-RR)       | Use case                        |
|------------------|------------------|-------------------|---------------------------------|
| `Weak`           | 0                | 10                | Data-dominated reference        |
| `Skeptical`      | 0                | 0.354             | ~95% mass: HR in [0.5, 2]       |
| `Enthusiastic`   | `log(mcid_rr)`   | 0.354             | Centred on the powered effect   |
| `EB`             | `post_mean · SE` | `post_sd · SE`    | Cochrane-calibrated shrinkage   |
| `DesignInformed` | `log(mcid_rr)`   | `SE[log-HR]`      | Trial-design prior (narrow)     |

If Stan is not installed, `run_pipeline.R` falls back to an analytic
Normal–Normal conjugate update under each prior — the `rct_pipeline.qmd`
Quarto document handles this per-chunk as well.

### Per-study summary (`summarise_posterior`)
Returns `median_RR`, `lo95`, `hi95`, `P_RR_lt_1`, `P_RR_le_mcid` for each
prior. `P_RR_le_mcid` is the probability of a clinically-meaningful effect
— the quantity most relevant to decision-making.

## Validation

The DECAF trial (`data/decaf.pdf`, JAMA 2025) serves as the end-to-end
validation case for the pipeline (see `eb2.qmd` for the lecture-style walk-through):

| Quantity                           | Value                      |
|------------------------------------|----------------------------|
| Observed HR (95% CI)               | 0.61 (0.42 – 0.89)         |
| z                                  | −2.58                      |
| p                                  | 0.0099                     |
| S-value                            | 6.66 bits                  |
| Type-S (at 80% power)              | 0.000                      |
| Type-M median                      | 1.09×                      |
| EB posterior mean SNR              | −1.80                      |
| EB overall shrinkage B̄             | 0.697                      |
| **EB posterior HR median**         | **0.659** (MCMC: 0.656)    |
| EB 95% CrI                         | 0.51 – 0.86                |
| EB `P(HR < 1)`                     | 0.999                      |
| EB `P(HR ≤ MCID = 0.59)`           | 0.206                      |

Two invariants are enforced as regression tests:

- `tests/testthat/test-decaf.R` — asserts the EB posterior HR median for
  DECAF lies in `[0.65, 0.90]`.
- `rct_pipeline.qmd` includes a `stopifnot()` with the same bounds inside the
  DECAF validation section.

### Run the test suite

Run from the **project root** so the path resolver can locate `R/`:

```sh
cd /path/to/projects
Rscript tests/testthat.R
# decaf: ..        (2)
# eb: ............ (12)
# priors: .......   (7)
# ══ DONE ═══════════════
```

## Repository layout

```
app_full.r                 # Shiny app — full feature (Bayesian + calibration)
app_bayes_alone.r          # Shiny app — Bayesian only (no calibration)
R/
  priors.R                 # sigma_from_ci, post_from_prior_normal,
                           # logit/invlogit, map_logRR_prior_to_alpha
  likelihoods.R            # obs_from_counts, laplace_binomial
  eb_mixture.R             # EB_DEFAULTS, eb_posterior_snr,
                           # eb_marginal_logRR_prior,
                           # eb_posterior_logRR_draws
  summaries.R              # summarise_posterior_rr, rr_plot_limits, resolve_mcid
  plots.R                  # prior/likelihood/posterior + forest helpers
                           # (Shiny side)
  pubmed.R                 # fetch_cv_rcts()  (rentrez + xml2)
  pipeline.R               # Pipeline-only: type_m_s, s_value,
                           # fit_brms_study, build_priors,
                           # summarise_posterior, plot_study_densities,
                           # plot_forest, plot_heatmap_mcid,
                           # plot_shrinkage_scatter
rct_pipeline.qmd           # Systematic Quarto pipeline
run_pipeline.R             # Stan-free / CLI driver for the pipeline
data/
  rct_cohort.csv           # v1 curated cohort (8 trials)
  decaf.pdf                # Source PDF for the DECAF validation case
results/                   # Rendered pipeline artifacts
eb1.qmd, eb2.qmd           # Lecture decks that motivate the methodology
pubmed.R                   # Thin wrapper over R/pubmed.R
tests/testthat/            # Unit + regression tests
```

## Citation

If you use this code or the EB prior, please cite:

> van Zwet E, Więcek W, Gelman A. **Meta-analysis with a single study.**
> *Statistical Methods in Medical Research.* 2025;34(12):2302–2312.
> [doi:10.1177/09622802251380628](https://doi.org/10.1177/09622802251380628)

Related background:
> Gelman A, Carlin J. **Beyond power calculations: assessing Type S (sign)
> and Type M (magnitude) errors.** *Perspect Psychol Sci.* 2014;9(6):641-651.


