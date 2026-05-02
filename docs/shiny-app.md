# Shiny app — parameter reference
Plain-text reference for every input and output in `app.r`. Kept in plain text
(no rendered math) so it's readable in any terminal, editor, or GitHub preview.

Launch with:

    Rscript -e 'shiny::runApp("app.r", launch.browser = TRUE)'

---

## Inputs (sidebar, top to bottom)

### Likelihood
Radio group. Selects how the likelihood is combined with each prior to form
the posterior.

- **Normal (analytic)** — default. Uses the observed log(RR) and its SE as a
  Normal likelihood; combines with each Normal prior in closed form. Fastest.
  Requires at least one event in each arm (otherwise the observed SE is not
  defined).
- **Binomial (Laplace, no normality)** — fits a Bayesian logistic regression
  on the 2x2 counts (intercept = control log-odds, slope = arm contrast) by
  Laplace approximation. Use this when events are sparse or the Normal
  approximation on log(RR) is suspect.
- **MCMC (brms + cmdstanr)** — opt-in. A full Bayesian binomial-logit fit via
  `brms`. Requires `cmdstanr` (or `rstan`) to be installed. The Stan model is
  pre-compiled once at app startup; subsequent fits reuse the binary via
  `brms::update()`, and results are memoised. If neither backend is available
  this radio choice is disabled with a visible note and the other two paths
  continue to work.

When MCMC is selected, two extra numeric inputs appear:

- **MCMC chains** — number of chains to run (default 2).
- **MCMC iterations** — total iterations per chain (half are warmup; default
  2000). Bump for publication-quality runs.

There is also a **Run MCMC** button: MCMC is only triggered on demand (not on
every input change), so the UI stays responsive.

### Observed RCT counts
Four numeric boxes:

- **Treatment n1** — total in the treatment arm.
- **Treatment events e1** — outcome events in the treatment arm (0 ≤ e1 ≤ n1).
- **Control n2** — total in the control arm.
- **Control events e2** — outcome events in the control arm (0 ≤ e2 ≤ n2).

From these the app computes:

- p1 = e1 / n1, p2 = e2 / n2
- Observed RR = p1 / p2
- Observed log(RR)
- SE[log(RR)] = sqrt( (1 - p1) / (p1 * n1) + (1 - p2) / (p2 * n2) )
- 95 % CI on RR = exp( log(RR) +/- 1.96 * SE[log(RR)] )
- z = log(RR) / SE;  two-sided p = 2 * pnorm(-|z|)

A **Frequentist summary** panel in the sidebar (immediately under the
count inputs) shows the resulting `RR (95 % CI)` along with `z` and `p`.
For the default inputs:

    RR = 0.986 (95% CI 0.852 - 1.142)
    z = -0.18, p = 0.855

The main-panel "Observed (from counts)" block shows the same numbers plus
log(RR) and SE[log(RR)]. If any count is degenerate the panels show NA and
the Normal path will refuse to compute the posterior (use Laplace instead).

### Prior entry mode
Radio toggle: **RR & 95% CI** vs **Mean log(RR) & SD**. Changes how the three
user-configurable priors below are entered. Internally they are always stored
as Normal(mu0, tau^2) on the log(RR) scale.

- `RR & 95% CI` mode: enter a central RR and the lower/upper bounds of a 95%
  CI on the RR scale. The app back-calculates
  tau = (log(upper) - log(lower)) / (2 * 1.96) and mu0 = log(RR).
- `Mean log(RR) & SD` mode: enter mu0 and tau directly.

### Weak / Enthusiastic / Skeptical priors
Three editable priors with sensible defaults:

- **Weak** — RR 1.00, 95% CI 0.50 – 2.00 (log-RR mu0 = 0, tau ≈ 0.354). Broad;
  lets the data dominate.
- **Enthusiastic** — RR 0.80, 95% CI 0.70 – 0.92 (centred on a modest
  treatment benefit).
- **Skeptical** — RR 1.10, 95% CI 1.00 – 1.21 (centred on a slight harm).

You can change any of these freely. Use them as "what would a reasonable
person who believes X think of this trial?" probes.

### Empirical Bayes prior
There is **no input row** for this prior — it has no user-tunable parameters.
It is automatically included alongside the three user-configurable priors.

It comes from van Zwet, Więcek & Gelman (2025, *Stat Methods Med Res*
34(12):2302–2312, doi:10.1177/09622802251380628). The prior is a 4-component
Normal mixture on the z-scale learned from ~23 551 Cochrane RCTs:

    component SDs (z-scale): 1.17, 1.74, 2.38, 5.73
    component weights     : 0.32, 0.31, 0.30, 0.07

The app uses these to form a posterior on the signal-to-noise scale
(SNR = log(RR) / SE), then back-transforms to log(RR) and RR. Interpretation:
this is an evidence-calibrated, data-driven shrinkage prior — it pulls the
raw estimate toward 0 by an amount learned from the distribution of all
published trials, not an amount chosen by the analyst.

The EB prior shows up in every downstream tab and table as a fourth
distribution / row labelled "Empirical Bayes".

### Clinically meaningful benefit (MCID)
The **Minimum Clinically Important Difference (MCID)** is the smallest
effect a clinician would consider worth acting on. It is often not
formally defined for a given trial. A common pragmatic estimate is **the
effect size the original trial was powered to detect** (its design
target) — i.e., assume that the investigators chose the trial size to
detect an effect they themselves regarded as clinically important.

Four related inputs:

- **MCID (RR)** — the minimum clinically important RR threshold you
  want to track. Default 0.80. Posterior tables report
  `P(RR <= MCID)` against this value, and posterior plots add a
  vertical reference line at this value with light shading of the
  region [0, MCID].
- **Power calc: n per arm / alpha / target power** — *optional*. If
  you leave `MCID (RR)` blank or set it <= 0 and fill these in, the
  app derives an implicit MCID from a Normal-approximation power
  calculation (control-arm event rate is taken from the observed
  counts if available). Use this when you know the trial was powered
  at e.g. n = 3500 per arm, alpha = 0.05, power = 0.80 and you want
  the resulting minimum-detectable effect as the MCID.

A footnote under the posterior tables echoes the active MCID value
*and* its source (sidebar input vs. derived from power calc vs.
default). The "MCID in use" line at the top of the main panel echoes
the resolved numeric value.

### Predictive settings
One numeric input:

- **SE_new (log RR)** — the standard error that a hypothetical future
  replication trial would have for its own estimate of log(RR). Used only by
  the *Predictive (RR)* tab. See the separate section below.

### Download summary CSV
Button that exports the combined posterior + predictive summary tables.

---

## The SE_new parameter, in depth
SE_new is the most frequently misunderstood input. Treat this section as
canonical.

### 1. What it represents
SE_new is a property of a *hypothetical future replication trial*, not of
the trial you just analysed. The current trial + your prior gives a
**posterior** over the true log(RR); a future trial of specified size/design
would estimate that true value with Normal noise of SD = SE_new.

Chain of reasoning:

    current data + prior   ->   posterior over the true log(RR)
    future trial           ->   observes log(RR)_new with Normal noise SE_new
    combine                ->   posterior predictive for log(RR)_new

### 2. Where it enters the math (plain text)
Let mu = the true log(RR).

- Posterior from current data:
      mu | data, prior  ~  Posterior(mu)

- Future replication likelihood:
      logRR_hat_new | mu  ~  Normal(mu, SE_new^2)

- Posterior predictive (what the app actually reports):
      logRR_hat_new | data  ~  integral over mu of
                               Normal(mu, SE_new^2) * p(mu | data)  dmu

Draw-by-draw implementation (what the app does internally):

    mu_pred <- rnorm(n_draws, mean = posterior_mu_draw, sd = SE_new)
    RR_pred <- exp(mu_pred)

Key identity — the **variance decomposition**:

    Var(logRR_hat_new | data)  =  Var(mu | data)  +  SE_new^2
                                  [posterior        [future sampling
                                   uncertainty]      noise]

SE_new adds variance. It does not move the centre.

### 3. How to choose a value
You're answering a design question: "how big do I imagine the replication
being?" For a balanced two-arm binary-outcome trial with per-arm sample
size n and control-arm event rate p, the large-sample SE is:

    SE[log(RR)]  ≈  sqrt( 2 * (1 - p) / (p * n) )

Reference values:

    Per-arm n   p = 0.05   p = 0.10   p = 0.20   p = 0.30
          50     0.87       0.60       0.40       0.31
         100     0.62       0.42       0.28       0.22
         500     0.28       0.19       0.13       0.10
       1 000     0.19       0.13       0.09       0.07
       4 000     0.10       0.07       0.04       0.04

Leaving the field **blank or ≤ 0** tells the app to reuse the *current*
trial's observed SE — i.e., predictive = "another trial exactly this size".

### 4. Effect on the Bayesian outputs — numerical demo
DECAF-shaped posterior: true log(RR) ~ Normal(mean = -0.417, sd = 0.187).
I.e., posterior HR median = 0.659, 95% CrI 0.45 – 0.95.

    Scenario                            SE_new    Median   95% interval   P(RR<1)  P(RR<=0.59)
    Current-sized replication           0.192      0.659   0.39 – 1.11      0.94       0.34
    n = 500 per arm,  p = 0.10          0.190      0.659   0.39 – 1.11      0.94       0.34
    n = 1 000 per arm, p = 0.10         0.130      0.659   0.42 – 1.03      0.97       0.31
    n = 4 000 per arm, p = 0.10         0.070      0.659   0.45 – 0.98      0.98       0.29
    n = 100 per arm,  p = 0.10          0.420      0.659   0.27 – 1.62      0.82       0.41
    "Infinite" replication (≈ 0)        ~ 0        0.659   0.46 – 0.95      0.99       0.28

### 5. How to read the table
- **Predictive median is invariant to SE_new** — anchored to the posterior
  centre (0.659). SE_new controls *spread*, not *location*.
- **95% predictive interval widens monotonically with SE_new.** It widens
  faster than the posterior because variances add: total Var = 0.187^2 +
  SE_new^2.
- **P(RR < 1) is NOT monotone in SE_new.** As SE_new -> 0, P(RR < 1)
  converges to the posterior probability (≈ 0.99 — "probability the *truth*
  is below 1"). As SE_new grows, mass is smeared symmetrically around the
  posterior median in log space; because the median is already left of 1,
  the extra mass crosses 1 more than it moves deeper into benefit, so
  P(RR < 1) falls (0.99 -> 0.94 -> 0.82).
- **P(RR <= MCID) usually rises with larger SE_new** in this example. Bigger
  noise spreads mass into both tails; with median (0.659) sitting just above
  the MCID (0.59), extra lower-tail mass crosses the threshold. A noisier
  future trial *looks* more likely to produce a clinically meaningful
  estimate, but this is sampling-noise assistance, not stronger evidence.

### 6. Two limiting cases worth memorising
- **SE_new -> 0 (imaginary infinite future trial):** predictive = posterior
  of the true effect. The Predictive tab's P(...) columns then coincide with
  the Posterior tab's numbers — "probability the truth is below X". Upper
  bound on how certain any empirical work could make you.
- **SE_new = observed SE (default when left blank):** predictive = "what
  would a same-size replication plausibly report?" The Bayesian analogue of
  Gelman–Carlin replication framing.

### 7. Relationship to Type-M / replication probability
The pipeline's type_m_s() function asks the frequentist question "if the
trial were run again at the same SE, how often would we see significance in
the correct direction?" The Predictive tab is the Bayesian analogue — it
gives you the full distribution of possible future estimates under your
posterior combined with a specified future trial design.

The Bayesian *replication probability* under a given SE_new is:

    P( |logRR_hat_new| > 1.96 * SE_new   AND   same sign as logRR_obs  |  data )

You can read it off the predictive draws by applying that threshold (not
currently surfaced as a stand-alone number in the app; a natural future
addition).

### 8. Common pitfalls
1. Confusing SE_new with the current trial's SE. SE_new describes a *future*
   trial. Smaller than current = asking about a bigger / more efficient
   design; larger = a smaller / less efficient one.
2. Treating SE_new = 0 as "no noise in the posterior". It only removes
   *future sampling* noise; posterior uncertainty itself still has full
   width = 0.187 in the demo.
3. Reading `P(RR < 1)` from the Predictive tab as the Bayesian probability
   of benefit. That quantity belongs on the Posteriors tab. Predictive
   `P(RR < 1)` is "probability a future point estimate lands below 1",
   which is systematically closer to 0.5 than the posterior probability.
4. Forgetting to match p. The SE formula assumes the future trial's
   control-arm event rate is similar. If the intervention has driven event
   rates down in contemporary practice, use a smaller p.

### 9. Practical recipe
1. Decide what decision the predictive is supporting.
2. Settle on a plausible replication design: n_per_arm and p.
3. Compute SE_new = sqrt( 2 * (1 - p) / (p * n_per_arm) ).
4. Enter it in the app; read off `Predictive 50%`, the 2.5% / 97.5%
   quantiles, and the `P(...)` columns.
5. Sanity-check with the two limiting cases: SE_new = 0 (posterior) and
   SE_new = observed SE (same-size replication).

### Bottom line for SE_new
SE_new is the dial that separates "what the *true* effect might be"
(predictive width = posterior width; set SE_new = 0) from "what a *future
trial's point estimate* might look like" (predictive width > posterior
width, with extra noise SE_new^2 added on top). Use the Posteriors tab for
scientific claims about the effect; use the Predictive tab with an honest
SE_new for design and replication questions.

---

## Outputs (main panel, top to bottom)

### Observed panel
Shows the observed RR, log(RR), and SE[log(RR)] as computed from the 2x2
counts. NA if any count is degenerate.

### MCID in use
Echoes the resolved MCID (either the value you entered or the one derived
from the power-calc inputs).

### Tabs
- **Priors (RR scale)** — overlay of the four prior densities (Weak,
  Enthusiastic, Skeptical, Empirical Bayes) on the RR scale. Dashed line at
  RR = 1; dotted line at the MCID.
- **Posteriors (RR)** — stat_slabinterval plot of posterior RR under each
  prior, with 95% CrI. Shaded region [0, MCID]. Below the plot: a small
  numeric table with `Mean [95% CrI]`, `P(RR <= MCID)`, `P(RR < 1)`,
  `P(RR < 0.80)`, `P(RR < 0.90)`, `P(RR < 1.15)` for each prior.
- **Predictive (RR)** — posterior-predictive distribution of RR in a future
  trial of SE_new. Below the plot: 2.5 / 50 / 97.5% quantiles for each
  prior.
- **Prior vs Likelihood vs Posterior** — faceted by prior, showing the
  three densities on a shared RR axis. Useful for "see where the prior
  pulled the posterior relative to the data".
- **Forest** — compact dot-whisker of posterior median + 95% CrI across the
  four priors. Vertical dashed line at 1, dotted line at MCID.
- **Tables** — the same numeric tables as under the plot tabs, gathered in
  one place.

### Download summary CSV
Combines the posterior and predictive numeric tables (one row per prior)
and downloads as CSV.

---

## Troubleshooting
- **"Normal path needs nonzero events in both groups"** — pick Binomial
  (Laplace) or MCMC instead, or adjust your counts so each arm has at
  least one event.
- **"MCID in use" shows 0.80 when you wanted something else** — MCID input
  is empty / zero / negative; fill in `MCID (RR)` directly, or fill all
  three power-calc inputs.
- **MCMC radio is greyed out** — `cmdstanr` (or `rstan`) is not installed
  in this R session. Install `cmdstanr` and run
  `cmdstanr::install_cmdstan()` once, then restart the app.
- **"RR & CI must be positive and ordered"** — one of your prior RR / CI
  boxes has a non-positive value or lower > upper. Fix the row.
- **MCMC "Run MCMC" button hangs** — first run compiles the Stan binary
  (10–60 s on a modern laptop). Subsequent runs reuse the compiled model
  and only re-sample; results are also memoised per input combination.
