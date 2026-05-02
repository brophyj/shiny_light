# R/likelihoods.R
#
# Likelihood / estimation helpers: observed summary from 2x2 counts,
# and a Laplace (quadratic) approximation to the joint posterior of the
# arm-effect on the logit scale under independent Normal priors.

#' Observed RR / logRR / SE from 2x2 counts
#'
#' @param e1,n1 Treatment events and total.
#' @param e2,n2 Control events and total.
#' @return list with RR (rr), log(RR) (y), SE[log(RR)] (sigma) and the
#'   arm event rates (p1, p2).  Returns NAs if counts are degenerate.
obs_from_counts <- function(e1, n1, e2, n2) {
  ok <- isTRUE(n1 > 0 & n2 > 0 & e1 >= 0 & e2 >= 0 & e1 <= n1 & e2 <= n2)
  if (!ok) return(list(rr = NA, y = NA, sigma = NA, p1 = NA, p2 = NA))
  p1 <- e1 / n1; p2 <- e2 / n2
  rr <- ifelse(p2 > 0, p1 / p2, NA_real_)
  y  <- ifelse(is.finite(rr) && rr > 0, log(rr), NA_real_)
  sigma <- if (p1 > 0 && p1 < 1 && p2 > 0 && p2 < 1)
    sqrt((1 - p1) / (p1 * n1) + (1 - p2) / (p2 * n2)) else NA_real_
  rd    <- p1 - p2
  se_rd <- if (p1 >= 0 && p1 <= 1 && p2 >= 0 && p2 <= 1 && n1 > 0 && n2 > 0)
    sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2) else NA_real_
  list(rr = rr, y = y, sigma = sigma, p1 = p1, p2 = p2,
       rd = rd, se_rd = se_rd)
}

#' Laplace posterior approximation for a 2x2 binomial trial
#'
#' Fits logit(p_t) = beta + alpha, logit(p_c) = beta with independent Normal
#' priors on (alpha, beta) and returns draws from a Normal approximation at
#' the posterior mode (Laplace approximation, Gelman et al. BDA3 ch. 4).
#'
#' @param e1,n1 Treatment events / total.
#' @param e2,n2 Control events / total.
#' @param m_alpha,s_alpha Normal prior on alpha (arm contrast).
#' @param m_beta,s_beta   Normal prior on beta (control logit rate).
#' @param ndraws Number of approximate posterior draws.
#' @return data.frame with columns mu (log-RR draws), RR, OR, p1, p2.
laplace_binomial <- function(e1, n1, e2, n2,
                             m_alpha = 0, s_alpha = 2,
                             m_beta  = 0, s_beta = 2,
                             ndraws = 40000) {
  stopifnot(e1 <= n1, e2 <= n2, e1 >= 0, e2 >= 0, n1 > 0, n2 > 0)
  logpost <- function(theta) {
    a <- theta[1]; b <- theta[2]
    p2 <- invlogit(b); p1 <- invlogit(b + a)
    if (!is.finite(p1) || !is.finite(p2) || p1 <= 0 || p1 >= 1 ||
        p2 <= 0 || p2 >= 1) return(-Inf)
    ll <- dbinom(e1, n1, p1, log = TRUE) + dbinom(e2, n2, p2, log = TRUE)
    lp <- dnorm(a, m_alpha, s_alpha, log = TRUE) +
          dnorm(b, m_beta,  s_beta,  log = TRUE)
    ll + lp
  }
  fit <- optim(c(0, m_beta), fn = function(th) -logpost(th),
               method = "BFGS", hessian = TRUE)
  if (fit$convergence != 0) warning("Laplace: optim did not fully converge")
  H <- fit$hessian
  Sigma <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(Sigma)) stop("Laplace: Hessian not invertible; try different priors")
  L <- chol(Sigma)
  Z <- matrix(rnorm(2 * ndraws), ncol = 2)
  TH <- sweep(Z %*% t(L), 2, fit$par, `+`)
  a <- TH[, 1]; b <- TH[, 2]
  p2 <- invlogit(b); p1 <- invlogit(b + a)
  RR <- p1 / p2
  mu <- log(RR)
  data.frame(mu = mu, RR = RR, OR = exp(a), p1 = p1, p2 = p2,
             check.names = FALSE)
}
