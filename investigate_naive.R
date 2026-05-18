# ══════════════════════════════════════════════════════════════════════════════
# GOAL
# Show that IPW recovers the true causal treatment effect while the naive model
# is clearly biased, in a longitudinal (non-survival) MSM setting.
#
# PROBLEM
# Getting the naive model to be substantially biased proved harder than expected.
# The root cause: in survivl, the copula only adds RESIDUAL L-Y dependence after
# the structural equation for Y. Since sum(A) and mean(L) are highly correlated
# (r = 0.87 with the feedback loop), most of the L signal is already captured by
# A in the structural model, leaving little residual for the copula to act on.
# Swapping copula sign (k_tau negative) and varying strength did not produce
# large consistent naive bias in the Monte Carlo.
#
# ATTEMPTED APPROACHES
# 1. Stronger L->A effect (beta c(30,-3,0.1)) + A->L feedback (A_l1 in L formula)
#    Result: cor(sum_A, mean_L) jumped to 0.87, but naive bias stayed small (~0.05)
#    because the copula residual was tiny once A explained most of Y's variance.
#
# 2. Negative copula (k_tau = -0.7) to flip L-Y direction
#    Result: cor(mean_L, Y) improved to -0.43, but naive bias was still ~0.018
#    and in the wrong direction — the two pathways (causal A->Y and copula L->Y)
#    partially cancel in complex ways due to the feedback loop.
#
# 3. Direct structural L effect on Y: added sum(L) to Y's formula
#    Result: sum(L) ~ 65 overwhelmed the logit (Y all zeros). Fixing with mean(L)
#    or small coefficient: the feedback loop (A->L->A) meant IPW was also biased
#    because the marginal causal effect (what IPW targets) is no longer -0.5 —
#    L's distribution changes under intervention when A->L feedback exists.
#    Removing the feedback loop: cor(sum_A, mean_L) collapsed to -0.16, so both
#    approaches gave similar (small) bias.
#
# POTENTIAL SOLUTIONS TO TRY
# A. Use a LINEAR (Gaussian) outcome — avoids logistic non-collapsibility, so
#    the marginal causal effect = conditional = -0.5 exactly, and OVB formula
#    is exact. Should make IPW vs naive comparison clean.
#
# B. Keep structural L effect on Y but engineer the DGP so that:
#    - L -> A is strong (low hemoglobin -> treatment)
#    - L -> Y is strong and direct (in the structural formula)
#    - A -> L feedback is ABSENT (so IPW marginal effect = structural effect)
#    - A_0 = A_1 = A_2 constraint is removed so there is more treatment variation
#      at all time points and stronger A-L confounding at each step
#
# C. Remove the A_0=A_1=A_2 constraint entirely — let all time points vary freely.
#    Currently only A_3 and A_4 are time-varying, limiting how much total sum(A)
#    co-varies with sum(L).
#
# D. Try a different family for L (e.g. binary L) so the logit-scale L->A
#    relationship is cleaner and OVB is easier to reason about.
# ══════════════════════════════════════════════════════════════════════════════

library(survivl)
options(digits = 4)

set.seed(42)
n <- 1e4

# ── Setup ──────────────────────────────────────────────────────────────────────

qtls <- data.frame(matrix(runif(n * 4), ncol = 4))
colnames(qtls) <- c("B", "C", "L_0", "A_0")
qtls <- cbind(qtls, qtls[["A_0"]], qtls[["A_0"]])
colnames(qtls) <- c("B", "C", "L_0", "A_0", "A_1", "A_2")

B   <- qbinom(qtls[["B"]], 1, 0.1)
C   <- qunif(qtls[["C"]], 25, 35)
L_0 <- qnorm(qtls[["L_0"]], 11 - 0.05 * B - 0.02 * C, 0.5)
A_0 <- A_1 <- A_2 <- qbinom(qtls[["A_0"]], 1, 0.5)
dat0 <- cbind(B, C, L_0, A_0, A_1, A_2)

forms <- list(
  list(B ~ 1, C ~ 1),
  L ~ B + C,                 # L does NOT depend on past A — no feedback loop
  A ~ L_l0 + A_l1,           # strong L->A confounding
  Y ~ B + C + I(A_0 + A_1 + A_2 + A_3 + A_4) + I(L_0 + L_1 + L_2 + L_3 + L_4),
  L ~ 1
)
fams <- list(list(5, 4), 1, 5, 5, 2)

# structural parameters: intercept, B, C, sum(A), sum(L)
# L has a direct negative effect on Y (high hemoglobin -> less anemia)
# analyst fits only Y ~ B + C + sum(A), omitting sum(L) -> OVB
# without feedback, L distribution is fixed under intervention, so
# IPW correctly recovers -0.5 (the marginal ~ conditional causal effect)
true_A_pars <- c(-0.5, 0.1, 0.02, -0.5)

# mean sum(L) ~ 5 * (11 - 0.05*0.1 - 0.02*30) ~ 50.5
# coef -0.05 gives contribution ~ -2.5; intercept calibrated for P(Y=1) ~ 0.3
causal_parameters <- c(1.9, 0.1, 0.02, -0.5, -0.05)

pars <- list(
  B   = list(beta = 0.1),
  C   = list(beta = 1),
  L   = list(beta = c(11, -0.05, -0.02), phi = 0.5),
  A   = list(beta = c(30, -3, 0.1)),    # strong L->A
  Y   = list(beta = causal_parameters),
  cop = list(Y = list(L = list(k_tau = 0.1, par2 = 5)))  # negligible residual copula
)

surv_model <- survivl_model(formulas = forms, family = fams,
                            pars = pars, T = 5, dat = dat0, qtls = qtls)
dat <- rmsm(n, surv_model)

# ── Estimators ─────────────────────────────────────────────────────────────────

# IPW
glm_A0     <- glm(A_0 ~ L_0,         family = binomial(), data = dat)
glm_A3     <- glm(A_3 ~ L_3 + A_2,   family = binomial(), data = dat)
glm_A4     <- glm(A_4 ~ L_4 + A_3,   family = binomial(), data = dat)
glm_A0_num <- glm(A_0 ~ 1,           family = binomial(), data = dat)
glm_A3_num <- glm(A_3 ~ A_0,         family = binomial(), data = dat)
glm_A4_num <- glm(A_4 ~ A_0 + A_3,   family = binomial(), data = dat)

ps_A0     <- predict(glm_A0,     type = "response")
ps_A3     <- predict(glm_A3,     type = "response")
ps_A4     <- predict(glm_A4,     type = "response")
ps_A0_num <- predict(glm_A0_num, type = "response")
ps_A3_num <- predict(glm_A3_num, type = "response")
ps_A4_num <- predict(glm_A4_num, type = "response")

w0 <- ifelse(dat$A_0 == 1, ps_A0, 1 - ps_A0)
w3 <- ifelse(dat$A_3 == 1, ps_A3, 1 - ps_A3)
w4 <- ifelse(dat$A_4 == 1, ps_A4, 1 - ps_A4)
w0_num <- ifelse(dat$A_0 == 1, ps_A0_num, 1 - ps_A0_num)
w3_num <- ifelse(dat$A_3 == 1, ps_A3_num, 1 - ps_A3_num)
w4_num <- ifelse(dat$A_4 == 1, ps_A4_num, 1 - ps_A4_num)
dat$iptw_weights <- w0_num * w3_num * w4_num / (w0 * w3 * w4)

glm_Y_ipw     <- glm(Y ~ B + C + I(A_0 + A_1 + A_2 + A_3 + A_4),
                     data = dat, weights = iptw_weights, family = binomial())
glm_Y_outcome <- glm(Y ~ B + C + I(L_0 + L_1 + L_2 + L_3 + L_4) +
                       I(A_0 + A_1 + A_2 + A_3 + A_4),
                     data = dat, family = binomial())
glm_Y_naive   <- glm(Y ~ B + C + I(A_0 + A_1 + A_2 + A_3 + A_4),
                     data = dat, family = binomial())

# ── Single-run results ─────────────────────────────────────────────────────────

sumY_ipw     <- summary(glm_Y_ipw)
sumY_outcome <- summary(glm_Y_outcome)
sumY_naive   <- summary(glm_Y_naive)
shared_nms   <- rownames(sumY_ipw$coefficients)

cat("True A parameters (what IPW and naive target):\n")
print(setNames(true_A_pars, shared_nms))

cat("\n── Raw coefficients ──────────────────────────────────────────────────\n")
coefs <- cbind(
  truth   = true_A_pars,
  ipw     = sumY_ipw$coefficients[, 1],
  outcome = sumY_outcome$coefficients[shared_nms, 1],
  naive   = sumY_naive$coefficients[, 1]
)
print(round(coefs, 4))

cat("\n── SEs from truth (< 2 is good) ──────────────────────────────────────\n")
ses <- cbind(
  ipw     = abs(sumY_ipw$coefficients[, 1]     - true_A_pars) / sumY_ipw$coefficients[, 2],
  outcome = abs(sumY_outcome$coefficients[shared_nms, 1] - true_A_pars) / sumY_outcome$coefficients[shared_nms, 2],
  naive   = abs(sumY_naive$coefficients[, 1]   - true_A_pars) / sumY_naive$coefficients[, 2]
)
print(round(ses, 3))

# ── Investigate confounding: does A predict L, and does L predict Y? ──────────

cat("\n── Confounding check ─────────────────────────────────────────────────\n")
cat("Correlation between sum(A) and mean(L):\n")
sum_A  <- with(dat, A_0 + A_1 + A_2 + A_3 + A_4)
mean_L <- with(dat, (L_0 + L_1 + L_2 + L_3 + L_4) / 5)
cat(" cor(sum_A, mean_L) =", round(cor(sum_A, mean_L), 3), "\n")
cat(" cor(mean_L, Y)     =", round(cor(mean_L, dat$Y), 3), "\n")
cat(" cor(sum_A, Y)      =", round(cor(sum_A, dat$Y), 3), "\n")

cat("\nTreatment prevalence by time point:\n")
trt_prev <- colMeans(dat[, paste0("A_", 0:4)])
print(round(trt_prev, 3))

cat("\nMarginal P(Y=1) by sum_A:\n")
breaks <- unique(quantile(sum_A, probs = 0:4/4))
if (length(breaks) > 2) {
  quartile_A <- cut(sum_A, breaks = breaks, include.lowest = TRUE)
  print(round(tapply(dat$Y, quartile_A, mean), 3))
} else {
  print(round(tapply(dat$Y, sum_A, mean), 3))
}

cat("\nIPW weight summary:\n")
print(summary(dat$iptw_weights))

# ── Monte Carlo: repeat over many seeds to check if naive is consistently good ─

cat("\n── Monte Carlo bias check (50 replications, n=1000) ──────────────────\n")

run_one <- function(seed) {
  set.seed(seed)
  n_mc <- 1e3

  qtls_mc <- data.frame(matrix(runif(n_mc * 4), ncol = 4))
  colnames(qtls_mc) <- c("B", "C", "L_0", "A_0")
  qtls_mc <- cbind(qtls_mc, qtls_mc[["A_0"]], qtls_mc[["A_0"]])
  colnames(qtls_mc) <- c("B", "C", "L_0", "A_0", "A_1", "A_2")

  B_mc   <- qbinom(qtls_mc[["B"]], 1, 0.1)
  C_mc   <- qunif(qtls_mc[["C"]], 25, 35)
  L0_mc  <- qnorm(qtls_mc[["L_0"]], 11 - 0.05 * B_mc - 0.02 * C_mc, 0.5)
  A0_mc  <- qbinom(qtls_mc[["A_0"]], 1, 0.5)
  dat0_mc <- cbind(B = B_mc, C = C_mc, L_0 = L0_mc,
                   A_0 = A0_mc, A_1 = A0_mc, A_2 = A0_mc)

  pars_mc <- pars  # picks up updated pars from outer scope
  sm <- suppressMessages(
    survivl_model(formulas = forms, family = fams,
                  pars = pars_mc, T = 5, dat = dat0_mc, qtls = qtls_mc)
  )
  d <- rmsm(n_mc, sm)

  # IPW weights
  a0d <- glm(A_0 ~ L_0,       family = binomial(), data = d)
  a3d <- glm(A_3 ~ L_3 + A_2, family = binomial(), data = d)
  a4d <- glm(A_4 ~ L_4 + A_3, family = binomial(), data = d)
  a0n <- glm(A_0 ~ 1,         family = binomial(), data = d)
  a3n <- glm(A_3 ~ A_0,       family = binomial(), data = d)
  a4n <- glm(A_4 ~ A_0 + A_3, family = binomial(), data = d)

  p0d <- predict(a0d, type = "response"); p3d <- predict(a3d, type = "response")
  p4d <- predict(a4d, type = "response"); p0n <- predict(a0n, type = "response")
  p3n <- predict(a3n, type = "response"); p4n <- predict(a4n, type = "response")

  w <- (ifelse(d$A_0==1, p0n, 1-p0n) * ifelse(d$A_3==1, p3n, 1-p3n) * ifelse(d$A_4==1, p4n, 1-p4n)) /
       (ifelse(d$A_0==1, p0d, 1-p0d) * ifelse(d$A_3==1, p3d, 1-p3d) * ifelse(d$A_4==1, p4d, 1-p4d))
  d$w <- w

  fit_ipw   <- glm(Y ~ B + C + I(A_0+A_1+A_2+A_3+A_4), data=d, weights=w,  family=binomial())
  fit_naive <- glm(Y ~ B + C + I(A_0+A_1+A_2+A_3+A_4), data=d, family=binomial())

  trt_nm <- grep("^I\\(", names(coef(fit_ipw)), value = TRUE)
  c(ipw_trt   = unname(coef(fit_ipw)[trt_nm]),
    naive_trt = unname(coef(fit_naive)[trt_nm]))
}

mc_res <- t(sapply(1:50, run_one))
truth_trt <- true_A_pars[4]  # -0.5

cat(sprintf("Treatment coef (truth = %.2f)\n", truth_trt))
cat(sprintf("  IPW   — mean: %6.3f  bias: %6.3f  SD: %.3f\n",
            mean(mc_res[,"ipw_trt"]),
            mean(mc_res[,"ipw_trt"]) - truth_trt,
            sd(mc_res[,"ipw_trt"])))
cat(sprintf("  Naive — mean: %6.3f  bias: %6.3f  SD: %.3f\n",
            mean(mc_res[,"naive_trt"]),
            mean(mc_res[,"naive_trt"]) - truth_trt,
            sd(mc_res[,"naive_trt"])))
