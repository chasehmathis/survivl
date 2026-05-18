# ══════════════════════════════════════════════════════════════════════════════
# GOAL
# Show that IPW recovers the true causal treatment effect while the naive model
# is clearly biased, in a longitudinal (non-survival) MSM setting.
#
# KEY INSIGHT FROM PREVIOUS ATTEMPTS
# The original DGP put sum(L) directly in Y's structural formula, leaving only a
# small residual for the copula to act on. Once A explained most of Y's variance
# structurally, the copula could not generate meaningful L–Y residual dependence,
# and naive bias stayed near zero across every copula setting we tried.
#
# NEW STRATEGY
# Strip sum(L) out of Y's structural formula entirely. Then ALL of the L→Y
# dependence must flow through the copula. This is a cleaner test of how copula
# parametrization controls naive bias, and it lets us crank the copula HARD.
#
# Below we sweep several "stronger / more complicated" copula parametrizations:
#   (1) baseline    — weak Gaussian copula      (k_tau = 0.1)              fam 1
#   (2) strong gauss — strong Gaussian copula   (k_tau = 0.85)              fam 1
#   (3) heavy tail t — strong t-copula, df=2    (k_tau = 0.85, tail dep.)   fam 2
#   (4) conditional — copula parameter depends on past treatment (A_l1)    fam 2
#   (5) clayton     — lower-tail-dependent Clayton                          fam 3
#   (6) gumbel      — upper-tail-dependent Gumbel                           fam 4
#
# We keep the rest of the DGP fixed so the only thing changing across runs is
# the copula spec, and we report IPW bias vs naive bias for each one.
# ══════════════════════════════════════════════════════════════════════════════

library(survivl)
options(digits = 4)

# rho_to_beta: invert the latent expit link used by the `beta`-form copula spec.
# The conditional correlation is rho_i = 2 * expit(beta' x_i) - 1, so to target
# a marginal correlation `rho` we set the intercept to logit((rho + 1) / 2).
rho_to_beta <- function(rho) {
  x <- (rho + 1) / 2
  log(x / (1 - x))
}

# ── Fixed pieces of the DGP ───────────────────────────────────────────────────
# Y's structural formula deliberately OMITS sum(L). All L→Y dependence comes
# from the copula. No A→L feedback (L formula has no A term), so the marginal
# causal effect targeted by IPW equals the structural coefficient on sum(A).

forms <- list(
  list(B ~ 1, C ~ 1),
  L ~ B + C,                    # NO A_l1 here — no feedback loop
  A ~ L_l0 + A_l1,              # strong L→A confounding
  Y ~ B + C + I(A_0 + A_1 + A_2 + A_3 + A_4),
  L ~ 1                         # default: constant copula (overridden in (4))
)
fams <- list(list(5, 4), 1, 5, 5, 1)  # last entry = copula family; overridden per scenario

# Truth: IPW and naive both target this vector
true_pars <- c(`(Intercept)` = -2, B = 0.1, C = 0.02, `I(sum_A)` = -0.5)
causal_parameters <- unname(true_pars)

# Strong L→A: high L (good hemoglobin) → A = 0; low L → A = 1
A_pars  <- c(10, -1, 0.1)

base_pars <- list(
  B = list(beta = 0.1),
  C = list(beta = 1),
  L = list(beta = c(11, -0.05, -0.02), phi = 0.5),
  A = list(beta = A_pars),
  Y = list(beta = causal_parameters)
)

# ── Helper to build dat0 / qtls for a given n ────────────────────────────────

make_init <- function(n) {
  qtls <- data.frame(matrix(runif(n * 4), ncol = 4))
  colnames(qtls) <- c("B", "C", "L_0", "A_0")
  qtls <- cbind(qtls, qtls[["A_0"]], qtls[["A_0"]])
  colnames(qtls) <- c("B", "C", "L_0", "A_0", "A_1", "A_2")

  B   <- qbinom(qtls[["B"]], 1, 0.1)
  C   <- qunif(qtls[["C"]], 25, 35)
  L_0 <- qnorm(qtls[["L_0"]], 11 - 0.05 * B - 0.02 * C, 0.5)
  A_0 <- A_1 <- A_2 <- qbinom(qtls[["A_0"]], 1, 0.5)
  dat0 <- cbind(B, C, L_0, A_0, A_1, A_2)
  list(dat0 = dat0, qtls = qtls)
}

# ── IPW + naive estimation, returns the treatment coefficient from each ──────

fit_estimators <- function(d) {
  a0d <- glm(A_0 ~ L_0,       family = binomial(), data = d)
  a3d <- glm(A_3 ~ L_3 + A_2, family = binomial(), data = d)
  a4d <- glm(A_4 ~ L_4 + A_3, family = binomial(), data = d)
  a0n <- glm(A_0 ~ 1,         family = binomial(), data = d)
  a3n <- glm(A_3 ~ A_0,       family = binomial(), data = d)
  a4n <- glm(A_4 ~ A_0 + A_3, family = binomial(), data = d)

  p0d <- predict(a0d, type = "response"); p3d <- predict(a3d, type = "response")
  p4d <- predict(a4d, type = "response"); p0n <- predict(a0n, type = "response")
  p3n <- predict(a3n, type = "response"); p4n <- predict(a4n, type = "response")

  w <- (ifelse(d$A_0 == 1, p0n, 1 - p0n) *
        ifelse(d$A_3 == 1, p3n, 1 - p3n) *
        ifelse(d$A_4 == 1, p4n, 1 - p4n)) /
       (ifelse(d$A_0 == 1, p0d, 1 - p0d) *
        ifelse(d$A_3 == 1, p3d, 1 - p3d) *
        ifelse(d$A_4 == 1, p4d, 1 - p4d))
  d$w <- w

  fit_ipw   <- glm(Y ~ B + C + I(A_0 + A_1 + A_2 + A_3 + A_4),
                   data = d, weights = w, family = binomial())
  fit_naive <- glm(Y ~ B + C + I(A_0 + A_1 + A_2 + A_3 + A_4),
                   data = d, family = binomial())

  trt_nm <- grep("^I\\(", names(coef(fit_ipw)), value = TRUE)
  c(ipw_trt   = unname(coef(fit_ipw)[trt_nm]),
    naive_trt = unname(coef(fit_naive)[trt_nm]),
    weight_max = max(w),
    weight_mean = mean(w))
}

# ── Scenario definitions ──────────────────────────────────────────────────────
# Each scenario overrides (a) the copula formula entry forms[[5]], (b) the
# copula family entry fams[[5]], and (c) the cop entry in pars.

scenarios <- list(

  weak_gaussian = list(
    descr = "Weak Gaussian copula (k_tau = 0.1). Establishes baseline.",
    forms5 = L ~ 1,
    fam5   = 1,
    cop    = list(Y = list(L = list(k_tau = 0.1, par2 = 5)))
  ),

  strong_gaussian = list(
    descr = "Strong Gaussian copula (k_tau = 0.85). Heavy residual L–Y dep.",
    forms5 = L ~ 1,
    fam5   = 1,
    cop    = list(Y = list(L = list(k_tau = 0.85, par2 = 5)))
  ),

  heavy_tail_t = list(
    descr = "Strong t-copula (k_tau = 0.85, df = 2). Joint tail dependence.",
    forms5 = L ~ 1,
    fam5   = 2,
    cop    = list(Y = list(L = list(k_tau = 0.85, df = 2, par2 = 2)))
  ),

  conditional_t = list(
    # The copula parameter is a regression rho_i = 2*expit(beta0 + beta1*A_l1) - 1.
    # At A_l1 = 0 we target rho ≈ 0.30, at A_l1 = 1 we target rho ≈ 0.90.
    # So the L-Y dependence is much stronger after a treatment period — a
    # genuinely conditional copula that the naive model can never reproduce.
    descr = "Conditional t-copula: L–Y dep depends on past treatment A_l1.",
    forms5 = list(Y = list(L ~ A_l1)),
    fam5   = 2,
    cop    = list(Y = list(L = list(
      beta = c(rho_to_beta(0.30), rho_to_beta(0.90) - rho_to_beta(0.30)),
      df   = 3,
      par2 = 3
    )))
  ),

  clayton = list(
    descr = "Clayton copula (k_tau = 0.7). Lower-tail dependence only.",
    forms5 = L ~ 1,
    fam5   = 3,
    cop    = list(Y = list(L = list(k_tau = 0.7, par2 = 5)))
  )
  # NOTE: Gumbel (family 4) hangs under this DGP — looks like a numerical
  # pathology in the Gumbel h-function inversion when combined with the
  # rejection-sampling loop, so it's omitted from the sweep.
)

# ── Build a survivl_model for a given scenario ────────────────────────────────

build_model <- function(scn, dat0, qtls) {
  f <- forms
  f[[5]] <- scn$forms5
  fa <- fams
  fa[[5]] <- scn$fam5
  p <- base_pars
  p$cop <- scn$cop
  survivl_model(formulas = f, family = fa, pars = p,
                T = 5, dat = dat0, qtls = qtls)
}

# ── Single-shot diagnostic (large n) for the most interesting scenario ──────

set.seed(42)
n <- 1e4
init <- make_init(n)

cat("══════════════════════════════════════════════════════════════════════\n")
cat("Single-run diagnostic on the conditional-t scenario, n =", n, "\n")
cat("══════════════════════════════════════════════════════════════════════\n")

sm  <- build_model(scenarios$conditional_t, init$dat0, init$qtls)
dat <- rmsm(n, sm)

est <- fit_estimators(dat)
cat(sprintf("truth on sum_A coef = %.3f\n", causal_parameters[4]))
cat(sprintf("  IPW   estimate     = %.3f  (bias %+0.3f)\n",
            est["ipw_trt"],   est["ipw_trt"]   - causal_parameters[4]))
cat(sprintf("  Naive estimate     = %.3f  (bias %+0.3f)\n",
            est["naive_trt"], est["naive_trt"] - causal_parameters[4]))
cat(sprintf("  max IPW weight = %.2f, mean weight = %.3f\n",
            est["weight_max"], est["weight_mean"]))

sum_A  <- with(dat, A_0 + A_1 + A_2 + A_3 + A_4)
mean_L <- with(dat, (L_0 + L_1 + L_2 + L_3 + L_4) / 5)
cat(sprintf("\ncor(sum_A, mean_L) = %.3f\n",  cor(sum_A, mean_L)))
cat(sprintf("cor(mean_L, Y)     = %.3f\n",    cor(mean_L, dat$Y)))
cat(sprintf("cor(sum_A, Y)      = %.3f\n",    cor(sum_A, dat$Y)))

# ── Monte Carlo sweep across scenarios ────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("Monte Carlo sweep: 50 reps, n = 1000, across copula parametrizations\n")
cat("══════════════════════════════════════════════════════════════════════\n")

run_one <- function(seed, scn) {
  set.seed(seed)
  init_mc <- make_init(1000)
  sm <- suppressMessages(build_model(scn, init_mc$dat0, init_mc$qtls))
  d  <- rmsm(1000, sm)
  fit_estimators(d)
}

results <- lapply(names(scenarios), function(nm) {
  scn <- scenarios[[nm]]
  cat(sprintf("\n--- %s ---\n%s\n", nm, scn$descr))
  mc <- t(sapply(1:50, function(s) {
    tryCatch(run_one(s, scn),
             error = function(e) c(ipw_trt = NA, naive_trt = NA,
                                   weight_max = NA, weight_mean = NA))
  }))
  truth <- causal_parameters[4]
  cat(sprintf("  IPW   mean = %6.3f   bias = %+0.3f   SD = %.3f\n",
              mean(mc[, "ipw_trt"],   na.rm = TRUE),
              mean(mc[, "ipw_trt"],   na.rm = TRUE) - truth,
              sd(mc[, "ipw_trt"],     na.rm = TRUE)))
  cat(sprintf("  Naive mean = %6.3f   bias = %+0.3f   SD = %.3f\n",
              mean(mc[, "naive_trt"], na.rm = TRUE),
              mean(mc[, "naive_trt"], na.rm = TRUE) - truth,
              sd(mc[, "naive_trt"],   na.rm = TRUE)))
  cat(sprintf("  weight: mean(max) = %.1f, mean(mean) = %.2f\n",
              mean(mc[, "weight_max"],  na.rm = TRUE),
              mean(mc[, "weight_mean"], na.rm = TRUE)))
  mc
})
names(results) <- names(scenarios)

# ── Summary table ─────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("Summary: bias relative to truth = ", causal_parameters[4], "\n")
cat("══════════════════════════════════════════════════════════════════════\n")
summary_tbl <- t(sapply(results, function(mc) {
  c(ipw_bias   = mean(mc[, "ipw_trt"],   na.rm = TRUE) - causal_parameters[4],
    naive_bias = mean(mc[, "naive_trt"], na.rm = TRUE) - causal_parameters[4],
    ipw_sd     = sd(mc[, "ipw_trt"],     na.rm = TRUE),
    naive_sd   = sd(mc[, "naive_trt"],   na.rm = TRUE),
    bias_gap   = abs(mean(mc[, "naive_trt"], na.rm = TRUE) - causal_parameters[4]) -
                 abs(mean(mc[, "ipw_trt"],   na.rm = TRUE) - causal_parameters[4]))
}))
print(round(summary_tbl, 3))
cat("\nThe `bias_gap` column = |naive bias| - |IPW bias|. Bigger = better\n")
cat("evidence that IPW separates from naive under this parametrization.\n")
