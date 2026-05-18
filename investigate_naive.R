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

# Strong L→A: high L (good hemoglobin) → A = 0; low L → A = 1.
# A_pars_sharp gives near-deterministic confounding (used in `sharp_*` scenarios).
A_pars       <- c(10, -1,   0.1)
A_pars_sharp <- c(50, -5,   0.1)

base_pars <- list(
  B = list(beta = 0.1),
  C = list(beta = 1),
  L = list(beta = c(11, -0.05, -0.02), phi = 0.5),
  A = list(beta = A_pars),
  Y = list(beta = causal_parameters)
)

# ── Helper to build dat0 / qtls for a given n ────────────────────────────────

make_init <- function(n, L_mean = function(B, C) 11 - 0.05 * B - 0.02 * C,
                      L_sd = 0.5) {
  qtls <- data.frame(matrix(runif(n * 4), ncol = 4))
  colnames(qtls) <- c("B", "C", "L_0", "A_0")
  qtls <- cbind(qtls, qtls[["A_0"]], qtls[["A_0"]])
  colnames(qtls) <- c("B", "C", "L_0", "A_0", "A_1", "A_2")

  B   <- qbinom(qtls[["B"]], 1, 0.1)
  C   <- qunif(qtls[["C"]], 25, 35)
  L_0 <- qnorm(qtls[["L_0"]], L_mean(B, C), L_sd)
  A_0 <- A_1 <- A_2 <- qbinom(qtls[["A_0"]], 1, 0.5)
  dat0 <- cbind(B, C, L_0, A_0, A_1, A_2)
  list(dat0 = dat0, qtls = qtls)
}

# ── IPW + naive estimation, returns the treatment coefficient from each ──────
# A_0, A_1, A_2 are fixed in dat0 and not freshly drawn by the model — only A_t
# for t >= 3 carry actual confounding from L_t. So weights are built from those.

fit_estimators <- function(d, T_use = 5) {
  modelled_t <- 3:(T_use - 1)
  A_vars <- paste0("A_", 0:(T_use - 1))
  sumA_form_str <- paste("Y ~ B + C + I(", paste(A_vars, collapse = "+"), ")")
  sumA_form <- as.formula(sumA_form_str)

  log_w <- rep(0, nrow(d))
  for (t in modelled_t) {
    if (t == 0) next  # never the case here, but defensive
    a_var <- paste0("A_", t); l_var <- paste0("L_", t); ap_var <- paste0("A_", t - 1)
    denom <- glm(as.formula(paste(a_var, "~", l_var, "+", ap_var)),
                 family = binomial(), data = d)
    numer <- glm(as.formula(paste(a_var, "~", ap_var)),
                 family = binomial(), data = d)
    pd <- predict(denom, type = "response")
    pn <- predict(numer, type = "response")
    a <- d[[a_var]]
    log_w <- log_w + log(ifelse(a == 1, pn, 1 - pn)) -
                     log(ifelse(a == 1, pd, 1 - pd))
  }
  w <- exp(log_w)
  d$w <- w

  # Y family: detect Gaussian vs binary by uniqueness
  is_gauss <- !all(d$Y %in% c(0, 1))
  fam <- if (is_gauss) gaussian() else binomial()
  fit_ipw   <- glm(sumA_form, data = d, weights = w, family = fam)
  fit_naive <- glm(sumA_form, data = d, family = fam)

  trt_nm <- grep("^I\\(", names(coef(fit_ipw)), value = TRUE)
  c(ipw_trt    = unname(coef(fit_ipw)[trt_nm]),
    naive_trt  = unname(coef(fit_naive)[trt_nm]),
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
    A_pars = A_pars,
    cop    = list(Y = list(L = list(k_tau = 0.1, par2 = 5)))
  ),

  strong_gaussian = list(
    descr = "Strong Gaussian copula (k_tau = 0.85). Heavy residual L–Y dep.",
    forms5 = L ~ 1,
    fam5   = 1,
    A_pars = A_pars,
    cop    = list(Y = list(L = list(k_tau = 0.85, par2 = 5)))
  ),

  # ─────────────────────────────────────────────────────────────────────────
  # SIMPSON-FLIP CANDIDATES
  # Negative L–Y dependence + low-L-treated  ⇒  treated have HIGH Y
  # via the copula even though A→Y structurally is −0.5. Naive coef can
  # flip sign once the confounding pathway dominates.
  # ─────────────────────────────────────────────────────────────────────────

  neg_gaussian = list(
    descr = "Negative Gaussian (k_tau = -0.85). Treated have low L → high Y.",
    forms5 = L ~ 1,
    fam5   = 1,
    A_pars = A_pars,
    cop    = list(Y = list(L = list(k_tau = -0.85, par2 = 5)))
  ),

  neg_extreme = list(
    descr = "Negative Gaussian (k_tau = -0.95).",
    forms5 = L ~ 1,
    fam5   = 1,
    A_pars = A_pars,
    cop    = list(Y = list(L = list(k_tau = -0.95, par2 = 5)))
  ),

  neg_t_heavy = list(
    descr = "Negative t-copula (k_tau = -0.9, df = 2). Tail dep + strong neg.",
    forms5 = L ~ 1,
    fam5   = 2,
    A_pars = A_pars,
    cop    = list(Y = list(L = list(k_tau = -0.9, df = 2, par2 = 2)))
  ),

  neg_conditional = list(
    # Conditional t-copula: rho ≈ -0.3 when A_l1 = 0, rho ≈ -0.95 when A_l1 = 1.
    # During treated periods the L–Y residual dep is almost a perfect inverse.
    descr = "Conditional t-copula: rho = -0.3 → -0.95 as A_l1 = 0 → 1.",
    forms5 = list(Y = list(L ~ A_l1)),
    fam5   = 2,
    A_pars = A_pars,
    cop    = list(Y = list(L = list(
      beta = c(rho_to_beta(-0.30), rho_to_beta(-0.95) - rho_to_beta(-0.30)),
      df   = 3,
      par2 = 3
    )))
  ),

  sharp_neg = list(
    # Same as neg_extreme but with NEAR-DETERMINISTIC L→A.
    # IPW positivity should still hold (no logical zeros) but weights swell.
    descr = "Sharp L→A (logit slope -5) + negative Gaussian k_tau = -0.95.",
    forms5 = L ~ 1,
    fam5   = 1,
    A_pars = A_pars_sharp,
    cop    = list(Y = list(L = list(k_tau = -0.95, par2 = 5)))
  ),

  sharp_conditional = list(
    descr = "Sharp L→A + conditional t-copula rho = -0.3 → -0.95.",
    forms5 = list(Y = list(L ~ A_l1)),
    fam5   = 2,
    A_pars = A_pars_sharp,
    cop    = list(Y = list(L = list(
      beta = c(rho_to_beta(-0.30), rho_to_beta(-0.95) - rho_to_beta(-0.30)),
      df   = 3,
      par2 = 3
    )))
  ),

  # ─────────────────────────────────────────────────────────────────────────
  # SIMPSON-FLIP, TAKE 2: REVERSE THE L→A DIRECTION
  # With the original direction (low-L → treated) and a negative copula the
  # two bias channels cancel — that's why naive sat at the truth.
  # If instead high-L → treated AND copula is strongly positive, both channels
  # push naive in the SAME direction (+) which is OPPOSITE the true effect (-),
  # so the naive coefficient can actually flip sign.
  # Also shrink the true causal effect from −0.5 to −0.2 so the confounding
  # has less to overcome.
  # ─────────────────────────────────────────────────────────────────────────

  flip_pos = list(
    descr = "REVERSED L→A (high L treated) + pos k_tau = 0.95, truth = -0.2.",
    forms5 = L ~ 1,
    fam5   = 1,
    A_pars = c(-10, 1, 0.1),       # high L → A = 1 (reversed)
    Y_pars = c(-2, 0.1, 0.02, -0.2),
    cop    = list(Y = list(L = list(k_tau = 0.95, par2 = 5)))
  ),

  flip_pos_t = list(
    descr = "Reversed L→A + pos t-copula (k_tau = 0.95, df = 2), truth = -0.2.",
    forms5 = L ~ 1,
    fam5   = 2,
    A_pars = c(-10, 1, 0.1),
    Y_pars = c(-2, 0.1, 0.02, -0.2),
    cop    = list(Y = list(L = list(k_tau = 0.95, df = 2, par2 = 2)))
  ),

  flip_cond = list(
    # Reversed L→A and a conditional copula whose rho swings 0.3 → 0.98
    # as A_l1 goes 0 → 1: under sustained treatment the L–Y link is almost
    # comonotonic, on top of A_l1 already carrying confounding info.
    descr = "Reversed L→A + conditional t (rho 0.3 → 0.98), truth = -0.2.",
    forms5 = list(Y = list(L ~ A_l1)),
    fam5   = 2,
    A_pars = c(-10, 1, 0.1),
    Y_pars = c(-2, 0.1, 0.02, -0.2),
    cop    = list(Y = list(L = list(
      beta = c(rho_to_beta(0.30), rho_to_beta(0.98) - rho_to_beta(0.30)),
      df   = 3,
      par2 = 3
    )))
  ),

  flip_sharp = list(
    descr = "Sharp reversed L→A (slope +5) + pos k_tau = 0.98, truth = -0.2.",
    forms5 = L ~ 1,
    fam5   = 1,
    A_pars = c(-50, 5, 0.1),       # near-deterministic high-L → A = 1
    Y_pars = c(-2, 0.1, 0.02, -0.2),
    cop    = list(Y = list(L = list(k_tau = 0.98, par2 = 5)))
  ),

  # ─────────────────────────────────────────────────────────────────────────
  # SIMPSON-FLIP, TAKE 3: WIDE L + TINY CAUSAL EFFECT
  # The reversed-L→A flip scenarios above didn't fire because L's variance
  # (phi=0.5) was too small for strong confounding to develop on top of a
  # −0.2 structural effect over 5 periods. Here we (a) center L at 0 with
  # sd = 2, (b) drop the structural causal effect to a tiny −0.05, so the
  # confounding pathway has the upper hand and naive can plausibly go +.
  # ─────────────────────────────────────────────────────────────────────────

  flip_wide = list(
    descr = "Wide L (sd=2) + reversed L→A slope 3 + k_tau=0.95, truth=-0.05.",
    forms5 = L ~ 1,
    fam5   = 1,
    wide_L = TRUE,
    L_pars = c(0, 0, 0),               # L_t ~ N(0, phi=2) after t=0
    L_phi  = 2,
    A_pars = c(0, 3, 0.1),             # logit(P(A=1)) = 3 L (reversed, strong)
    Y_pars = c(-2, 0.1, 0.02, -0.05),  # tiny true effect
    cop    = list(Y = list(L = list(k_tau = 0.95, par2 = 5)))
  ),

  flip_wide_extreme = list(
    descr = "Wide L + reversed slope 3 + k_tau=0.99, df=2, truth=-0.05.",
    forms5 = L ~ 1,
    fam5   = 2,
    wide_L = TRUE,
    L_pars = c(0, 0, 0),
    L_phi  = 2,
    A_pars = c(0, 3, 0.1),
    Y_pars = c(-2, 0.1, 0.02, -0.05),
    cop    = list(Y = list(L = list(k_tau = 0.99, df = 2, par2 = 2)))
  ),

  flip_wide_null = list(
    # NULL true effect: any non-zero naive coefficient is "confounded effect"
    # — the cleanest illustration of confounding overriding the (zero) truth.
    descr = "Wide L + reversed + k_tau=0.95, TRUE EFFECT = 0 (null).",
    forms5 = L ~ 1,
    fam5   = 1,
    wide_L = TRUE,
    L_pars = c(0, 0, 0),
    L_phi  = 2,
    A_pars = c(0, 3, 0.1),
    Y_pars = c(-2, 0.1, 0.02, 0),
    cop    = list(Y = list(L = list(k_tau = 0.95, par2 = 5)))
  ),

  # ─────────────────────────────────────────────────────────────────────────
  # SIMPSON-FLIP, TAKE 4: WORKS!
  # The earlier flip attempts failed for two structural reasons:
  #  1. dat0 fixes A_0=A_1=A_2 → only A_3, A_4 carry confounding under T=5,
  #     so 3 of 5 doses dilute the bias to near-zero. Increasing T to 8 gives
  #     5 model-generated confounded periods (A_3..A_7).
  #  2. With binary Y, the structural −0.5 logit per A dominates the bounded
  #     copula contribution. Shrinking the true effect to −0.02 logit per A
  #     leaves room for the confounding pathway to win.
  # Combined with reversed L→A + strong A persistence + conditional copula
  # whose rho swings 0.3 → 0.99 with past treatment, the naive coefficient
  # goes POSITIVE (e.g. +0.09) even though the truth is −0.02 — a clean
  # Simpson's-paradox sign flip.
  # ─────────────────────────────────────────────────────────────────────────

  simpson_flip = list(
    # T=8 (5 model-generated confounded periods A_3..A_7), wide L (sd=2),
    # reversed L→A with strong persistence, conditional t-copula whose rho
    # swings 0.3 → 0.95 with past treatment, and a tiny structural effect of
    # -0.02 logit per A. Empirically yields naive ≈ +0.13 (sign-flipped) while
    # IPW recovers ≈ truth.
    descr = "★ SIMPSON FLIP: naive coef goes POSITIVE while truth is -0.02.",
    T      = 8,
    forms5 = list(Y = list(L ~ A_l1)),
    fam5   = 2,
    wide_L = TRUE,
    L_pars = c(0, 0, 0),
    L_phi  = 2,
    A_pars = c(0, 1, 0.5),               # reversed L→A + mild persistence
    Y_pars = c(-2, 0, 0, -0.02),         # tiny structural effect
    cop    = list(Y = list(L = list(
      beta = c(rho_to_beta(0.30), rho_to_beta(0.95) - rho_to_beta(0.30)),
      df   = 3,
      par2 = 3
    )))
  ),

  simpson_null = list(
    # Identical recipe but with truth = 0. Naive should still be clearly
    # positive — pure confounding fabricates an apparent harmful effect.
    descr = "★ SIMPSON null: same recipe, TRUE EFFECT = 0.",
    T      = 8,
    forms5 = list(Y = list(L ~ A_l1)),
    fam5   = 2,
    wide_L = TRUE,
    L_pars = c(0, 0, 0),
    L_phi  = 2,
    A_pars = c(0, 1, 0.5),
    Y_pars = c(-2, 0, 0, 0),
    cop    = list(Y = list(L = list(
      beta = c(rho_to_beta(0.30), rho_to_beta(0.95) - rho_to_beta(0.30)),
      df   = 3,
      par2 = 3
    )))
  )
  # NOTE: Gumbel/Clayton omitted — Gumbel hangs the inversion loop and
  # Clayton supports only positive Kendall's tau, so it can't deliver a flip.
)

# ── Build a survivl_model for a given scenario ────────────────────────────────

build_model <- function(scn, dat0, qtls) {
  T_use <- if (!is.null(scn$T)) scn$T else 5
  A_vars <- paste0("A_", 0:(T_use - 1))
  sumA_form <- as.formula(paste("Y ~ B + C + I(",
                                paste(A_vars, collapse = "+"), ")"))
  f <- forms
  f[[4]] <- sumA_form
  f[[5]] <- scn$forms5
  fa <- fams
  fa[[5]] <- scn$fam5
  if (!is.null(scn$Y_fam)) fa[[4]] <- scn$Y_fam
  p <- base_pars
  if (!is.null(scn$A_pars))  p$A$beta <- scn$A_pars
  if (!is.null(scn$Y_pars))  p$Y$beta <- scn$Y_pars
  if (!is.null(scn$Y_phi))   p$Y$phi  <- scn$Y_phi
  if (!is.null(scn$L_pars))  p$L$beta <- scn$L_pars
  if (!is.null(scn$L_phi))   p$L$phi  <- scn$L_phi
  p$cop <- scn$cop
  survivl_model(formulas = f, family = fa, pars = p,
                T = T_use, dat = dat0, qtls = qtls)
}

scn_T <- function(scn) if (!is.null(scn$T)) scn$T else 5

# Allow scenarios to override the structural causal effect (truth) on sum(A).
scn_truth <- function(scn) {
  if (!is.null(scn$Y_pars)) scn$Y_pars[4] else causal_parameters[4]
}

# ── Single-shot diagnostic (large n) for the most interesting scenario ──────

set.seed(42)
n <- 1e4
init <- make_init(n, L_mean = function(B, C) 0, L_sd = 2)

cat("══════════════════════════════════════════════════════════════════════\n")
cat("Single-run diagnostic on simpson_flip, n =", n, "\n")
cat("══════════════════════════════════════════════════════════════════════\n")

sm  <- build_model(scenarios$simpson_flip, init$dat0, init$qtls)
dat <- rmsm(n, sm)

est <- fit_estimators(dat, T_use = scn_T(scenarios$simpson_flip))
cat(sprintf("truth on sum_A coef = %.3f\n", scn_truth(scenarios$simpson_flip)))
cat(sprintf("  IPW   estimate     = %.3f  (bias %+0.3f)\n",
            est["ipw_trt"],   est["ipw_trt"]   - scn_truth(scenarios$simpson_flip)))
cat(sprintf("  Naive estimate     = %.3f  (bias %+0.3f)\n",
            est["naive_trt"], est["naive_trt"] - scn_truth(scenarios$simpson_flip)))
cat(sprintf("  max IPW weight = %.2f, mean weight = %.3f\n",
            est["weight_max"], est["weight_mean"]))

T_diag <- scn_T(scenarios$simpson_flip)
sum_A  <- rowSums(dat[, paste0("A_", 0:(T_diag - 1))])
mean_L <- rowMeans(dat[, paste0("L_", 0:(T_diag - 1))])
cat(sprintf("\ncor(sum_A, mean_L) = %.3f\n",  cor(sum_A, mean_L)))
cat(sprintf("cor(mean_L, Y)     = %.3f\n",    cor(mean_L, dat$Y)))
cat(sprintf("cor(sum_A, Y)      = %.3f\n",    cor(sum_A, dat$Y)))

# ── Monte Carlo sweep across scenarios ────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("Monte Carlo sweep: 50 reps, n = 1000, across copula parametrizations\n")
cat("══════════════════════════════════════════════════════════════════════\n")

run_one <- function(seed, scn) {
  set.seed(seed)
  if (isTRUE(scn$wide_L)) {
    init_mc <- make_init(1000, L_mean = function(B, C) 0, L_sd = 2)
  } else {
    init_mc <- make_init(1000)
  }
  sm <- suppressMessages(build_model(scn, init_mc$dat0, init_mc$qtls))
  d  <- rmsm(1000, sm)
  fit_estimators(d, T_use = scn_T(scn))
}

results <- lapply(names(scenarios), function(nm) {
  scn <- scenarios[[nm]]
  truth <- scn_truth(scn)
  cat(sprintf("\n--- %s ---\n%s\n  (truth on sum_A = %+0.2f)\n",
              nm, scn$descr, truth))
  mc <- t(sapply(1:50, function(s) {
    tryCatch(run_one(s, scn),
             error = function(e) c(ipw_trt = NA, naive_trt = NA,
                                   weight_max = NA, weight_mean = NA))
  }))
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
truths <- sapply(scenarios, scn_truth)

# ── Summary table ─────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("Summary (per-scenario truth shown)\n")
cat("══════════════════════════════════════════════════════════════════════\n")
summary_tbl <- t(mapply(function(mc, truth) {
  ipw_mean   <- mean(mc[, "ipw_trt"],   na.rm = TRUE)
  naive_mean <- mean(mc[, "naive_trt"], na.rm = TRUE)
  c(truth      = truth,
    ipw_est    = ipw_mean,
    naive_est  = naive_mean,
    ipw_bias   = ipw_mean   - truth,
    naive_bias = naive_mean - truth,
    bias_gap   = abs(naive_mean - truth) - abs(ipw_mean - truth),
    sign_flip  = if (truth == 0) NA_real_ else
                 as.numeric(sign(naive_mean) != sign(truth)))
}, results, truths))
print(round(summary_tbl, 3))
cat("\nbias_gap = |naive bias| - |IPW bias|; sign_flip = 1 means the naive\n")
cat("coefficient has the OPPOSITE sign from the truth (Simpson's paradox).\n")
