# deterministic dose-escalation treatment rule (dose_escalation())

set.seed(1)
n <- 1e4
T <- 7
L <- 3   # dose levels 0 (placebo), 1, ..., L (highest)

formulas <- list(C ~ 1,
                 Z ~ X_l1 + C,                      # time-varying side effect (binary)
                 X ~ 1,                              # treatment: deterministic, RHS ignored
                 Y ~ I(X_0 + X_1 + X_2 + X_3 + X_4 + X_5 + X_6) + C,
                 cop ~ 1)
family <- list(5,
               5,                                    # binary side-effect indicator
               dose_escalation(max_level = L, side_effect = "Z", init_prob = 0.5),
               3,
               1)
pars <- list(C = list(beta = 0),
             Z = list(beta = c(-1/2, 0.4, 0.25)),    # intercept + X_{t-1} + C
             X = list(),                             # deterministic: no parameters
             Y = list(beta = c(0.05, 0.45, 0.05), phi = 1),
             cop = list(Y = list(Z = list(beta = 0.8472979))))

surv_model <- survivl_model(T = T, formulas = formulas, family = family, pars = pars)
dat <- rmsm(n, surv_model)

A <- as.matrix(dat[, paste0("X_", 0:(T - 1))])
Z <- as.matrix(dat[, paste0("Z_", 0:(T - 1))])

test_that("dose_escalation() produces a valid deterministic regime", {
  ## doses stay within the allowed range
  expect_true(all(A >= 0 & A <= L))

  ## time-0 dose is a coin flip restricted to {placebo, lowest active dose}
  expect_setequal(unique(A[, 1]), c(0, 1))
  expect_equal(mean(A[, 1]), 0.5, tolerance = 0.05)

  ## the escalate / de-escalate rule holds exactly at every transition
  for (k in 2:T) {
    prev     <- A[, k - 1]
    bad      <- Z[, k] > 0
    expected <- ifelse(bad, pmax(prev - 1, 0), pmin(prev + 1, L))
    expect_identical(A[, k], expected)
  }
})

test_that("dose_escalation() input validation works", {
  expect_error(dose_escalation(side_effect = "Z"), "max_level")
  expect_error(dose_escalation(max_level = 3), "side_effect")
  expect_error(dose_escalation(max_level = 3, side_effect = "Z", slip = 0.1), "slip")
})

test_that("is_deterministic() recognises rule objects", {
  fam <- dose_escalation(max_level = 3, side_effect = "Z")
  expect_true(is_deterministic(fam))
  expect_true(is_deterministic(list(fam)))
  expect_false(is_deterministic(5))
  expect_false(is_deterministic(causl::gaussian_causl_fam()))
})

test_that("g-computation recovers the structural (MSM) parameters", {
  ## Under a deterministic regime the treatment propensity is degenerate (0/1),
  ## so IPW (as used in the other longitudinal tests) is undefined.  The correct
  ## estimator is the g-formula: intervene on the whole dose path, do(X_k = c),
  ## let the time-varying confounder Z follow its model, and read off Y.  Fitting
  ## the MSM to this interventional pseudo-data must recover pars$Y$beta.
  ##
  ## A constant-dose intervention is itself a deterministic family: a rule that
  ## ignores Z and returns the fixed level c.
  static_dose <- function (c) {
    structure(list(name = "deterministic",
                   rule = function (dat, k, vnm) rep(c, nrow(dat)),
                   pars = character(0), link = "identity", custom_links = NULL),
              class = c("causl_deterministic", "causl_family"))
  }

  do_pars <- pars$Y$beta   # c(0.05, 0.45, 0.05): intercept, cumulative dose, C
  n_int   <- 2e4

  gcomp <- do.call(rbind, lapply(0:L, function(c) {
    fam_c <- family
    fam_c[[3]] <- static_dose(c)
    sm  <- survivl_model(T = T, formulas = formulas, family = fam_c, pars = pars)
    d   <- rmsm(n_int, sm)
    data.frame(Y = d$Y, cumX = rowSums(d[, paste0("X_", 0:(T - 1))]), C = d$C)
  }))

  msm_fit <- glm(Y ~ cumX + C, data = gcomp, family = Gamma(link = "log"))
  sm_fit  <- summary(msm_fit)$coefficients

  expect_lt(max(abs(sm_fit[, 1] - do_pars) / sm_fit[, 2]), 2.5)
})
