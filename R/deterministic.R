##' Deterministic treatment-assignment rules
##'
##' Helpers for specifying a treatment that is a *deterministic* function of the
##' past, rather than a draw from a GLM.  The returned object is dropped into the
##' treatment slot of the `family` list passed to [msm_samp()] / [rmsm()], in the
##' same position where one would otherwise put a numeric family code (e.g. `5`
##' for binomial).
##'
##' @details
##' [dose_escalation()] encodes a dose-escalation regime over ordered dose
##' levels `0` (placebo), `1` (lowest active dose), ..., `max_level` (highest):
##'
##' * **time 0:** the dose is a coin flip between placebo (`0`) and the lowest
##'   active dose (`1`), with `P(dose = 1) = init_prob`;
##' * **no side effect** (`trigger(Z_k)` is `FALSE`): escalate one level, capped
##'   at `max_level`, i.e. `A_k = min(A_{k-1} + 1, max_level)`;
##' * **bad side effect** (`trigger(Z_k)` is `TRUE`): de-escalate one level,
##'   floored at placebo, i.e. `A_k = max(A_{k-1} - 1, 0)`.
##'
##' At the top level the same rule applies: with no side effect the patient
##' stays at `max_level`, and a side effect knocks them down a level.
##'
##' The rule is evaluated inside the per-time-step simulation block *after* the
##' time-varying covariates (including the side-effect variable) have been drawn,
##' so `Z_k` is available when `A_k` is assigned.  The treatment formula in
##' `formulas` is only used to register the treatment variable name and its place
##' in the variable ordering; the right-hand side is ignored (use `A ~ 1`).  No
##' `beta`/parameter entry is required for a deterministic treatment.
##'
##' Currently the regime is fully deterministic given the past and the
##' side-effect process (the only stochastic input is the time-0 coin flip), so
##' the conditional treatment "propensity" is degenerate (0/1).  The `slip`
##' argument is reserved for a future extension that adds a small per-transition
##' slip probability, restoring positivity for inverse-probability weighting.
##'
##' @param max_level highest dose level `L` (levels run `0, 1, ..., max_level`)
##' @param side_effect name (stem) of the time-varying side-effect variable
##'   (e.g. `"Z"`); the column `paste0(side_effect, "_", k)` is read at time `k`
##' @param init_prob probability the time-0 coin flip lands on the lowest active
##'   dose (`1`) rather than placebo (`0`)
##' @param trigger function mapping the side-effect variable to a logical vector
##'   that is `TRUE` where a bad side effect occurred (default treats any
##'   positive value as a side effect)
##' @param lock_placebo if `TRUE`, subjects whose time-0 assignment is placebo
##'   (`0`) are held at `0` for the entire follow-up: the placebo arm never
##'   escalates.  This distinguishes a randomized-placebo subject from an active
##'   subject who happens to have de-escalated down to `0` (the latter can still
##'   escalate again).  Defaults to `TRUE`.
##' @param slip reserved for future use; must be `0` (fully deterministic)
##'
##' @return An object of class `c("causl_deterministic", "causl_family")`.
##'
##' @examples
##' \dontrun{
##' forms <- list(list(B ~ 1),
##'               Z ~ B + A_l1,                    # side-effect covariate
##'               A ~ 1,                            # treatment: RHS ignored
##'               Y ~ B + I(A_0 + A_1 + A_2),
##'               Z ~ 1)
##' fams <- list(list(1),
##'              5,                                 # binary side effect
##'              dose_escalation(max_level = 3, side_effect = "Z"),
##'              1,
##'              1)
##' }
##'
##' @export
dose_escalation <- function (max_level, side_effect, init_prob = 0.5,
                             trigger = function (z) z > 0, lock_placebo = TRUE,
                             slip = 0) {
  if (missing(max_level)) stop("Must specify 'max_level' (highest dose level)")
  if (missing(side_effect)) stop("Must specify 'side_effect' (name of side-effect variable)")
  if (!isTRUE(all.equal(slip, 0))) {
    stop("Non-zero 'slip' is not supported yet; the regime is currently fully deterministic")
  }

  rule <- function (dat, k, vnm) {
    n <- nrow(dat)
    ## time 0: coin flip between placebo (0) and lowest active dose (1)
    if (k <= 0) return(stats::rbinom(n, size = 1, prob = init_prob))

    stem <- sub("_[0-9]+$", "", vnm)
    prev_col <- paste0(stem, "_", k - 1)
    if (is.null(dat[[prev_col]])) {
      stop(paste0("Previous treatment '", prev_col, "' not found when applying dose escalation"))
    }
    se_col <- paste0(side_effect, "_", k)
    if (is.null(dat[[se_col]])) {
      stop(paste0("Side-effect variable '", se_col, "' not found; it must be a ",
                  "time-varying covariate simulated before the treatment at each time"))
    }

    prev <- dat[[prev_col]]
    bad  <- as.logical(trigger(dat[[se_col]]))
    out  <- ifelse(bad, pmax(prev - 1, 0), pmin(prev + 1, max_level))

    ## lock the placebo arm: subjects whose time-0 assignment was placebo (0)
    ## are held at 0 for the entire follow-up.
    if (lock_placebo) {
      init_col <- paste0(stem, "_0")
      if (is.null(dat[[init_col]])) {
        stop(paste0("Time-0 treatment '", init_col, "' not found; required when ",
                    "lock_placebo = TRUE"))
      }
      out[dat[[init_col]] == 0] <- 0
    }

    out
  }

  structure(list(name = "deterministic",
                 rule = rule,
                 max_level = max_level,
                 side_effect = side_effect,
                 init_prob = init_prob,
                 trigger = trigger,
                 lock_placebo = lock_placebo,
                 slip = slip,
                 pars = character(0),
                 link = "identity",
                 custom_links = NULL),
            class = c("causl_deterministic", "causl_family"))
}

##' @describeIn dose_escalation Test whether a family entry is a deterministic rule
##' @param x an object to test
##' @export
is_deterministic <- function (x) {
  if (methods::is(x, "causl_deterministic")) return(TRUE)
  is.list(x) && length(x) > 0 && methods::is(x[[1]], "causl_deterministic")
}
