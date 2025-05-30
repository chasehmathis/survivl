% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msm_samp.R
\name{msm_samp}
\alias{msm_samp}
\alias{coxSamp}
\alias{cox_samp}
\title{Simulate from Cox MSM Model}
\usage{
msm_samp(
  n,
  dat = NULL,
  qtls = NULL,
  T,
  formulas,
  family,
  pars,
  link = NULL,
  method = "inversion",
  control = list()
)

coxSamp(n, T, formulas, family, pars, link = NULL, control = list())

cox_samp(n, T, formulas, family, pars, link = NULL, control = list())
}
\arguments{
\item{n}{number of samples}

\item{dat}{optional data frame for plasmode simulation}

\item{qtls}{associated quantiles from dat for plasmode simulation}

\item{T}{number of time points}

\item{formulas}{list of formulas to use}

\item{family}{list of families to use}

\item{pars}{list of parameter settings to use}

\item{link}{link functions for GLM-like models}

\item{method}{sampling method (defaults to \code{"inversion"})}

\item{control}{list of control parameters}
}
\value{
An object of class \code{survivl_dat} containing the simulated data.
}
\description{
Obtain samples from a specified Cox Marginal Structural Model
using the frugal parameterization.
}
\details{
Samples from a Marginal Structural Model specified
by a frugal parameterization; that is, one must specify the
marginal model (i.e. the dependence of the outcome on the
treatments), the distribution of the covariates, observed
confounders and treatments, and then a copula to join the distribution
of the outcome to that of the confounders.

Among the left-hand sides of outcome variables, the variable 'Cen' has a
special meaning as censoring.  This keyword can be changed to something
else by using the argument \code{censor} in the \code{control} list.
}
\section{Functions}{
\itemize{
\item \code{coxSamp()}: old name

\item \code{cox_samp()}: old name

}}
