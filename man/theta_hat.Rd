% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_calc_myreg_helpers_coef.R
\name{theta_hat}
\alias{theta_hat}
\title{Create a vector of coefficients from the outcome model (yreg)}
\usage{
theta_hat(yreg, yreg_fit, avar, mvar, cvar, interaction)
}
\arguments{
\item{yreg}{A character vector of length 1. Outcome regression type: \code{"linear"}, \code{"logistic"}, \code{"loglinear"}, \code{"poisson"}, \code{"negbin"}, \code{"survCox"}, \code{"survAFT_exp"}, or \code{"survAFT_weibull"}.}

\item{yreg_fit}{Model fit object for yreg (outcome model).}

\item{avar}{A character vector of length 1. Treatment variable name.}

\item{mvar}{A character vector of length 1. Mediator variable name.}

\item{cvar}{A character vector of length > 0. Covariate names. Use \code{NULL} if there is no covariate. However, this is a highly suspicious situation. Even if \code{avar} is randomized, \code{mvar} is not. Thus, there should usually be some confounder(s) to account for the common cause structure (confounding) between \code{mvar} and \code{yvar}.}

\item{interaction}{A logical vector of length 1. Default to TRUE. Whether to include a mediator-treatment interaction term in the outcome regression model.}
}
\value{
A named numeric vector of coefficients.
}
\description{
This function extracts \code{\link{coef}} from \code{yreg_fit} and pads with zeros appropriately to create a named vector consistently having the following elements:
\code{(Intercept)}: A zero element is added for \code{yreg = "survCox"} for which no intercept is estimated (the baseline hazard is left unspecified).
\code{avar}
\code{mvar}
\code{avar:mvar}: A zero element is added when \code{interaction = FALSE}.
\code{cvar}: This part is eliminated when \code{cvar = NULL}.
}
