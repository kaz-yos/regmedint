################################################################################
### For package documentation only
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################

##' regmedint: A package for regression-based causal mediation analysis
##'
##' The package is an R implementation of regression-based closed-form causal mediation as originally described in Valeri & VanderWeele 2013 and Valeri & VanderWeele 2015 \url{https://hsph.harvard.edu/research/vanderweele-lab/tools-and-tutorials/}. The earlier version is a sister program of the SAS macro. The current extended version (version 1.0 and later) supports effect modification by covariates (treatment-covariate and mediator-covariate product terms) in mediator and outcome models.
##'
##' @section Fitting models:
##' Use the regmedint function to fit models and set up regression-based causal mediation analysis.
##'
##' @section Examining results:
##' Several methods are available to examine the regmedint object.
##' print
##' summary
##' coef
##' confint
##'
##' @importFrom stats coef confint pnorm qnorm sigma vcov pt residuals as.formula na.omit
##' @importFrom survival Surv
##'
##' @keywords internal
"_PACKAGE"



## http://r-pkgs.had.co.nz/data.html#data-data

##' Example dataset from Valeri and VanderWeele 2015.
##'
##' An example dataset from Valeri and VanderWeele (2015) <doi:10.1097/EDE.0000000000000253>.
##'
##' @format A tibble with 100 rows and 7 variables:
##' \describe{
##'   \item{id}{Positive integer id.}
##'   \item{x}{Binary treatment assignment variable.}
##'   \item{m}{Binary mediator variable.}
##'   \item{y}{Time to event outcome variable.}
##'   \item{cens}{Binary censoring indicator. Censored is 1.}
##'   \item{c}{Continuous confounder variable.}
##'   \item{event}{Binary event indicator. Event is 1.}
##' }
##'
##' @source \url{https://hsph.harvard.edu/research/vanderweele-lab/tools-and-tutorials/}
"vv2015"
