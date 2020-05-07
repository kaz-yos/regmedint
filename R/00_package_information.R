################################################################################
### For package documentation only
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################

##' regmedint: A package for regression-based causal mediation analysis
##'
##' The package is a simple R implementation of the SAS macro as described in Valeri & VanderWeele 2013 and Valeri & VanderWeele 2015 \url{https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/}.
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
##' FIXME: Document once implemented.
##'
##' @importFrom stats coef confint pnorm qnorm sigma vcov pt residuals
##' @importFrom survival Surv
##'
##' @docType package
##' @name regmedint
NULL


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
##' @source \url{https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/}
"vv2015"
