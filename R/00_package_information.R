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
##' @importFrom stats coef confint pnorm qnorm sigma vcov
##' @importFrom survival Surv
##'
##' @docType package
##' @name regmedint
NULL
