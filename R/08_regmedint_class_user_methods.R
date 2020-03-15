################################################################################
### User interface methods
##
## Created on: 2020-03-14
## Author: Kazuki Yoshida
################################################################################


###
### print method
################################################################################

##' print method for regmedint object
##'
##' Print the \code{mreg_fit}, \code{yreg_fit}, and the mediation analysis effect estimates.
##'
##' @param x An object of the \code{regmedint} class.
##' @param m_cde A numeric vector of length one. A mediator value at which the controlled direct effect (CDE) conditional on the adjustment covariates is evaluated. If not provided, the default value supplied to the call to \code{\link{regmedint}} will be used. Only the CDE is affected.
##' @param c_cond A numeric vector as long as the number of adjustment covariates. A set of covariate values at which the conditional natural effects are evaluated.
##' @param ...
##'
##' @return Invisibly return the \code{regmedint} class object as is.
print.regmedint <- function(x, m_cde = NULL, c_cond = NULL, , ...) {

    stop("print.regmedint not implemented")

}

###
### summary method
################################################################################

##' summary method for regmedint object
##'
##'
##'
##' .. content for \details{} ..
##'
##'
##' @param x
##' @param ...
##'
##' @return
summary.regmedint <- function(x, ...) {

    stop("summary.regmedint not implemented")

}


###
### Others
################################################################################

coef.regmedint <- function(x, ...) {

    stop("coef.regmedint not implemented")

}

confint.regmedint <- function(x, ...) {

    stop("confint.regmedint not implemented")

}
