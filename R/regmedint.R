################################################################################
### User interface functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################


###
### Main user interface function
################################################################################

##
## https://adv-r.hadley.nz/s3.html#s3-classes
## Follow the recommendation for S3 classes in Advanced R.
## Provide a low-level constructor, validator, and user-friendly helper.
##
## https://adv-r.hadley.nz/s3.html#s3-constructor
## https://adv-r.hadley.nz/s3.html#validators
## https://adv-r.hadley.nz/s3.html#helpers


##' Conduct regression-based causal mediation analysis
##'
##' This is a user-interface for regression-based causal mediation analysis as described in Valeri & VanderWeele 2013 and Valeri & VanderWeele 2015.
##'
##' @param data Data frame containing the relevant variables.
##' @param yvar A character vector of length 1. Outcome variable name. It should be the time variable for survival outcomes.
##' @param avar A character vector of length 1. Treatment variable name.
##' @param mvar A character vector of length 1. Mediator variable name.
##' @param cvar A character vector of length > 0. Covariate names.
##' @param a0 A numeric vector of length 1. Reference level of treatment variable that is considered "untreated" or "unexposed".
##' @param a1 A numeric vector of length 1.
##' @param m_cde A numeric vector of length 1. Mediator level at which controlled direct effect is evaluated at.
##' @param yreg A character vector of length 1. Outcome regression type. One of
##' @param mreg A character vector of length 1. Outcome regression type. One of "linear"
##' @param interaction A logical vector of length 1. Default to TRUE. Whether to include a mediator-treatment interaction term in the outcome regression model.
##' @param casecontrol A logical vector of length 1. Default to FALSE. Whether data comes from a case-control study.
##' @param full_output A logical vector of length 1. Default to FALSE. Whether to give a full output containing both pure and total natural effects. When FALSE, the natural direct effect is the pure natural direct effect
##' @param c_cond A numeric vector of the same length as cvar. Required for the ful output. The conditional effect
##' @param boot A numeric vector of length 1. Default to 0. Number of bootstrap iterations used to give bootstrap confidence intervals.
##' @param eventvar An character vector of length 1. Only required for survival outcome regression models. Note that the coding is 1 for event and 0 for censoring, following the R survival package convention.
##'
##' @return regmedint object, which is a list containing the mediator regression object, the outcome regression object, and the regression-based mediation results.
##'
##' @export
regmedint <- function(data,
                      yvar,
                      avar,
                      mvar,
                      cvar,
                      a0,
                      a1,
                      m_cde,
                      yreg,
                      mreg,
                      interaction = TRUE,
                      casecontrol = FALSE,
                      full_output = FALSE,
                      c_cond = NULL,
                      boot = 0L,
                      eventvar = NULL) {
    ## This is the user-friendly helper function with a name that is the class name.
    ## https://adv-r.hadley.nz/s3.html#helpers

    ## Check data contains yvar, avar, mvar, cvar, eventvar (if provided)
    validate_args(data,
                  yvar,
                  avar,
                  mvar,
                  cvar,
                  a0,
                  a1,
                  m_cde,
                  yreg,
                  mreg,
                  interaction,
                  casecontrol,
                  full_output,
                  c_cond,
                  boot,
                  eventvar)

    ## Construct a regmedint object after argument validation.
    ## This is the low-level constructor function.
    ## https://adv-r.hadley.nz/s3.html#s3-constructor
    res <- new_regmedint(data,
                         yvar,
                         avar,
                         mvar,
                         cvar,
                         a0,
                         a1,
                         m_cde,
                         yreg,
                         mreg,
                         interaction,
                         casecontrol,
                         full_output,
                         c_cond,
                         boot,
                         eventvar)

    ## Check the resulting object for anomalies.
    ## This is the class validator function.
    ## https://adv-r.hadley.nz/s3.html#validators
    validate_regmedint(res)

    ## Return results
    res
}


###
### User interface methods for the regmedint class
################################################################################

## print method
print.regmedint <- function(x, ...) {

}

## summary method
summary.regmedint <- function(x, ...) {

}
