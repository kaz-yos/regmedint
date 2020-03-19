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
##' @param cvar A character vector of length > 0. Covariate names. Use \code{NULL} if there is no covariate. However, this is a highly suspicious situation. Even if \code{avar} is randomized, \code{mvar} is not. Thus, there should usually be some confounder(s) to account for the common cause structure (confounding) between \code{avar} and \code{yvar}.
##' @param a0 A numeric vector of length 1. Reference level of treatment variable that is considered "untreated" or "unexposed".
##' @param a1 A numeric vector of length 1.
##' @param m_cde A numeric vector of length 1. Mediator level at which controlled direct effect is evaluated at.
##' @param c_cond A numeric vector of the same length as cvar. Required for the ful output. The conditional effect
##' @param mreg A character vector of length 1. Outcome regression type. One of "linear"
##' @param yreg A character vector of length 1. Outcome regression type. One of
##' @param interaction A logical vector of length 1. Default to TRUE. Whether to include a mediator-treatment interaction term in the outcome regression model.
##' @param casecontrol A logical vector of length 1. Default to FALSE. Whether data comes from a case-control study.
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
                      c_cond,
                      mreg,
                      yreg,
                      interaction = TRUE,
                      casecontrol = FALSE,
                      eventvar = NULL) {
    ## This is the user-friendly helper function with a name that is the class name.
    ## https://adv-r.hadley.nz/s3.html#helpers

    ## Check data contains yvar, avar, mvar, cvar, eventvar (if provided)
    validate_args(data = data,
                  yvar = yvar,
                  avar = avar,
                  mvar = mvar,
                  cvar = cvar,
                  a0 = a0,
                  a1 = a1,
                  m_cde = m_cde,
                  c_cond = c_cond,
                  mreg = mreg,
                  yreg = yreg,
                  interaction = interaction,
                  casecontrol = casecontrol,
                  eventvar = eventvar)

    ## Construct a regmedint object after argument validation.
    ## This is the low-level constructor function.
    ## https://adv-r.hadley.nz/s3.html#s3-constructor
    res <- new_regmedint(data = data,
                         yvar = yvar,
                         avar = avar,
                         mvar = mvar,
                         cvar = cvar,
                         a0 = a0,
                         a1 = a1,
                         m_cde = m_cde,
                         c_cond = c_cond,
                         mreg = mreg,
                         yreg = yreg,
                         interaction = interaction,
                         casecontrol = casecontrol,
                         eventvar = eventvar)

    ## Check the resulting object for anomalies.
    ## This is the class validator function.
    ## https://adv-r.hadley.nz/s3.html#validators
    validate_regmedint(res)

    ## Return results
    res
}


###
### Argument validation function
################################################################################

##' Validate arguments to regmedint before passing to other functions
##'
##' Internal functions (usually) do not validate arguments, thus, we need to make sure informative errors are raised when the arguments are not safe for subsequent computation.
##'
##' @inheritParams regmedint
##'
##' @return NULL
validate_args <- function(data,
                          yvar,
                          avar,
                          mvar,
                          cvar,
                          a0,
                          a1,
                          m_cde,
                          c_cond,
                          mreg,
                          yreg,
                          interaction,
                          casecontrol,
                          eventvar) {

    assertthat::assert_that(is.data.frame(data))
    ##
    assertthat::assert_that(is.character(yvar))
    assertthat::assert_that(length(yvar) == 1)
    ##
    assertthat::assert_that(is.character(avar))
    assertthat::assert_that(length(avar) == 1)
    ##
    assertthat::assert_that(is.character(mvar))
    assertthat::assert_that(length(mvar) == 1)
    ##
    assertthat::assert_that(is.null(cvar) | is.character(cvar))
    ##
    assertthat::assert_that(is.numeric(a0))
    assertthat::assert_that(length(a0) == 1)
    ##
    assertthat::assert_that(is.numeric(a1))
    assertthat::assert_that(length(a1) == 1)
    ##
    assertthat::assert_that(is.numeric(m_cde))
    assertthat::assert_that(length(m_cde) == 1)
    ##
    assertthat::assert_that(is.null(c_cond) | is.numeric(c_cond))
    assertthat::assert_that(length(c_cond) == length(cvar))
    ##
    assertthat::assert_that(is.character(mreg))
    assertthat::assert_that(length(mreg) == 1)
    assertthat::assert_that(mreg %in% c("linear","logistic"))
    ##
    assertthat::assert_that(is.character(yreg))
    assertthat::assert_that(length(yreg) == 1)
    assertthat::assert_that(yreg %in% c("linear",
                                        "logistic","loglinear","poisson","negbin",
                                        "survCox","survAFT_exp","survAFT_weibull"))
    ##
    assertthat::assert_that(is.logical(interaction))
    assertthat::assert_that(length(interaction) == 1)
    ##
    assertthat::assert_that(is.logical(casecontrol))
    assertthat::assert_that(length(casecontrol) == 1)
    ##
    assertthat::assert_that(is.null(eventvar) | is.character(eventvar))
    assertthat::assert_that(length(eventvar) <= 1)

    NULL
}


###
### Argument validation function
################################################################################

##' Validate soundness of a regmedint object.
##'
##' Check the structure of a proposed regmedint object for soundness.
##'
##' @param x A \code{regmedint} object.
##'
##' @return NULL
validate_regmedint <- function(x) {

    assertthat::assert_that(class(x)[[1]] == "regmedint")
    ##
    assertthat::assert_that("myreg" %in% names(x))
    assertthat::assert_that("args" %in% names(x))
    ##
    assertthat::assert_that(is.list(x$args))

    NULL
}
