################################################################################
### Internal helper functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################

##' Low level constructor for a regmedint S3 class object.
##'
##' This is not a user function and meant to be executed within the regmedint function after validatingthe arguments.
##'
##' @inheritParams regmedint
##'
##' @return A regmedint object.
new_regmedint <- function(data,
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
                          eventvar) {

    ## Perform mreg
    mreg_fit <- fit_mreg(mreg, data, avar, mvar, cvar)

    ## Perform yreg
    yreg_fit <- fit_yreg(yreg, data, yvar, avar, mvar, cvar, interaction, eventvar)

    ## Construct the result object
    res <- list(mreg = mreg_fit,
                yreg = yreg_fit,
                ## FIXME: made up numbers
                med = 42L,
                boot = 24L)
    ## The main class is regmedint.
    class(res) <- c("regmedint", class(res))

    ##
    res
}


##' Fit a model for the mediator given the treatment and covariates.
##'
##' \code{\link{lm}} is called if \code{mreg = "linear"}. \code{\link{glm}} is called with \code{family = binomial()} if \code{mreg = "logistic"}.
##'
##' @inheritParams regmedint
##'
##' @return A regression object of class lm (linear) or glm (logistic)
fit_mreg <- function(mreg,
                     data,
                     avar,
                     mvar,
                     cvar) {

}


## The third and subsequent paragraphs go into details.
## http://r-pkgs.had.co.nz/man.html#roxygen-comments

##' Fit a model for the outcome given the treatment, mediator, and covariates.
##'
##' The outcome model type \code{yreg} can be one of the following \code{"linear"}, \code{"logistic"}, \code{"loglinear"}, \code{"poisson"}, \code{"negbin"}, \code{"survCox"}, \code{"survAFT_exp"}, or \code{"survAFT_weibull"}.
##'
##' The outcome regression functions to be called are the following:
##' \itemize{
##'   \item \code{"linear"} \code{\link{lm}}
##'   \item \code{"logistic"} \code{\link{glm}}
##'   \item \code{"loglinear"} \code{\link{glm}}
##'   \item \code{"poisson"} \code{\link{glm}}
##'   \item \code{"negbin"} \code{\link[MASS]{glm.nb}}
##'   \item \code{"survCox"} \code{\link[survival]{coxph}}
##'   \item \code{"survAFT_exp"} \code{\link[survival]{survreg}}
##'   \item \code{"survAFT_weibull"} \code{\link[survival]{survreg}}
##' }
##'
##' @inheritParams regmedint
##'
##' @return
fit_yreg <- function(yreg,
                     data,
                     yvar,
                     avar,
                     mvar,
                     cvar,
                     interaction,
                     eventvar) {

}
