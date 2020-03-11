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

    ## Perform mediation analysis
    myreg_fit <- calc_myreg(mreg,
                            mreg_fit,
                            yreg,
                            yreg_fit,
                            interaction,
                            a0,
                            a1,
                            m_cde,
                            c_cond)

    ## Construct the result object
    res <- list(mreg = mreg_fit,
                yreg = yreg_fit,
                ## FIXME: made up numbers
                myreg = myreg_fit,
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

    ## Create a string representation of the formula
    string_formula <- string_mreg_formula(mvar,
                                          avar,
                                          cvar)

    ## Quasi-quoting to make the formula readable.
    ## bquote suppresses evaluation except within .(...).
    ## Evaluate restart the evaluation with the .() part
    ## already expanded.
    if (mreg == "linear") {

        eval(
            bquote(
                lm(formula = .(as.formula(string_formula)),
                   data = data)
            )
        )

    } else if (mreg == "logistic") {

        eval(
            bquote(
                glm(formula = .(as.formula(string_formula)),
                    family = binomial(link = "logit"),
                    data = data)
            )
        )

    } else {

        stop("Unsupported model type in yreg")

    }
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

    ## Create a string representation of the formula
    string_formula <- string_yreg_formula(yvar,
                                          avar,
                                          mvar,
                                          cvar,
                                          interaction,
                                          eventvar)

    ## Quasi-quoting to make the formula readable.
    ## bquote suppresses evaluation except within .(...).
    ## Evaluate restart the evaluation with the .() part
    ## already expanded.
    if (yreg == "linear") {

        eval(
            bquote(
                lm(formula = .(as.formula(string_formula)),
                   data = data)
            )
        )

    } else if (yreg == "logistic") {

        eval(
            bquote(
                glm(formula = .(as.formula(string_formula)),
                    family = binomial(link = "logit"),
                    data = data)
            )
        )

    } else if (yreg == "loglinear") {

    } else if (yreg == "poisson") {

        eval(
            bquote(
                glm(formula = .(as.formula(string_formula)),
                    family = poisson(link = "log"),
                    data = data)
            )
        )

    } else if (yreg == "negbin") {

    } else if (yreg == "survCox") {

    } else if (yreg == "survAFT_exp") {

        eval(
            bquote(
                survival::survreg(formula = .(as.formula(string_formula)),
                                  data = data,
                                  dist = "exponential")
            )
        )

    } else if (yreg == "survAFT_weibull") {

        eval(
            bquote(
                survival::survreg(formula = .(as.formula(string_formula)),
                                  data = data,
                                  dist = "weibull")
            )
        )

    } else {

        stop("Unsupported model type in yreg")

    }
}


###
### Formula string creators
################################################################################
string_mreg_formula <- function(mvar,
                                avar,
                                cvar) {

    if (is.null(cvar)) {
        acvar_string <- avar
    } else {
        cvar_string <- paste0(cvar, collapse = " + ")
        acvar_string <- paste0(c(avar, cvar_string), collapse = " + ")
    }

    sprintf("%s ~ %s", mvar, acvar_string)
}

string_yreg_formula <- function(yvar,
                                avar,
                                mvar,
                                cvar,
                                interaction,
                                eventvar) {

    ## Create A*M or A + M depending on interaction.
    if (interaction) {
        amvar_string <- sprintf("%s*%s", avar, mvar)
    } else {
        amvar_string <- sprintf("%s + %s", avar, mvar)
    }

    ## Add covariates if they exist.
    if (is.null(cvar)) {
        amcvar_string <- amvar_string
    } else {
        cvar_string <- paste0(cvar, collapse = " + ")
        amcvar_string <- sprintf("%s + %s", amvar_string, cvar_string)
    }

    ## eventvar must be NULL for a non-survival outcome model.
    if (is.null(eventvar)) {

        return(sprintf("%s ~ %s", yvar, amcvar_string))

    } else {

        ## Survival outcome
        surv_string <- sprintf("Surv(%s, %s)", yvar, eventvar)
        return(sprintf("%s ~ %s", surv_string, amcvar_string))

    }
}


###
### Functions for regression-based causal mediation analysis given two models
################################################################################

##' Calculate effect measures based on two regression fits.
##'
##' Causal effect measures are calculated given the mediator model fit (\code{mreg_fit}) and the outcome model fit (\code{yreg_fit}).
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return \code{myreg} object containing effect measures.
calc_myreg <- function(mreg,
                       mreg_fit,
                       yreg,
                       yreg_fit,
                       interaction,
                       a0,
                       a1,
                       m_cde,
                       c_cond) {

    ## Only four patterns as the non-linear yreg cases are the
    ## same as the logistic yreg case.
    ## See VanderWeele 2015 Appendix.
    ##     Valeri & VanderWeele 2013 Appendix.
    ##     README for this repo.
    if (mreg == "linear" & yreg == "linear") {

        ## VanderWeele 2015 p466 Proposition 2.3
        calc_myreg_mreg_linear_yreg_linear(mreg,
                                           mreg_fit,
                                           yreg,
                                           yreg_fit,
                                           interaction,
                                           a0,
                                           a1,
                                           m_cde,
                                           c_cond)

    } else if (mreg == "linear" & yreg != "linear") {

        ## VanderWeele 2015 p468 Proposition 2.4
        calc_myreg_mreg_linear_yreg_logistic(mreg,
                                             mreg_fit,
                                             yreg,
                                             yreg_fit,
                                             interaction,
                                             a0,
                                             a1,
                                             m_cde,
                                             c_cond)

    } else if (mreg == "logistic" & yreg == "linear") {

        ## VanderWeele 2015 p471 Proposition 2.5
        calc_myreg_mreg_logistic_yreg_linear(mreg,
                                             mreg_fit,
                                             yreg,
                                             yreg_fit,
                                             interaction,
                                             a0,
                                             a1,
                                             m_cde,
                                             c_cond)


    } else if (mreg == "logistic" & yreg != "linear") {

        ## VanderWeele 2015 p473 Proposition 2.6
        calc_myreg_mreg_logistic_yreg_logistic(mreg,
                                               mreg_fit,
                                               yreg,
                                               yreg_fit,
                                               interaction,
                                               a0,
                                               a1,
                                               m_cde,
                                               c_cond)

    } else  {

        stop("Unsupported mreg or yreg.")

    }

}
