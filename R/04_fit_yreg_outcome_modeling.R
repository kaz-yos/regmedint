################################################################################
### Internal helper functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################


###
### Model fitters for yreg
################################################################################

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

        ## https://github.com/mdonoghoe/logbin
        eval(
            bquote(
                logbin::logbin(formula = .(as.formula(string_formula)),
                               data = data)
            )
        )

    } else if (yreg == "poisson") {

        eval(
            bquote(
                glm(formula = .(as.formula(string_formula)),
                    family = poisson(link = "log"),
                    data = data)
            )
        )

    } else if (yreg == "negbin") {

        eval(
            bquote(
                MASS::glm.nb(formula = .(as.formula(string_formula)),
                             data = data)
            )
        )

    } else if (yreg == "survCox") {

        eval(
            bquote(
                survival::coxph(formula = .(as.formula(string_formula)),
                                data = data,
                                ties = "efron")
            )
        )

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

string_yreg_formula <- function(yvar,
                                avar,
                                mvar,
                                cvar,
                                interaction,
                                eventvar) {

    assertthat::assert_that(!is.null(mvar))
    assertthat::assert_that(!is.null(avar))
    assertthat::assert_that(!is.null(yvar))

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
