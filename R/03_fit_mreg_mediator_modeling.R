################################################################################
### Internal helper functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################


###
### Model fitters for mreg
################################################################################

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


###
### Formula string creators
################################################################################
string_mreg_formula <- function(mvar,
                                avar,
                                cvar) {

    assertthat::assert_that(!is.null(mvar))
    assertthat::assert_that(!is.null(avar))

    if (is.null(cvar)) {
        acvar_string <- avar
    } else {
        cvar_string <- paste0(cvar, collapse = " + ")
        acvar_string <- paste0(c(avar, cvar_string), collapse = " + ")
    }

    sprintf("%s ~ %s", mvar, acvar_string)
}
