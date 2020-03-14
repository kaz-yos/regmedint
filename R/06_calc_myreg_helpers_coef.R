################################################################################
### Helper functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


###
### coef extractors
################################################################################

##' Create a vector of coefficients from the mediator model (mreg)
##'
##' This function extracts \code{\link{coef}} from \code{mreg_fit} and pads with zeros appropriately to create a named vector consistently having the following elements:
##' \code{(Intercept)}:
##' \code{avar}
##' \code{cvar}: A zero element is added when \code{cvar = NULL}.
##'
##' @inheritParams regmedint
##'
##' @return A named numeric vector of coefficients.
beta_hat <- function(mreg, mreg_fit, avar, cvar) {
    vars <- c("(Intercept)", avar, cvar)
    coef(mreg_fit)[vars]
}

##' Create a vector of coefficients from the outcome model (yreg)
##'
##' This function extracts \code{\link{coef}} from \code{yreg_fit} and pads with zeros appropriately to create a named vector consistently having the following elements:
##' \code{(Intercept)}: A zero element is added for \code{yreg = "survCox"} for which no intercept is estimated (the baseline hazard is left unspecified).
##' \code{avar}
##' \code{mvar}
##' \code{avar:mvar}: A zero element is added when \code{interaction = FALSE}.
##' \code{cvar}: A zero element is added when \code{cvar = NULL}.
##'
##' @inheritParams regmedint
##'
##' @return A named numeric vector of coefficients.
theta_hat <- function(yreg, yreg_fit, avar, mvar, cvar, interaction) {

    coef_raw <- coef(yreg_fit)

    if (yreg == "survCox") {

        ## Pad with a zero for the missing (Intercept)
        coef_ready <- c(0, coef_raw)
        names(coef_ready) <- c("(Intercept)", names(coef_raw))

    } else {

        coef_ready <- coef_raw

    }


    if (!interaction) {

        ## Always have a position for an interaction term to ease subsequent manipulation.
        coef_ready <- c(coef_ready[c("(Intercept)",avar,mvar)],
                        0,
                        coef_ready[cvar])
        names(coef_ready) <- c("(Intercept)",
                               avar,mvar,
                               paste0(avar,":",mvar),
                               cvar)
    }

    ## Subset to ensure the ordering and error on non-existent element.
    vars <- c("(Intercept)", avar, mvar, paste0(avar,":",mvar), cvar)
    coef_ready[vars]
}

sigma_hat_sq <- function(mreg_fit) {
    sigma(mreg_fit)^2
}
