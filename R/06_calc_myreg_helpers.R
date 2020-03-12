################################################################################
### Helper functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


###
### coef extractors
################################################################################

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
##' \code{cvar}
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


###
### vcov extractors
################################################################################

Sigma_beta_hat <- function(mreg, mreg_fit, avar, cvar) {
    vars <- c("(Intercept)", avar, cvar)
    vcov(mreg_fit)[vars,vars, drop = FALSE]
}

Sigma_theta_hat <- function(yreg, yreg_fit, avar, mvar, cvar, interaction) {

    vcov_raw <- vcov(yreg_fit)

    ## Older versions of survival:::vcov.survreg() did not give dimension names.
    ## https://github.com/therneau/survival/commit/ed1c71b3817d4bfced43ed374e5e598e5f229bb8

    if (yreg == "survAFT_exp") {

        vcov_ready <- vcov_raw
        if (is.null(dimnames(vcov_ready))) {
            ## Older vcov.survreg gives an unnamed vcov matrix
            dimnames(vcov_ready) <- list(names(coef(yreg_fit)),
                                         names(coef(yreg_fit)))
        }

    } else if (yreg == "survAFT_weibull") {

        ## vcov.survreg(weibull_fit) has an extra row and column corresponding
        ## to the log(scale) parameter (See the above commit).
        coef_ind <- seq_along(coef(yreg_fit))
        vcov_ready <- vcov_raw[coef_ind,coef_ind]
        if (is.null(dimnames(vcov_ready))) {
            ## Older vcov.survreg gives an unnamed vcov matrix
            dimnames(vcov_ready) <- list(names(coef(yreg_fit)),
                                         names(coef(yreg_fit)))
        }

    } else if (yreg == "survCox") {

        ## Pad the left and upper edges with zeros by creating a block diagonal.
        vcov_ready <- Matrix::bdiag(matrix(0),
                                    vcov_raw)
        vars_cox <- c("(Intercept)", names(coef(yreg_fit)))
        dimnames(vcov_ready) <- list(vars_cox,
                                     vars_cox)

    } else {

        vcov_ready <- vcov_raw

    }

    if (!interaction) {
        ## Always have a position for an interaction term to ease subsequent manipulation.
        vcov_ready <- Matrix::bdiag(vcov_ready[c("(Intercept)",avar,mvar),
                                               c("(Intercept)",avar,mvar),
                                               drop = FALSE],
                                    ## Padding for the avar:mvar position.
                                    matrix(0),
                                    vcov_ready[cvar,
                                               cvar,
                                               drop = FALSE])
        vars_add_int <- c("(Intercept)",
                          avar,mvar,
                          paste0(avar,":",mvar),
                          cvar)
        dimnames(vcov_ready) <- list(vars_add_int,
                                     vars_add_int)
    }

    ## Subset to ensure the ordering and error on non-existent element.
    vars <- c("(Intercept)", avar, mvar, paste0(avar,":",mvar), cvar)
    vcov_ready[vars,vars, drop = FALSE]
}

Sigma_sigma_hat_sq <- function(mreg_fit) {

    ## FIXME: This is tentative and suspicious.
    ## Chapter 2
    ## Quadratic Forms of Random Variables
    ## http://pages.stat.wisc.edu/~st849-1/lectures/Ch02.pdf
    ## Corollary 6. In a full-rank Gaussian model for M with p covariates X.
    ## ||M - X beta_hat||^2 / sigma^2 ~ (central chi-squared DF = (n - p))
    ## Var(central chi-squared DF = (n - p)) = (n - p)
    ## Var(sigma_hat^2) = Var(1/(n-p) * ||M - X beta_hat||^2)
    ##                  = 1/(n-p)^2 * Var(||M - X beta_hat||^2 * sigma^2/sigma^2)
    ##                  = 1/(n-p)^2 * (sigma^2)^2 * Var(||M - X beta_hat||^2 /sigma^2)
    ##                  = 1/(n-p)^2 * (sigma^2)^2 * (n-p)
    ##                  = (sigma^2)^2 / (n-p)
    ##
    ## (sigma_hat^2)^2 / (n-p)
    ## VanderWeele 2015. p470 states??
    ## 2 * (sigma_hat^2)^2 / (n-p)
    matrix(((sigma(mreg_fit))^2)^2 / mreg_fit$df.residual)
}
