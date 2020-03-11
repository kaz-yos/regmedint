################################################################################
### Main functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
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
    ## Each helper function is defined in a dedicated file.
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


###
### calc_myreg helpers
################################################################################


Sigma_beta_hat <- function(mreg, mreg_fit, avar, cvar) {
    vars <- c("(Intercept)", avar, cvar)
    vcov(mreg_fit)[vars,vars, drop = FALSE]
}

Sigma_theta_hat <- function(yreg, yreg_fit, avar, mvar, cvar, interaction) {

    vcov_raw <- vcov(yreg_fit)

    if (yreg == "survAFT_exp") {

        ## vcov.survreg gives an unnamed vcov matrix
        vcov_ready <- vcov_raw
        dimnames(vcov_ready) <- list(coef(yreg_fit),
                                     coef(yreg_fit))

    } else if (yreg == "survAFT_weibull") {

        ## vcov.survreg(weibull_fit) has a scale parameter variance
        ## as the right lower corner.
        coef_ind <- seq_along(coef(yreg_fit))
        vcov_ready <- vcov_raw[coef_ind,coef_ind]
        dimnames(vcov_ready) <- list(coef(yreg_fit),
                                     coef(yreg_fit))

    } else if (yreg == "survCox") {

        ## Pad the left and upper edges with zeros by creating a block diagonal.
        vcov_ready <- Matrix::bdiag(matrix(0), vcov_raw)
        dimnames(vcov_ready) <- list(coef(yreg_fit),
                                     coef(yreg_fit))

    } else {

        vcov_ready <- vcov_raw

    }

    if (interaction) {
        vars <- c("(Intercept)", avar, mvar, paste0(avar,":",mvar), cvar)
    } else {
        vars <- c("(Intercept)", avar, mvar,                        cvar)
    }
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
    matrix(((sigma(mreg_fit))^2)^2 / mreg_fit$df.residual)
}
