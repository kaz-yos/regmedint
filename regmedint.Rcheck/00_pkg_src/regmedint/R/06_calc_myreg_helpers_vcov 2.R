################################################################################
### Helper functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


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


    if (interaction) {

        ## Interaction case
        ## No data manipulation is necessary.
        ## Technically, the first case can be used in both because NULL
        ## drops out in c(..., NULL). Here it is made explicit.
        if (!is.null(cvar)) {
            vars <- c("(Intercept)", avar, mvar, paste0(avar,":",mvar), cvar)
        } else {
            vars <- c("(Intercept)", avar, mvar, paste0(avar,":",mvar))
        }

    } else {

        if (!is.null(cvar)) {
            vars <- c("(Intercept)", avar, mvar, paste0(avar,":",mvar), cvar)
            ## Always have a position for an interaction term
            ## to ease subsequent manipulation.
            vcov_ready <-
                Matrix::bdiag(vcov_ready[c("(Intercept)",avar,mvar),
                                         c("(Intercept)",avar,mvar),
                                         drop = FALSE],
                              ## Padding for the avar:mvar position.
                              matrix(0),
                              vcov_ready[cvar,
                                         cvar,
                                         drop = FALSE])
            dimnames(vcov_ready) <- list(vars,
                                         vars)
        } else {
            vars <- c("(Intercept)", avar, mvar, paste0(avar,":",mvar))
            ## Always have a position for an interaction term
            ## to ease subsequent manipulation.
            vcov_ready <-
                Matrix::bdiag(vcov_ready[c("(Intercept)",avar,mvar),
                                         c("(Intercept)",avar,mvar),
                                         drop = FALSE],
                              ## Padding for the avar:mvar position.
                              matrix(0))
            dimnames(vcov_ready) <- list(vars,
                                         vars)
        }

    }

    vcov_ready[vars,vars, drop = FALSE]
}

Sigma_sigma_sq_hat <- function(mreg_fit) {

    ## VanderWeele 2015. p470 states
    ## 2 * (sigma_hat^2)^2 / (n-p)
    ##
    ## Chapter 2 Quadratic Forms of Random Variables
    ## http://pages.stat.wisc.edu/~st849-1/lectures/Ch02.pdf
    ## Corollary 6. In a full-rank Gaussian model for M with p covariates X.
    ## (||M - X beta_hat||^2 / sigma^2) ~ (central chi-squared DF = (n - p))
    ##
    ## https://en.wikipedia.org/wiki/Chi-squared_distribution
    ## Var(central chi-squared with DF = (n - p)) = 2 * (n - p)
    ##
    ## Var(sigma_hat^2) = Var(1/(n-p) * ||M - X beta_hat||^2)
    ##                  = 1/(n-p)^2 * Var(||M - X beta_hat||^2 * sigma^2/sigma^2)
    ##                  = 1/(n-p)^2 * (sigma^2)^2 * Var(||M - X beta_hat||^2 /sigma^2)
    ##                  = 1/(n-p)^2 * (sigma^2)^2 * (2 * (n-p))
    ##                  = 2 * (sigma^2)^2 / (n-p)
    ##
    matrix(2 * ((sigma(mreg_fit))^2)^2 / mreg_fit$df.residual)
}
