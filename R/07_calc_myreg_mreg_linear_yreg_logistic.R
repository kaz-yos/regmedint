################################################################################
### Main functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


## VanderWeele 2015 p468 Proposition 2.4
##' Create calculators for effects and se (mreg linear / yreg logistic)
##'
##' Construct functions for the conditional effect estimates and their standard errors in the mreg linear / yreg logistic setting. Internally, this function deconstruct model objects and feed parameter estiamtes to the internal worker functions \code{\link{calc_myreg_mreg_linear_yreg_logistic_est}} and \code{\link{calc_myreg_mreg_linear_yreg_logistic_se}}.
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return A list contraining a function for effect estimates and a function for corresponding standard errors.
calc_myreg_mreg_linear_yreg_logistic <- function(mreg,
                                                 mreg_fit,
                                                 yreg,
                                                 yreg_fit,
                                                 avar,
                                                 mvar,
                                                 cvar, # This can be NULL.
                                                 interaction) {

    ## mreg coefficients
    beta_hat <- beta_hat(mreg = mreg,
                         mreg_fit = mreg_fit,
                         avar = avar,
                         cvar = cvar)
    beta0 <- beta_hat["(Intercept)"]
    beta1 <- beta_hat[avar]
    if (is.null(cvar)) {
        ## beta_hat does not contain the cvar part in this case.
        beta2 <- NULL
    } else {
        beta2 <- beta_hat[cvar]
    }
    ## This mreg linear yreg logistic is the only case that uses sigma^2.
    sigma_sq <- sigma_hat_sq(mreg_fit = mreg_fit)
    ## yreg coefficients
    theta_hat <- theta_hat(yreg = yreg,
                           yreg_fit = yreg_fit,
                           avar = avar,
                           mvar = mvar,
                           cvar = cvar,
                           interaction = interaction)
    theta0 <- theta_hat["(Intercept)"]
    theta1 <- theta_hat[avar]
    theta2 <- theta_hat[mvar]
    theta3 <- theta_hat[paste0(avar,":",mvar)]
    if (is.null(cvar)) {
        ## theta_hat does not contain the cvar part in this case.
        theta4 <- NULL
    } else {
        theta4 <- theta_hat[cvar]
    }
    ## Construct a function of (a1, a0, m_cde, c_cond) that returns
    ## a vector of point estimates for quantities of interest.
    myreg_est_fun <-
        calc_myreg_mreg_linear_yreg_logistic_est(beta0 = beta0,
                                                 beta1 = beta1,
                                                 beta2 = beta2,
                                                 theta0 = theta0,
                                                 theta1 = theta1,
                                                 theta2 = theta2,
                                                 theta3 = theta3,
                                                 theta4 = theta4,
                                                 sigma_sq = sigma_sq)

    ## vcovs
    Sigma_beta <- Sigma_beta_hat(mreg = mreg,
                                 mreg_fit = mreg_fit,
                                 avar = avar,
                                 cvar = cvar)
    Sigma_theta <- Sigma_theta_hat(yreg = yreg,
                                   yreg_fit = yreg_fit,
                                   avar = avar,
                                   mvar = mvar,
                                   cvar = cvar,
                                   interaction = interaction)
    Sigma_sigma_sq <- Sigma_sigma_sq_hat(mreg_fit = mreg_fit)
    ## Construct a function of (a0, a1, m_cde, c_cond) that returns
    ## a vector of estimates.
    myreg_se_fun <-
        calc_myreg_mreg_linear_yreg_logistic_se(beta0 = beta0,
                                                beta1 = beta1,
                                                beta2 = beta2,
                                                theta0 = theta0,
                                                theta1 = theta1,
                                                theta2 = theta2,
                                                theta3 = theta3,
                                                theta4 = theta4,
                                                sigma_sq = sigma_sq,
                                                Sigma_beta = Sigma_beta,
                                                Sigma_theta = Sigma_theta,
                                                Sigma_sigma_sq = Sigma_sigma_sq)

    ## Return a list of functions.
    list(
        ## args (a0, a1, m_cde, c_cond)
        ## -> vector c(cde, pnde, tnie, tnde, pnie, te, pm)
        myreg_est_fun = myreg_est_fun,
        ## args (a0, a1, m_cde, c_cond)
        ## -> vector c(se_cde, se_pnde, se_tnie, se_tnde, se_pnie, se_te, se_pm)
        myreg_se_fun = myreg_se_fun)
}


calc_myreg_mreg_linear_yreg_logistic_est <- function(beta0,
                                                     beta1,
                                                     beta2,
                                                     theta0,
                                                     theta1,
                                                     theta2,
                                                     theta3,
                                                     theta4,
                                                     sigma_sq) {

    assertthat::assert_that(length(beta0) == 1,
                            length(beta1) == 1,
                            length(beta2) == length(theta4),
                            length(theta1) == 1,
                            length(theta2) == 1,
                            length(theta3) == 1,
                            length(sigma_sq) == 1)

    ## Construct a function for point estimates given (a0, a1, m_cde, c_cond).
    fun_est <- function(a0, a1, m_cde, c_cond) {

        ## Term involving an inner product of beta2 and c_cond
        ## matrix operation to error on non-conformant structure.
        if (is.null(beta2)) {
            assertthat::assert_that(is.null(c_cond))
            beta2_c <- 0
        } else {
            assertthat::assert_that(!is.null(c_cond))
            assertthat::assert_that(length(c_cond) == length(beta2))
            beta2_c <- sum(t(matrix(beta2)) %*% matrix(c_cond))
        }

        ## VanderWeele 2015 p468
        ## Adopted from mediation.sas and modified.
        ## Look up the third occurrence of the following:
        ## %if &yreg^=linear & &mreg=linear & &interaction=true %then %do;
        cde <- (theta1 + (theta3 * m_cde)) * (a1 - a0)
        ## Pearl decomposition (Regular NDE and NIE)
        ## Note the a0 in the first line.                      vv
        pnde <- (theta1 + (theta3 * beta0) + (theta3 * beta1 * a0) +
                 (theta3 * beta2_c) + (theta3 * theta2 * sigma_sq)) * (a1 - a0) +
            ((1/2) * theta3^2 * sigma_sq) * (a1^2 - a0^2)
        ## Note the a1.                              vv
        tnie <- ((theta2 * beta1) + theta3 * beta1 * a1) * (a1 - a0)
        ## Another decomposition
        ## Note the a0 -> a1 change in the first line.         vv
        tnde <- (theta1 + (theta3 * beta0) + (theta3 * beta1 * a1) +
                 (theta3 * beta2_c) + (theta3 * theta2 * sigma_sq)) * (a1 - a0) +
            ((1/2) * theta3^2 * sigma_sq) * (a1^2 - a0^2)
        ## Note the a1 -> a0 change.                 vv
        pnie <- ((theta2 * beta1) + theta3 * beta1 * a0) * (a1 - a0)
        ## It is the sum of NDE and NIE on the log scale.
        te <- pnde + tnie
        ## VanderWeele 2015 p48.
        pm <- (exp(pnde) * (exp(tnie) - 1)) / (exp(te) - 1)

        ## Return a named vector
        c(cde  = unname(cde),
          pnde = unname(pnde),
          tnie = unname(tnie),
          tnde = unname(tnde),
          pnie = unname(pnie),
          te   = unname(te),
          pm   = unname(pm))
    }

    return(fun_est)
}


calc_myreg_mreg_linear_yreg_logistic_se <- function(beta0,
                                                    beta1,
                                                    beta2,
                                                    theta0,
                                                    theta1,
                                                    theta2,
                                                    theta3,
                                                    theta4,
                                                    sigma_sq,
                                                    Sigma_beta,
                                                    Sigma_theta,
                                                    Sigma_sigma_sq) {

    Sigma <- Matrix::bdiag(Sigma_beta,
                           Sigma_theta,
                           Sigma_sigma_sq)

    ## The dimension
    size_expected <- sum(1, # beta0 (Intercept)
                         1, # beta1 for avar
                         ## This can be 0 = length(NULL) when cvar = NULL
                         length(beta2), # beta2 vector cvar
                         ##
                         1, # theta0 (Intercept). Never used so not in args.
                         1, # theta1 for avar
                         1, # theta2 for mvar
                         1, # theta3 for avar:mvar
                         ## This can be 0 = length(NULL) when cvar = NULL
                         length(theta4), # theta4 for cvar
                         ##
                         1) # sigma_sq
    assertthat::assert_that(dim(Sigma)[1] == size_expected)
    assertthat::assert_that(dim(Sigma)[2] == size_expected)

    ## Construct a function for SE estimates given (a0, a1, m_cde, c_cond)
    fun_se <- function(a0, a1, m_cde, c_cond) {

        ## Term involving an inner product of beta2 and c_cond
        ## matrix operation to error on non-conformant structure.
        if (is.null(beta2)) {
            assertthat::assert_that(is.null(c_cond))
            beta2_c <- 0
        } else {
            assertthat::assert_that(!is.null(c_cond))
            assertthat::assert_that(length(c_cond) == length(beta2))
            beta2_c <- sum(t(matrix(beta2)) %*% matrix(c_cond))
        }

        ## VanderWeele 2015. p468
        ## Valeri & VanderWeele 2013. Appendix p6-9
        ## These are the gradient vector of each scalar quantity of interest.
        ## Obtain the first partial derivative wrt to each parameter.
        Gamma_cde <-
            matrix(c(0,                       # beta0
                     0,                       # beta1
                     rep(0, length(beta2)),   # beta2 vector
                     ##
                     0,                       # theta0
                     1,                       # theta1
                     0,                       # theta2
                     m_cde,                   # theta3
                     rep(0, length(theta4)),  # theta4 vector
                     ##
                     0))                      # sigma^2
        ##
        Gamma_pnde <-
            matrix(c(
                theta3,                                             # beta0
                (theta3 * a0),                                      # beta1
                (theta3 * c_cond),                                  # beta2 vector
                ##
                0,                                                  # theta0
                1,                                                  # theta1
                (theta3 * sigma_sq),                                # theta2
                (beta0 + beta1 * a0 + beta2_c + theta2 *
                 sigma_sq + theta3 * sigma_sq * (a1 + a0)),         # theta3
                rep(0, length(theta4)),                             # theta4 vector
                ##
                (theta3 * theta2 + (1/2) * theta3^2 * (a1 + a0))))  # sigma^2
        ##
        Gamma_tnie <-
            matrix(c(
                0,                       # beta0
                (theta2 + theta3 * a1),  # beta1
                rep(0, length(beta2)),   # beta2 vector
                ##
                0,                       # theta0
                0,                       # theta1
                beta1,                   # theta2
                (beta1 * a1),            # theta3
                rep(0, length(theta4)),  # theta4 vector
                ##
                0))                      # sigma^2
        ##
        Gamma_tnde <-
            matrix(c(
                theta3,                                             # beta0
                (theta3 * a1),                                      # beta1 a0 -> a1
                (theta3 * c_cond),                                  # beta2 vector
                ##
                0,                                                  # theta0
                1,                                                  # theta1
                (theta3 * sigma_sq),                                # theta2
                (beta0 + beta1 * a1 + beta2_c + theta2 *
                 sigma_sq + theta3 * sigma_sq * (a1 + a0)),         # theta3 a0 -> a1
                rep(0, length(theta4)),                             # theta4 vector
                ##
                (theta3 * theta2 + (1/2) * theta3^2 * (a1 + a0))))  # sigma^2
        ##
        Gamma_pnie <-
            matrix(c(
                0,                       # beta0
                (theta2 + theta3 * a0),  # beta1 a1 -> a0
                rep(0, length(beta2)),   # beta2 vector
                ##
                0,                       # theta0
                0,                       # theta1
                beta1,                   # theta2
                (beta1 * a0),            # theta3 a1 -> a0
                rep(0, length(theta4)),  # theta4 vector
                ##
                0))                      # sigma^2
        ## Note VV2013 Appendix p9 seems to have an error in the gradient
        ## vector element correspoding to sigma^2. Otherwise, it is a
        ## sum of gradients for pnde and tnie.
        Gamma_te <-
            Gamma_pnde + Gamma_tnie # By linearity of differentiation
        ##
        ## Not implemented in mediation.sas.
        ## Not mentioned in VV2013, VV2015, or VanderWeele 2015.
        Gamma_pm <- matrix(rep(as.numeric(NA), nrow(Sigma)))


        ## SEs
        a1_sub_a0 <- (a1 - a0)
        se_cde  <- sqrt(as.numeric(t(Gamma_cde)  %*% Sigma %*% Gamma_cde)) * a1_sub_a0
        se_pnde <- sqrt(as.numeric(t(Gamma_pnde) %*% Sigma %*% Gamma_pnde)) * a1_sub_a0
        se_tnie <- sqrt(as.numeric(t(Gamma_tnie) %*% Sigma %*% Gamma_tnie)) * a1_sub_a0
        se_tnde <- sqrt(as.numeric(t(Gamma_tnde) %*% Sigma %*% Gamma_tnde)) * a1_sub_a0
        se_pnie <- sqrt(as.numeric(t(Gamma_pnie) %*% Sigma %*% Gamma_pnie)) * a1_sub_a0
        se_te   <- sqrt(as.numeric(t(Gamma_te)   %*% Sigma %*% Gamma_te)) * a1_sub_a0
        se_pm   <- sqrt(as.numeric(t(Gamma_pm)   %*% Sigma %*% Gamma_pm)) * a1_sub_a0

        ## Return a vector
        c(se_cde  = unname(se_cde),
          se_pnde = unname(se_pnde),
          se_tnie = unname(se_tnie),
          se_tnde = unname(se_tnde),
          se_pnie = unname(se_pnie),
          se_te   = unname(se_te),
          se_pm   = unname(se_pm))
    }

    return(fun_se)
}
