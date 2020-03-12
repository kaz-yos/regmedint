################################################################################
### Main functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


## VanderWeele 2015 p468 Proposition 2.4
##' Calculate effect measures based on two regression fits (linear/linear)
##'
##' Causal effect parameters are calculated.
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
                                                 cvar,
                                                 interaction) {
    ## mreg coefficients
    beta_hat <- beta_hat(mreg = mreg,
                         mreg_fit = mreg_fit,
                         avar = avar,
                         cvar = cvar)
    beta0 <- beta_hat["(Intercept)"]
    beta1 <- beta_hat[avar]
    beta2 <- beta_hat[cvar]
    ## This mreg linear yreg logistic is the only case that uses sigma^2.
    sigma_sq <- sigma_hat_sq(mreg_fit = mreg_fit)
    ## yreg coefficients
    theta_hat <- theta_hat(yreg = yreg,
                           yreg_fit = yreg_fit,
                           avar = avar,
                           mvar = mvar,
                           cvar = cvar,
                           interaction = interaction)
    theta1 <- theta_hat[avar]
    theta2 <- theta_hat[mvar]
    theta3 <- theta_hat[paste0(avar,":",mvar)]
    theta4 <- theta_hat[cvar]
    ## Construct a function of (a1, a0, m_cde, c_cond) that returns
    ## a vector of point estimates for quantities of interest.
    myreg_est_fun <-
        calc_myreg_mreg_linear_yreg_logistic_est(beta0 = beta0,
                                                 beta1 = beta1,
                                                 beta2 = beta2,
                                                 theta1 = theta1,
                                                 theta2 = theta2,
                                                 theta3 = theta3,
                                                 theta4 = theta4,
                                                 sigma_sq = sigma_sq)

    ## vcovs
    Sigma_beta_hat <- Sigma_beta_hat(mreg = mreg,
                                     mreg_fit = mreg_fit,
                                     avar = avar,
                                     cvar = cvar)
    Sigma_theta_hat <- Sigma_theta_hat(yreg = yreg,
                                       yreg_fit = yreg_fit,
                                       avar = avar,
                                       mvar = mvar,
                                       cvar = cvar,
                                       interaction = interaction)
    Sigma_sigma_hat_sq <- Sigma_sigma_hat_sq(mreg_fit = mreg_fit)
    ## Construct a function of (a0, a1, m_cde, c_cond) that returns
    ## a vector of estimates.
    myreg_se_fun <-
        calc_myreg_mreg_linear_yreg_logistic_se(beta0 = beta0,
                                                beta1 = beta1,
                                                beta2 = beta2,
                                                theta1 = theta1,
                                                theta2 = theta2,
                                                theta3 = theta3,
                                                theta4 = theta4,
                                                sigma_sq = sigma_sq,
                                                Sigma_beta = Sigma_beta,
                                                Sigma_theta = Sigma_theta,
                                                Sigma_sigma = Sigma_sigma)

    ## Return a list of functions.
    list(myreg_est_fun = myreg_est_fun,
         myreg_se_fun = myreg_est_fun)
}


calc_myreg_mreg_linear_yreg_logistic_est <- function(beta0,
                                                     beta1,
                                                     beta2,
                                                     theta1,
                                                     theta2,
                                                     theta3,
                                                     theta4,
                                                     sigma_sq) {
    ## Construct a function for point estimates given (a0, a1, m_cde, c_cond).
    fun_est <- function(a0, a1, m_cde, c_cond) {

        ## Term involving an inner product of beta2 and c_cond
        ## matrix operation to error on non-conformant structure.
        beta2_c <- sum(t(matrix(beta2)) %*% matrix(c_cond))

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

        ## Return a vector
        c(cde = cde,
          pnde = pnde,
          tnie = tnie,
          tnde = tnde,
          pnie = pnie,
          te = te,
          pm = pm)
    }

    return(fun_est)
}


calc_myreg_mreg_linear_yreg_logistic_se <- function(beta0,
                                                    beta1,
                                                    beta2,
                                                    theta1,
                                                    theta2,
                                                    theta3,
                                                    theta4,
                                                    sigma_sq,
                                                    Sigma_beta,
                                                    Sigma_theta,
                                                    Sigma_sigma) {
    ## Construct a function for SE estimates given (a0, a1, m_cde, c_cond)
    fun_se <- function(a0, a1, m_cde, c_cond) {
        return(NULL)
    }

    return(fun_se)
}
