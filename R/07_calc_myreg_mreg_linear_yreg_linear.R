################################################################################
### Functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


## VanderWeele 2015 p466 Proposition 2.3
##' Create calculators for effects and se (mreg linear / yreg linear)
##'
##' Construct functions for the conditional effect estimates and their standard errors in the mreg linear / yreg linear setting.
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return A list contraining a function for effect estimates and a function for corresponding standard errors.
calc_myreg_mreg_linear_yreg_linear <- function(mreg,
                                               mreg_fit,
                                               yreg,
                                               yreg_fit,
                                               avar,
                                               mvar,
                                               cvar, # This can be NULL.
                                               interaction) {

    ## FIXME: This needs to be allowed and handled nicely!
    assertthat::assert_that(!is.null(cvar))

    ## mreg coefficients
    beta_hat <- beta_hat(mreg = mreg,
                         mreg_fit = mreg_fit,
                         avar = avar,
                         cvar = cvar)
    beta0 <- beta_hat["(Intercept)"]
    beta1 <- beta_hat[avar]
    beta2 <- beta_hat[cvar] # This can be length 0.
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
        calc_myreg_mreg_linear_yreg_linear_est(beta0 = beta0,
                                               beta1 = beta1,
                                               beta2 = beta2,
                                               theta1 = theta1,
                                               theta2 = theta2,
                                               theta3 = theta3,
                                               theta4 = theta4)

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
    ## Construct a function of (a0, a1, m_cde, c_cond) that returns
    ## a vector of estimates.
    myreg_se_fun <-
        calc_myreg_mreg_linear_yreg_linear_se(beta0 = beta0,
                                              beta1 = beta1,
                                              beta2 = beta2,
                                              theta1 = theta1,
                                              theta2 = theta2,
                                              theta3 = theta3,
                                              theta4 = theta4,
                                              Sigma_beta = Sigma_beta_hat,
                                              Sigma_theta = Sigma_theta_hat)

    ## Return a list of functions.
    list(myreg_est_fun = myreg_est_fun,
         myreg_se_fun = myreg_est_fun)
}


calc_myreg_mreg_linear_yreg_linear_est <- function(beta0,
                                                   beta1,
                                                   beta2,
                                                   theta1,
                                                   theta2,
                                                   theta3,
                                                   theta4) {
    ## Construct a function for point estimates given (a0, a1, m_cde, c_cond).
    fun_est <- function(a0, a1, m_cde, c_cond) {

        ## Term involving an inner product of beta2 and c_cond
        ## matrix operation to error on non-conformant structure.
        beta2_c <- sum(t(matrix(beta2)) %*% matrix(c_cond))

        ## VanderWeele 2015 p466
        ## Adopted from mediation.sas and modified.
        ## Look up the third occurrence of the following:
        ## %if &yreg^=linear & &mreg=linear & &interaction=true %then %do;
        cde <- (theta1 + (theta3 * m_cde)) * (a1 - a0)
        ## Pearl decomposition (Regular NDE and NIE)
        ## Note the a0 in the first line.                      vv
        pnde <- (theta1 + (theta3 * beta0) + (theta3 * beta1 * a0) +
                 (theta3 * beta2_c)) * (a1 - a0)
        ## Note the a1.                               vv
        tnie <- ((theta2 * beta1) + (theta3 * beta1 * a1)) * (a1 - a0)
        ## Another decomposition
        ## Note the a0 -> a1 change in the first line.         vv
        tnde <- (theta1 + (theta3 * beta0) + (theta3 * beta1 * a1) +
                 (theta3 * beta2_c)) * (a1 - a0)
        ## Note the a1 -> a0 change.                  vv
        pnie <- ((theta2 * beta1) + (theta3 * beta1 * a0)) * (a1 - a0)
        ## It is the sum of NDE and NIE on the log scale.
        te <- pnde + tnie
        ## VanderWeele 2015 p47
        pm <- tnie / te

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


calc_myreg_mreg_linear_yreg_linear_se <- function(beta0,
                                                  beta1,
                                                  beta2,
                                                  theta1,
                                                  theta2,
                                                  theta3,
                                                  theta4,
                                                  Sigma_beta,
                                                  Sigma_theta) {

    Sigma <- Matrix::bdiag(Sigma_beta,
                           Sigma_theta)

    ## Construct a function for SE estimates given (a0, a1, m_cde, c_cond)
    fun_se <- function(a0, a1, m_cde, c_cond) {

        ## Term involving an inner product of beta2 and c_cond
        ## matrix operation to error on non-conformant structure.
        beta2_c <- sum(t(matrix(beta2)) %*% matrix(c_cond))

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
                     rep(0, length(theta4)))) # theta4 vector
        ##
        Gamma_pnde <-
            matrix(c(
                theta3,                            # beta0
                (theta3 * a0),                     # beta1
                (theta3 * c_cond),                 # beta2 vector
                ##
                0,                                 # theta0
                1,                                 # theta1
                0,                                 # theta2
                (beta0 + (beta1 * a0) + beta2_c),  # theta3
                rep(0, length(theta4))))           # theta4 vector
        ##
        Gamma_tnie <-
            matrix(c(
                0,                         # beta0
                (theta2 + (theta3 * a1)),  # beta1
                rep(0, length(beta2)),     # beta2 vector
                ##
                0,                         # theta0
                0,                         # theta1
                beta1,                     # theta2
                (beta1 * a1),              # theta3
                rep(0, length(theta4))))   # theta4 vector
        ##
        Gamma_tnde <-
            matrix(c(
                theta3,                            # beta0
                (theta3 * a1),                     # beta1 a0 -> a1
                (theta3 * c_cond),                 # beta2 vector
                ##
                0,                                 # theta0
                1,                                 # theta1
                (theta3 * sigma_sq),               # theta2
                (beta0 + (beta1 * a1) + beta2_c),  # theta3 a0 -> a1
                rep(0, length(theta4))))           # theta4 vector
        ##
        Gamma_pnie <-
            matrix(c(
                0,                         # beta0
                (theta2 + (theta3 * a0)),  # beta1 a1 -> a0
                rep(0, length(beta2)),     # beta2 vector
                ##
                0,                         # theta0
                0,                         # theta1
                beta1,                     # theta2
                (beta1 * a0),              # theta3 a1 -> a0
                rep(0, length(theta4))))   # theta4 vector
        ##
        Gamma_te <-
            Gamma_pnde + Gamma_tnie # By linearity of differentiation
        ##
        ## PM part. VV2013 Appendix p5-6.
        d1_numer <- -theta3 * ((theta2 * beta1) + (theta3 * beta1 * a1))
        d1_denom_sqrt <- (theta1 + (theta3 * beta0) + (theta3 * beta1 * a0) + (theta3 * beta2_c) + (theta2 * beta1) + (theta3 * beta1 * a1))
        ##
        d1 <- d1_numer / d1_denom_sqrt^2
        ##
        d2 <- ((theta2 + (theta3 * a1) * (-1 * ((theta2 * beta1) + (theta3 * best1 * a1)) + d1_denom_sqrt)) - (theta3 * a0)) / d1_denom_sqrt^2
        ## Vector valued
        d3 <- c_cond * (d1_numer / d1_denom_sqrt^2)
        ##
        d4 <- 0
        ##
        d5 <- ((theta2 * beta1) + (theta3 * beta1 * a1)) / d1_denom_sqrt^2
        ##
        d6 <- beta1 * (-1 * ((theta2 * beta1) + (theta3 * beta1 * a1)) + d1_denom_sqrt) / d1_denom_sqrt^2
        ##
        d7 <- ((beta1 * a1 * d1_denom_sqrt) - ((beta0 + (beta1 * (a1 + a0)) + beta2_c) * ((theta2 * beta1) + (theta3 * beta1 * a1)))) / d1_denom_sqrt^2
        ##
        d8 <- rep(0, length(theta4))
        ##
        Gamma_pm <- c(d1, # beta0
                      d2, # beta1
                      d3, # beta2 vector
                      ##
                      d4, # theta0
                      d5, # theta1
                      d6, # theta2
                      d7, # theta3
                      d8) # theta4 vector

        ## SE calcuation via multivariate delta method
        ## https://en.wikipedia.org/wiki/Delta_method# Multivariate_delta_method
        a1_sub_a0 <- (a1 - a0)
        se_cde <- sqrt(as.numeric(t(Gamma_cde) %*% Sigma %*% Gamma_cde)) * a1_sub_a0
        se_pnde <- sqrt(as.numeric(t(Gamma_pnde) %*% Sigma %*% Gamma_pnde)) * a1_sub_a0
        se_tnie <- sqrt(as.numeric(t(Gamma_tnie) %*% Sigma %*% Gamma_tnie)) * a1_sub_a0
        se_tnde <- sqrt(as.numeric(t(Gamma_tnde) %*% Sigma %*% Gamma_tnde)) * a1_sub_a0
        se_pnie <- sqrt(as.numeric(t(Gamma_pnie) %*% Sigma %*% Gamma_pnie)) * a1_sub_a0
        se_te <- sqrt(as.numeric(t(Gamma_te) %*% Sigma %*% Gamma_te)) * a1_sub_a0
        se_pm <- sqrt(as.numeric(t(Gamma_pm) %*% Sigma %*% Gamma_pm)) * a1_sub_a0

        ## Return a vector
        c(se_cde = se_cde,
          se_pnde = pnde,
          se_tnie = tnie,
          se_tnde = tnde,
          se_pnie = pnie,
          se_te = te,
          se_pm = pm)
    }

    return(fun_se)
}
