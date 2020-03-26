################################################################################
### Main functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


## VanderWeele 2015 p471 Proposition 2.5
##' Create calculators for effects and se (mreg linear / yreg linear)
##'
##' Construct functions for the conditional effect estimates and their standard errors in the mreg linear / yreg linear setting.
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return A list contraining a function for effect estimates and a function for corresponding standard errors.
calc_myreg_mreg_logistic_yreg_linear <- function(mreg,
                                                 mreg_fit,
                                                 yreg,
                                                 yreg_fit,
                                                 avar,
                                                 mvar,
                                                 cvar, # This can be NULL.
                                                 interaction) {

    ## mreg coefficients
    beta_hat <- beta_hat_helper(mreg = mreg,
                                mreg_fit = mreg_fit,
                                avar = avar,
                                cvar = cvar)
    ## yreg coefficients
    theta_hat <- theta_hat_helper(yreg = yreg,
                                  yreg_fit = yreg_fit,
                                  avar = avar,
                                  mvar = mvar,
                                  cvar = cvar,
                                  interaction = interaction)
    ## Construct a function of (a1, a0, m_cde, c_cond) that returns
    ## a vector of point estimates for quantities of interest.
    est_fun <-
        calc_myreg_mreg_logistic_yreg_linear_est(beta0 = beta_hat$beta0,
                                                 beta1 = beta_hat$beta1,
                                                 beta2 = beta_hat$beta2,
                                                 theta0 = theta_hat$theta0,
                                                 theta1 = theta_hat$theta1,
                                                 theta2 = theta_hat$theta2,
                                                 theta3 = theta_hat$theta3,
                                                 theta4 = theta_hat$theta4)

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
    se_fun <-
        calc_myreg_mreg_logistic_yreg_linear_se(beta0 = beta_hat$beta0,
                                                beta1 = beta_hat$beta1,
                                                beta2 = beta_hat$beta2,
                                                theta0 = theta_hat$theta0,
                                                theta1 = theta_hat$theta1,
                                                theta2 = theta_hat$theta2,
                                                theta3 = theta_hat$theta3,
                                                theta4 = theta_hat$theta4,
                                                Sigma_beta = Sigma_beta_hat,
                                                Sigma_theta = Sigma_theta_hat)

    ## Return a list of functions.
    list(
        ## args (a0, a1, m_cde, c_cond)
        ## -> vector c(cde, pnde, tnie, tnde, pnie, te, pm)
        est_fun = est_fun,
        ## args (a0, a1, m_cde, c_cond)
        ## -> vector c(se_cde, se_pnde, se_tnie, se_tnde, se_pnie, se_te, se_pm)
        se_fun = se_fun)
}


calc_myreg_mreg_logistic_yreg_linear_est <- function(beta0,
                                                     beta1,
                                                     beta2,
                                                     theta0,
                                                     theta1,
                                                     theta2,
                                                     theta3,
                                                     theta4) {

    validate_myreg_coefs(beta0 = beta0,
                         beta1 = beta1,
                         beta2 = beta2,
                         theta0 = theta0,
                         theta1 = theta1,
                         theta2 = theta2,
                         theta3 = theta3,
                         theta4 = theta4)

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

        ## VanderWeele 2015 p471
        cde <- (theta1 + (theta3 * m_cde)) * (a1 - a0)
        ## Pearl decomposition (Regular NDE and NIE)
        ## Note the a0 in the second term.
        pnde <- (theta1 * (a1 - a0)) + (theta3 * (a1 - a0)) *
            (exp(beta0 + (beta1 * a0) + beta2_c) /
             (1 + exp(beta0 + (beta1 * a0) + beta2_c)))
        ## Note the a1 in the first term.
        tnie <- (theta2 + (theta3 * a1)) *
            ((exp(beta0 + (beta1 * a1) + beta2_c) /
              (1 + exp(beta0 + (beta1 * a1) + beta2_c)))
                - (exp(beta0 + (beta1 * a0) + beta2_c) /
                   (1 + exp(beta0 + (beta1 * a0) + beta2_c))))
        ## Another decomposition
        ## Note the a0 -> a1 change in the second term.
        tnde <- (theta1 * (a1 - a0)) + (theta3 * (a1 - a0)) *
            (exp(beta0 + (beta1 * a1) + beta2_c) /
             (1 + exp(beta0 + (beta1 * a1) + beta2_c)))
        ## Note the a1 -> a0 change in the first term.
        pnie <- (theta2 + (theta3 * a0)) *
            ((exp(beta0 + (beta1 * a0) + beta2_c) /
              (1 + exp(beta0 + (beta1 * a1) + beta2_c)))
                - (exp(beta0 + (beta1 * a0) + beta2_c) /
                   (1 + exp(beta0 + (beta1 * a0) + beta2_c))))
        ## It is the sum of NDE and NIE on the log scale.
        te <- pnde + tnie
        ## VanderWeele 2015 p47
        pm <- tnie / te

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


calc_myreg_mreg_logistic_yreg_linear_se <- function(beta0,
                                                    beta1,
                                                    beta2,
                                                    theta0,
                                                    theta1,
                                                    theta2,
                                                    theta3,
                                                    theta4,
                                                    Sigma_beta,
                                                    Sigma_theta) {

    validate_myreg_coefs(beta0 = beta0,
                         beta1 = beta1,
                         beta2 = beta2,
                         theta0 = theta0,
                         theta1 = theta1,
                         theta2 = theta2,
                         theta3 = theta3,
                         theta4 = theta4)

    validate_myreg_vcovs(beta0 = beta0,
                         beta1 = beta1,
                         beta2 = beta2,
                         theta0 = theta0,
                         theta1 = theta1,
                         theta2 = theta2,
                         theta3 = theta3,
                         theta4 = theta4,
                         Sigma_beta = Sigma_beta,
                         Sigma_theta = Sigma_theta)

    Sigma <- Matrix::bdiag(Sigma_beta,
                           Sigma_theta)

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
                     rep(0, length(theta4)))) # theta4 vector
        ##
        pnde_d1 <- theta3 * (exp(beta0 + (beta1 * a0) + beta2_c) / (1 + exp(beta0 + (beta1 * a0) + beta2_c))^2)
        ## d2 and d3 in VanderWeele 2015 p471 and VV 2013 Appendix p12 have typos.
        ## a0 and c_cond and theta3 are outside the fraction.
        pnde_d2 <- a0 * pnde_d1
        pnde_d3 <- c_cond * pnde_d1
        pnde_d4 <- 0
        pnde_d5 <- 1 # (a1 - a0) is factored out.
        pnde_d6 <- 0
        pnde_d7 <- exp(beta0 + (beta1 * a0) + beta2_c) /
            (1 + exp(beta0 + (beta1 * a0) + beta2_c)) #
        pnde_d8 <- rep(0, length(theta4))
        Gamma_pnde <-
            matrix(c(
                pnde_d1,   # beta0
                pnde_d2,   # beta1
                pnde_d3,   # beta2 vector
                ##
                pnde_d4,   # theta0
                pnde_d5,   # theta1
                pnde_d6,   # theta2
                pnde_d7,   # theta3
                pnde_d8))  # theta4 vector
        ##
        tnie_A <- NA
        tnie_B <- NA
        tnie_K <- NA
        tnie_D <- NA
        tnie_d1 <- (theta2 + (theta3 * a1)) * (tnie_A - tnie_B)
        tnie_d2 <- (theta2 + (theta3 * a1)) * ((a1 * tnie_A) - (a0 * tnie_B))
        tnie_d3 <- (theta2 + (theta3 * a1)) * c_cond * (tnie_A - tnie_B)
        tnie_d4 <- 0
        tnie_d5 <- 0
        tnie_d6 <- tnie_K - tnie_D
        tnie_d7 <- a1 * (tnie_K - tnie_D)
        tnie_d8 <- rep(0, length(theta4))
        Gamma_tnie <-
            matrix(c(
                tnie_d1,   # beta0
                tnie_d2,   # beta1
                tnie_d3,   # beta2 vector
                ##
                tnie_d4,   # theta0
                tnie_d5,   # theta1
                tnie_d6,   # theta2
                tnie_d7,   # theta3
                tnie_d8))  # theta4 vector
        ##
        Gamma_tnde <-
            matrix(c(
                theta3,                            # beta0
                (theta3 * a1),                     # beta1 a0 -> a1
                (theta3 * c_cond),                 # beta2 vector
                ##
                0,                                 # theta0
                1,                                 # theta1
                0,                                 # theta2
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
        d2 <- ((theta2 + (theta3 * a1) * (-1 * ((theta2 * beta1) + (theta3 * beta1 * a1)) + d1_denom_sqrt)) - (theta3 * a0)) / d1_denom_sqrt^2
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
