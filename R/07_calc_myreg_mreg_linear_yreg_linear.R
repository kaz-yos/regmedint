################################################################################
### Functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


## VanderWeele 2015 p466 Proposition 2.3
##' Create calculators for effects and se (mreg linear / yreg linear)
##'
##' Construct functions for the conditional effect estimates and their standard errors in the mreg linear / yreg linear setting. Internally, this function deconstructs model objects and feeds parameter estiamtes to the internal worker functions \code{calc_myreg_mreg_linear_yreg_linear_est} and \code{calc_myreg_mreg_linear_yreg_linear_se}.
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return A list containing a function for effect estimates and a function for corresponding standard errors.
calc_myreg_mreg_linear_yreg_linear <- function(mreg,
                                               mreg_fit,
                                               yreg,
                                               yreg_fit,
                                               avar,
                                               mvar,
                                               cvar, 
                                               emm_ac_mreg,
                                               emm_ac_yreg,
                                               emm_mc_yreg,
                                               interaction) {

    ## mreg coefficients
    beta_hat <- beta_hat_helper(mreg = mreg,
                                mreg_fit = mreg_fit,
                                avar = avar,
                                cvar = cvar,
                                emm_ac_mreg = emm_ac_mreg)
    ## yreg coefficients
    theta_hat <- theta_hat_helper(yreg = yreg,
                                  yreg_fit = yreg_fit,
                                  avar = avar,
                                  mvar = mvar,
                                  cvar = cvar,
                                  emm_ac_yreg = emm_ac_yreg,
                                  emm_mc_yreg = emm_mc_yreg,
                                  interaction = interaction)
    ## Construct a function of (a1, a0, m_cde, c_cond) that returns
    ## a vector of point estimates for quantities of interest.
    est_fun <-
        calc_myreg_mreg_linear_yreg_linear_est(beta0 = beta_hat$beta0,
                                               beta1 = beta_hat$beta1,
                                               beta2 = beta_hat$beta2,
                                               beta3 = beta_hat$beta3,
                                               theta0 = theta_hat$theta0,
                                               theta1 = theta_hat$theta1,
                                               theta2 = theta_hat$theta2,
                                               theta3 = theta_hat$theta3,
                                               theta4 = theta_hat$theta4,
                                               theta5 = theta_hat$theta5,
                                               theta6 = theta_hat$theta6)

    ## vcovs
    Sigma_beta_hat <- Sigma_beta_hat(mreg = mreg,
                                     mreg_fit = mreg_fit,
                                     avar = avar,
                                     cvar = cvar,
                                     emm_ac_mreg = emm_ac_mreg)
    Sigma_theta_hat <- Sigma_theta_hat(yreg = yreg,
                                       yreg_fit = yreg_fit,
                                       avar = avar,
                                       mvar = mvar,
                                       cvar = cvar,
                                       emm_ac_yreg = emm_ac_yreg,
                                       emm_mc_yreg = emm_mc_yreg,
                                       interaction = interaction)
    ## Construct a function of (a0, a1, m_cde, c_cond) that returns
    ## a vector of estimates.
    se_fun <-
        calc_myreg_mreg_linear_yreg_linear_se(beta0 = beta_hat$beta0,
                                              beta1 = beta_hat$beta1,
                                              beta2 = beta_hat$beta2,
                                              beta3 = beta_hat$beta3,
                                              theta0 = theta_hat$theta0,
                                              theta1 = theta_hat$theta1,
                                              theta2 = theta_hat$theta2,
                                              theta3 = theta_hat$theta3,
                                              theta4 = theta_hat$theta4,
                                              theta5 = theta_hat$theta5,
                                              theta6 = theta_hat$theta6,
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


calc_myreg_mreg_linear_yreg_linear_est <- function(beta0,
                                                   beta1,
                                                   beta2,
                                                   beta3,
                                                   theta0,
                                                   theta1,
                                                   theta2,
                                                   theta3,
                                                   theta4,
                                                   theta5,
                                                   theta6) {

    validate_myreg_coefs(beta0 = beta0,
                         beta1 = beta1,
                         beta2 = beta2,
                         beta3 = beta3,
                         theta0 = theta0,
                         theta1 = theta1,
                         theta2 = theta2,
                         theta3 = theta3,
                         theta4 = theta4,
                         theta5 = theta5,
                         theta6 = theta6)

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
        
        if (is.null(beta3)) {
            beta3_c <- 0
        } else {
            assertthat::assert_that(length(c_cond) == length(beta3))
            beta3_c <- sum(t(matrix(beta3)) %*% matrix(c_cond))
        }
        
        if (is.null(theta4)) {
            assertthat::assert_that(is.null(c_cond))
            theta4_c <- 0
        } else {
            assertthat::assert_that(!is.null(c_cond))
            assertthat::assert_that(length(c_cond) == length(theta4))
            theta4_c <- sum(t(matrix(theta4)) %*% matrix(c_cond))
        }
        
        if (is.null(theta5)) {
            theta5_c <- 0
        } else {
            assertthat::assert_that(length(c_cond) == length(theta5))
            theta5_c <- sum(t(matrix(theta5)) %*% matrix(c_cond))
        }
        
        if (is.null(theta6)) {
            theta6_c <- 0
        } else {
            assertthat::assert_that(length(c_cond) == length(theta6))
            theta6_c <- sum(t(matrix(theta6)) %*% matrix(c_cond))
        }
        
        ## Extension of VanderWeele 2015 p466
        ## Adopted from mediation.sas and modified.
        ## Look up the third occurrence of the following:
        ## %if &yreg=linear & &mreg=linear & &interaction=true %then %do;
        cde <- (theta1 + theta3*m_cde + theta5_c) * (a1 - a0)
        ## Pearl decomposition (Regular NDE and NIE)
        ## Note the a0 in the first line.                      
        pnde <- (theta1 + theta3*(beta0 + beta1*a0 + beta2_c + beta3_c*a0) + theta5_c) * (a1 - a0)
        ## Note the a1.                               
        tnie <- (theta2 + theta3*a1 + theta6_c) * (beta1 + beta3_c) * (a1 - a0)
        ## Another decomposition
        ## Note the a0 -> a1 change in the first line.         
        tnde <- (theta1 + theta3*(beta0 + beta1*a1 + beta2_c + beta3_c*a1) + theta5_c) * (a1 - a0)
        ## Note the a1 -> a0 change.                  
        pnie <- (theta2 + theta3*a0 + theta6_c) * (beta1 + beta3_c) * (a1 - a0)
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


calc_myreg_mreg_linear_yreg_linear_se <- function(beta0,
                                                  beta1,
                                                  beta2,
                                                  beta3,
                                                  theta0,
                                                  theta1,
                                                  theta2,
                                                  theta3,
                                                  theta4,
                                                  theta5,
                                                  theta6,
                                                  Sigma_beta,
                                                  Sigma_theta) {

    validate_myreg_coefs(beta0 = beta0,
                         beta1 = beta1,
                         beta2 = beta2,
                         beta3 = beta3,
                         theta0 = theta0,
                         theta1 = theta1,
                         theta2 = theta2,
                         theta3 = theta3,
                         theta4 = theta4,
                         theta5 = theta5,
                         theta6 = theta6)

    validate_myreg_vcovs(beta0 = beta0,
                         beta1 = beta1,
                         beta2 = beta2,
                         beta3 = beta3,
                         theta0 = theta0,
                         theta1 = theta1,
                         theta2 = theta2,
                         theta3 = theta3,
                         theta4 = theta4,
                         theta5 = theta5,
                         theta6 = theta6,
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
        
        if (is.null(beta3)) {
            beta3_c <- 0
        } else {
            assertthat::assert_that(length(c_cond) == length(beta3))
            beta3_c <- sum(t(matrix(beta3)) %*% matrix(c_cond))
        }
        
        if (is.null(theta4)) {
            assertthat::assert_that(is.null(c_cond))
            theta4_c <- 0
        } else {
            assertthat::assert_that(!is.null(c_cond))
            assertthat::assert_that(length(c_cond) == length(theta4))
            theta4_c <- sum(t(matrix(theta4)) %*% matrix(c_cond))
        }
        
        if (is.null(theta5)) {
            theta5_c <- 0
        } else {
            assertthat::assert_that(length(c_cond) == length(theta5))
            theta5_c <- sum(t(matrix(theta5)) %*% matrix(c_cond))
        }
        
        if (is.null(theta6)) {
            theta6_c <- 0
        } else {
            assertthat::assert_that(length(c_cond) == length(theta6))
            theta6_c <- sum(t(matrix(theta6)) %*% matrix(c_cond))
        }

        ## Extension of VanderWeele 2015. p468
        ## Valeri & VanderWeele 2013. Appendix p6-9
        ## These are the gradient vector of each scalar quantity of interest.
        ## Obtain the first partial derivative wrt to each parameter.
       if(is.null(theta5)){
           pd_cde_theta5 <- rep(0, length(theta5))
       }else{
           pd_cde_theta5 <- c_cond
       }
        Gamma_cde <-
            matrix(c(0,                       # beta0
                     0,                       # beta1
                     rep(0, length(beta2)),   # beta2 vector
                     rep(0, length(beta3)),   # beta3 vector
                     ##
                     0,                       # theta0
                     1,                       # theta1
                     0,                       # theta2
                     m_cde,                   # theta3
                     rep(0, length(theta4)),  # theta4 vector
                     pd_cde_theta5,           # theta5 vector
                     rep(0, length(theta6))   # theta6 vector
                     )) 
        
        ##
        if(is.null(beta3)){
            pd_pnde_beta3 <- rep(0, length(beta3))
        }else{
            pd_pnde_beta3 <- theta3*c_cond
        }
        if(is.null(theta5)){
            pd_pnde_theta5 <- rep(0, length(theta5))
        }else{
            pd_pnde_theta5 <- c_cond
        }
        Gamma_pnde <-
            matrix(c(
                theta3,                            # beta0
                theta3*a0,                         # beta1
                theta3*c_cond,                     # beta2 vector
                pd_pnde_beta3,                     # beta3 vector
                ##
                0,                                 # theta0
                1,                                 # theta1
                0,                                 # theta2
                beta0 + beta1*a0 + beta2_c + beta3_c*a0,  # theta3
                rep(0, length(theta4)),            # theta4 vector
                pd_pnde_theta5,                    # theta5 vector
                rep(0, length(theta6))             # theta6 vector
                ))  
        
        ##
        if(is.null(beta3)){
            pd_tnie_beta3 <- rep(0, length(beta3))
        }else{
            pd_tnie_beta3 <- c_cond * (theta2 + theta3*a1 + theta6_c)
        }
        if(is.null(theta6)){
            pd_tnie_theta6 <- rep(0, length(theta6))
        }else{
            pd_tnie_theta6 <- c_cond * (beta1 + beta3_c)
        }
        Gamma_tnie <-
            matrix(c(
                0,                         # beta0
                theta2 + theta3*a1 + theta6_c,  # beta1
                rep(0, length(beta2)),     # beta2 vector
                pd_tnie_beta3,             # beta3 vector
                ##
                0,                         # theta0
                0,                         # theta1
                beta1 + beta3_c,           # theta2
                a1 * (beta1 + beta3_c),    # theta3
                rep(0, length(theta4)),    # theta4 vector
                rep(0, length(theta5)),    # theta5 vector
                pd_tnie_theta6             # theta6
                ))  
        
        ##
        if(is.null(beta3)){
            pd_tnde_beta3 <- rep(0, length(beta3))
        }else{
            pd_tnde_beta3 <- theta3*a1*c_cond
        }
        if(is.null(theta5)){
            pd_tnde_theta5 <- rep(0, length(theta5))
        }else{
            pd_tnde_theta5 <- c_cond
        }
        Gamma_tnde <-
            matrix(c(
                theta3,                            # beta0
                theta3*a1,                         # beta1 a0 -> a1
                theta3*c_cond,                     # beta2 vector
                pd_tnde_beta3,                     # beta3 vector
                ##
                0,                                 # theta0
                1,                                 # theta1
                0,                                 # theta2
                beta0 + beta1*a1 + beta2_c + beta3_c*a1,  # theta3 a0 -> a1
                rep(0, length(theta4)),            # theta4 vector
                pd_tnde_theta5,                    # theta5 vector
                rep(0, length(theta6))             # theta6 vector
                ))     
        
        ##
        if(is.null(beta3)){
            pd_pnie_beta3 <- rep(0, length(beta3))
        }else{
            pd_pnie_beta3 <- c_cond * (theta2 + theta3*a0 + theta6_c)
        }
        if(is.null(theta6)){
            pd_pnie_theta6 <- rep(0, length(theta6))
        }else{
            pd_pnie_theta6 <- c_cond * (beta1 + beta3_c)
        }
        Gamma_pnie <-
            matrix(c(
                0,                         # beta0
                theta2 + theta3*a0 + theta6_c,        # beta1 a1 -> a0
                rep(0, length(beta2)),     # beta2 vector
                pd_pnie_beta3,             # beta3 vector
                ##
                0,                         # theta0
                0,                         # theta1
                beta1 + beta3_c,           # theta2
                a0*(beta1 + beta3_c),      # theta3 a1 -> a0
                rep(0, length(theta4)),    # theta4 vector
                rep(0, length(theta5)),    # theta5 vector
                pd_pnie_theta6             # theta6 vector
                ))   
        ##
        Gamma_te <-
            Gamma_pnde + Gamma_tnie # By linearity of differentiation
        ##
        ## PM
        ## Copied from calc_myreg_mreg_linear_yreg_linear_est
        ## Note the a0 in the first line.                      
        pnde <- (theta1 + theta3*(beta0 + beta1*a0 + beta2_c + beta3_c*a0) + theta5_c) * (a1 - a0)
        ## Note the a1.                               
        tnie <- (theta2 + theta3*a1 + theta6_c) * (beta1 + beta3_c) * (a1 - a0)
        ##
        ## Need to unname argument vectors to get c(pnde = , tnie = ).
        d_pm <- grad_prop_med_yreg_linear(pnde = unname(pnde), tnie = unname(tnie))
        ## Multivariate chain rule.
        ## https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Calculus_(OpenStax)/14%3A_Differentiation_of_Functions_of_Several_Variables/14.5%3A_The_Chain_Rule_for_Multivariable_Functions)
        ## d_pm / d_params = (d_pm / d_(pnde, tnie)) %*% (d_(pnde, tnie) / d_params)
        ##                 = (d_pm / d_pnde) * (d_pnde / d_params) +
        ##                   (d_pm / d_tnie) * (d_tnie / d_params)
        ## where (d_pnde / d_params) is (a1 - a0) * Gamma_pnde and
        ##       (d_tnie / d_params) is (a1 - a0) * Gamma_tnie.
        ## Factor out (a1 - a0)
        ## FIXME: This is not tested aginst a reference standard.
        Gamma_pm <-
            (d_pm[["pnde"]] * Gamma_pnde) +
            (d_pm[["tnie"]] * Gamma_tnie)

        ## SE calcuation via multivariate delta method
        ## https://en.wikipedia.org/wiki/Delta_method# Multivariate_delta_method
        a1_sub_a0 <- abs(a1 - a0)
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
