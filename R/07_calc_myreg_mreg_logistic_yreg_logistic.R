################################################################################
### Main functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


## Extension of VanderWeele 2015 p473 Proposition 2.6
##' Create calculators for effects and se (mreg logistic / yreg logistic)
##'
##' Construct functions for the conditional effect estimates and their standard errors in the mreg logistic / yreg logistic setting. Internally, this function deconstructs model objects and feeds parameter estimates to the internal worker functions \code{calc_myreg_mreg_logistic_yreg_logistic_est} and \code{calc_myreg_mreg_logistic_yreg_logistic_se}.
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return A list containing a function for effect estimates and a function for corresponding standard errors.

calc_myreg_mreg_logistic_yreg_logistic <- function(mreg,
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
        calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = beta_hat$beta0,
                                                   beta1 = beta_hat$beta1,
                                                   beta2 = beta_hat$beta2,
                                                   beta3 = beta_hat$beta3,
                                                   theta0 = theta_hat$theta0,
                                                   theta1 = theta_hat$theta1,
                                                   theta2 = theta_hat$theta2,
                                                   theta3 = theta_hat$theta3,
                                                   theta4 = theta_hat$theta4,
                                                   theta5 = theta_hat$theta5,
                                                   theta6 = theta_hat$theta6
                                                   )

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
        calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = beta_hat$beta0,
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


calc_myreg_mreg_logistic_yreg_logistic_est <- function(beta0,
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
        
        expit <- function(x){exp(x)/(1+exp(x))}

        ## Extension of VanderWeele 2015 p473. Estimates are all on log OR scale.
        ## Note that only cde is presented on the log scale
        cde <- (theta1 + theta3*m_cde + theta5_c) * (a1 - a0)
        ## Pearl decomposition (Regular NDE and NIE)
        ## Note the a0 in the second term.
        pnde <- ((theta1 + theta5_c) * (a1 - a0)) +
            log(1 + exp(beta0 + beta1*a0 + beta2_c + beta3_c*a0 + theta2 + theta3*a1 + theta6_c)) -
            log(1 + exp(beta0 + beta1*a0 + beta2_c + beta3_c*a0 + theta2 + theta3*a0 + theta6_c))
        ## Note the a1 in the first term.
        tnie <- log(1 + exp(beta0 + beta1*a1 + beta2_c + beta3_c*a1 + theta2 + theta3*a1 + theta6_c)) - 
            log(1 + exp(beta0 + beta1*a0 + beta2_c + beta3_c*a0 + theta2 + theta3*a1 + theta6_c)) +
            log(1 + exp(beta0 + beta1*a0 + beta2_c + beta3_c*a0)) - 
            log(1 + exp(beta0 + beta1*a1 + beta2_c + beta3_c*a1))
        ## Another decomposition
        ## Note the a0 -> a1 changes associated with beta1.
        tnde <- ((theta1 + theta5_c) * (a1 - a0)) +
            log(1 + exp(beta0 + beta1*a1 + beta2_c + beta3_c*a1 + theta2 + theta3*a1 + theta6_c)) -
            log(1 + exp(beta0 + beta1*a1 + beta2_c + beta3_c*a1 + theta2 + theta3*a0 + theta6_c))
        ## Note the a1 -> a0 changes associated with theta3.
        pnie <-
            log(1 + exp(beta0 + beta1*a1 + beta2_c + beta3_c*a1 + theta2 + theta3*a0 + theta6_c)) - 
            log(1 + exp(beta0 + beta1*a0 + beta2_c + beta3_c*a0 + theta2 + theta3*a0 + theta6_c)) +
            log(1 + exp(beta0 + beta1*a0 + beta2_c + beta3_c*a0)) - 
            log(1 + exp(beta0 + beta1*a1 + beta2_c + beta3_c*a1))
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


calc_myreg_mreg_logistic_yreg_logistic_se <- function(beta0,
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
        
        expit <- function(x){exp(x)/(1+exp(x))}

        ## Extension of VanderWeele 2015. p473
        ## Note that Gammas are all for log OR effects.
        ## Valeri & VanderWeele 2013. Appendix p6-9
        ## These are the gradient vector of each scalar quantity of interest.
        ## Obtain the first partial derivative wrt to each parameter.
        if(is.null(theta5)){
            pd_cde_theta5 <- rep(0, length(theta5))
        }else{
            pd_cde_theta5 <- c_cond
        }
        Gamma_cde <-
            matrix(c(0,                        # beta0
                     0,                        # beta1
                     rep(0, length(beta2)),    # beta2 vector
                     rep(0, length(beta3)),    # beta3 vector
                     ##
                     0,                        # theta0
                     (a1 - a0),                # theta1
                     0,                        # theta2
                     (a1 - a0) * m_cde,        # theta3
                     rep(0, length(theta4)),   # theta4 vector
                     pd_cde_theta5,            # theta5 vector
                     rep(0, length(theta6))    # theta6 vector
                     ))  
        ##
        pnde_expit_a1 <- expit(a0*(beta1 + beta3_c) + a1*theta3 + beta0 + beta2_c + theta6_c + theta2)
        pnde_expit_a0 <- expit(a0*(beta1 + beta3_c) + a0*theta3 + beta0 + beta2_c + theta6_c + theta2)
        
        pnde_d1 <- pnde_expit_a1 - pnde_expit_a0
        pnde_d2 <- a0*pnde_d1
        pnde_d3 <- c_cond*pnde_d1
        # pnde_d4 <- c_cond*a0*pnde_d1
        if(is.null(beta3)){
            pnde_d4 <- rep(0, length(beta3))
        }else{
            pnde_d4 <- c_cond*a0*pnde_d1
        }
        pnde_d5 <- 0
        pnde_d6 <- a1 - a0
        pnde_d7 <- pnde_d1
        pnde_d8 <- a1*pnde_expit_a1 - a0*pnde_expit_a0
        pnde_d9 <- rep(0, length(theta4))
        # pnde_d10 <- c_cond * (a1 - a0)
        if(is.null(theta5)){
            pnde_d10 <- rep(0, length(theta5))
        }else{
            pnde_d10 <- c_cond*(a1 - a0)
        }
        # pnde_d11 <- c_cond*pnde_d1
        if(is.null(theta6)){
            pnde_d11 <- rep(0, length(theta6))
        }else{
            pnde_d11 <- c_cond*pnde_d1
        }
        
        Gamma_pnde <-
            matrix(c(
                pnde_d1,   # beta0
                pnde_d2,   # beta1
                pnde_d3,   # beta2 vector
                pnde_d4,   # beta3 vector
                ##
                pnde_d5,   # theta0
                pnde_d6,   # theta1
                pnde_d7,   # theta2
                pnde_d8,   # theta3
                pnde_d9,   # theta4 vector
                pnde_d10,  # theta5 vector
                pnde_d11   # theta6 vector
                ))  
        ##
        tnie_expit_q1 <- expit(a0*(beta1 + beta3_c) + beta0 + beta2_c)
        tnie_expit_q2 <- expit(a0*(beta1 + beta3_c) + a1*theta3 + beta0 + beta2_c + theta6_c + theta2)
        tnie_expit_q3 <- expit(a1*(beta1 + beta3_c) + beta0 + beta2_c)
        tnie_expit_q4 <- expit(a1*(beta1 + beta3_c) + a1*theta3 + beta0 + beta2_c + theta6_c + theta2)
        
        tnie_d1 <- tnie_expit_q1 - tnie_expit_q2 - tnie_expit_q3 + tnie_expit_q4
        tnie_d2 <- a0 * (tnie_expit_q1-tnie_expit_q2) + a1 * (-tnie_expit_q3+tnie_expit_q4)
        tnie_d3 <- c_cond*tnie_d1
        # tnie_d4 <- c_cond*tnie_d2
        if(is.null(beta3)){
            tnie_d4 <- rep(0, length(beta3))
        }else{
            tnie_d4 <- c_cond*tnie_d2
        }
        tnie_d5 <- 0
        tnie_d6 <- 0
        tnie_d7 <- tnie_expit_q4 - tnie_expit_q2  # watch out the order
        tnie_d8 <- a1*tnie_d7
        tnie_d9 <- rep(0, length(theta4))
        tnie_d10 <- rep(0, length(theta5))
        # tnie_d11 <- c_cond*tnie_d7 
        if(is.null(theta6)){
            tnie_d11 <- rep(0, length(theta6))
        }else{
            tnie_d11 <- c_cond*tnie_d7
        }
        
        Gamma_tnie <-
            matrix(c(
                tnie_d1,   # beta0
                tnie_d2,   # beta1
                tnie_d3,   # beta2 vector
                tnie_d4,   # beta3 vector
                ##
                tnie_d5,   # theta0
                tnie_d6,   # theta1
                tnie_d7,   # theta2
                tnie_d8,   # theta3
                tnie_d9,   # theta4 vector
                tnie_d10,  # theta5 vector
                tnie_d11   # theta6 vector
            ))
        ##
        ## a's from mreg beta1 should be a1.
        tnde_expit_a1 <- expit(a1*(beta1+beta3_c) + a1*theta3 + beta0 + beta2_c + theta6_c + theta2)
        tnde_expit_a0 <- expit(a1*(beta1+beta3_c) + a0*theta3 + beta0 + beta2_c + theta6_c + theta2)
        
        tnde_d1 <- tnde_expit_a1 - tnde_expit_a0
        tnde_d2 <- a1*tnde_d1
        tnde_d3 <- c_cond*tnde_d1
        # tnde_d4 <- c_cond*a1*tnde_d1
        if(is.null(beta3)){
            tnde_d4 <- rep(0, length(beta3))
        }else{
            tnde_d4 <- c_cond*a1*tnde_d1
        }
        tnde_d5 <- 0
        tnde_d6 <- a1 - a0
        tnde_d7 <- tnde_d1
        tnde_d8 <- a1*tnde_expit_a1 - a0*tnde_expit_a0
        tnde_d9 <- rep(0, length(theta4))
        # tnde_d10 <- c_cond*(a1 - a0)
        if(is.null(theta5)){
            tnde_d10 <- rep(0, length(theta5))
        }else{
            tnde_d10 <- c_cond*(a1 - a0)
        }
        # tnde_d11 <- tnde_d3
        if(is.null(theta6)){
            tnde_d11 <- rep(0, length(theta6))
        }else{
            tnde_d11 <- tnde_d3
        }
        
        Gamma_tnde <-
            matrix(c(
                tnde_d1,   # beta0
                tnde_d2,   # beta1
                tnde_d3,   # beta2 vector
                tnde_d4,   # beta3 vector
                ##
                tnde_d5,   # theta0
                tnde_d6,   # theta1
                tnde_d7,   # theta2
                tnde_d8,   # theta3
                tnde_d9,   # theta4 vector
                tnde_d10,  # theta5 vector
                tnde_d11   # theta6 vector
                ))  
        ##
        ## a's from yreg theta3 should be a0.
        pnie_expit_q1 <- expit(a0*(beta1+beta3_c) + beta0 + beta2_c)
        pnie_expit_q2 <- expit(a0*(beta1+beta3_c) + a0*theta3 + beta0 + beta2_c + theta6_c + theta2)
        pnie_expit_q3 <- expit(a1*(beta1+beta3_c) + beta0 + beta2_c)
        pnie_expit_q4 <- expit(a1*(beta1+beta3_c) + a0*theta3 + beta0 + beta2_c + theta6_c + theta2)
        
        pnie_d1 <- pnie_expit_q1 - pnie_expit_q2 - pnie_expit_q3 + pnie_expit_q4
        pnie_d2 <- a0 * (pnie_expit_q1 - pnie_expit_q2) + a1 * (- pnie_expit_q3 + pnie_expit_q4)
        pnie_d3 <- c_cond * pnie_d1

        if(is.null(beta3)){
            pnie_d4 <- rep(0, length(beta3))
        }else{
            pnie_d4 <- c_cond * pnie_d2
        }
        pnie_d5 <- 0
        pnie_d6 <- 0
        pnie_d7 <- pnie_expit_q4 - pnie_expit_q2
        pnie_d8 <- a0 * pnie_d7 # change from a1
        pnie_d9 <- rep(0, length(theta4))
        pnie_d10 <- rep(0, length(theta5))
        
        if(is.null(theta6)){
            pnie_d11 <- rep(0, length(theta6))
        }else{
            pnie_d11 <- c_cond * pnie_d7
        }
        
        Gamma_pnie <-
            matrix(c(
                pnie_d1,   # beta0
                pnie_d2,   # beta1
                pnie_d3,   # beta2 vector
                pnie_d4,   # beta3 vector
                ##
                pnie_d5,   # theta0
                pnie_d6,   # theta1
                pnie_d7,   # theta2
                pnie_d8,   # theta3
                pnie_d9,   # theta4 vector
                pnie_d10,  # theta5 vector
                pnie_d11   # theta6 vector
                ))  
        ## because Gamma_pnie does not have a common factor.
        Gamma_te <-
            Gamma_pnde + Gamma_tnie # By linearity of differentiation
        ##
        ## Not implemented in mediation.sas.
        ## Not mentioned in VV2013, VV2015, or VanderWeele 2015.
        ## Gradient of pm wrt pnde and tnie. A vector of two.
        ## Copied from calc_myreg_mreg_logistic_yreg_logistic_est
        pnde <-
            ((theta1 + theta5_c) * (a1 - a0)) +
            log(1 + exp(beta0 + beta1*a0 + beta2_c + beta3_c*a0 + theta2 + theta3*a1 + theta6_c)) -
            log(1 + exp(beta0 + beta1*a0 + beta2_c + beta3_c*a0 + theta2 + theta3*a0 + theta6_c))
        tnie <-
            log(1 + exp(beta0 + beta1*a1 + beta2_c + beta3_c*a1 + theta2 + theta3*a1 + theta6_c)) - 
            log(1 + exp(beta0 + beta1*a0 + beta2_c + beta3_c*a0 + theta2 + theta3*a1 + theta6_c)) +
            log(1 + exp(beta0 + beta1*a0 + beta2_c + beta3_c*a0)) - 
            log(1 + exp(beta0 + beta1*a1 + beta2_c + beta3_c*a1))
        ## Need to unname argument vectors to get c(pnde = , tnie = ).
        d_pm <- grad_prop_med_yreg_logistic(pnde = unname(pnde), tnie = unname(tnie))
        ## Multivariate chain rule.
        ## https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Calculus_(OpenStax)/14%3A_Differentiation_of_Functions_of_Several_Variables/14.5%3A_The_Chain_Rule_for_Multivariable_Functions)
        ## d_pm / d_params = (d_pm / d_(pnde, tnie)) %*% (d_(pnde, tnie) / d_params)
        ##                 = (d_pm / d_pnde) * (d_pnde / d_params) +
        ##                   (d_pm / d_tnie) * (d_tnie / d_params)
        ## where (d_pnde / d_params) is Gamma_pnde and
        ##       (d_tnie / d_params) is Gamma_tnie.
        ## FIXME: This is not tested aginst a reference standard.
        Gamma_pm <-
            (d_pm[["pnde"]] * Gamma_pnde) + (d_pm[["tnie"]] * Gamma_tnie)

        ## SE calcuation via multivariate delta method
        ## https://en.wikipedia.org/wiki/Delta_method# Multivariate_delta_method
        se_cde <- sqrt(as.numeric(t(Gamma_cde) %*% Sigma %*% Gamma_cde))
        se_pnde <- sqrt(as.numeric(t(Gamma_pnde) %*% Sigma %*% Gamma_pnde))
        se_tnie <- sqrt(as.numeric(t(Gamma_tnie) %*% Sigma %*% Gamma_tnie))
        se_tnde <- sqrt(as.numeric(t(Gamma_tnde) %*% Sigma %*% Gamma_tnde))
        se_pnie <- sqrt(as.numeric(t(Gamma_pnie) %*% Sigma %*% Gamma_pnie))
        se_te <- sqrt(as.numeric(t(Gamma_te) %*% Sigma %*% Gamma_te))
        se_pm <- sqrt(as.numeric(t(Gamma_pm) %*% Sigma %*% Gamma_pm))

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
