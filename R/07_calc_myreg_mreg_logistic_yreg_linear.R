################################################################################
### Main functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


## Extension of VanderWeele 2015 p471 Proposition 2.5
##' Create calculators for effects and se (mreg logistic / yreg linear)
##'
##' Construct functions for the conditional effect estimates and their standard errors in the mreg logistic / yreg linear setting. Internally, this function deconstructs model objects and feeds parameter estimates to the internal worker functions \code{calc_myreg_mreg_logistic_yreg_linear_est} and \code{calc_myreg_mreg_logistic_yreg_linear_se}.
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return A list containing a function for effect estimates and a function for corresponding standard errors.

calc_myreg_mreg_logistic_yreg_linear <- function(mreg,
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
        calc_myreg_mreg_logistic_yreg_linear_est(beta0 = beta_hat$beta0,
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
        calc_myreg_mreg_logistic_yreg_linear_se(beta0 = beta_hat$beta0,
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


calc_myreg_mreg_logistic_yreg_linear_est <- function(beta0,
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

        ## Extension of VanderWeele 2015 p471. Estimates are all on log OR scale.
        cde <- (theta1 + theta3*m_cde + theta5_c) * (a1 - a0)
        ## Pearl decomposition (Regular NDE and NIE)
        ## Note the a0 in the second term.
        pnde <- (theta1 + theta3*expit(beta0 + beta1*a0 + beta2_c + beta3_c * a0) + theta5_c) * (a1 - a0)
        ## Note the a1 in the first term.
        tnie <- (theta2 + theta3*a1 + theta6_c) * (expit(beta0 + beta1*a1 + beta2_c + beta3_c*a1) - expit(beta0 + beta1*a0 + beta2_c + beta3_c*a0))
        ## Another decomposition
        ## Note the a0 -> a1 change in the second term.
        tnde <- (theta1 + theta3*expit(beta0 + beta1*a1 + beta2_c + beta3_c * a1) + theta5_c) * (a1 - a0)
        ## Note the a1 -> a0 change in the first term.
        pnie <- (theta2 + theta3*a0 + theta6_c) * (expit(beta0 + beta1*a1 + beta2_c + beta3_c*a1) - expit(beta0 + beta1*a0 + beta2_c + beta3_c*a0))
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
        
        # intermediate quantities: from Deriv
        .e3 <- a0 * (beta1 + beta3_c) + beta0 + beta2_c
        .e4 <- exp(.e3)
        .e5 <- 1 + .e4
        .e6 <- 1 - .e4/.e5
        
        pnde_d1 <- theta3 * .e6 * .e4/.e5
        pnde_d2 <- a0 * theta3 * .e6 * .e4/.e5
        pnde_d3 <- c_cond* theta3 * .e6 * .e4/.e5
        if(is.null(beta3)){
          pnde_d4 <- rep(0, length(beta3))
        }else{
          pnde_d4 <- a0 * c_cond * theta3 * .e6 * .e4/.e5
        }
        pnde_d5 <- 0
        pnde_d6 <- 1
        pnde_d7 <- 0
        pnde_d8 <- expit(.e3)
        pnde_d9 <- rep(0, length(theta4))
        if(is.null(theta5)){
          pnde_d10 <- rep(0, length(theta5))
        }else{
          pnde_d10 <- c_cond
        }
        pnde_d11 <- rep(0, length(theta6))
        
        ## 20211212 Yl: no typo in PNDE
        ## d2 and d3 in VanderWeele 2015 p471 and VV 2013 Appendix p12 have typos.
        ## a0 and c_cond and theta3 are outside the fraction.
        # pnde_d1 <- theta3 * expit(a0*(beta1+beta3_c) + beta0 + beta2_c) * (1 - expit(a0*(beta1+beta3_c) + beta0 + beta2_c))
        # pnde_d2 <- a0*pnde_d1
        # pnde_d3 <- c_cond * pnde_d1
        # if(is.null(beta3)){
        #   pnde_d4 <- rep(0, length(beta3))
        # }else{
        #   pnde_d4 <- c_cond*a0*pnde_d1
        # }
        # pnde_d5 <- 0
        # pnde_d6 <- 1
        # pnde_d7 <- 0
        # pnde_d8 <- expit(a0*(beta1+beta3_c) + beta0 + beta2_c)
        # pnde_d9 <- rep(0, length(theta4))
        # if(is.null(theta5)){
        #   pnde_d10 <- rep(0, length(theta5))
        # }else{
        #   pnde_d10 <- c_cond
        # }
        # pnde_d11 <- rep(0, length(theta6))
        
        
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
        
        # intermediate quantities: from Deriv
        .e1 <- beta1 + beta3_c
        .e2 <- beta2_c
        .e5 <- a0 * .e1 + beta0 + .e2
        .e8 <- a1 * .e1 + beta0 + .e2
        .e9 <- exp(.e5)
        .e10 <- exp(.e8)
        .e11 <- 1 + .e9
        .e12 <- 1 + .e10
        .e13 <- 1 - .e9/.e11
        .e14 <- 1 - .e10/.e12
        .e17 <- a1 * theta3 + theta6_c + theta2
        .e20 <- expit(.e8) - expit(.e5)
        .e25 <- .e14 * .e10/.e12 - .e13 * .e9/.e11
        .e32 <- a1 * .e14 * .e10/.e12 - a0 * .e13 * .e9/.e11

        tnie_d1 <- .e25 * .e17
        tnie_d2 <- .e32 * .e17
        tnie_d3 <- c_cond * .e25 * .e17
        if(is.null(beta3)){
          tnie_d4 <- rep(0, length(beta3))
        }else{
          tnie_d4 <- c_cond * .e32 * .e17
        }
        tnie_d5 <- 0
        tnie_d6 <- 0
        tnie_d7 <- .e20
        tnie_d8 <- a1 * .e20
        tnie_d9 <- rep(0, length(theta4))
        tnie_d10 <- rep(0, length(theta5))
        if(is.null(theta6)){
          tnie_d11 <- rep(0, length(theta6))
        }else{
          tnie_d11 <- c_cond * .e20
        }

        ## 20211212 YL: FIXED typo in tnide_d2. can use either below code or from Deriv
        # tnie_expit_a1 <- expit(a1*(beta1+beta3_c) + beta0 + beta2_c)
        # tnie_expit_a0 <- expit(a0*(beta1+beta3_c) + beta0 + beta2_c)
        # tnie_d1 <- (theta2 + theta3*a1 + theta6_c) * (tnie_expit_a1*(1 - tnie_expit_a1) - tnie_expit_a0*(1 - tnie_expit_a0))
        # tnie_d2 <- (theta2 + theta3*a1 + theta6_c) * (a1*tnie_expit_a1*(1 - tnie_expit_a1) - a0*tnie_expit_a0*(1 - tnie_expit_a0))
        # tnie_d3 <- c_cond * tnie_d1
        # 
        # if(is.null(beta3)){
        #   tnie_d4 <- rep(0, length(beta3))
        # }else{
        #   tnie_d4 <- c_cond * (theta2 + theta3*a1 + theta6_c) * (a1*tnie_expit_a1*(1 - tnie_expit_a1) - a0*tnie_expit_a0*(1 - tnie_expit_a0))
        # }
        # tnie_d5 <- 0
        # tnie_d6 <- 0
        # tnie_d7 <- tnie_expit_a1 - tnie_expit_a0
        # tnie_d8 <- a1 * tnie_d7
        # tnie_d9 <- rep(0, length(theta4))
        # tnie_d10 <- rep(0, length(theta5))
        # if(is.null(theta6)){
        #   tnie_d11 <- rep(0, length(theta6))
        # }else{
        #   tnie_d11 <- c_cond * tnie_d7
        # }
        
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
        
        # intermediate quantities: from Deriv       
        .e3 <- a1 * (beta1 + beta3_c) + beta0 + beta2_c
        .e4 <- exp(.e3)
        .e5 <- 1 + .e4
        .e6 <- 1 - .e4/.e5
        
        tnde_d1 <- theta3 * .e6 * .e4/.e5
        tnde_d2 <- a1 * theta3 * .e6 * .e4/.e5
        tnde_d3 <- c_cond* theta3 * .e6 * .e4/.e5
        if(is.null(beta3)){
          tnde_d4 <- rep(0, length(beta3))
        }else{
          tnde_d4 <- a1 * c_cond * theta3 * .e6 * .e4/.e5
        }
        tnde_d5 <- 0
        tnde_d6 <- 1
        tnde_d7 <- 0
        tnde_d8 <- expit(.e3)
        tnde_d9 <- rep(0, length(theta4))
        if(is.null(theta5)){
          tnde_d10 <- rep(0, length(theta5))
        }else{
          tnde_d10 <- c_cond
        }
        tnde_d11 <- rep(0, length(theta6))
        
        ## 20211212 YL: no typo in TNDE
        # tnde_d1 <- theta3 * (tnie_expit_a1*(1-tnie_expit_a1) - tnie_expit_a0*(1-tnie_expit_a0))
        # tnde_d2 <- a1 * tnde_d1
        # tnde_d3 <- c_cond * tnde_d1
        # if(is.null(beta3)){
        #   tnde_d4 <- rep(0, length(beta3))
        # }else{
        #   tnde_d4 <- c_cond * a1 * tnde_d1
        # }
        # tnde_d5 <- 0 
        # tnde_d6 <- 1
        # tnde_d7 <- 0
        # tnde_d8 <- tnie_expit_a1
        # tnde_d9 <- rep(0, length(theta4))
        # if(is.null(theta5)){
        #   tnde_d10 <- rep(0, length(theta5))
        # }else{
        #   tnde_d10 <- c_cond
        # }
        # tnde_d11 <- rep(0, length(theta6))
        
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
        
        # intermediate quantities: from Deriv
        .e1 <- beta1 + beta3_c
        .e2 <- beta2_c
        .e5 <- a0 * .e1 + beta0 + .e2
        .e8 <- a1 * .e1 + beta0 + .e2
        .e9 <- exp(.e5)
        .e10 <- exp(.e8)
        .e11 <- 1 + .e9
        .e12 <- 1 + .e10
        .e13 <- 1 - .e9/.e11
        .e14 <- 1 - .e10/.e12
        .e17 <- a0 * theta3 + theta6_c + theta2
        .e20 <- expit(.e8) - expit(.e5)
        .e25 <- .e14 * .e10/.e12 - .e13 * .e9/.e11
        .e32 <- a1 * .e14 * .e10/.e12 - a0 * .e13 * .e9/.e11

        pnie_d1 <- .e25 * .e17
        pnie_d2 <- .e32 * .e17
        pnie_d3 <- c_cond * .e25 * .e17
        if(is.null(beta3)){
          pnie_d4 <- rep(0, length(beta3))
        }else{
          pnie_d4 <- c_cond * .e32 * .e17
        }
        pnie_d5 <- 0
        pnie_d6 <- 0
        pnie_d7 <- .e20
        pnie_d8 <- a0 * .e20
        pnie_d9 <- rep(0, length(theta4))
        pnie_d10 <- rep(0, length(theta5))
        if(is.null(theta6)){
          pnie_d11 <- rep(0, length(theta6))
        }else{
          pnie_d11 <- c_cond * .e20
        }
        
        ## 20211212 YL: FIXED typo in pnide_d2. can use either below code or from Deriv
        # pnie_d1 <- (theta2 + theta3*a0 + theta6_c) * (tnie_expit_a1*(1-tnie_expit_a1) - tnie_expit_a0*(1-tnie_expit_a0))
        # pnie_d2 <- (theta2 + theta3*a0 + theta6_c) * (a1*tnie_expit_a1*(1 - tnie_expit_a1) - a0*tnie_expit_a0*(1 - tnie_expit_a0))
        # pnie_d3 <- c_cond * pnie_d1
        # if(is.null(beta3)){
        #   pnie_d4 <- rep(0, length(beta3))
        # }else{
        #   pnie_d4 <- c_cond * pnie_d2
        # }
        # pnie_d5 <- 0
        # pnie_d6 <- 0
        # pnie_d7 <- tnie_expit_a1 - tnie_expit_a0
        # pnie_d8 <- a0 * pnie_d7
        # pnie_d9 <- rep(0, length(theta4))
        # pnie_d10 <- rep(0, length(theta5))
        # # pnie_d11 <- c_cond*pnie_d7
        # if(is.null(theta6)){
        #   pnie_d11 <- rep(0, length(theta6))
        # }else{
        #   pnie_d11 <- c_cond * pnie_d7
        # }
        # 
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
        
        ## (a1 - a0) without abs must enter here for pnde
        ## because Gamma_pnie does not have a common factor.
        ## d_te/d_params = d_(pnde + tnie)/d_params
        ##               = d_pnde/d_params + d_tnie/d_params
        ##               = (a1-a0) * Gamma_pnde + Gamma_tnie
        ## Do not use abs(a1 - a0 in the se_te derivation below.
        ## VV2013 Appendix p14 for Gamma_te has an error in d3 where
        ## Gamma_tnie's contribution misses c'. Otherwise, it has
        ## the following linear combination form.
        Gamma_te <-
            ((a1 - a0) * Gamma_pnde) + Gamma_tnie # By linearity of differentiation
        ##
        ## Not implemented in mediation.sas.
        ## Not mentioned in VV2013, VV2015, or VanderWeele 2015.
        ## Gradient of pm wrt pnde and tnie. A vector of two.
        ## Copied from calc_myreg_mreg_logistic_yreg_linear_est
        pnde <- (theta1 + theta3*expit(beta0 + beta1*a0 + beta2_c + beta3_c * a0) + theta5_c) * (a1 - a0)
        tnie <- (theta2 + theta3*a1 + theta6_c) * (expit(beta0 + beta1*a1 + beta2_c + beta3_c*a1) - expit(beta0 + beta1*a0 + beta2_c + beta3_c*a0))
        ## Need to unname argument vectors to get c(pnde = , tnie = ).
        d_pm <- grad_prop_med_yreg_linear(pnde = unname(pnde), tnie = unname(tnie))
        ## Multivariate chain rule.
        ## https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Calculus_(OpenStax)/14%3A_Differentiation_of_Functions_of_Several_Variables/14.5%3A_The_Chain_Rule_for_Multivariable_Functions)
        ## d_pm / d_params = (d_pm / d_(pnde, tnie)) %*% (d_(pnde, tnie) / d_params)
        ##                 = (d_pm / d_pnde) * (d_pnde / d_params) +
        ##                   (d_pm / d_tnie) * (d_tnie / d_params)
        ## where (d_pnde / d_params) is (a1 - a0) * Gamma_pnde and
        ##       (d_tnie / d_params) is Gamma_tnie.
        ## FIXME: This is not tested aginst a reference standard.
        Gamma_pm <-
            (d_pm[["pnde"]] * (a1 - a0) * Gamma_pnde) + (d_pm[["tnie"]] * Gamma_tnie)

        ## SE calcuation via multivariate delta method
        ## https://en.wikipedia.org/wiki/Delta_method# Multivariate_delta_method
        ## NIEs do not have common factor abs(a1 - a0), thus, it does not show up
        ## in se_te and se_pm.
        a1_sub_a0 <- abs(a1 - a0)
        se_cde <- sqrt(as.numeric(t(Gamma_cde) %*% Sigma %*% Gamma_cde)) * a1_sub_a0
        se_pnde <- sqrt(as.numeric(t(Gamma_pnde) %*% Sigma %*% Gamma_pnde)) * a1_sub_a0
        se_tnie <- sqrt(as.numeric(t(Gamma_tnie) %*% Sigma %*% Gamma_tnie))
        se_tnde <- sqrt(as.numeric(t(Gamma_tnde) %*% Sigma %*% Gamma_tnde)) * a1_sub_a0
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
