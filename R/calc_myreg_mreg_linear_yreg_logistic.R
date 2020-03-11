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
##' @return
calc_myreg_mreg_linear_yreg_logistic <- function(mreg,
                                                 mreg_fit,
                                                 yreg,
                                                 yreg_fit,
                                                 interaction,
                                                 a0,
                                                 a1,
                                                 m_cde,
                                                 c_cond) {

    list(est = calc_myreg_mreg_linear_yreg_logistic_est(beta0,
                                                        beta1,
                                                        beta2,
                                                        theta1,
                                                        theta2,
                                                        theta3,
                                                        theta4,
                                                        sigma),
         se = calc_myreg_mreg_linear_yreg_logistic_se(beta0,
                                                      beta1,
                                                      beta2,
                                                      theta1,
                                                      theta2,
                                                      theta3,
                                                      theta4,
                                                      sigma,
                                                      Sigma_beta,
                                                      Sigma_theta,
                                                      Sigma_sigma))

}
