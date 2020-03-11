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
                                                        sigma,
                                                        a0,
                                                        a1,
                                                        m_cde,
                                                        c_cond),
         se = calc_myreg_mreg_linear_yreg_logistic_se(beta0,
                                                      beta1,
                                                      beta2,
                                                      theta1,
                                                      theta2,
                                                      theta3,
                                                      theta4,
                                                      sigma,
                                                      ## vcov
                                                      Sigma_beta,
                                                      Sigma_theta,
                                                      Sigma_sigma,
                                                      ## Values at which to evaluate effects
                                                      a0,
                                                      a1,
                                                      m_cde,
                                                      c_cond))
}


calc_myreg_mreg_linear_yreg_logistic_est <- function(beta0,
                                                     beta1,
                                                     beta2,
                                                     theta1,
                                                     theta2,
                                                     theta3,
                                                     theta4,
                                                     sigma,
                                                     a0,
                                                     a1,
                                                     m_cde,
                                                     c_cond) {

    m <- m_cde
    rm <- sigma^2 # FIXME
    asq <- a1^2
    a1sq <- a0^2 #
    tsq <- theta3^2 # FIXME

    ## Adopted from mediation.sas
    ## Look up
    ## %if &yreg^=linear & &mreg=linear & &interaction=true %then %do;
    ## */MARGINAL=CONDITIONAL CDE*/;
    x1=(theta1+theta3*m)*(a1-a0);
    ## */MARGINAL=CONDITIONAL NDE*/;
    x2=(theta1+theta3*beta0+theta3*beta1*a0+theta3*theta2*rm)*(a1-a0)+(1/2)*tsq*rm*(asq-a1sq);
    ## */MARGINAL=CONDITIONAL NIE*/;
    x3=(theta2*beta1+theta3*beta1*a0)*(a1-a0);
    ## */MARGINAL=CONDITIONAL TNDE*/;
    x4=(theta1+theta3*beta0+theta3*beta1*a1+theta3*theta2*rm)*(a1-a0)+(1/2)*tsq*rm*(asq-a1sq);
    ## */ MARGINAL=CONDITIONAL TNIE*/;
    x5=(theta2*beta1+theta3*beta1*a1)*(a1-a0);

}
