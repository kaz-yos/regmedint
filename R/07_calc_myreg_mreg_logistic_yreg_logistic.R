################################################################################
### Main functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


## VanderWeele 2015 p473 Proposition 2.6
##' Calculate effect measures based on two regression fits (linear/linear)
##'
##' Causal effect parameters are calculated.
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return
calc_myreg_mreg_logistic_yreg_logistic <- function(mreg,
                                                   mreg_fit,
                                                   yreg,
                                                   yreg_fit,
                                                   interaction,
                                                   a0,
                                                   a1,
                                                   m_cde,
                                                   c_cond) {

    ## FIXME: This needs to be allowed and handled nicely!
    assertthat::assert_that(!is.null(cvar))
}
