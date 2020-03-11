################################################################################
### Main functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


##' Calculate effect measures based on two regression fits.
##'
##' Causal effect measures are calculated given the mediator model fit (\code{mreg_fit}) and the outcome model fit (\code{yreg_fit}).
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return \code{myreg} object containing effect measures.
calc_myreg <- function(mreg,
                       mreg_fit,
                       yreg,
                       yreg_fit,
                       interaction,
                       a0,
                       a1,
                       m_cde,
                       c_cond) {

    ## Only four patterns as the non-linear yreg cases are the
    ## same as the logistic yreg case.
    ## See VanderWeele 2015 Appendix.
    ##     Valeri & VanderWeele 2013 Appendix.
    ##     README for this repo.
    ## Each helper function is defined in a dedicated file.
    if (mreg == "linear" & yreg == "linear") {

        ## VanderWeele 2015 p466 Proposition 2.3
        calc_myreg_mreg_linear_yreg_linear(mreg,
                                           mreg_fit,
                                           yreg,
                                           yreg_fit,
                                           interaction,
                                           a0,
                                           a1,
                                           m_cde,
                                           c_cond)

    } else if (mreg == "linear" & yreg != "linear") {

        ## VanderWeele 2015 p468 Proposition 2.4
        calc_myreg_mreg_linear_yreg_logistic(mreg,
                                             mreg_fit,
                                             yreg,
                                             yreg_fit,
                                             interaction,
                                             a0,
                                             a1,
                                             m_cde,
                                             c_cond)

    } else if (mreg == "logistic" & yreg == "linear") {

        ## VanderWeele 2015 p471 Proposition 2.5
        calc_myreg_mreg_logistic_yreg_linear(mreg,
                                             mreg_fit,
                                             yreg,
                                             yreg_fit,
                                             interaction,
                                             a0,
                                             a1,
                                             m_cde,
                                             c_cond)


    } else if (mreg == "logistic" & yreg != "linear") {

        ## VanderWeele 2015 p473 Proposition 2.6
        calc_myreg_mreg_logistic_yreg_logistic(mreg,
                                               mreg_fit,
                                               yreg,
                                               yreg_fit,
                                               interaction,
                                               a0,
                                               a1,
                                               m_cde,
                                               c_cond)

    } else  {

        stop("Unsupported mreg or yreg.")

    }
}


## VanderWeele 2015 p466 Proposition 2.3
##' Calculate effect measures based on two regression fits (linear/linear)
##'
##' Causal effect parameters are calculated.
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return
calc_myreg_mreg_linear_yreg_linear <- function(mreg,
                                               mreg_fit,
                                               yreg,
                                               yreg_fit,
                                               interaction,
                                               a0,
                                               a1,
                                               m_cde,
                                               c_cond) {

    calc_myreg_mreg_linear_yreg_linear_est
    calc_myreg_mreg_linear_yreg_linear_se

    ## class myreg_mreg_linear_yreg_linear, myreg
}


## VanderWeele 2015 p468 Proposition 2.4
calc_myreg_mreg_linear_yreg_logistic <- function(mreg,
                                                 mreg_fit,
                                                 yreg,
                                                 yreg_fit,
                                                 interaction,
                                                 a0,
                                                 a1,
                                                 m_cde,
                                                 c_cond) {
}


## VanderWeele 2015 p471 Proposition 2.5
calc_myreg_mreg_logistic_yreg_linear <- function(mreg,
                                                 mreg_fit,
                                                 yreg,
                                                 yreg_fit,
                                                 interaction,
                                                 a0,
                                                 a1,
                                                 m_cde,
                                                 c_cond) {
}


## VanderWeele 2015 p473 Proposition 2.6
calc_myreg_mreg_logistic_yreg_logistic <- function(mreg,
                                                   mreg_fit,
                                                   yreg,
                                                   yreg_fit,
                                                   interaction,
                                                   a0,
                                                   a1,
                                                   m_cde,
                                                   c_cond) {
}
