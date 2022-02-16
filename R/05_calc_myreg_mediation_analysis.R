################################################################################
### Main functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


##' Return mediation analysis functions given mediator and outcome models specifications.
##'
##' This function returns functions that can be used to calculate the causal effect measures, given the mediator model fit (\code{mreg_fit}) and the outcome model fit (\code{yreg_fit}).
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit from \code{\link{fit_mreg}}
##' @param yreg_fit Model fit from \code{\link{fit_yreg}}
##'
##' @return A list containing two functions. The first is for calculating point estimates. The second is for calculating the correspoding
calc_myreg <- function(mreg,
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

    supported_nonlinear_yreg <- c("logistic","loglinear","poisson","negbin",
                                  "survCox","survAFT_exp","survAFT_weibull")

    ## FIXME: Use this to do.call()
    args <- list(mreg = mreg,
                 mreg_fit = mreg_fit,
                 yreg = yreg,
                 yreg_fit = yreg_fit,
                 avar = avar,
                 mvar = mvar,
                 cvar = cvar,
                 emm_ac_mreg = emm_ac_mreg,
                 emm_ac_yreg = emm_ac_yreg,
                 emm_mc_yreg = emm_mc_yreg,
                 interaction = interaction)

    ## Only four patterns as the non-linear yreg cases are the
    ## same as the logistic yreg case.
    ## See VanderWeele 2015 Appendix.
    ##     Valeri & VanderWeele 2013 Appendix.
    ##     README for this repo.
    ## Each helper function is defined in a dedicated file.
    if (mreg == "linear" & yreg == "linear") {

        ## Extended VanderWeele 2015 p466 Proposition 2.3
        list_est_fun_se_fun <-
            calc_myreg_mreg_linear_yreg_linear(mreg = mreg,
                                               mreg_fit = mreg_fit,
                                               yreg = yreg,
                                               yreg_fit = yreg_fit,
                                               avar = avar,
                                               mvar = mvar,
                                               cvar = cvar,
                                               emm_ac_mreg = emm_ac_mreg,
                                               emm_ac_yreg = emm_ac_yreg,
                                               emm_mc_yreg = emm_mc_yreg,
                                               interaction = interaction)

    } else if (mreg == "linear" & yreg %in% supported_nonlinear_yreg) {

        ## Extension of VanderWeele 2015 p468 Proposition 2.4
        list_est_fun_se_fun <-
            calc_myreg_mreg_linear_yreg_logistic(mreg = mreg,
                                                 mreg_fit = mreg_fit,
                                                 yreg = yreg,
                                                 yreg_fit = yreg_fit,
                                                 avar = avar,
                                                 mvar = mvar,
                                                 cvar = cvar,
                                                 emm_ac_mreg = emm_ac_mreg,
                                                 emm_ac_yreg = emm_ac_yreg,
                                                 emm_mc_yreg = emm_mc_yreg,
                                                 interaction = interaction)

    } else if (mreg == "logistic" & yreg == "linear") {

        ## Extension of VanderWeele 2015 p471 Proposition 2.5
        list_est_fun_se_fun <-
            calc_myreg_mreg_logistic_yreg_linear(mreg = mreg,
                                                 mreg_fit = mreg_fit,
                                                 yreg = yreg,
                                                 yreg_fit = yreg_fit,
                                                 avar = avar,
                                                 mvar = mvar,
                                                 cvar = cvar,
                                                 emm_ac_mreg = emm_ac_mreg,
                                                 emm_ac_yreg = emm_ac_yreg,
                                                 emm_mc_yreg = emm_mc_yreg,
                                                 interaction = interaction)


    } else if (mreg == "logistic" & yreg %in% supported_nonlinear_yreg) {

        ## Extended VanderWeele 2015 p473 Proposition 2.6
        list_est_fun_se_fun <-
            calc_myreg_mreg_logistic_yreg_logistic(mreg = mreg,
                                                   mreg_fit = mreg_fit,
                                                   yreg = yreg,
                                                   yreg_fit = yreg_fit,
                                                   avar = avar,
                                                   mvar = mvar,
                                                   cvar = cvar,
                                                   emm_ac_mreg = emm_ac_mreg,
                                                   emm_ac_yreg = emm_ac_yreg,
                                                   emm_mc_yreg = emm_mc_yreg,
                                                   interaction = interaction)

    } else  {

        stop("Unsupported mreg or yreg!")

    }

    ## Return a list of the form:
    ##  list(est_fun = est_fun,
    ##       se_fun = est_fun)
    return(list_est_fun_se_fun)
}
