################################################################################
### Internal helper functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################


construct_regmedint <- function(data,
                                yvar,
                                avar,
                                mvar,
                                cvar,
                                a0,
                                a1,
                                m_cde,
                                yreg,
                                mreg,
                                interaction,
                                casecontrol,
                                full_output,
                                c_cond,
                                boot,
                                eventvar) {

    ## Perform mreg
    mreg_fit <- fit_mreg(mreg, data, avar, mvar, cvar)

    ## Perform yreg
    yreg_fit <- fit_yreg(yreg, data, yvar, avar, mvar, cvar, eventvar)

    ## Construct result object
    res <- list(mreg = mreg_fit,
                yreg = yreg_fit,
                ## FIXME: made up numbers
                med = 42L,
                boot = 24L)
    class(res) <- c("regmedint", class(res))

    ##
    res
}
