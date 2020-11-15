################################################################################
### Internal helper functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################

###
### Constructor for internal use
################################################################################

##' Low level constructor for a regmedint S3 class object.
##'
##' This is not a user function and meant to be executed within the regmedint function after validatingthe arguments.
##'
##' @inheritParams regmedint
##'
##' @return A regmedint object.
new_regmedint <- function(data,
                          yvar,
                          avar,
                          mvar,
                          cvar,
                          eventvar,
                          a0,
                          a1,
                          m_cde,
                          c_cond,
                          yreg,
                          mreg,
                          interaction,
                          casecontrol) {

    ## If the data is from a case-control study the mediator model
    ## should be fit in the controls. This requires that the outcome
    ## be rare in the population from which the controls were sampled.
    ## Valeri & VanderWeele 2013. p144
    if (casecontrol) {
        ## Single index to extract as an unnamed vector.
        ## This works on both data.frame and tbl_df.
        mreg_data <- data[data[[yvar]] == 0,]
    } else {
        mreg_data <- data
    }

    ## Perform mreg
    mreg_fit <- fit_mreg(mreg = mreg,
                         data = mreg_data,
                         avar = avar,
                         mvar = mvar,
                         cvar = cvar)

    ## Perform yreg
    yreg_fit <- fit_yreg(yreg = yreg,
                         data = data,
                         yvar = yvar,
                         avar = avar,
                         mvar = mvar,
                         cvar = cvar,
                         interaction = interaction,
                         eventvar = eventvar)

    ## Return a list of functions
    myreg_funs <- calc_myreg(mreg = mreg,
                             mreg_fit = mreg_fit,
                             yreg = yreg,
                             yreg_fit = yreg_fit,
                             avar = avar,
                             mvar = mvar,
                             cvar = cvar,
                             interaction = interaction)

    ## Construct the result object
    res <- list(mreg_fit = mreg_fit,
                yreg_fit = yreg_fit,
                myreg = myreg_funs,
                ## Remember arguments
                args = list(yvar = yvar,
                            avar = avar,
                            mvar = mvar,
                            cvar = cvar,
                            a0 = a0,
                            a1 = a1,
                            m_cde = m_cde,
                            yreg = yreg,
                            mreg = mreg,
                            interaction = interaction,
                            casecontrol = casecontrol,
                            c_cond = c_cond,
                            eventvar = eventvar))
    ## The main class is regmedint.
    class(res) <- c("regmedint", class(res))

    ##
    res
}
