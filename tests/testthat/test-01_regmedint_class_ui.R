################################################################################
### Tests for regmedint
##
## Created on: 2020-03-16
## Author: Kazuki Yoshida
################################################################################

library(testthat)
library(survival)


###
### regmedint user interface
################################################################################

describe("regmedint", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               ## Combine transplant and death for testing purpose
               status = if_else(status == 0, 0L, 1L))

    describe("regmedint mreg linear yreg linear", {
        it("runs witmh zero cvar with no interaction", {
            fit_regmedint <- regmedint(data = pbc_cc,
                                       yvar = "alk.phos",
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       a0 = 1,
                                       a1 = 2,
                                       m_cde = 0,
                                       c_cond = NULL,
                                       mreg = "linear",
                                       yreg = "linear",
                                       interaction = FALSE,
                                       casecontrol = FALSE,
                                       eventvar = NULL)
        })
    })
    ##
    describe("regmedint mreg linear yreg logistic", {
        it("runs witmh zero cvar with no interaction", {
            fit_regmedint <- regmedint(data = pbc_cc,
                                       yvar = "spiders",
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       a0 = 1,
                                       a1 = 2,
                                       m_cde = 0,
                                       c_cond = NULL,
                                       mreg = "linear",
                                       yreg = "logistic",
                                       interaction = FALSE,
                                       casecontrol = FALSE,
                                       eventvar = NULL)
        })
    })
})
