################################################################################
### Specifications for the main mediation analysis function
##
## Created on: 2020-03-12
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)
library(survival)
library(tidyverse)


###
### calc_myreg
################################################################################

## The only job of calc_myreg is to delegate the subsequent work to the correct
## specialized functions like calc_myreg_mreg_linear_yreg_logistic

describe("calc_myreg", {
    describe("calc_myreg mreg linear", {
        it("calls calc_myreg_mreg_linear_yreg_linear when mreg linear/yreg linear", {
            with_mock(
                ## Mock
                calc_myreg_mreg_linear_yreg_linear =
                    function(...) {
                        message("calc_myreg_mreg_linear_yreg_linear was called!")
                    },
                ## Body
                {
                    expect_message(calc_myreg(mreg,
                                              mreg_fit,
                                              yreg,
                                              yreg_fit,
                                              avar,
                                              mvar,
                                              cvar,
                                              interaction),
                                   "calc_myreg_mreg_linear_yreg_logistic was called!")
                })
        })
        ##
        it("calls calc_myreg_mreg_linear_yreg_logistic when mreg linear/yreg logistic", {
            with_mock(
                ## Mock
                calc_myreg_mreg_linear_yreg_logistic =
                    function(...) {
                        message("calc_myreg_mreg_linear_yreg_logistic was called!")
                    },
                ## Body
                {
                    expect_message(calc_myreg(mreg,
                                              mreg_fit,
                                              yreg,
                                              yreg_fit,
                                              avar,
                                              mvar,
                                              cvar,
                                              interaction),
                                   "calc_myreg_mreg_linear_yreg_logistic was called!")
                })
        })
    })
})
