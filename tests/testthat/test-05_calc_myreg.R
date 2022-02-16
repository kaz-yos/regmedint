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

data(pbc)
## Missing data should be warned in validate_args()
pbc_cc <- pbc[complete.cases(pbc),] %>%
    mutate(male = if_else(sex == "m", 1L, 0L),
           ## Combine transplant and death for testing purpose
           status = if_else(status == 0, 0L, 1L),
           ## Binary mvar
           bili_bin = if_else(bili > median(bili), 1L, 0L),
           alk_phos = alk.phos)

describe("calc_myreg", {
    ##
    it("errors informatively when mreg is unsupported", {
        mreg_fit <- fit_mreg(mreg = "linear",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili",
                             cvar = NULL,
                             emm_ac_mreg = NULL)
        yreg_fit <- fit_yreg(yreg = "linear",
                             data = pbc_cc,
                             yvar = "alk_phos",
                             avar = "trt",
                             mvar = "bili",
                             cvar = NULL,
                             emm_ac_yreg = NULL,
                             emm_mc_yreg = NULL,
                             interaction = TRUE,
                             eventvar = NULL)
        expect_error(calc_myreg(mreg = "unsupported",
                                mreg_fit = mreg_fit,
                                yreg = "linear",
                                yreg_fit = yreg_fit,
                                avar = "trt",
                                mvar = "bili",
                                cvar = NULL,
                                emm_ac_mreg = NULL,
                                emm_ac_yreg = NULL,
                                emm_mc_yreg = NULL,
                                interaction = TRUE),
                     regexp = "Unsupported mreg or yreg!")
    })
    ##
    it("errors informatively when yreg is unsupported", {
        mreg_fit <- fit_mreg(mreg = "linear",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili",
                             cvar = NULL,
                             emm_ac_mreg = NULL)
        yreg_fit <- fit_yreg(yreg = "linear",
                             data = pbc_cc,
                             yvar = "alk_phos",
                             avar = "trt",
                             mvar = "bili",
                             cvar = NULL,
                             interaction = TRUE,
                             emm_ac_yreg = NULL,
                             emm_mc_yreg = NULL,
                             eventvar = NULL)
        expect_error(calc_myreg(mreg = "linear",
                                mreg_fit = mreg_fit,
                                yreg = "unsupported",
                                yreg_fit = yreg_fit,
                                avar = "trt",
                                mvar = "bili",
                                cvar = NULL,
                                emm_ac_mreg = NULL,
                                emm_ac_yreg = NULL,
                                emm_mc_yreg = NULL,
                                interaction = TRUE),
                     regexp = "Unsupported mreg or yreg!")
    })
    ##
    it("calls calc_myreg_mreg_linear_yreg_linear when mreg linear / yreg linear", {
        mreg_fit <- fit_mreg(mreg = "linear",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili",
                             cvar = NULL,
                             emm_ac_mreg = NULL)
        yreg_fit <- fit_yreg(yreg = "linear",
                             data = pbc_cc,
                             yvar = "alk_phos",
                             avar = "trt",
                             mvar = "bili",
                             cvar = NULL,
                             emm_ac_yreg = NULL,
                             emm_mc_yreg = NULL,
                             interaction = TRUE,
                             eventvar = NULL)
        with_mock(
            ## Mock
            ## https://github.com/r-lib/testthat/issues/734
            "regmedint:::calc_myreg_mreg_linear_yreg_linear" =
                function(...) {
                    message("calc_myreg_mreg_linear_yreg_linear was called!")
                },
            ## Body
            {
                expect_message(calc_myreg(mreg = "linear",
                                          mreg_fit = mreg_fit,
                                          yreg = "linear",
                                          yreg_fit = yreg_fit,
                                          avar = "trt",
                                          mvar = "bili",
                                          cvar = NULL,
                                          emm_ac_mreg = NULL,
                                          emm_ac_yreg = NULL,
                                          emm_mc_yreg = NULL,
                                          interaction = TRUE),
                               regexp = "calc_myreg_mreg_linear_yreg_linear was called!")
            })
    })
    ##
    it("calls calc_myreg_mreg_linear_yreg_logistic when mreg linear / yreg logistic", {
        mreg_fit <- fit_mreg(mreg = "linear",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili",
                             cvar = NULL,
                             emm_ac_mreg = NULL)
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili",
                             cvar = NULL,
                             emm_ac_yreg = NULL,
                             emm_mc_yreg = NULL,
                             interaction = TRUE,
                             eventvar = NULL)
        with_mock(
            ## Mock
            ## https://github.com/r-lib/testthat/issues/734
            "regmedint:::calc_myreg_mreg_linear_yreg_logistic" =
                function(...) {
                    message("calc_myreg_mreg_linear_yreg_logistic was called!")
                },
            ## Body
            {
                expect_message(calc_myreg(mreg = "linear",
                                          mreg_fit = mreg_fit,
                                          yreg = "logistic",
                                          yreg_fit = yreg_fit,
                                          avar = "trt",
                                          mvar = "bili",
                                          cvar = NULL,
                                          emm_ac_mreg = NULL,
                                          emm_ac_yreg = NULL,
                                          emm_mc_yreg = NULL,
                                          interaction = TRUE),
                               regexp = "calc_myreg_mreg_linear_yreg_logistic was called!")
            })
    })
    ##
    it("calls calc_myreg_mreg_logistic_yreg_linear when mreg logistic / yreg linear", {
        mreg_fit <- fit_mreg(mreg = "logistic",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = NULL,
                             emm_ac_mreg = NULL)
        yreg_fit <- fit_yreg(yreg = "linear",
                             data = pbc_cc,
                             yvar = "alk_phos",
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = NULL,
                             emm_ac_yreg = NULL,
                             emm_mc_yreg = NULL,
                             interaction = TRUE,
                             eventvar = NULL)
        with_mock(
            ## Mock
            ## https://github.com/r-lib/testthat/issues/734
            "regmedint:::calc_myreg_mreg_logistic_yreg_linear" =
                function(...) {
                    message("calc_myreg_mreg_logistic_yreg_linear was called!")
                },
            ## Body
            {
                expect_message(calc_myreg(mreg = "logistic",
                                          mreg_fit = mreg_fit,
                                          yreg = "linear",
                                          yreg_fit = yreg_fit,
                                          avar = "trt",
                                          mvar = "bili",
                                          cvar = NULL,
                                          emm_ac_mreg = NULL,
                                          emm_ac_yreg = NULL,
                                          emm_mc_yreg = NULL,
                                          interaction = TRUE),
                               regexp = "calc_myreg_mreg_logistic_yreg_linear was called!")
            })
    })
    ##
    it("calls calc_myreg_mreg_logistic_yreg_logistic when mreg logistic / yreg logistic", {
        mreg_fit <- fit_mreg(mreg = "logistic",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = NULL,
                             emm_ac_mreg = NULL)
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = NULL,
                             emm_ac_yreg = NULL,
                             emm_mc_yreg = NULL,
                             interaction = TRUE,
                             eventvar = NULL)
        with_mock(
            ## Mock
            ## https://github.com/r-lib/testthat/issues/734
            "regmedint:::calc_myreg_mreg_logistic_yreg_logistic" =
                function(...) {
                    message("calc_myreg_mreg_logistic_yreg_logistic was called!")
                },
            ## Body
            {
                expect_message(calc_myreg(mreg = "logistic",
                                          mreg_fit = mreg_fit,
                                          yreg = "logistic",
                                          yreg_fit = yreg_fit,
                                          avar = "trt",
                                          mvar = "bili",
                                          cvar = NULL,
                                          emm_ac_mreg = NULL,
                                          emm_ac_yreg = NULL,
                                          emm_mc_yreg = NULL,
                                          interaction = TRUE),
                               regexp = "calc_myreg_mreg_logistic_yreg_logistic was called!")
            })
    })
})
