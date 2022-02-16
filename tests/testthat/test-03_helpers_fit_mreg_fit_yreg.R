################################################################################
### Tests for internal functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)


###
### Internal function for mreg string formula creation
################################################################################

describe("string_mreg_formula", {
    describe("string_mreg_formula (good args)", {
        it("handles NULL cvar by omitting", {
            expect_equal(string_mreg_formula(mvar = "M",
                                             avar = "A",
                                             cvar = NULL,
                                             emm_ac_mreg = NULL),
                         "M ~ A")
        })
        it("handles one cvar by adding", {
            expect_equal(string_mreg_formula(mvar = "M",
                                             avar = "A",
                                             cvar = c("C"),
                                             emm_ac_mreg = NULL),
                         "M ~ A + C")
        })
        it("handles three cvar by adding all", {
            expect_equal(string_mreg_formula(mvar = "M",
                                             avar = "A",
                                             cvar = c("C1","C2","C3"),
                                             emm_ac_mreg = NULL),
                         "M ~ A + C1 + C2 + C3")
        })
    })
    ##
    describe("string_mreg_formula (bad args)", {
        it("throws an error on NULL mvar", {
            expect_error(string_mreg_formula(mvar = NULL,
                                             avar = "A",
                                             cvar = NULL,
                                             emm_ac_mreg = NULL))
        })
        it("throws an error on NULL avar", {
            expect_error(string_mreg_formula(mvar = "M",
                                             avar = NULL,
                                             cvar = NULL,
                                             emm_ac_mreg = NULL))
        })
    })
})


###
### Internal function for yreg string formula creation
################################################################################

describe("string_yreg_formula", {
    describe("string_yreg_formula (non-survival yvar)", {
        describe("string_yreg_formula (non-survival yvar; good args)", {
            it("handles NULL cvar by omitting", {
                ## Zero covariates
                expect_equal(string_yreg_formula(yvar = "Y",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = NULL,
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = NULL),
                             "Y ~ A + M")
                expect_equal(string_yreg_formula(yvar = "Y",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = NULL,
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = TRUE,
                                                 eventvar = NULL),
                             "Y ~ A + M + A:M")
            })
            it("handles one cvar by adding", {
                ## One covariate
                expect_equal(string_yreg_formula(yvar = "Y",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = c("C"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = NULL),
                             "Y ~ A + M + C")
                expect_equal(string_yreg_formula(yvar = "Y",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = c("C"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = TRUE,
                                                 eventvar = NULL),
                             "Y ~ A + M + A:M + C")
            })
            it("handles three cvar by adding all", {
                ## Three covariates
                expect_equal(string_yreg_formula(yvar = "Y",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = c("C1","C2","C3"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = NULL),
                             "Y ~ A + M + C1 + C2 + C3")
                expect_equal(string_yreg_formula(yvar = "Y",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = c("C1","C2","C3"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = TRUE,
                                                 eventvar = NULL),
                             "Y ~ A + M + A:M + C1 + C2 + C3")
            })
        })
        describe("string_yreg_formula (non-survival yvar; bad args)", {
            it("throws an error on NULL yvar", {
                expect_error(string_yreg_formula(yvar = NULL,
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = c("C1","C2","C3"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = NULL))
            })
            it("throws an error on NULL mvar", {
                expect_error(string_yreg_formula(yvar = "Y",
                                                 avar = NULL,
                                                 mvar = "M",
                                                 cvar = c("C1","C2","C3"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = NULL))
            })
            it("throws an error on NULL avar", {
                expect_error(string_yreg_formula(yvar = "Y",
                                                 avar = "A",
                                                 mvar = NULL,
                                                 cvar = c("C1","C2","C3"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = NULL))
            })
        })
    })
    ##
    describe("string_yreg_formula (survival yvar)", {
        describe("string_yreg_formula (survival yvar; good args)", {
            it("handles NULL cvar by omitting", {
                ## Zero covariates
                expect_equal(string_yreg_formula(yvar = "time",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = NULL,
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = "event"),
                             "Surv(time, event) ~ A + M")
                expect_equal(string_yreg_formula(yvar = "time",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = NULL,
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = TRUE,
                                                 eventvar = "event"),
                             "Surv(time, event) ~ A + M + A:M")
            })
            it("handles one cvar by adding", {
                ## One covariate
                expect_equal(string_yreg_formula(yvar = "time",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = c("C"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = "event"),
                             "Surv(time, event) ~ A + M + C")
                expect_equal(string_yreg_formula(yvar = "time",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = c("C"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = TRUE,
                                                 eventvar = "event"),
                             "Surv(time, event) ~ A + M + A:M + C")
            })
            it("handles three cvar by adding all", {
                ## Three covariates
                expect_equal(string_yreg_formula(yvar = "time",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = c("C1","C2","C3"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = "event"),
                             "Surv(time, event) ~ A + M + C1 + C2 + C3")
                expect_equal(string_yreg_formula(yvar = "time",
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = c("C1","C2","C3"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = TRUE,
                                                 eventvar = "event"),
                             "Surv(time, event) ~ A + M + A:M + C1 + C2 + C3")
            })
        })
        describe("string_yreg_formula (survival yvar; bad args)", {
            it("throws an error on NULL yvar", {
                expect_error(string_yreg_formula(yvar = NULL,
                                                 avar = "A",
                                                 mvar = "M",
                                                 cvar = c("C1","C2","C3"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = "event"))
            })
            it("throws an error on NULL mvar", {
                expect_error(string_yreg_formula(yvar = "time",
                                                 avar = NULL,
                                                 mvar = "M",
                                                 cvar = c("C1","C2","C3"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = "event"))
            })
            it("throws an error on NULL avar", {
                expect_error(string_yreg_formula(yvar = "time",
                                                 avar = "A",
                                                 mvar = NULL,
                                                 cvar = c("C1","C2","C3"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL,
                                                 interaction = FALSE,
                                                 eventvar = "event"))
            })
        })
    })
})
