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
                                             cvar = NULL),
                         "M ~ A")
        })
        it("handles one cvar by adding", {
            expect_equal(string_mreg_formula(mvar = "M",
                                             avar = "A",
                                             cvar = c("C")),
                         "M ~ A + C")
        })
        it("handles three cvar by adding all", {
            expect_equal(string_mreg_formula(mvar = "M",
                                             avar = "A",
                                             cvar = c("C1","C2","C3")),
                         "M ~ A + C1 + C2 + C3")
        })
    })
    ##
    describe("string_mreg_formula (bad args)", {
        it("throws an error on NULL mvar", {
            expect_error(string_mreg_formula(mvar = NULL,
                                             avar = "A",
                                             cvar = NULL))
        })
        it("throws an error on NULL avar", {
            expect_error(string_mreg_formula(mvar = "M",
                                             avar = NULL,
                                             cvar = NULL))
        })
    })
})


###
### Internal function for yreg string formula creation
################################################################################

test_that("string_yreg_formula create sound formula strings for non-survival outcomes", {

    ## Zero covariates
    expect_equal(string_yreg_formula("Y","A","M",NULL, interaction = FALSE,
                                     eventvar = NULL),
                 "Y ~ A + M")
    expect_equal(string_yreg_formula("Y","A","M",NULL, interaction = TRUE,
                                     eventvar = NULL),
                 "Y ~ A*M")

    ## One covariate
    expect_equal(string_yreg_formula("Y","A","M",c("C"), interaction = FALSE,
                                     eventvar = NULL),
                 "Y ~ A + M + C")
    expect_equal(string_yreg_formula("Y","A","M",c("C"), interaction = TRUE,
                                     eventvar = NULL),
                 "Y ~ A*M + C")

    ## Three covariates
    expect_equal(string_yreg_formula("Y","A","M",c("C1","C2","C3"), interaction = FALSE,
                                     eventvar = NULL),
                 "Y ~ A + M + C1 + C2 + C3")
    expect_equal(string_yreg_formula("Y","A","M",c("C1","C2","C3"), interaction = TRUE,
                                     eventvar = NULL),
                 "Y ~ A*M + C1 + C2 + C3")
})


test_that("string_yreg_formula create sound formula strings for survival outcomes", {

    ## Zero covariates
    expect_equal(string_yreg_formula("time","A","M",NULL, interaction = FALSE,
                                     eventvar = "event"),
                 "Surv(time, event) ~ A + M")
    expect_equal(string_yreg_formula("time","A","M",NULL, interaction = TRUE,
                                     eventvar = "event"),
                 "Surv(time, event) ~ A*M")

    ## One covariate
    expect_equal(string_yreg_formula("time","A","M",c("C"), interaction = FALSE,
                                     eventvar = "event"),
                 "Surv(time, event) ~ A + M + C")
    expect_equal(string_yreg_formula("time","A","M",c("C"), interaction = TRUE,
                                     eventvar = "event"),
                 "Surv(time, event) ~ A*M + C")

    ## Three covariates
    expect_equal(string_yreg_formula("time","A","M",c("C1","C2","C3"), interaction = FALSE,
                                     eventvar = "event"),
                 "Surv(time, event) ~ A + M + C1 + C2 + C3")
    expect_equal(string_yreg_formula("time","A","M",c("C1","C2","C3"), interaction = TRUE,
                                     eventvar = "event"),
                 "Surv(time, event) ~ A*M + C1 + C2 + C3")
})
