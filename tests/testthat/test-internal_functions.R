################################################################################
### Tests for internal functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)


###
### Internal functions for string formula creation
################################################################################

test_that("string_mreg_formula create sound formula strings", {

    expect_equal(string_mreg_formula("M","A",NULL),
                 "M ~ A")

    expect_equal(string_mreg_formula("M","A",c("C")),
                 "M ~ A + C")

    expect_equal(string_mreg_formula("M","A",c("C1","C2","C3")),
                 "M ~ A + C1 + C2 + C3")

})


test_that("string_yreg_formula create sound formula strings for non-survival outcomes", {

    ## Zero covariates
    expect_equal(string_yreg_formula("Y","A","M",NULL, interaction = FALSE),
                 "Y ~ A + M")
    expect_equal(string_yreg_formula("Y","A","M",NULL, interaction = TRUE),
                 "Y ~ A*M")

    ## One covariate
    expect_equal(string_yreg_formula("Y","A","M",c("C"), interaction = FALSE),
                 "Y ~ A + M + C")
    expect_equal(string_yreg_formula("Y","A","M",c("C"), interaction = TRUE),
                 "Y ~ A*M + C")

    ## Three covariates
    expect_equal(string_yreg_formula("Y","A","M",c("C1","C2","C3"),
                                     interaction = FALSE),
                 "Y ~ A + M + C1 + C2 + C3")
    expect_equal(string_yreg_formula("Y","A","M",c("C1","C2","C3"),
                                     interaction = TRUE),
                 "Y ~ A*M + C1 + C2 + C3")
})

test_that("string_yreg_formula create sound formula strings for survival outcomes", {

    ## Zero covariates
    expect_equal(string_yreg_formula("time","A","M",NULL,"event", interaction = FALSE),
                 "Surv(time, event) ~ A + M")
    expect_equal(string_yreg_formula("time","A","M",NULL,"event", interaction = TRUE),
                 "Surv(time, event) ~ A*M")

    ## One covariate
    expect_equal(string_yreg_formula("time","A","M",c("C"),"event", interaction = FALSE),
                 "Surv(time, event) ~ A + M + C")
    expect_equal(string_yreg_formula("time","A","M",c("C"),"event", interaction = TRUE),
                 "Surv(time, event) ~ A*M + C")

    ## Three covariates
    expect_equal(string_yreg_formula("time","A","M",c("C1","C2","C3"),"event",
                                     interaction = FALSE),
                 "Surv(time, event) ~ A + M + C1 + C2 + C3")
    expect_equal(string_yreg_formula("time","A","M",c("C1","C2","C3"),"event",
                                     interaction = TRUE),
                 "Surv(time, event) ~ A*M + C1 + C2 + C3")
})
