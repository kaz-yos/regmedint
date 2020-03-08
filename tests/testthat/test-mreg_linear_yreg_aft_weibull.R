################################################################################
### Tests
##
## Created on: 2020-03-08
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)
library(survival)


###
### Set up SAS reference results
################################################################################


###
### Fit separate models to obtain reference R results
################################################################################

mreg_fit <- NA
## lm(formula = ,
##    data = )

yreg_fit <- NA
## survreg(formula = ,
##         data = ,
##         dist = "weibull")


###
### Run mediation analysis of interest
################################################################################

med_fit <- NULL


###
### Test against reference results
################################################################################

test_that("the regmedint object structure is sound", {

    expect_equal(class(med_fit$mreg), class(mreg_fit))
    expect_equal(class(med_fit$yreg), class(yreg_fit))

    ## FIXME: Drop once tests are written.
    expect_true(FALSE)

})


test_that("point estimates are compatible with SAS results", {

    expect_equal(med_fit$med$te, "total effect point estimate from SAS")


    ## FIXME: Drop once tests are written.
    expect_true(FALSE)

})


test_that("variance estimates are compatible with SAS results", {

    expect_equal(med_fit$med$te, "total effect variance estimate from SAS")

    ## FIXME: Drop once tests are written.
    expect_true(FALSE)

})
