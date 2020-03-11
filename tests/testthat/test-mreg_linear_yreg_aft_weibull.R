################################################################################
### Tests
##
## Created on: 2020-03-08
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)
library(survival)
library(tidyverse)

source("./utilities_for_tests.R")

###
### Set up SAS reference results
################################################################################

sas_res <- read_parsed_sas_mediation_output(
    "../reference_results/sas-mreg_linear_yreg_aft_weibull.txt")


###
### Fit separate models to obtain reference R results
################################################################################

data1 <- read.delim(file = "../reference_results/data-valeri-vanderweele-2015.txt",
                    header = TRUE,
                    sep = " ")

mreg_fit <-
    lm(formula = m ~ x + c,
       data = data1)

yreg_fit <-
    survreg(formula = Surv(y, cens) ~ x*m + c,
            data = data1,
            dist = "weibull")


###
### Run mediation analysis of interest
################################################################################

med_fit <- regmedint(data = data1,
                     yvar = "y",
                     avar = "x",
                     mvar = "m",
                     cvar = "c",
                     a0 = 0,
                     a1 = 1,
                     m_cde = 0,
                     yreg = "survAFT_weibull",
                     mreg = "logistic",
                     interaction = TRUE,
                     casecontrol = FALSE,
                     full_output = FALSE,
                     c_cond = NULL,
                     boot = FALSE,
                     eventvar = "cens")

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
