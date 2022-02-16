################################################################################
### Tests for internal functions
##
## Created on: 2020-03-10
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)
library(survival)
library(tidyverse)


###
### Internal function for yreg model fitting (logistic)
################################################################################

describe("fit_yreg Cox (no interaction)", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               ## Combine transplant and death for testing purpose
               status = if_else(status == 0, 0L, 1L))

    it("fits a correct model with no covariates", {
        ## No covariates
        yreg_fit0 <- fit_yreg(yreg = "survCox",
                              data = pbc_cc,
                              yvar = "time",
                              avar = "trt",
                              mvar = "bili",
                              cvar = NULL,
                              emm_ac_yreg = NULL,
                              emm_mc_yreg = NULL,
                              interaction = FALSE,
                              eventvar = "status")
        ref_fit0 <- coxph(formula = Surv(time,status) ~ trt + bili,
                          data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit0),
                     class(ref_fit0))
        ## Same formula
        expect_equal(as.character(yreg_fit0$call$formula),
                     as.character(ref_fit0$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit0),
                     coef(ref_fit0))
        ## Same vcov
        expect_equal(vcov(yreg_fit0),
                     vcov(ref_fit0))
    })

    it("fits a correct model with one covariate", {
        ## One covariates
        yreg_fit1 <- fit_yreg(yreg = "survCox",
                              data = pbc_cc,
                              yvar = "time",
                              avar = "trt",
                              mvar = "bili",
                              cvar = c("age"),
                              emm_ac_yreg = NULL,
                              emm_mc_yreg = NULL,
                              interaction = FALSE,
                              eventvar = "status")
        ref_fit1 <- coxph(formula = Surv(time,status) ~ trt + bili + age,
                          data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit1),
                     class(ref_fit1))
        ## Same formula
        expect_equal(as.character(yreg_fit1$call$formula),
                     as.character(ref_fit1$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit1),
                     coef(ref_fit1))
        ## Same vcov
        expect_equal(vcov(yreg_fit1),
                     vcov(ref_fit1))
    })

    it("fits a correct model with three covariates", {
        ## Three covariates
        yreg_fit3 <- fit_yreg(yreg = "survCox",
                              data = pbc_cc,
                              yvar = "time",
                              avar = "trt",
                              mvar = "bili",
                              cvar = c("age","male","stage"),
                              emm_ac_yreg = NULL,
                              emm_mc_yreg = NULL,
                              interaction = FALSE,
                              eventvar = "status")
        ref_fit3 <- coxph(formula = Surv(time,status) ~ trt + bili + age + male + stage,
                          data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit3),
                     class(ref_fit3))
        ## Same formula
        expect_equal(as.character(yreg_fit3$call$formula),
                     as.character(ref_fit3$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit3),
                     coef(ref_fit3))
        ## Same vcov
        expect_equal(vcov(yreg_fit3),
                     vcov(ref_fit3))
    })
    
    # only test when emm_ac_yreg and emm_mc_yreg are both not null:
    it("fits a correct model with three covariates", {
        ## Three covariates
        yreg_fit6 <- fit_yreg(yreg = "survCox",
                              data = pbc_cc,
                              yvar = "time",
                              avar = "trt",
                              mvar = "bili",
                              cvar = c("age","male","stage"),
                              emm_ac_yreg = c("age"),
                              emm_mc_yreg = c("male", "stage"),
                              interaction = FALSE,
                              eventvar = "status")
        ref_fit6 <- coxph(formula = Surv(time,status) ~ trt + bili + age + male + stage +
                              trt:age + bili:male + bili:stage,
                          data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit6),
                     class(ref_fit6))
        ## Same formula
        expect_equal(as.character(yreg_fit6$call$formula),
                     as.character(ref_fit6$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit6),
                     coef(ref_fit6))
        ## Same vcov
        expect_equal(vcov(yreg_fit6),
                     vcov(ref_fit6))
    })

})


describe("fit_yreg Cox (interaction)", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               ## Combine transplant and death for testing purpose
               status = if_else(status == 0, 0L, 1L))

    it("fits a correct model with no covariates", {
        ## No covariates
        yreg_fit0 <- fit_yreg(yreg = "survCox",
                              data = pbc_cc,
                              yvar = "time",
                              avar = "trt",
                              mvar = "bili",
                              cvar = NULL,
                              emm_ac_yreg = NULL,
                              emm_mc_yreg = NULL,
                              interaction = TRUE,
                              eventvar = "status")
        ref_fit0 <- coxph(formula = Surv(time,status) ~ trt + bili + trt:bili,
                          data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit0),
                     class(ref_fit0))
        ## Same formula
        expect_equal(as.character(yreg_fit0$call$formula),
                     as.character(ref_fit0$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit0),
                     coef(ref_fit0))
        ## Same vcov
        expect_equal(vcov(yreg_fit0),
                     vcov(ref_fit0))
    })

    it("fits a correct model with one covariate", {
        ## One covariates
        yreg_fit1 <- fit_yreg(yreg = "survCox",
                              data = pbc_cc,
                              yvar = "time",
                              avar = "trt",
                              mvar = "bili",
                              cvar = c("age"),
                              emm_ac_yreg = NULL,
                              emm_mc_yreg = NULL,
                              interaction = TRUE,
                              eventvar = "status")
        ref_fit1 <- coxph(formula = Surv(time,status) ~ trt + bili + trt:bili + age,
                          data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit1),
                     class(ref_fit1))
        ## Same formula
        expect_equal(as.character(yreg_fit1$call$formula),
                     as.character(ref_fit1$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit1),
                     coef(ref_fit1))
        ## Same vcov
        expect_equal(vcov(yreg_fit1),
                     vcov(ref_fit1))
    })

    it("fits a correct model with three covariates", {
        ## Three covariates
        yreg_fit3 <- fit_yreg(yreg = "survCox",
                              data = pbc_cc,
                              yvar = "time",
                              avar = "trt",
                              mvar = "bili",
                              cvar = c("age","male","stage"),
                              emm_ac_yreg = NULL,
                              emm_mc_yreg = NULL,
                              interaction = TRUE,
                              eventvar = "status")
        ref_fit3 <- coxph(formula = Surv(time,status) ~ trt + bili + trt:bili + age + male + stage,
                          data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit3),
                     class(ref_fit3))
        ## Same formula
        expect_equal(as.character(yreg_fit3$call$formula),
                     as.character(ref_fit3$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit3),
                     coef(ref_fit3))
        ## Same vcov
        expect_equal(vcov(yreg_fit3),
                     vcov(ref_fit3))
    })
    
    # only test when emm_ac_yreg and emm_mc_yreg are both not null:
    it("fits a correct model with three covariates", {
        ## Three covariates
        yreg_fit6 <- fit_yreg(yreg = "survCox",
                              data = pbc_cc,
                              yvar = "time",
                              avar = "trt",
                              mvar = "bili",
                              cvar = c("age","male","stage"),
                              emm_ac_yreg = c("age"),
                              emm_mc_yreg = c("male", "stage"),
                              interaction = TRUE,
                              eventvar = "status")
        ref_fit6 <- coxph(formula = Surv(time,status) ~ trt + bili + trt:bili + age + male + stage +
                              trt:age + bili:male + bili:stage,
                          data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit6),
                     class(ref_fit6))
        ## Same formula
        expect_equal(as.character(yreg_fit6$call$formula),
                     as.character(ref_fit6$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit6),
                     coef(ref_fit6))
        ## Same vcov
        expect_equal(vcov(yreg_fit6),
                     vcov(ref_fit6))
    })

})
