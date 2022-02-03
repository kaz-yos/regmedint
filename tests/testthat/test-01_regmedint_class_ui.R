################################################################################
### Tests for regmedint
##
## Created on: 2020-03-16
## Author: Kazuki Yoshida
################################################################################

library(testthat)
library(survival)
library(tidyverse)


###
### regmedint user interface
################################################################################

describe("regmedint", {

    data(pbc)

    describe("na_omit", {
        ## Create a numeric dataset with NAs
        pbc_na <- pbc %>%
            as_tibble() %>%
            mutate(male = if_else(sex == "m", 1L, 0L),
                   ## Combine transplant and death for testing purpose
                   status = if_else(status == 0, 0L, 1L),
                   ##
                   event = if_else(status == 1L, 1L, 0L),
                   ## Binary mvar
                   bili_bin = if_else(bili > median(bili), 1L, 0L),
                   alk_phos = alk.phos)

        describe("na_omit (mreg linear, yreg linear)", {
            ## Specify variables of interest
            yvar <- "alk.phos"
            avar <- "trt"
            mvar <- "bili_bin"
            cvar <- c("age","male","stage")
            ## Gold standard sample size
            expected_sample_size <- nrow(pbc_na) -
                sum(rowSums(is.na(pbc_na[,c(yvar, avar, mvar, cvar)])) > 0)

            fit_regmedint_obj <- regmedint(data = pbc_na,
                                           yvar = yvar,
                                           avar = avar,
                                           mvar = mvar,
                                           cvar = cvar,
                                           a0 = 1,
                                           a1 = 2,
                                           m_cde = 0,
                                           c_cond = c(50,2,4),
                                           mreg = "linear",
                                           yreg = "linear",
                                           interaction = FALSE,
                                           casecontrol = FALSE,
                                           ## No need for eventvar
                                           na_omit = TRUE)
            out_regmedint <- summary(fit_regmedint_obj)

            ## Sample size extracted from the models
            mreg_sample_size <- out_regmedint$summary_mreg_fit$df[2] +
                nrow(out_regmedint$summary_mreg_fit$coefficients)
            yreg_sample_size <- out_regmedint$summary_yreg_fit$df[2] +
                nrow(out_regmedint$summary_yreg_fit$coefficients)

            ## Check sample sizes of both mediator and outcome models
            it("drops rows containing NAs in variables of interest when TRUE", {
                expect_equal(mreg_sample_size, expected_sample_size)
                expect_equal(yreg_sample_size, expected_sample_size)
            })
        })

        describe("na_omit (mreg linear, yreg survCox)", {
            ## Specify variables of interest
            yvar = "time"
            avar = "trt"
            mvar = "bili_bin"
            cvar = c("age")
            eventvar = "status"

            ## Gold standard sample size
            expected_sample_size <- nrow(pbc_na) -
                sum(rowSums(is.na(pbc_na[,c(yvar, avar, mvar, cvar, eventvar)])) > 0)

            fit_regmedint_obj <- regmedint(data = pbc_na,
                                           yvar = yvar,
                                           avar = avar,
                                           mvar = mvar,
                                           cvar = cvar,
                                           eventvar = eventvar,
                                           a0 = 1,
                                           a1 = 2,
                                           m_cde = 0,
                                           c_cond = 40,
                                           mreg = "linear",
                                           yreg = "survCox",
                                           interaction = FALSE,
                                           casecontrol = FALSE,
                                           na_omit = TRUE)
            out_regmedint <- summary(fit_regmedint_obj)

            ## Sample size extracted from the models
            mreg_sample_size <- out_regmedint$summary_mreg_fit$df[2] +
                nrow(out_regmedint$summary_mreg_fit$coefficients)
            ## summary.coxph object does not have df.
            yreg_sample_size <- out_regmedint$summary_yreg_fit$n

            ## Check sample sizes of both mediator and outcome models
            it("drops rows containing NAs in variables of interest when NA", {
                expect_equal(mreg_sample_size, expected_sample_size)
                expect_equal(yreg_sample_size, expected_sample_size)
            })
        })
    })


    pbc_cc <- pbc %>%
        as_tibble() %>%
        ## Missing data should be warned in validate_args()
        drop_na() %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               ## Combine transplant and death for testing purpose
               status = if_else(status == 0, 0L, 1L),
               ##
               event = if_else(status == 1L, 1L, 0L),
               ## Binary mvar
               bili_bin = if_else(bili > median(bili), 1L, 0L), # changed from bili_bin
               alk_phos = alk.phos)

    describe("validate_args (regmedint argument validation)", {
        it("rejects missing data in the variales of interest", {
            ## This is a regexp. () needs escapes.
            msg_missing <- "Missingness is not allowed in variables of interest!
For multiple imputation, see the multiple imputation vignette: vignette\\(\"vig_04_mi\"\\)
To perform complete case analysis, use na_omit = TRUE.
Variables with missingness: "
            expect_error(pbc_cc %>%
                         mutate(alk_phos = NA) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         paste0(msg_missing, "alk_phos"))
            expect_error(pbc_cc %>%
                         mutate(trt = NA) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         paste0(msg_missing, "trt"))
            expect_error(pbc_cc %>%
                         mutate(bili = NA) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         paste0(msg_missing, "bili"))
            expect_error(pbc_cc %>%
                         mutate(age = NA) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         paste0(msg_missing, "age"))
            expect_error(pbc_cc %>%
                         mutate(age = NA,
                                male = NA) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         paste0(msg_missing, "age, male"))
        })
        ##
        it("allows missing data in unused variables", {
            expect_equal(pbc_cc %>%
                         mutate(sex = NA) %>%
                         regmedint(data = .,
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
                                   eventvar = NULL),
                         regmedint(data = pbc_cc,
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
                                   eventvar = NULL))
        })
        ##
        it("rejects factor variables", {
            msg_factor <- "Factors are not allowed! Use numeric variables only."
            expect_error(pbc_cc %>%
                         mutate(alk_phos = factor(alk_phos)) %>%
                         regmedint(data = ,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_factor)
            expect_error(pbc_cc %>%
                         mutate(trt = factor(trt)) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_factor)
            expect_error(pbc_cc %>%
                         mutate(bili = factor(bili)) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_factor)
            expect_error(pbc_cc %>%
                         mutate(status = factor(status)) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = "status"),
                         msg_factor)
        })
        ##
        it("rejects missing in evaluation arguments", {
            msg_na_assertion_c_cond <- "are not true"
            msg_na_assertion <- "is not TRUE"
            expect_error(pbc_cc %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = as.numeric(NA),
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_na_assertion)
            expect_error(pbc_cc %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = as.numeric(NA),
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_na_assertion)
            expect_error(pbc_cc %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = as.numeric(NA),
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_na_assertion)
            expect_error(pbc_cc %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,NA),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_na_assertion_c_cond)
        })
    })

    describe("regmedint mreg linear yreg linear", {
        it("runs with zero cvar with no interaction", {
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
        it("runs with zero cvar with no interaction", {
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
