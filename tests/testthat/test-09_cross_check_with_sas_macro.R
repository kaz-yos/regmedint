################################################################################
### Test against SAS macro results
##
## Created on: 2020-03-24
## Author: Kazuki Yoshida
################################################################################

library(survival)
library(testthat)
library(tidyverse)

source("./utilities_for_tests.R")


###
### Placeholder test to turn on the display of errors
################################################################################

test_that("placeholder", {
    expect_true(TRUE)
})


###
### Prepare data
################################################################################

data(pbc)

## For the purpose of this cross testing, complete case analysis is fine.
pbc_cc <- pbc %>%
    as_tibble() %>%
    ## Missing data should be warned in validate_args()
    drop_na() %>%
    mutate(male = if_else(sex == "m", 1L, 0L),
           ## Combine transplant and death for testing purpose
           status = if_else(status == 0, 0L, 1L),
           ## censoring status reverse coded for SAS macro
           cens = if_else(status == 1L, 0L, 1L),
           ## Binary mvar
           bili_bin = if_else(bili > median(bili), 1L, 0L),
           ## Fake count yvar
           edema = 2 * edema,
           alk_phos = alk.phos) %>%
    select(
        ## avar
        trt,
        ##
        ## mvar (continuous)
        bili,
        ## mvar (binary)
        bili_bin,
        ##
        ## yvar (continuous)
        alk_phos,
        ## yvar (binary)
        spiders,
        ## yvar (fake count)
        edema,
        ## yvar (survival)
        time,
        ## eventvar (survival)
        status,
        ## censvar (survival)
        cens,
        ##
        ## cvar (continuous/binary/handled continuous)
        age, male, stage
    )


###
### Prepare to cover all available patterns
################################################################################

## Copied from tests/reference_results/02_generate_sas_macro_calls.R and modified
macro_args <-
    crossing(mreg = c("linear",
                      "logistic"),
             yreg = c("linear",
                      "logistic",
                      "loglinear",
                      "poisson",
                      "negbin",
                      "survCox",
                      "survAFT_exp",
                      "survAFT_weibull"),
             interaction = c("false","true"),
             casecontrol = c("false","true"),
             ## empty cvar may not work with output = full
             ncvar = c(3)) %>%
    ## casecontrol is only relevant for binary yreg
    filter(casecontrol == "false" |
           (casecontrol == "true" & (yreg == "logistic" | yreg == "loglinear"))) %>%
    ## Choose appropriate variable names
    mutate(avar = "trt",
           ##
           mvar = case_when(mreg == "linear" ~ "bili",
                            mreg == "logistic" ~ "bili_bin"),
           ##
           yvar = case_when(yreg == "linear" ~ "alk_phos",
                            ##
                            yreg == "logistic" ~ "spiders",
                            yreg == "loglinear" ~ "spiders",
                            ##
                            yreg == "poisson" ~ "edema",
                            yreg == "negbin" ~ "edema",
                            ##
                            yreg == "survCox" ~ "time",
                            yreg == "survAFT_exp" ~ "time",
                            yreg == "survAFT_weibull" ~ "time"),
           ##
           eventvar = case_when(yreg == "survCox" ~ "status",
                                yreg == "survAFT_exp" ~ "status",
                                yreg == "survAFT_weibull" ~ "status",
                                TRUE ~ ""),
           cens = "cens",
           ##
           cvar = case_when(ncvar == 0 ~ "",
                            ncvar == 3 ~ "age male stage"),
           ##
           c_cond = case_when(ncvar == 0 ~ "",
                              ncvar == 3 ~ "50 1 2"),
           a0 = "1",
           a1 = "2",
           m_cde = "1") %>%
    ## File name
    mutate(filename = sprintf("sas-mreg_%s_yreg_%s_int_%s_caco_%s_ncvar%s.sas",
                              ## All variables in crossing() should appear hear.
                              mreg,
                              yreg,
                              str_sub(interaction, start = 1, end = 1),
                              str_sub(casecontrol, start = 1, end = 1),
                              ncvar))


###
### Add SAS results
################################################################################

macro_args_sas <- macro_args %>%
        ## Load SAS results in comparable form
    mutate(sas = map(filename, function(filename) {
        res_filename <- stringr::str_replace_all(filename, "sas$", "txt")
        ## The working directory here is the enclosing folder ./tests/testthat/.
        ## print(getwd())
        relpath <- paste0("../reference_results/", res_filename)
        read_parsed_sas_mediation_output(relpath)
    }))

junk <- macro_args_sas %>%
    mutate(junk = map2(filename, sas, function(filename, sas) {
        res_filename <- stringr::str_replace_all(filename, "sas$", "txt")
        test_that(paste0(res_filename, " was extracted successfully"), {
            expect_equal(class(sas)[1], "tbl_df")
            expect_true(ncol(sas) %in% c(5,6))
        })
    }))


###
### Add R results
################################################################################

macro_args_sas_r_prelim <- macro_args_sas %>%
    ## Run R analyses
    mutate(res = pmap(
               list(yvar,
                    avar,
                    mvar,
                    cvar,
                    a0,
                    a1,
                    m_cde,
                    c_cond,
                    mreg,
                    yreg,
                    interaction,
                    casecontrol,
                    eventvar),
               function(yvar,
                        avar,
                        mvar,
                        cvar,
                        a0,
                        a1,
                        m_cde,
                        c_cond,
                        mreg,
                        yreg,
                        interaction,
                        casecontrol,
                        eventvar) {
                   ##
                   if (cvar == "") {
                       if (eventvar == "") {
                           try(regmedint(data = pbc_cc,
                                         yvar = yvar,
                                         avar = avar,
                                         mvar = mvar,
                                         cvar = NULL,
                                         a0 = as.numeric(a0),
                                         a1 = as.numeric(a1),
                                         m_cde = as.numeric(m_cde),
                                         c_cond = NULL,
                                         mreg = mreg,
                                         yreg = yreg,
                                         interaction = if_else(interaction == "true", TRUE, FALSE),
                                         casecontrol = if_else(casecontrol == "true", TRUE, FALSE),
                                         eventvar = NULL))
                       } else {
                           try(regmedint(data = pbc_cc,
                                         yvar = yvar,
                                         avar = avar,
                                         mvar = mvar,
                                         cvar = NULL,
                                         a0 = as.numeric(a0),
                                         a1 = as.numeric(a1),
                                         m_cde = as.numeric(m_cde),
                                         c_cond = NULL,
                                         mreg = mreg,
                                         yreg = yreg,
                                         interaction = if_else(interaction == "true", TRUE, FALSE),
                                         casecontrol = if_else(casecontrol == "true", TRUE, FALSE),
                                         eventvar = eventvar))
                       }
                   } else {
                       cvar <- str_split(cvar, " ")[[1]]
                       c_cond <- str_split(c_cond, " ")[[1]]
                       if (eventvar == "") {
                           try(regmedint(data = pbc_cc,
                                         yvar = yvar,
                                         avar = avar,
                                         mvar = mvar,
                                         cvar = cvar,
                                         a0 = as.numeric(a0),
                                         a1 = as.numeric(a1),
                                         m_cde = as.numeric(m_cde),
                                         c_cond = as.numeric(c_cond),
                                         mreg = mreg,
                                         yreg = yreg,
                                         interaction = if_else(interaction == "true", TRUE, FALSE),
                                         casecontrol = if_else(casecontrol == "true", TRUE, FALSE),
                                         eventvar = NULL))
                       } else {
                           try(regmedint(data = pbc_cc,
                                         yvar = yvar,
                                         avar = avar,
                                         mvar = mvar,
                                         cvar = cvar,
                                         a0 = as.numeric(a0),
                                         a1 = as.numeric(a1),
                                         m_cde = as.numeric(m_cde),
                                         c_cond = as.numeric(c_cond),
                                         mreg = mreg,
                                         yreg = yreg,
                                         interaction = if_else(interaction == "true", TRUE, FALSE),
                                         casecontrol = if_else(casecontrol == "true", TRUE, FALSE),
                                         eventvar = eventvar))
                       }
                   }
                   ##
               }))


###
### Extract R results for testing
################################################################################

macro_args_sas_r <- macro_args_sas_r_prelim %>%
    ## Extract useful elements when available
    mutate(
        coef = map(res, function(res) {
            if(class(res) == "try-error") {
                return(NULL)
            } else {
                return(coef(res))
            }
        }),
        p = map(res, function(res) {
            if(class(res) == "try-error") {
                return(NULL)
            } else {
                capture_output(summary_mat <- summary(res))
                return(summary_mat[,"p"])
            }
        }),
        ## 2 * (1 - pnorm(1.96)) to get confint corresponding to 1.96 * se
        lower = map(res, function(res) {
            if(class(res) == "try-error") {
                return(NULL)
            } else {
                return(confint(res, alpha = 2 * (1 - pnorm(1.96)))[,"lower"])
            }
        }),
        upper = map(res, function(res) {
            if(class(res) == "try-error") {
                return(NULL)
            } else {
                return(confint(res, alpha = 2 * (1 - pnorm(1.96)))[,"upper"])
            }
        })
    )


###
### Tests
################################################################################

## This stops on the first error when run interactively.
## It continues when running
for (i in seq_len(nrow(macro_args_sas_r))) {
    ##
    title <- sprintf("Index %d: File %s equivalent fits ok", i, macro_args_sas_r$filename[i])
    res <- macro_args_sas_r$res[[i]]
    ##
    test_that(title, {
        expect_equal(class(res)[[1]], "regmedint")
    })
}


## Result comparison
junk <- macro_args_sas_r %>%
    mutate(
        junk = pmap(
            list(filename, sas, res, coef, p, lower, upper),
            function(filename, sas, res, coef, p, lower, upper) {

                if (class(res) == "try-error") {
                    stop(paste0("R fit for ", filename, " gave an try-error object!"))
                } else {

                    ## Create named vectors
                    sas_estimate <- sas$estimate %>% `names<-`(sas$effect)
                    sas_p <- sas$p %>% `names<-`(sas$effect)
                    sas_lower <- sas$lower %>% `names<-`(sas$effect)
                    sas_upper <- sas$upper %>% `names<-`(sas$effect)

                    ## [["cde"]] to unname it.
                    if (nrow(sas) == 12) {

                        ## Full result case with interaction
                        if (ncol(sas) == 6) {
                            ## yreg linear with se column
                            test_that(paste0("coef match with ", filename), {
                                expect_equal(coef[["cde"]], sas_estimate[["conditional_cde"]])
                                expect_equal(coef[["pnde"]], sas_estimate[["conditional_pnde"]])
                                expect_equal(coef[["tnie"]], sas_estimate[["conditional_tnie"]])
                                expect_equal(coef[["tnde"]], sas_estimate[["conditional_tnde"]])
                                expect_equal(coef[["pnie"]], sas_estimate[["conditional_pnie"]])
                                expect_equal(coef[["te"]], sas_estimate[["conditional_te"]])
                            })

                            test_that(paste0("p match with ", filename), {
                                expect_equal(p[["cde"]], sas_p[["conditional_cde"]])
                                expect_equal(p[["pnde"]], sas_p[["conditional_pnde"]])
                                expect_equal(p[["tnie"]], sas_p[["conditional_tnie"]])
                                expect_equal(p[["tnde"]], sas_p[["conditional_tnde"]])
                                expect_equal(p[["pnie"]], sas_p[["conditional_pnie"]])
                                expect_equal(p[["te"]], sas_p[["conditional_te"]])
                            })

                            test_that(paste0("lower confint match with ", filename), {
                                expect_equal(lower[["cde"]], sas_lower[["conditional_cde"]])
                                expect_equal(lower[["pnde"]], sas_lower[["conditional_pnde"]])
                                expect_equal(lower[["tnie"]], sas_lower[["conditional_tnie"]])
                                expect_equal(lower[["tnde"]], sas_lower[["conditional_tnde"]])
                                expect_equal(lower[["pnie"]], sas_lower[["conditional_pnie"]])
                                expect_equal(lower[["te"]], sas_lower[["conditional_te"]])
                            })

                            test_that(paste0("upper confint match with ", filename), {
                                expect_equal(upper[["cde"]], sas_upper[["conditional_cde"]])
                                expect_equal(upper[["pnde"]], sas_upper[["conditional_pnde"]])
                                expect_equal(upper[["tnie"]], sas_upper[["conditional_tnie"]])
                                expect_equal(upper[["tnde"]], sas_upper[["conditional_tnde"]])
                                expect_equal(upper[["pnie"]], sas_upper[["conditional_pnie"]])
                                expect_equal(upper[["te"]], sas_upper[["conditional_te"]])
                            })
                        } else {
                            ## yreg non-linear with exp(coef), p, exp(lower), exp(upper)
                            test_that(paste0("coef match with ", filename), {
                                expect_equal(exp[(coef["cde"])], sas_estimate[["conditional_cde"]])
                                expect_equal(exp[(coef["pnde"])], sas_estimate[["conditional_pnde"]])
                                expect_equal(exp[(coef["tnie"])], sas_estimate[["conditional_tnie"]])
                                expect_equal(exp[(coef["tnde"])], sas_estimate[["conditional_tnde"]])
                                expect_equal(exp[(coef["pnie"])], sas_estimate[["conditional_pnie"]])
                                expect_equal(exp[(coef["te"])], sas_estimate[["conditional_te"]])
                            })

                            test_that(paste0("p match with ", filename), {
                                expect_equal(p[["cde"]], sas_p[["conditional_cde"]])
                                expect_equal(p[["pnde"]], sas_p[["conditional_pnde"]])
                                expect_equal(p[["tnie"]], sas_p[["conditional_tnie"]])
                                expect_equal(p[["tnde"]], sas_p[["conditional_tnde"]])
                                expect_equal(p[["pnie"]], sas_p[["conditional_pnie"]])
                                expect_equal(p[["te"]], sas_p[["conditional_te"]])
                            })

                            test_that(paste0("lower confint match with ", filename), {
                                expect_equal(exp[(lower["cde"])], sas_lower[["conditional_cde"]])
                                expect_equal(exp[(lower["pnde"])], sas_lower[["conditional_pnde"]])
                                expect_equal(exp[(lower["tnie"])], sas_lower[["conditional_tnie"]])
                                expect_equal(exp[(lower["tnde"])], sas_lower[["conditional_tnde"]])
                                expect_equal(exp[(lower["pnie"])], sas_lower[["conditional_pnie"]])
                                expect_equal(exp[(lower["te"])], sas_lower[["conditional_te"]])
                            })

                            test_that(paste0("upper confint match with ", filename), {
                                expect_equal(exp[(upper["cde"])], sas_upper[["conditional_cde"]])
                                expect_equal(exp[(upper["pnde"])], sas_upper[["conditional_pnde"]])
                                expect_equal(exp[(upper["tnie"])], sas_upper[["conditional_tnie"]])
                                expect_equal(exp[(upper["tnde"])], sas_upper[["conditional_tnde"]])
                                expect_equal(exp[(upper["pnie"])], sas_upper[["conditional_pnie"]])
                                expect_equal(exp[(upper["te"])], sas_upper[["conditional_te"]])
                            })
                        }

                    } else {

                        ## Reduced result case with no interaction
                        if (ncol(sas) == 6) {
                            ## yreg linear with se column
                            test_that(paste0("coef match with ", filename), {
                                expect_equal(coef[["cde"]], sas_estimate[["cde=nde"]])
                                expect_equal(coef[["pnde"]], sas_estimate[["cde=nde"]])
                                expect_equal(coef[["tnie"]], sas_estimate[["nie"]])
                                expect_equal(coef[["tnde"]], sas_estimate[["cde=nde"]])
                                expect_equal(coef[["pnie"]], sas_estimate[["nie"]])
                                expect_equal(coef[["te"]], sas_estimate[["te"]])
                            })

                            test_that(paste0("p match with ", filename), {
                                expect_equal(p[["cde"]], sas_p[["cde=nde"]])
                                expect_equal(p[["pnde"]], sas_p[["cde=nde"]])
                                expect_equal(p[["tnie"]], sas_p[["nie"]])
                                expect_equal(p[["tnde"]], sas_p[["cde=nde"]])
                                expect_equal(p[["pnie"]], sas_p[["nie"]])
                                expect_equal(p[["te"]], sas_p[["te"]])
                            })

                            test_that(paste0("lower confint match with ", filename), {
                                expect_equal(lower[["cde"]], sas_lower[["cde=nde"]])
                                expect_equal(lower[["pnde"]], sas_lower[["cde=nde"]])
                                expect_equal(lower[["tnie"]], sas_lower[["nie"]])
                                expect_equal(lower[["tnde"]], sas_lower[["cde=nde"]])
                                expect_equal(lower[["pnie"]], sas_lower[["nie"]])
                                expect_equal(lower[["te"]], sas_lower[["te"]])
                            })

                            test_that(paste0("upper confint match with ", filename), {
                                expect_equal(upper[["cde"]], sas_upper[["cde=nde"]])
                                expect_equal(upper[["pnde"]], sas_upper[["cde=nde"]])
                                expect_equal(upper[["tnie"]], sas_upper[["nie"]])
                                expect_equal(upper[["tnde"]], sas_upper[["cde=nde"]])
                                expect_equal(upper[["pnie"]], sas_upper[["nie"]])
                                expect_equal(upper[["te"]], sas_upper[["te"]])
                            })
                        } else {
                            ## yreg non-linear with exp(coef), p, exp(lower), exp(upper)
                            test_that(paste0("coef match with ", filename), {
                                expect_equal(exp[(coef["cde"])], sas_estimate[["cde=nde"]])
                                expect_equal(exp[(coef["pnde"])], sas_estimate[["cde=nde"]])
                                expect_equal(exp[(coef["tnie"])], sas_estimate[["nie"]])
                                expect_equal(exp[(coef["tnde"])], sas_estimate[["cde=nde"]])
                                expect_equal(exp[(coef["pnie"])], sas_estimate[["nie"]])
                                expect_equal(exp[(coef["te"])], sas_estimate[["te"]])
                            })

                            test_that(paste0("p match with ", filename), {
                                expect_equal(p[["cde"]], sas_p[["cde=nde"]])
                                expect_equal(p[["pnde"]], sas_p[["cde=nde"]])
                                expect_equal(p[["tnie"]], sas_p[["nie"]])
                                expect_equal(p[["tnde"]], sas_p[["cde=nde"]])
                                expect_equal(p[["pnie"]], sas_p[["nie"]])
                                expect_equal(p[["te"]], sas_p[["te"]])
                            })

                            test_that(paste0("lower confint match with ", filename), {
                                expect_equal(exp[(lower["cde"])], sas_lower[["cde=nde"]])
                                expect_equal(exp[(lower["pnde"])], sas_lower[["cde=nde"]])
                                expect_equal(exp[(lower["tnie"])], sas_lower[["nie"]])
                                expect_equal(exp[(lower["tnde"])], sas_lower[["cde=nde"]])
                                expect_equal(exp[(lower["pnie"])], sas_lower[["nie"]])
                                expect_equal(exp[(lower["te"])], sas_lower[["te"]])
                            })

                            test_that(paste0("upper confint match with ", filename), {
                                expect_equal(exp[(upper["cde"])], sas_upper[["cde=nde"]])
                                expect_equal(exp[(upper["pnde"])], sas_upper[["cde=nde"]])
                                expect_equal(exp[(upper["tnie"])], sas_upper[["nie"]])
                                expect_equal(exp[(upper["tnde"])], sas_upper[["cde=nde"]])
                                expect_equal(exp[(upper["pnie"])], sas_upper[["nie"]])
                                expect_equal(exp[(upper["te"])], sas_upper[["te"]])
                            })
                        }

                    }
                }
            }))
