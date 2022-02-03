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
### Adjust tolerance by masking the original expect_equal
################################################################################

## Need to remove from global_env() to make sure we are not overwriting.
if (rlang::env_has(rlang::global_env(), "expect_equal")) {
    rlang::env_unbind(rlang::global_env(), "expect_equal")
}

## scale = 1 to test on the absolute scale
expect_equal <- purrr::partial(expect_equal, tolerance = 0.009, scale = 1)


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
        ## yvar (count)
        platelet,
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

## Factored out to avoid inconsistency between SAS call and R call.
source("../reference_results/02_generate_sas_macro_calls_helpers.R")

macro_args <- generate_sas_macro_args()


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
            expect_equal(class(sas)[[1]], "tbl_df")
            expect_true(ncol(sas) %in% c(5,6))
        })
    }))


###
### Add R results
################################################################################

macro_args_sas_r_prelim <- macro_args_sas %>%
    ## Add cmean vectors (R re-evaluation results can be validated agains "marginals").
    mutate(cmean = map(cvar, function(cvar) {
        if (cvar == "") {
            return(NULL)
        } else {
            ## We need a character vector.
            cvar_vec <- str_split(cvar, " ")[[1]]
            return(colMeans(pbc_cc[,cvar_vec]))
        }
    })) %>%
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
### Write R results to files
################################################################################

junk <- macro_args_sas_r_prelim %>%
    mutate(filename_r_res = stringr::str_replace_all(filename, "sas$", "r.txt")) %>%
    mutate(junk = pmap(list(res, filename_r_res, cmean), function(res, filename_r_res, cmean) {
        ## Create a textual representation of the results
        res_char_vec <- capture.output(summary(res, exponentiate = TRUE))
        res_char <- paste0(paste0(res_char_vec, collapse = "\n"), "\n")
        ## Create a textual representation of the results at cmean when avaialble
        if (!is.null(cmean)) {
            res_char_vec_cmean <- capture.output(coef(summary(res, c_cond = cmean, exponentiate = TRUE)))
            res_char_cmean <- paste0(paste0(res_char_vec_cmean, collapse = "\n"), "\n")
            res_char <- paste0(res_char,
                               "\n### Re-evaluation at c_cond = cmean\n",
                               res_char_cmean)
        }
        ## The working directory here is the enclosing folder ./tests/testthat/.
        ## print(getwd())
        relpath <- paste0("../reference_results/", filename_r_res)
        cat(res_char, file = relpath)
    }))


###
### Extract R results for testing against SAS macro results
################################################################################


macro_args_sas_r <- macro_args_sas_r_prelim %>%
    ## Extract useful elements when available
    mutate(
        ## Evaluated at the specified c_cond values
        ## Check against SAS macro's "conditional" results
        coef = map(res, function(res) {
            if(is.error(res)) {
                return(NULL)
            } else {
                return(coef(res))
            }
        }),
        se = map(res, function(res) {
            if(is.error(res)) {
                return(NULL)
            } else {
                return(sqrt(diag(vcov(res))))
            }
        }),
        p = map(res, function(res) {
            if(is.error(res)) {
                return(NULL)
            } else {
                return(coef(summary(res))[,"p"])
            }
        }),
        ## 2 * (1 - pnorm(1.96)) to get confint corresponding to 1.96 * se
        lower = map(res, function(res) {
            if(is.error(res)) {
                return(NULL)
            } else {
                ## confidence level 0.9500042 corresponding to +/- 1.96 * se
                return(confint(res, level = 1 - (2 * (1 - pnorm(1.96))))[,"lower"])
            }
        }),
        upper = map(res, function(res) {
            if(is.error(res)) {
                return(NULL)
            } else {
                ## confidence level 0.9500042 corresponding to +/- 1.96 * se
                return(confint(res, level = 1 - (2 * (1 - pnorm(1.96))))[,"upper"])
            }
        }),
        ##
        ## Re-evaluated at mean cvar vector values.
        ## Check against SAS macro's "marginal" results
        coef_cmean = pmap(list(res, cmean), function(res, cmean) {
            if(is.error(res) | is.null(cmean)) {
                return(NULL)
            } else {
                return(coef(res, c_cond = cmean))
            }
        }),
        se_cmean = pmap(list(res, cmean), function(res, cmean) {
            if(is.error(res)  | is.null(cmean)) {
                return(NULL)
            } else {
                return(sqrt(diag(vcov(res, c_cond = cmean))))
            }
        }),
        p_cmean = pmap(list(res, cmean), function(res, cmean) {
            if(is.error(res)  | is.null(cmean)) {
                return(NULL)
            } else {
                return(coef(summary(res, c_cond = cmean))[,"p"])
            }
        }),
        ## 2 * (1 - pnorm(1.96)) to get confint corresponding to 1.96 * se
        lower_cmean = pmap(list(res, cmean), function(res, cmean) {
            if(is.error(res)  | is.null(cmean)) {
                return(NULL)
            } else {
                ## confidence level 0.9500042 corresponding to +/- 1.96 * se
                return(confint(res, level = 1 - (2 * (1 - pnorm(1.96))), c_cond = cmean)[,"lower"])
            }
        }),
        upper_cmean = pmap(list(res, cmean), function(res, cmean) {
            if(is.error(res)  | is.null(cmean)) {
                return(NULL)
            } else {
                ## confidence level 0.9500042 corresponding to +/- 1.96 * se
                return(confint(res, level = 1 - (2 * (1 - pnorm(1.96))), c_cond = cmean)[,"upper"])
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
    ##
    ## FIXME: yreg loglinear must be implemented or dropped.
    filter(yreg != "loglinear") %>%
    ##
    mutate(
        junk = pmap(
            list(filename, sas, res, coef, se, p, lower, upper,
                 coef_cmean, se_cmean, p_cmean, lower_cmean, upper_cmean),
            function(filename, sas, res, coef, se, p, lower, upper,
                     coef_cmean, se_cmean, p_cmean, lower_cmean, upper_cmean) {

                ## First rule out error.
                if (is.error(res)) {
                    stop(paste0("R fit for ", filename, " gave an try-error object!"))
                } else {

                    ## Create named vectors
                    sas_estimate <- sas$estimate %>% `names<-`(sas$effect)
                    sas_p <- sas$p %>% `names<-`(sas$effect)
                    sas_lower <- sas$lower %>% `names<-`(sas$effect)
                    sas_upper <- sas$upper %>% `names<-`(sas$effect)

                    ## Use [["cde"]] etc to unname it.
                    if (nrow(sas) == 12) {

                        ## Full result case with interaction
                        if (ncol(sas) == 6) {
                            ## yreg linear with se column
                            sas_se <- sas$se %>% `names<-`(sas$effect)

                            test_that(paste0("coef match with ", filename), {
                                expect_equal(coef[["cde"]], sas_estimate[["conditional_cde"]])
                                expect_equal(coef[["pnde"]], sas_estimate[["conditional_pnde"]])
                                expect_equal(coef[["tnie"]], sas_estimate[["conditional_tnie"]])
                                expect_equal(coef[["tnde"]], sas_estimate[["conditional_tnde"]])
                                expect_equal(coef[["pnie"]], sas_estimate[["conditional_pnie"]])
                                expect_equal(coef[["te"]], sas_estimate[["conditional_te"]])
                            })

                            test_that(paste0("se match with ", filename), {
                                expect_equal(se[["cde"]], sas_se[["conditional_cde"]])
                                expect_equal(se[["pnde"]], sas_se[["conditional_pnde"]])
                                expect_equal(se[["tnie"]], sas_se[["conditional_tnie"]])
                                expect_equal(se[["tnde"]], sas_se[["conditional_tnde"]])
                                expect_equal(se[["pnie"]], sas_se[["conditional_pnie"]])
                                expect_equal(se[["te"]], sas_se[["conditional_te"]])
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

                            test_that(paste0("coef_cmean match with ", filename), {
                                expect_equal(coef_cmean[["cde"]], sas_estimate[["marginal_cde"]])
                                expect_equal(coef_cmean[["pnde"]], sas_estimate[["marginal_pnde"]])
                                expect_equal(coef_cmean[["tnie"]], sas_estimate[["marginal_tnie"]])
                                expect_equal(coef_cmean[["tnde"]], sas_estimate[["marginal_tnde"]])
                                expect_equal(coef_cmean[["pnie"]], sas_estimate[["marginal_pnie"]])
                                expect_equal(coef_cmean[["te"]], sas_estimate[["marginal_te"]])
                            })

                            test_that(paste0("se_cmean match with ", filename), {
                                expect_equal(se_cmean[["cde"]], sas_se[["marginal_cde"]])
                                expect_equal(se_cmean[["pnde"]], sas_se[["marginal_pnde"]])
                                expect_equal(se_cmean[["tnie"]], sas_se[["marginal_tnie"]])
                                expect_equal(se_cmean[["tnde"]], sas_se[["marginal_tnde"]])
                                expect_equal(se_cmean[["pnie"]], sas_se[["marginal_pnie"]])
                                expect_equal(se_cmean[["te"]], sas_se[["marginal_te"]])
                            })

                            test_that(paste0("p_cmean match with ", filename), {
                                expect_equal(p_cmean[["cde"]], sas_p[["marginal_cde"]])
                                expect_equal(p_cmean[["pnde"]], sas_p[["marginal_pnde"]])
                                expect_equal(p_cmean[["tnie"]], sas_p[["marginal_tnie"]])
                                expect_equal(p_cmean[["tnde"]], sas_p[["marginal_tnde"]])
                                expect_equal(p_cmean[["pnie"]], sas_p[["marginal_pnie"]])
                                expect_equal(p_cmean[["te"]], sas_p[["marginal_te"]])
                            })

                            test_that(paste0("lower_cmean confint match with ", filename), {
                                expect_equal(lower_cmean[["cde"]], sas_lower[["marginal_cde"]])
                                expect_equal(lower_cmean[["pnde"]], sas_lower[["marginal_pnde"]])
                                expect_equal(lower_cmean[["tnie"]], sas_lower[["marginal_tnie"]])
                                expect_equal(lower_cmean[["tnde"]], sas_lower[["marginal_tnde"]])
                                expect_equal(lower_cmean[["pnie"]], sas_lower[["marginal_pnie"]])
                                expect_equal(lower_cmean[["te"]], sas_lower[["marginal_te"]])
                            })

                            test_that(paste0("upper_cmean confint match with ", filename), {
                                expect_equal(upper_cmean[["cde"]], sas_upper[["marginal_cde"]])
                                expect_equal(upper_cmean[["pnde"]], sas_upper[["marginal_pnde"]])
                                expect_equal(upper_cmean[["tnie"]], sas_upper[["marginal_tnie"]])
                                expect_equal(upper_cmean[["tnde"]], sas_upper[["marginal_tnde"]])
                                expect_equal(upper_cmean[["pnie"]], sas_upper[["marginal_pnie"]])
                                expect_equal(upper_cmean[["te"]], sas_upper[["marginal_te"]])
                            })

                        } else {
                            ## yreg non-linear with exp(coef), p, exp(lower), exp(upper) in SAS results.
                            test_that(paste0("coef match with ", filename), {
                                expect_equal(coef[["cde"]], log(sas_estimate[["conditional_cde"]]))
                                expect_equal(coef[["pnde"]], log(sas_estimate[["conditional_pnde"]]))
                                expect_equal(coef[["tnie"]], log(sas_estimate[["conditional_tnie"]]))
                                expect_equal(coef[["tnde"]], log(sas_estimate[["conditional_tnde"]]))
                                expect_equal(coef[["pnie"]], log(sas_estimate[["conditional_pnie"]]))
                                expect_equal(coef[["te"]], log(sas_estimate[["conditional_te"]]))
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
                                expect_equal(lower[["cde"]], log(sas_lower[["conditional_cde"]]))
                                expect_equal(lower[["pnde"]], log(sas_lower[["conditional_pnde"]]))
                                expect_equal(lower[["tnie"]], log(sas_lower[["conditional_tnie"]]))
                                expect_equal(lower[["tnde"]], log(sas_lower[["conditional_tnde"]]))
                                expect_equal(lower[["pnie"]], log(sas_lower[["conditional_pnie"]]))
                                expect_equal(lower[["te"]], log(sas_lower[["conditional_te"]]))
                            })

                            test_that(paste0("upper confint match with ", filename), {
                                expect_equal(upper[["cde"]], log(sas_upper[["conditional_cde"]]))
                                expect_equal(upper[["pnde"]], log(sas_upper[["conditional_pnde"]]))
                                expect_equal(upper[["tnie"]], log(sas_upper[["conditional_tnie"]]))
                                expect_equal(upper[["tnde"]], log(sas_upper[["conditional_tnde"]]))
                                expect_equal(upper[["pnie"]], log(sas_upper[["conditional_pnie"]]))
                                expect_equal(upper[["te"]], log(sas_upper[["conditional_te"]]))
                            })

                            test_that(paste0("coef_cmean match with ", filename), {
                                expect_equal(coef_cmean[["cde"]], log(sas_estimate[["marginal_cde"]]))
                                expect_equal(coef_cmean[["pnde"]], log(sas_estimate[["marginal_pnde"]]))
                                expect_equal(coef_cmean[["tnie"]], log(sas_estimate[["marginal_tnie"]]))
                                expect_equal(coef_cmean[["tnde"]], log(sas_estimate[["marginal_tnde"]]))
                                expect_equal(coef_cmean[["pnie"]], log(sas_estimate[["marginal_pnie"]]))
                                expect_equal(coef_cmean[["te"]], log(sas_estimate[["marginal_te"]]))
                            })

                            test_that(paste0("p_cmean match with ", filename), {
                                expect_equal(p_cmean[["cde"]], sas_p[["marginal_cde"]])
                                expect_equal(p_cmean[["pnde"]], sas_p[["marginal_pnde"]])
                                expect_equal(p_cmean[["tnie"]], sas_p[["marginal_tnie"]])
                                expect_equal(p_cmean[["tnde"]], sas_p[["marginal_tnde"]])
                                expect_equal(p_cmean[["pnie"]], sas_p[["marginal_pnie"]])
                                expect_equal(p_cmean[["te"]], sas_p[["marginal_te"]])
                            })

                            test_that(paste0("lower_cmean confint match with ", filename), {
                                expect_equal(lower_cmean[["cde"]], log(sas_lower[["marginal_cde"]]))
                                expect_equal(lower_cmean[["pnde"]], log(sas_lower[["marginal_pnde"]]))
                                expect_equal(lower_cmean[["tnie"]], log(sas_lower[["marginal_tnie"]]))
                                expect_equal(lower_cmean[["tnde"]], log(sas_lower[["marginal_tnde"]]))
                                expect_equal(lower_cmean[["pnie"]], log(sas_lower[["marginal_pnie"]]))
                                expect_equal(lower_cmean[["te"]], log(sas_lower[["marginal_te"]]))
                            })

                            test_that(paste0("upper_cmean confint match with ", filename), {
                                expect_equal(upper_cmean[["cde"]], log(sas_upper[["marginal_cde"]]))
                                expect_equal(upper_cmean[["pnde"]], log(sas_upper[["marginal_pnde"]]))
                                expect_equal(upper_cmean[["tnie"]], log(sas_upper[["marginal_tnie"]]))
                                expect_equal(upper_cmean[["tnde"]], log(sas_upper[["marginal_tnde"]]))
                                expect_equal(upper_cmean[["pnie"]], log(sas_upper[["marginal_pnie"]]))
                                expect_equal(upper_cmean[["te"]], log(sas_upper[["marginal_te"]]))
                            })
                        }

                    } else if (nrow(sas) == 6) {
                        ## output = full gives 6 rows if cvar is empty.

                        ## Full result case with interaction
                        if (ncol(sas) == 6) {
                            ## yreg linear with se column
                            sas_se <- sas$se %>% `names<-`(sas$effect)

                            test_that(paste0("coef match with ", filename), {
                                expect_equal(coef[["cde"]], sas_estimate[["cde"]])
                                expect_equal(coef[["pnde"]], sas_estimate[["pnde"]])
                                expect_equal(coef[["tnie"]], sas_estimate[["tnie"]])
                                expect_equal(coef[["tnde"]], sas_estimate[["tnde"]])
                                expect_equal(coef[["pnie"]], sas_estimate[["pnie"]])
                                expect_equal(coef[["te"]], sas_estimate[["te"]])
                            })

                            test_that(paste0("se match with ", filename), {
                                expect_equal(se[["cde"]], sas_se[["cde"]])
                                expect_equal(se[["pnde"]], sas_se[["pnde"]])
                                expect_equal(se[["tnie"]], sas_se[["tnie"]])
                                expect_equal(se[["tnde"]], sas_se[["tnde"]])
                                expect_equal(se[["pnie"]], sas_se[["pnie"]])
                                expect_equal(se[["te"]], sas_se[["te"]])
                            })

                            test_that(paste0("p match with ", filename), {
                                expect_equal(p[["cde"]], sas_p[["cde"]])
                                expect_equal(p[["pnde"]], sas_p[["pnde"]])
                                expect_equal(p[["tnie"]], sas_p[["tnie"]])
                                expect_equal(p[["tnde"]], sas_p[["tnde"]])
                                expect_equal(p[["pnie"]], sas_p[["pnie"]])
                                expect_equal(p[["te"]], sas_p[["te"]])
                            })

                            test_that(paste0("lower confint match with ", filename), {
                                expect_equal(lower[["cde"]], sas_lower[["cde"]])
                                expect_equal(lower[["pnde"]], sas_lower[["pnde"]])
                                expect_equal(lower[["tnie"]], sas_lower[["tnie"]])
                                expect_equal(lower[["tnde"]], sas_lower[["tnde"]])
                                expect_equal(lower[["pnie"]], sas_lower[["pnie"]])
                                expect_equal(lower[["te"]], sas_lower[["te"]])
                            })

                            test_that(paste0("upper confint match with ", filename), {
                                expect_equal(upper[["cde"]], sas_upper[["cde"]])
                                expect_equal(upper[["pnde"]], sas_upper[["pnde"]])
                                expect_equal(upper[["tnie"]], sas_upper[["tnie"]])
                                expect_equal(upper[["tnde"]], sas_upper[["tnde"]])
                                expect_equal(upper[["pnie"]], sas_upper[["pnie"]])
                                expect_equal(upper[["te"]], sas_upper[["te"]])
                            })
                        } else {
                            ## yreg non-linear with exp(coef), p, exp(lower), exp(upper) in SAS results.
                            test_that(paste0("coef match with ", filename), {
                                expect_equal(coef[["cde"]], log(sas_estimate[["cde"]]))
                                expect_equal(coef[["pnde"]], log(sas_estimate[["pnde"]]))
                                expect_equal(coef[["tnie"]], log(sas_estimate[["tnie"]]))
                                expect_equal(coef[["tnde"]], log(sas_estimate[["tnde"]]))
                                expect_equal(coef[["pnie"]], log(sas_estimate[["pnie"]]))
                                expect_equal(coef[["te"]], log(sas_estimate[["te"]]))
                            })

                            test_that(paste0("p match with ", filename), {
                                expect_equal(p[["cde"]], sas_p[["cde"]])
                                expect_equal(p[["pnde"]], sas_p[["pnde"]])
                                expect_equal(p[["tnie"]], sas_p[["tnie"]])
                                expect_equal(p[["tnde"]], sas_p[["tnde"]])
                                expect_equal(p[["pnie"]], sas_p[["pnie"]])
                                expect_equal(p[["te"]], sas_p[["te"]])
                            })

                            test_that(paste0("lower confint match with ", filename), {
                                expect_equal(lower[["cde"]], log(sas_lower[["cde"]]))
                                expect_equal(lower[["pnde"]], log(sas_lower[["pnde"]]))
                                expect_equal(lower[["tnie"]], log(sas_lower[["tnie"]]))
                                expect_equal(lower[["tnde"]], log(sas_lower[["tnde"]]))
                                expect_equal(lower[["pnie"]], log(sas_lower[["pnie"]]))
                                expect_equal(lower[["te"]], log(sas_lower[["te"]]))
                            })

                            test_that(paste0("upper confint match with ", filename), {
                                expect_equal(upper[["cde"]], log(sas_upper[["cde"]]))
                                expect_equal(upper[["pnde"]], log(sas_upper[["pnde"]]))
                                expect_equal(upper[["tnie"]], log(sas_upper[["tnie"]]))
                                expect_equal(upper[["tnde"]], log(sas_upper[["tnde"]]))
                                expect_equal(upper[["pnie"]], log(sas_upper[["pnie"]]))
                                expect_equal(upper[["te"]], log(sas_upper[["te"]]))
                            })
                        }

                    } else {

                        ## Reduced result case with no interaction
                        if (ncol(sas) == 6) {
                            ## yreg linear with se column
                            sas_se <- sas$se %>% `names<-`(sas$effect)

                            test_that(paste0("coef match with ", filename), {
                                expect_equal(coef[["cde"]], sas_estimate[["cde=nde"]])
                                expect_equal(coef[["pnde"]], sas_estimate[["cde=nde"]])
                                expect_equal(coef[["tnie"]], sas_estimate[["nie"]])
                                expect_equal(coef[["tnde"]], sas_estimate[["cde=nde"]])
                                expect_equal(coef[["pnie"]], sas_estimate[["nie"]])
                                expect_equal(coef[["te"]], sas_estimate[["te"]])
                            })

                            test_that(paste0("se match with ", filename), {
                                expect_equal(se[["cde"]], sas_se[["cde=nde"]])
                                expect_equal(se[["pnde"]], sas_se[["cde=nde"]])
                                expect_equal(se[["tnie"]], sas_se[["nie"]])
                                expect_equal(se[["tnde"]], sas_se[["cde=nde"]])
                                expect_equal(se[["pnie"]], sas_se[["nie"]])
                                expect_equal(se[["te"]], sas_se[["te"]])
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
                            ## yreg non-linear with exp(coef), p, exp(lower), exp(upper) in SAS results.
                            test_that(paste0("coef match with ", filename), {
                                expect_equal(coef[["cde"]], log(sas_estimate[["cde=nde"]]))
                                expect_equal(coef[["pnde"]], log(sas_estimate[["cde=nde"]]))
                                expect_equal(coef[["tnie"]], log(sas_estimate[["nie"]]))
                                expect_equal(coef[["tnde"]], log(sas_estimate[["cde=nde"]]))
                                expect_equal(coef[["pnie"]], log(sas_estimate[["nie"]]))
                                expect_equal(coef[["te"]], log(sas_estimate[["te"]]))
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
                                expect_equal(lower[["cde"]], log(sas_lower[["cde=nde"]]))
                                expect_equal(lower[["pnde"]], log(sas_lower[["cde=nde"]]))
                                expect_equal(lower[["tnie"]], log(sas_lower[["nie"]]))
                                expect_equal(lower[["tnde"]], log(sas_lower[["cde=nde"]]))
                                expect_equal(lower[["pnie"]], log(sas_lower[["nie"]]))
                                expect_equal(lower[["te"]], log(sas_lower[["te"]]))
                            })

                            test_that(paste0("upper confint match with ", filename), {
                                expect_equal(upper[["cde"]], log(sas_upper[["cde=nde"]]))
                                expect_equal(upper[["pnde"]], log(sas_upper[["cde=nde"]]))
                                expect_equal(upper[["tnie"]], log(sas_upper[["nie"]]))
                                expect_equal(upper[["tnde"]], log(sas_upper[["cde=nde"]]))
                                expect_equal(upper[["pnie"]], log(sas_upper[["nie"]]))
                                expect_equal(upper[["te"]], log(sas_upper[["te"]]))
                            })
                        }

                    }
                }
            }))


###
### Remove the modified version of expect_equal.
################################################################################

## Need to remove from global_env() to make sure we are not overwriting.
if (rlang::env_has(rlang::global_env(), "expect_equal")) {
    rlang::env_unbind(rlang::global_env(), "expect_equal")
}
