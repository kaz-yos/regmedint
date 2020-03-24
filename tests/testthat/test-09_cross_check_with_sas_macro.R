################################################################################
### Test against SAS macro results
##
## Created on: 2020-03-24
## Author: Kazuki Yoshida
################################################################################

library(testthat)
library(tidyverse)


test_that("placeholder fails", {
    expect_equal(TRUE,
                 FALSE)
})


###
### Read data from csv
################################################################################

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
### Conduct R analysis
################################################################################

macro_args_r_res <- macro_args %>%
    mutate(res = pmap(list(yvar,
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
                                  regmedint(data = pbc_cc,
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
                                            eventvar = NULL)
                              } else {
                                  regmedint(data = pbc_cc,
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
                                            eventvar = eventvar)
                              }
                          } else {
                              if (eventvar == "") {
                                  regmedint(data = pbc_cc,
                                            yvar = yvar,
                                            avar = avar,
                                            mvar = mvar,
                                            cvar = cvar,
                                            a0 = as.numeric(a0),
                                            a1 = as.numeric(a1),
                                            m_cde = as.numeric(m_cde),
                                            c_cond = c_cond,
                                            mreg = mreg,
                                            yreg = yreg,
                                            interaction = if_else(interaction == "true", TRUE, FALSE),
                                            casecontrol = if_else(casecontrol == "true", TRUE, FALSE),
                                            eventvar = NULL)
                              } else {
                                  regmedint(data = pbc_cc,
                                            yvar = yvar,
                                            avar = avar,
                                            mvar = mvar,
                                            cvar = cvar,
                                            a0 = as.numeric(a0),
                                            a1 = as.numeric(a1),
                                            m_cde = as.numeric(m_cde),
                                            c_cond = c_cond,
                                            mreg = mreg,
                                            yreg = yreg,
                                            interaction = if_else(interaction == "true", TRUE, FALSE),
                                            casecontrol = if_else(casecontrol == "true", TRUE, FALSE),
                                            eventvar = eventvar)
                              }
                          }
                          ##
                      }))
