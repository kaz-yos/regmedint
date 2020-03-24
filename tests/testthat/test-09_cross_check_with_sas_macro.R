################################################################################
### Test against SAS macro results
##
## Created on: 2020-03-24
## Author: Kazuki Yoshida
################################################################################

library(testthat)
library(tidyverse)


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
