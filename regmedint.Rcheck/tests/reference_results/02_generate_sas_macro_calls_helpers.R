library(tidyverse)


## Defined here as it will be called by
## tests/reference_results/02_generate_sas_macro_calls_helpers.R
## and
## tests/testthat/test-09_cross_check_with_sas_macro.R.
generate_sas_macro_args <- function() {

    sas_common_string <- "
/** Set libname */
libname w './';

/* Load SAS macro */
%include './mediation.sas';

/* Load data */
proc import datafile = './data-pbc_cc.csv'
    out = pbc_cc
    dbms = csv
    replace;
run;
"

    sas_macro_call_string_fmt <- "
%%mediation(
    data = pbc_cc,
    yvar = %s,
    avar = %s,
    mvar = %s,
    cvar = %s,
    a0 = %s,
    a1 = %s,
    m = %s,
    yreg = %s,
    mreg = %s,
    interaction = %s,
    casecontrol = %s,
    output = %s,
    c = %s,
    boot = ,
    cens = %s);
run;"

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
             ncvar = c(0,3)) %>%
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
                                yreg == "poisson" ~ "platelet",
                                yreg == "negbin" ~ "platelet",
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
               output = case_when(ncvar == 0 ~ "full",
                                  ncvar == 3 ~ "full"),
               ##
               c_cond = case_when(ncvar == 0 ~ "",
                                  ncvar == 3 ~ "50 1 2"),
               ## Stress test with nonsensical tricky values
               a0 = "1.1",
               a1 = "2.3",
               m_cde = "1.4") %>%
        ## Contruct macro text
        mutate(macro_call = sprintf(sas_macro_call_string_fmt,
                                    yvar,
                                    avar,
                                    mvar,
                                    cvar,
                                    a0,
                                    a1,
                                    m_cde,
                                    yreg,
                                    mreg,
                                    interaction,
                                    casecontrol,
                                    output,
                                    c_cond,
                                    cens)) %>%
        ## Construct the entire sas script
        mutate(sas_script = paste0(sas_common_string,
                                   macro_call)) %>%
        ## File name
        mutate(filename = sprintf("sas-mreg_%s_yreg_%s_int_%s_caco_%s_ncvar%s.sas",
                                  ## All variables in crossing() should appear hear.
                                  mreg,
                                  yreg,
                                  str_sub(interaction, start = 1, end = 1),
                                  str_sub(casecontrol, start = 1, end = 1),
                                  ncvar))

}
