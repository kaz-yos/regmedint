################################################################################
### Generate all relevant SAS macro calls
##
## Created on: 2020-03-24
## Author: Kazuki Yoshida
################################################################################

## When running non-interactively
.script_name. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(.script_name.) == 1) {
    cat("### Running:", paste(commandArgs()), "\n")
    options(width = 100)
    options(echo = TRUE)
}

cat("
###
### Prepare environment
################################################################################\n")

## Record start time
start_time <- Sys.time()
cat("### Started ", as.character(start_time), "\n")


## Load packages
library(tidyverse)


cat("
###
### Construct macro call templates
################################################################################\n")

common_string <- "
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

macro_call_string_fmt <- "
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
    output = full,
    c = %s,
    boot = ,
    cens = %s);
run;"


cat("
###
### Generate all relevant patterns of %mediation() calls
################################################################################\n")

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
    ## Contruct macro text
    mutate(macro_call = sprintf(macro_call_string_fmt,
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
                                c_cond,
                                cens)) %>%
    ## Construct the entire sas script
    mutate(sas_script = paste0(common_string,
                               macro_call)) %>%
    ## File name
    mutate(filename = sprintf("sas-mreg_%s_yreg_%s_int_%s_caco_%s_ncvar%s.sas",
                              ## All variables in crossing() should appear hear.
                              mreg,
                              yreg,
                              str_sub(interaction, start = 1, end = 1),
                              str_sub(casecontrol, start = 1, end = 1),
                              ncvar))
macro_args %>%
    print(n = Inf)

macro_args %>%
    select(filename) %>%
    print(n = Inf)


cat("
###
### Write to sas scripts
################################################################################\n")

junk <- macro_args %>%
    mutate(junk = map2(sas_script, filename,
                       function(sas_script, filename) {
                           write_lines(x = sas_script, path = filename)
                       }))


################################################################################
cat("
###
### Record package versions etc
################################################################################\n")
print(sessionInfo())
## Record execution time and multicore use
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units = "auto")
