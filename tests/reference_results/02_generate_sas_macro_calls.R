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
### Generate all relevant patterns of %mediation() calls
################################################################################\n")

macro_args <- crossing(mreg = c("linear",
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
                       casecontrol = c("false","true")) %>%
    ## casecontrol is only relevant for binary yreg
    filter(casecontrol == "false" |
           (casecontrol == "true" & (yreg == "logistic" | yreg == "loglinear"))) %>%
    ## File name
    mutate(filename = sprintf("sas-mreg_%s_yreg_%s_int_%s_caco_%s.sas",
                              mreg,
                              yreg,
                              str_sub(interaction, start = 1, end = 1),
                              str_sub(casecontrol, start = 1, end = 1)))
macro_args %>%
    print(n = Inf)


cat("
###
### Construct macro calls
################################################################################\n")




################################################################################
cat("
###
### Record package versions etc
################################################################################\n")
print(sessionInfo())
## Record execution time and multicore use
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units = "auto")
