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

source("./02_generate_sas_macro_calls_helpers.R")

macro_args <- generate_sas_macro_args()

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
