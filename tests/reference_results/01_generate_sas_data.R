################################################################################
### Generate SAS ready dataset from survival::pbc dataset
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
library(survival)
library(tableone)
library(tidyverse)


cat("
###
### Load data
################################################################################\n")

data(pbc)


cat("
###
### Manipulate
################################################################################\n")

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


cat("
###
### Show resulting data
################################################################################\n")

tab1 <- CreateTableOne(data = pbc_cc,
                       vars = c("bili","bili_bin",
                                "alk_phos","spiders","platelet","time","status","cens",
                                "age","male","stage"),
                       strata = c("trt"),
                       test = FALSE)
print(tab1, smd = TRUE)


cat("
###
### Write to a CSV file for SAS
################################################################################\n")

write_csv(pbc_cc,
          path = "./data-pbc_cc.csv")


################################################################################
cat("
###
### Record package versions etc
################################################################################\n")
print(sessionInfo())
## Record execution time and multicore use
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units = "auto")
