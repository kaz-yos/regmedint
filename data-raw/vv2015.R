################################################################################
###
##
## Created on: 2020-05-06
## Author: Kazuki Yoshida
################################################################################

## Load packages
library(tidyverse)
library(usethis)


###
### Load data
################################################################################

## The file is from:
##
## https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/
## 
## A tutorial on mediation with SAS, Stata, and SPSS macros
## Valeri, L. and VanderWeele, T.J. (2013).
## Mediation analysis allowing for exposure-mediator interactions
## and causal interpretation: theoretical assumptions and implementation
## with SAS and SPSS macros. Psychological Methods, 18:137-150.

vv2015 <- read_delim(file = "./vv2015.txt",
                     delim = " ") %>%
    ## Following R convention, an event indicator is used.
    mutate(event = if_else(cens == 0, 1L, 0L))


###
### Save data
################################################################################

usethis::use_data(vv2015)
