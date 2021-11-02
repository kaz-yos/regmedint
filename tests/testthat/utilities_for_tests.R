################################################################################
### Utilitiy functions for testing purpose only
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################

## Required packages
library(tidyverse)
library(stringr)

read_parsed_sas_mediation_output <- function(file) {
    ## Open a connection and save as a file/connection object.
    conn <- file(file, open = "r")
    ## Parse
    res <- conn %>%
        readLines() %>%
        stringr::str_replace_all("^ *", "") %>%
        stringr::str_replace_all("^[0-9]+ *", "") %>%
        stringr::str_replace_all("total effect", "te") %>%
        stringr::str_replace_all("nal ", "nal_") %>%
        stringr::str_replace_all(" +", " ") %>%
        read.delim(text = .,
                   header = FALSE,
                   sep = " ",
                   stringsAsFactors = FALSE) %>%
        as_tibble()
    ## Make sure the connection is closed
    close(conn)
    if (ncol(res) == 6) {
        names(res) <- c("effect","estimate","se","p","lower","upper")
    } else if (ncol(res) == 5) {
        names(res) <- c("effect","estimate","p","lower","upper")
    } else {
        print(ncol(res))
        stop("The number of columns is unexpected.")
    }
    ## Finally return the parsed result
    return(res)
}

## Functions to skip describe() and it()
xdescribe <- function(description, ...) {
    cat("Skipping", description, "\n")
}
xit <- function(description, ...) {
    cat("Skipping", description, "\n")
}

## http://adv-r.had.co.nz/Exceptions-Debugging.html
is.error <- function(x) {
    inherits(x, "try-error")
}
