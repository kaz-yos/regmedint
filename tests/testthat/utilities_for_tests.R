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
                   stringsAsFactors = FALSE,
                   col.names = c("effect","estimate","p","lower","upper")) %>%
        as_tibble()
    ## Make sure the connection is closed
    close(conn)
    ## Finally return the parsed result
    return(res)
}
