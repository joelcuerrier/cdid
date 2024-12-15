library(devtools)
devtools::build() 

load_all()


library(cdid)

document()

devtools::check()

# devtools::document()

rm(list = ls())  # Clear the environment


source("R/pre_process_cdid.R")




