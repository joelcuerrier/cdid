library(devtools)
devtools::build() 

load_all()






devtools::check()

devtools::document()

rm(list = ls())  # Clear the environment




# remotes::install_github("joelcuerrier/cdid", ref = "august26", build_vignettes = FALSE, force = TRUE)


devtools::build_vignettes()
knitr::knit("/Users/joelcuerrier/Desktop/cdid-august26/vignettes/vignettes.Rmd")
getwd()