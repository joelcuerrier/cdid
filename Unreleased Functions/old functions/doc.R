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



# archives

  # # Generate the dummy variable "P" with a randomized fraction
  # data$P_Y1_longDID <- as.numeric(runif(dim(data)[1], min = 0, max = 1)>=0)
  
  # # Filter indivs never observed across all years
  # data0 <- data[data$P_Y1_longDID == 1,]
  # data0 <- data0 %>%
  #   group_by(id) %>%
  #   filter(!(all(P_Y1_longDID == 0 & annee %in% 1:8))) %>%
  #   ungroup()
  
  
  # data0 <- data0[,c(1,2,3,6,7,8)]
  # data0 <- rename(data0,Y = Y1_chaine)
  # data0 <- rename(data0,date = annee)
  # data0 <- rename(data0,date_G = annee_G)
  # data0$P_Y1_chaine <- NULL
  
  # # Sort data0 by id and date
  # data0 <- data0 %>%
  #   arrange(id, date)
  
  # ### Correct id's
  # list_id <-as.data.frame(unique(data0[,"id"]))
  # list_id$iden_num<-1:dim(list_id)[1] 
  # data0 <- merge(data0, list_id, by.x = "id", by.y = "id", all.x = TRUE)
  # data0$id <- data0$iden_num
  # data0$iden_num <- NULL
  
  # # Add a treatment var
  # data0$treated <- 0
  # data0$treated[data0$date_G!=0] <- 1
  
  # # Add a new variable S
  # data0 <- data0 %>%
  #   group_by(id) %>%                   # Group data by id
  #   mutate(S = ifelse(date == min(date), 1, 1)) %>% # S == 1 for lowest year, 0 otherwise
  #   ungroup()       







# <!-- source("R/fonction_simu_attrition.R")
# source("R/gg.R")
# source("R/pre_process_cdid.R")
# source("gmm_compute_delta_att.R")
# source("gmm_convert_delta_to_att.R")
# source("gmm_convert_result.R")
# source('process_attgt_gmm.R') -->