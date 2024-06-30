#' @title Process Results from [compute.att_gt()]
#'
#' @param attgt.list list of results from [compute.att_gt()]
#'
#' @return list with elements:
#' \item{group}{which group a set of results belongs to}
#' \item{tt}{which time period a set of results belongs to}
#' \item{att}{the group time average treatment effect}
#'
#' @export
process_attgt_gmm <- function(attgt.list) {
  # nG <- length(unique(unlist(BMisc::getListElement(attgt.list, "group"))))
  # nT <- length(unique(unlist(BMisc::getListElement(attgt.list, "year"))))
    
  library(purrr)
  ntbis=(length(unique(map_dfr(attgt.list, ~as.data.frame(.x))$tbis)))
  nG=(length(unique(map_dfr(attgt.list, ~as.data.frame(.x))$group)))
  nT=length(unique(map_dfr(attgt.list, ~as.data.frame(.x))$year))
  
  # create vectors to hold the results
  group <- c()
  att <- c()
  tt <- c()
  i <- 1

  # counter=1
  # populate result vectors and matrices
  for (f in nG) {
    for (s in nT) { #added the dimension k. 
      for (t in 1:dim(map_dfr(attgt.list, ~as.data.frame(.x)))[1]/(nG*nT)) {
      
      
      group[i] <- attgt.list[[i]]$group
      tt[i] <- attgt.list[[i]]$year
      att[i] <- attgt.list[[i]]$att
      i <- i+1
      }
    }
  }

  
  

  list(group=group, att=att, tt=tt)
}