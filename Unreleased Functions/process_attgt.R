# Callaway, B. (2024). The did Library. Department of Economics, University of Georgia. Available at: https://github.com/bcallaway11/did

#' @title Process Results from [compute.att_gt()]
#'
#' @param attgt.list list of results from [compute.att_gt()]
#'
#' @return list with elements:
#' \item{group}{which group a set of results belongs to}
#' \item{tt}{which time period a set of results belongs to}
#' \item{att}{the group time average treatment effect}
#' @references Bellego, Benatia,Dortet-Bernardet (2022)  \"The Chained Difference-in-Differences\" , \url{https://www.davidbenatia.com/publication/chaineddid/}
#' @export
process_attgt <- function(attgt.list) {
  nG <- length(unique(unlist(BMisc::getListElement(attgt.list, "group"))))
  nT <- length(unique(unlist(BMisc::getListElement(attgt.list, "year"))))
  
  
  # create vectors to hold the results
  group <- c()
  att <- c()
  tt <- c()
  i <- 1

  # populate result vectors and matrices
  for (f in 1:nG) {
    for (s in 1:nT) {
      group[i] <- attgt.list[[i]]$group
      tt[i] <- attgt.list[[i]]$year
      att[i] <- attgt.list[[i]]$att
      i <- i+1
    }
  }

  list(group=group, att=att, tt=tt)
}