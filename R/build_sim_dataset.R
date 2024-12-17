#' Load the DGP1 Dataset
#'
#' This function loads the DGP1 dataset from the package's inst/extdata directory.
#'
#' @return A data frame containing the DGP1 dataset.
#' @export


build_sim_dataset <- function(dgp=1) {
  # Define the path to the file
  file_name <- sprintf("DGP%d.xlsx", dgp)  # or use paste0("DGP", dgp, ".xlsx")
  file_path <- system.file("extdata", file_name, package = "cdid")
  
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the Excel file
    library(readxl)
    data <- read_xlsx(file_path)
    return(data)
  } else {
    stop("DGP1 dataset not found.")
  }
}
