#' @title gg
#'
#' @description Function to simplify weight computations. For more details on the methodology, see:
#' Bellego, Benatia, and Dortet-Bernadet (2024), "The Chained Difference-in-Differences",
#' Journal of Econometrics, https://doi.org/10.1016/j.jeconom.2023.11.002.
#' @param x predictors
#' @param thet parameters
#' @export

gg<-function (x, thet) {
  x <- as.matrix(x)
  thet <- as.matrix(thet)
  gval <- 1/((1 + exp(x %*% thet))^2)
  as.numeric(gval)}


