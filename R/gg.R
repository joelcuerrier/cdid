gg<-function (x, thet) {
  x <- as.matrix(x)
  thet <- as.matrix(thet)
  gval <- 1/((1 + exp(x %*% thet))^2)
  as.numeric(gval)}


