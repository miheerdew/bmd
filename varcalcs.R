varcalc1 <- function (Xmat, Y) {
  
  Xmat <- scale(Xmat)
  Y <- scale(Y)
  
  m <- ncol(Xmat)
  n <- nrow(Xmat)
  
  if (n != length(Y))
    stop('X and Y variables must be of same dimension\n')
  
  # General calcs
  xyCors <- as.vector(cor(Y, X))
  xyDoubleCors <- as.vector(cor(Y^2, Xmat^2))
  xxDoubleCors <- cor(Xmat^2)
  y4 <- sum(Y^4)
  
  # Calc for star 1
  xy2mat <- Xmat * Y^2
  star1 <- sum(cov(xy2mat))
  
  # Calc for star 2
  star2 <- y4 * sum(xyCors)^2
  
  # Calc for star 3
  star3 <- 2 * sum(xyCors) * sum(xyCors + xyDoubleCors)
  
  # Calc for star 4
  star4 <- t(xyCors) %*% xxDoubleCors %*% xyCors
  
  # Calc for dagger 1
  y3x <- as.vector(cor(Y^3, X))
  dagger1 <- sum(xyCors) * sum(y3x)
  
  # Calc for dagger 2
  yxjjkmat <- cor(Xmat * Y, Xmat^2)
  dagger2 <- sum(t(xyCors) %*% yxjjkmat)
  
  
}