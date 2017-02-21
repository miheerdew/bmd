varcalc1 <- function (Xmat, Yvec) {
  
  Xmat <- scale(Xmat)
  Yvec <- as.vector(scale(Yvec))
  
  n <- nrow(Xmat)
  
  if (n != length(Yvec))
    stop('X and Yvec variables must be of same dimension\n')
  
  # General calcs
  xyCors <- as.vector(cor(Yvec, Xmat))
  xyDoubleCors <- as.vector(crossprod(Yvec^2, Xmat^2))
  xxDoubleCors <- crossprod(Xmat^2)
  y4 <- sum(Yvec^4)
  
  # Calc for star 1
  xyMat <- Xmat * Yvec
  star1 <- sum(crossprod(xyMat))
  
  # Calc for star 2
  star2 <- y4 * sum(xyCors)^2
  
  # Calc for star 3
  star3 <- 2 * sum(xyCors) * sum(xyCors * xyDoubleCors)
  
  # Calc for star 4
  star4 <- t(xyCors) %*% xxDoubleCors %*% xyCors
  
  # Calc for dagger 1
  y3x <- as.vector(crossprod(Yvec^3, Xmat))
  dagger1 <- sum(xyCors) * sum(y3x)
  
  # Calc for dagger 2
  yxjjkmat <- crossprod(Xmat * Yvec, Xmat^2)
  dagger2 <- sum(yxjjkmat %*% xyCors)
  
  return((star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2) / (n - 1))
  
  
}