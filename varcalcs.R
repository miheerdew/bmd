varcalc1 <- function (Yvec, Xmat) {
  
  Xmat <- scale(Xmat)
  Yvec <- as.vector(scale(Yvec))
  
  n <- nrow(Xmat)
  
  if (n != length(Yvec))
    stop('X and Y variables must be of same dimension\n')
  
  # General calcs
  xyCors <- as.vector(cor(Yvec, Xmat))
  y4 <- sum(Yvec^4)
  xRowSum <- rowSums(Xmat)
  xRowSum2 <- tcrossprod(xyCors, Xmat^2)
  
  # Calc for star 1
  star1 <- crossprod(Yvec^2, xRowSum^2)
  
  # Calc for star 2
  star2 <- y4 * sum(xyCors)^2
  
  # Calc for star 3
  star3 <- 2 * sum(xyCors) * crossprod(Yvec^2, t(xRowSum2))
  
  # Calc for star 4
  star4 <- tcrossprod(xRowSum2)
  
  # Calc for dagger 1
  dagger1 <- sum(xyCors) * crossprod(Yvec^3, xRowSum)
  
  # Calc for dagger 2
  dagger2 <- crossprod(t(xRowSum * xRowSum2), Yvec)
  
  
  
  return((star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2) / (n - 1))
  
  
}

varcalc1_multi <- function (Ymat, Xmat) {
  
  Xmat <- scale(Xmat)
  Ymat <- scale(Ymat)
  
  n <- nrow(Xmat)
  
  if (n != nrow(Ymat))
    stop('X and Y variables must be of same dimension\n')
  
  # General calcs
  xyCors <- cor(Ymat, Xmat)
  y4 <- colSums(Ymat^4)
  xRowSum <- rowSums(Xmat)
  xRowSum2 <- tcrossprod(xyCors, Xmat^2)
  
  # Calc for star 1
  star1 <- crossprod(Ymat^2, xRowSum^2)
  
  # Calc for star 2
  star2 <- y4 * rowSums(xyCors)^2
  
  # Calc for star 3
  star3 <- 2 * rowSums(xyCors) * colSums(Ymat^2 * t(xRowSum2))
  
  # Calc for star 4
  star4 <- rowSums(xRowSum2^2)
  
  # Calc for dagger 1
  dagger1 <- rowSums(xyCors) * crossprod(Ymat^3, xRowSum)
  
  # Calc for dagger 2
  dagger2 <- colSums(xRowSum * t(xRowSum2) * Ymat)
  
  return((star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2) / (n - 1))
  
  
}