makeVars <- function(Uis, M_A){
  
  # stdized to ss = 1
  
  n1 = ncol(Uis)
  
  r1s = M_A %*% t(Uis)
  k = nrow(r1s)
  
  Uis = Uis*sqrt(n1)
  M_A = M_A*sqrt(n1)
  W = t(array(colMeans(M_A), c(ncol(Uis), nrow(Uis))))
  Y = (t(r1s) %*% M_A^2)/k
  rA = colMeans(r1s)
  
  mat = 1/4*rA^2*rowMeans(Uis^4) + rA*rowMeans(Y/2*Uis^2) + rowMeans(W^2*Uis^2 + Y^2/4 - W*Y*Uis) - rA*rowMeans(W*Uis^3)
  
  return(mat/n1)
  
}

makeVar <- function(Ui, M_A){
  
  # stdized to ss = 1
  
  n1 = length(Ui)
  k = dim(M_A)[1]
  
  r1s = M_A %*% Ui
  
  Ui = Ui*sqrt(n1)
  M_A = M_A*sqrt(n1)
  W = colMeans(M_A)
  Y = (t(r1s) %*% M_A^2)/k
  rA = mean(r1s)
  
  mat = 1/4*rA^2*Ui^4 + rA*Y/2*Ui^2 + W^2*Ui^2 + Y^2/4 - W*Y*Ui - rA*W*Ui^3
  
  
  return(mean(mat)/n1)
  
}