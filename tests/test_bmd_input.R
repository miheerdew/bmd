library(Rcpp)
library(testthat)
library(microbenchmark)

sourceCpp("../bmd_input.cpp")
source("../varcalcs.R")

context("BmdInput")

XY <- MASS::mvrnorm(25, rep(0,20), diag(20))
X <- scale(XY[,1:10])
Y <- scale(XY[,10:20])
n <- nrow(X)

inp <- new(BmdInput, X,Y)

test_that("BmdInput object meta data", {
  expect_equal(inp$n, n)
  expect_equal(inp$dx, ncol(X))
  expect_equal(inp$dy, ncol(Y))
  expect_equivalent(inp$X,X)
  expect_equivalent(inp$Y,Y)
 })

test_that("BmdInput computes cross_cors correctly", {
   expect_equal(inp$cross_cors(3:5, 5:10), cor(X[,3:5],Y[,5:10]))
   expect_equal(inp$cross_cors(5, 2), cor(X[,5,drop=FALSE],Y[,2,drop=FALSE]))
})

expect_variance_correctness <- function(A, test_x){
  if(test_x){
    fixdMat <- Y[,A]
    thisMat <- X
  } else {
    fixdMat <- X[,A]
    thisMat <- Y
  }
  v1 <- inp$vars(A,test_x)
  v2 <- t(varcalc1_multi(thisMat, fixdMat))
  expect_equal(v1,v2)
}

test_that("BmdInput computes variance correctly", {
  expect_variance_correctness(1:3, FALSE)
  expect_variance_correctness(1:3, TRUE)
  expect_variance_correctness(1, FALSE)
  expect_variance_correctness(4, TRUE)
})

microbenchmark()