# Useful functions
d <- function (r) {1 - 3 * r^2 + 2 * r^3}
pbr <- function (r) {Beta * (1 - r + m * r)}
pxy <- function (r) {pbr(r) / sqrt(varY(r))}

# Function for single variance calculation
varY <- function (r) {m * Beta * (1 - r + r * m) + s2}
mu1111 <- function (r) {3}
mu1112 <- function (r) {3 * r}
mu1123 <- function (r) {2 * r^2 + r}
mu1122 <- function (r) {1 + 2 * r^2}
mu1234 <- function (r) {3 * r * (r - 1)^2 * (2 * r^2 + r) / d(r)}

# Functions for total variance calculation
pyyjk <- function (r) {(2 * pbr(r)^2 * (2 * r + 1) / (r + 1) + 
                          r * Beta * pbr(r) * (m + s2 / (Beta * pbr(r)) - 
                                                 2 * (1 - r + m * r) / (1 + r))) / varY(r)}
pyyyy <- function (r) {3}
pyyjj <- function (r) {(2 * pbr(r)^2 + varY(r)) / varY(r)}
pyyyj <- function (r) {3 * varY(r) * pbr(r) / varY(r)^(3 / 2)}
pyjjk <- function (r) {pbr(r) * (2 * r + 1) / sqrt(varY(r))}
pyjjj <- function (r) {pbr(r) * 3 / sqrt(varY(r))}

# Variance calculations
popvar1234 <- function (r) {mu1234(r) + r^2 * mu1122(r) - r * 2 * mu1123(r)}

popvar1213 <- function (r) {mu1123(r) + 0.25 * r^2 * 
    (mu1111(r) + 3 * mu1122(r)) - r * (mu1112(r) + mu1123(r))}

popvar <- function (r) {
  (
    pyyjj(r) + 0.25 * pxy(r)^2 * (pyyyy(r) + 2 * pyyjj(r) +
                                    mu1111(r)) - pxy(r) * (pyyyj(r) + pyjjj(r)) + 
      (
        pyyjk(r) + 0.25 * pxy(r)^2 * (pyyyy(r) + 2 * pyyjj(r) +
                                        mu1122(r)) - pxy(r) * (pyyyj(r) + pyjjk(r))
      ) * (m - 1)
    
  ) * m}