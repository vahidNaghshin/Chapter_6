# If the sample autocorrelation function of a time series appears to cut off
# after lag q (i.e., autocorrelations at lags higher than q are not significantly
# different from 0 and do not follow any clear patterns), then an MA(q)
# model might be suitable. An AR(p) model is indicated when the partial
# autocorrelation function cuts off after lag p. If there are no convincing
# cutoff points for either function, an ARMA model may provide the best
# fit.
rho <- function(k, beta) {
  q <- length(beta) - 1
  if (k > q) ACF <- 0 else {
    s1 <- 0; s2 <- 0
    for (i in 1:(q-k+1)) s1 <- s1 + beta[i] * beta[i+k]
    for (i in 1:(q+1)) s2 <- s2 + beta[i]^2
    ACF <- s1 / s2}
  ACF}
beta <- c(1, 0.5)
rho.k <- rep(1, 10)
for (k in 1:10) rho.k[k] <- rho(k, beta)
plot(0:10, c(1, rho.k), pch = 4, ylab = expression(rho[k]))
abline(0, 0)
beta <- c(1, 2)
rho.k <- rep(1, 10)
for (k in 1:10) rho.k[k] <- rho(k, beta)
plot(0:10, c(1, rho.k), pch = 4, ylab = expression(rho[k]))
abline(0, 0)
set.seed(1)
b <- c(1, 0.5)
x1 <- w <- rnorm(1000)
for (t in 3:1000) {
  for (j in 1:2) x1[t] <- x1[t] + b[j] * w[t - j]
}
plot(x1, type = "l")
acf(x1)
set.seed(1)
b <- c(1, 2)
x2 <- w <- rnorm(1000)
for (t in 3:1000) {
  for (j in 1:2) x2[t] <- x2[t] + b[j] * w[t - j]
}
layout(1:2)
plot(x1, type = "l")
plot(x2, type = "l")
acf(x2)
print(var(x1))
print(var(x2))



data(AirPassengers)
AP <- AirPassengers
plot(AP)
plot(log(AP))
SIN <- COS <- matrix(nr = length(AP), nc = 6)
for (i in 1:6) {
  SIN[, i] <- sin(2 * pi * i * time(AP))
  COS[, i] <- cos(2 * pi * i * time(AP))
}
TIME <- (time(AP) - mean(time(AP)))/sd(time(AP))
library(nlme)
AP.gls <- lm(log(AP) ~ TIME + I(TIME^2) + SIN[,1] + COS[,1] +
               SIN[,2] + COS[,2] + SIN[,3] + SIN[,4] + COS[,4] + SIN[,5])
coef(AP.gls)/sqrt(diag(vcov(AP.gls)))
pacf(resid(AP.gls))
# deriving the best ARMA order based on AIC
best.order <- c(0, 0, 0)
best.aic <- Inf
for (i in 0:2) for (j in 0:2) {fit.aic <- AIC(arima(resid(AP.gls), order = c(i, 0,j)))}
if (fit.aic < best.aic) {
  best.order <- c(i, 0, j)
  best.arma <- arima(resid(AP.gls), order = best.order)
  best.aic <- fit.aic
}
acf(resid(best.arma))
new.t <- time(ts(start = 1961, end = c(1961, 12), fr = 12))
TIME <- (new.t - mean(time(AP)))/sd(time(AP))
SIN <- COS <- matrix(nr = length(new.t), nc = 6)
for (i in 1:6) {
  COS[, i] <- cos(2 * pi * i * new.t)
  SIN[, i] <- sin(2 * pi * i * new.t)
}
SIN <- SIN[, -6]
new.dat <- data.frame(TIME = as.vector(TIME), SIN = SIN,
                      COS = COS)
AP.pred.ts <- exp(ts(predict(AP.gls, new.dat), st = 1961,
                     fr = 12)+ts(predict(best.arma, n.ahead = 12)$pred, st = 1961,
                                 fr = 12))
ts.plot(log(AP), log(AP.pred.ts), lty = 1:2)
ts.plot(AP, AP.pred.ts, lty = 1:2)


alpha <- 0.7
beta <- -0.5
rho <- function(k, alpha, beta)
{
  alpha^(k-1)*(alpha+beta)*(1+(alpha*beta))/(1+(alpha*beta)+beta^2)
}
rho.k <- rep(1, 20)
for (k in 1:20) rho.k[k] <- rho(k, alpha, beta)
plot(0:20, c(1, rho.k), pch = 4, ylab = expression(rho[k]))
abline(0, 0)
set.seed(1)
a <- 0.7
b <- -0.5
x <- w <- rnorm(100)
for (t in 2:100) {
  x[t] <- x[t] + a*x[t-1] + b * w[t-1]
}
plot(x, type = "l")
acf(x)
set.seed(1)
a <- 0.7
b <- -0.5
x <- w <- rnorm(1000)
for (t in 2:1000) {
  x[t] <- x[t] + a*x[t-1] + b * w[t-1]
}
plot(x, type = "l")
acf(x)



set.seed(1)
m <- rep(0, 100)
for (i in 1:100) m[i] <- mean(rnorm(20))
var(m)
set.seed(1)
m <- rep(0, 100)
for (i in 1:100)
{
  b <- 0.5
  x <- w <- rnorm(20)
  for (t in 2:20) x[t] <- x[t] + b * w[t - 1]
  m[i] <- mean(x)
}
var(m)
set.seed(1)
m <- rep(0, 100)
for (i in 1:100)
{
  b <- -0.5
  x <- w <- rnorm(20)
  for (t in 2:20) x[t] <- x[t] + b * w[t - 1]
  m[i] <- mean(x)
}
var(m)

set.seed(1)
b <- c(0.8, 0.6, 0.4)
x <- w <- rnorm(1000)
for (t in 4:1000) {
  for (j in 1:3) x[t] <- x[t] + b[j] * w[t - j]
}
plot(x, type = "l")
acf(x)
pacf(x)
x.ma <- arima(x, order = c(0, 0, 3))
x.ma$aic
x.ar <- ar(x, order.max = 3)
x.ar$aic
x.arma <- arima(x, order = c(1, 0, 1))
x.arma$aic

