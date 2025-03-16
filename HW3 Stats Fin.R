# Load necessary libraries
library(tidyquant)
library(PerformanceAnalytics)
library(ggplot2)
library(quadprog)
library(dplyr)
library(lubridate)

# Read the CSV file, treating 'Date' column as character
data = read.csv("C:/Users/aersh/Downloads/HW3 (1).csv", header = T)
price = cbind(data$CAT_AC, data$IBM_AC, data$MSFT_AC)
n = dim(price)[1]
return = price[2:n,]/price[1:(n-1),] - 1
mu = colMeans(return)
sigma = cov(return)


library(quadprog) # library for solve.QP
m = 200 # no. of points to evaluate
muP = seq(.0001,.001,length=m) # target portfolio return
sdP = rep(0, length(muP)) # sd of portfolio return
weight = matrix(0,nrow=m,ncol=3) # storage for portfolio weights
for (i in 1:length(muP)) { # find the optimal portfolios
  result = solve.QP(Dmat=2*sigma,dvec=rep(0,3),
                    Amat = cbind(rep(1,3),mu),bvec=c(1,muP[i]),meq=2)
  sdP[i] = sqrt(result$value)
  weight[i,] = result$solution
} 

# Define padding factor (e.g., 10% extra space)
x_padding <- diff(range(c(sdP, sqrt(diag(sigma))))) * 0.1
y_padding <- diff(range(c(muP, mu))) * 0.1

# Adjust xlim and ylim to ensure nothing gets cut off
plot(sdP[GMP:m], muP[GMP:m], type = "l", 
     xlim = range(c(sdP, sqrt(diag(sigma)))) + c(-x_padding, x_padding), 
     ylim = range(c(muP, mu)) + c(-y_padding, y_padding), 
     lwd = 3, col = "red", 
     xlab = "SD of portfolio return", 
     ylab = "Mean of portfolio return")

points(sdP[1:(GMP-1)], muP[1:(GMP-1)], type = "l",
       lty = 2, lwd = 3, col = "red") 

points(sqrt(diag(sigma)), mu, pch = 4) 

text(sqrt(diag(sigma)), mu, c("CAT", "IBM", "MSFT"), pos = 4) 


muP_noSS = seq(min(mu),max(mu),length=m) # target portfolio return
sdP_noSS = rep(0, length(muP_noSS))
for (i in 1:length(muP_noSS)) { # find the optimal portfolios
  result = solve.QP(Dmat=2*sigma,dvec=rep(0,3),
                    Amat=cbind(rep(1,3),mu,diag(1,3)),
                    bvec=c(1,muP_noSS[i],rep(0,3)),meq=2)
  sdP_noSS[i] = sqrt(result$value)
  weight[i,] = result$solution
}

# Define padding factor (e.g., 10% extra space)
x_padding <- diff(range(c(sdP_noSS, sqrt(diag(sigma))))) * 0.1
y_padding <- diff(range(c(muP_noSS, mu))) * 0.1

# Adjust xlim and ylim to ensure full visibility
plot(sdP_noSS[GMP_noSS:m], muP_noSS[GMP_noSS:m], type = "l", 
     xlim = range(c(sdP_noSS, sqrt(diag(sigma)))) + c(-x_padding, x_padding), 
     ylim = range(c(muP_noSS, mu)) + c(-y_padding, y_padding), 
     lwd = 3, col = "red", 
     xlab = "SD of Portfolio Return", 
     ylab = "Mean of Portfolio Return", 
     main = "Efficient Frontier without Short-Selling")

print(mu)

# Add the inefficient part of the frontier
points(sdP_noSS[1:(GMP_noSS - 1)], muP_noSS[1:(GMP_noSS - 1)], type = "l", 
       lty = 2, lwd = 3, col = "red")

# Plot individual assets
points(sqrt(diag(sigma)), mu, pch = 4, col = "blue")

# Label individual assets
text(sqrt(diag(sigma)), mu, labels = c("CAT", "IBM", "MSFT"), col = "blue", pos = 4)


annual_rf_rate <- 0.05
trading_days <- 253
daily_rf_rate <- annual_rf_rate / trading_days

excess_returns <- returns[, -1] - daily_rf_rate 
mu_excess <- colMeans(excess_returns)
sigma_excess <- cov(excess_returns)

one_vec <- rep(1, length(mu_excess))
inv_sigma_excess <- solve(sigma_excess)
weights_tangency <- inv_sigma_excess %*% mu_excess
weights_tangency <- weights_tangency / sum(weights_tangency)

mu_tangency <- sum(weights_tangency * mu)
sigma_tangency <- sqrt(t(weights_tangency) %*% sigma %*% weights_tangency)

# Define padding factors (e.g., 20% extra space)
x_padding <- diff(range(c(sdP, sigma_tangency))) * 0.2
y_padding <- diff(range(c(muP, cml_return))) * 0.2

# Adjust xlim and ylim dynamically
plot(sdP, muP, type = "l", col = "blue", lwd = 2,
     xlim = range(c(sdP, sigma_tangency)) + c(-x_padding, x_padding),
     ylim = range(c(muP, cml_return, daily_rf_rate)) + c(-y_padding, y_padding),
     xlab = "Portfolio Standard Deviation (Risk)",
     ylab = "Portfolio Expected Return",
     main = "Efficient Frontier and Capital Market Line")

# Add the CML
lines(cml_sd, cml_return, col = "green", lwd = 2)

# Mark the risk-free asset
points(0, daily_rf_rate, col = "red", pch = 4, cex = 1.5)

# Mark the tangency portfolio
points(sigma_tangency, mu_tangency, col = "red", pch = 4, cex = 1.5)

# Add a legend
legend("topright", legend = c("Efficient Frontier", "Capital Market Line", "Risk-Free Asset", "Tangency Portfolio"),
       col = c("blue", "green", "red", "red"), lty = c(1, 1, NA, NA), pch = c(NA, NA, 4, 4), cex = 0.8)
