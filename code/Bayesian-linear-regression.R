# Linear Regression
# R coding for Gibbs sampling

# Example: QTL affecting a specific quantitative trait

# y = Xb + e

y<-matrix(c(95.9, 108.0, 96.5, 92.9, 101.0, 94.5, 93.7, 89.8,
            101.2, 103.9, 85.9, 109.4, 105.7, 98.4, 84.1, 103.1,
            117.1, 95.2, 106.4, 104.7, 92.5, 123.9, 97.8, 100.5),nrow=24)

X<-matrix(c(1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1,
            -1, -1, -1, -1, -1, -1, -1, -1,
            0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1),nrow=24)

#solution
b <- solve(t(X) %*% X) %*% t(X) %*% y

# Using package lm

qtl <- X[,2]
plot(qtl,y,main="Scatter Plot",xlab="Genotype",ylab="Phenotypic Score",col="red", pch = 16)

reg <- lm(y ~ qtl)
reg
summary(reg)

# Gibbs sampling

library(MASS)
#library(asbio)
library(geoR)

# number of iterations
numit <- 5000

# Setting up a storage matrix for MCMC output
# 4 columns: iteration, intercept, slope, residual variance
MCMC <- array(0,c(numit,4))

# Starting value for residual variance
vare <- 500

# Some calculations
XX <- solve(t(X) %*% X)
m <- XX %*% t(X) %*% y
n <- nrow(X)

for (l in 1:numit) {

  # Sampling the location parameters (intercept and slope)
  b <- mvrnorm(1,m,XX * vare)

  # Sampling vare
  e <- y - X%*%b
  s <- (t(e)%*%e)/n
  vare <- as.numeric(rinvchisq(1,  df=n, scale=s))

  # store the sampled values

  MCMC[l,1] <- l
  MCMC[l,2:3] <- t(b)
  MCMC[l,4] <- vare

}

# For visual inspection of convergence of MCMC for each parameter
plot(MCMC[,1],MCMC[,4],type='l', main="Trace Plot (Residual Variance)", xlab="Iteration", ylab="Residual Variance", col = "blue")

plot(MCMC[,1],MCMC[,3],type='l', main="Trace Plot (Slope)", xlab="Iteration", ylab="Slope", col = "blue")

# To look at the posterior density of each parameter
# Here we are discarding the first 100 values as burn-in
plot(density(MCMC[101:numit,4]),main="Marginal Posterior Distribution", xlab="Residual Variance", col = "red")

# To get the posterior means of each parameter
mean(MCMC[101:numit,4])
