# Packages ----------------------------------------------------------------
require(zeallot)
require(tidyverse)
require(magrittr)
require(RcppZiggurat)
require(doParallel)
require(foreach)

Rcpp::sourceCpp("Simulation/Experiment1/Rcpp/MVNGen.cpp")

# Generation -----------------------------------------------------
SimData <- function(J = 2, TJ = c(10, 20), beta0 = c(1.5,1.5),
                    sigma20 = 1, phi0 = 1, tau20 = 1, large=F){
  # Time rate
  thetaT <- 5
  
  # Auxiliary
  ncov <- length(beta0)-J
  Dat <- data.frame()
  
  for(ind in 1:J)
  {
    # Time points
    tPoints <- cumsum(rexp(TJ[ind], thetaT))
    
    # Covariates and coefficients
    IndIntercepts <- matrix(rep(0, TJ[ind]*J), nrow=TJ[ind], ncol=J)
    IndIntercepts[, ind] <- 1
    X <- IndIntercepts
    colnames(X) <- paste("Xbeta0", 1:J, sep="")
    if(ncov >= 1){
      Covariates <- zrnorm(TJ[ind], ncov)
      X <- cbind(X, Covariates)
      colnames(X)[(J+1):(ncov+J)] <- paste("X", 1:ncov, sep = "")
    }
    
    # Mean
    mu <- X%*%beta0
    
    # Latent gaussian
    w <- GPSimul(tPoints, sigma20, phi0)
    
    # Error term
    eps <-  (sqrt(tau20)*zrnorm(TJ[[ind]], 1))[,1]
    
    # Observed process
    Z <- mu + w +  eps
    
    # Appending datasets
    Dat = rbind(Dat, data.frame(IndIdx = rep(ind, times = TJ[ind]), t = tPoints, Z = Z, X, w=w))
  }
  return(as_tibble(Dat))
}

