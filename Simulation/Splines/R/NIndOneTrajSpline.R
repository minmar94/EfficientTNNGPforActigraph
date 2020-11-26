# Packages ----------------------------------------------------------------
require(tidyverse)
require(magrittr)
require(RcppZiggurat)

# Source Cpp file
Rcpp::sourceCpp("Simulation/Splines/Rcpp/MVNGenSpline.cpp")

#########################################################################################################################
#
# Simulation of a temporal Gaussian Process with exponential covariance function, sampled over unidimensional spatial trajectories 
# evolving on a squared region. 
# The mean term includes independent covariates and a spline surface effect.
#
# Inputs: 
#   - J = number of individuals
#   - TJ = number of observations for each individual j = 1,..,J
#   - thetaT = time rate (we assume exponential waiting time)
#   - beta0 = individual intercepts default values
#   - smin, smax = boundaries of the space domain  
#   - kns = number of knots
#   - deg = degree of the spline basis
#   - lambda20 = shrinkage parameter
#   - sigma20 = sill 
#   - phi0 = decay parameter 
#   - tau20 = nugget
# 
# Output:
#   - dat = simulated dataset
#   - pars = parameters used for the simulation
# 
# 
#########################################################################################################################


# Generation -----------------------------------------------------
SimData <- function(J = 2, TJ = c(10, 20), thetaT = 5, 
                    beta0 = c(1.5,1.5),
                    smin = 0, smax = 10, kns = 10, deg = 2, lambda20 = 10,
                     sigma20 = 1, phi0 = 1, tau20 = 1){
  
  
  # Auxiliary
  ncov <- length(beta0)-J # number of covariates
  Dat <- data.frame()
  etas <- rnorm(((kns-2)+deg)*((kns-2)+deg))*sqrt(lambda20^(-1)) # shrinkage splines coefficients
  
  # For each individual 
  for(ind in 1:J)
  {
    # Time points
    tDiffs <- rexp(TJ[ind], thetaT)
    tPoints <- cumsum(tDiffs)
    
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
    
    # Splines
    sx <- rep(runif(1, smin, smax), (TJ[ind]))
    sy <- rep(runif(1, smin, smax), (TJ[ind]))
    for (j in 2:(TJ[ind]))
    {
      a <- sx[j-1]+rnorm(1, 0, .5*sqrt(tDiffs[j-1]))
      b <- sy[j-1]+rnorm(1, 0, .5*sqrt(tDiffs[j-1]))
      
      sx[j] <- ifelse(smin<=a & a<=smax, a, ifelse(a<smin, smin, smax))
      sy[j] <- ifelse(smin<=b & b<=smax, b, ifelse(b<smin, smin, smax))
    }

    x_spline <- splines2::bSpline(x = sx, degree = deg, Boundary.knots = c(smin, smax),
                                  knots = seq(smin, smax, length.out = kns)[-c(1, kns)])
    y_spline <- splines2::bSpline(x = sy, degree = deg, Boundary.knots = c(smin, smax), 
                                  knots = seq(smin, smax, length.out = kns)[-c(1, kns)])
    
    S <- model.matrix(~-1+x_spline:y_spline)
    colnames(S) <- paste("S", 1:(ncol(S)), sep="")
    rm(x_spline)
    rm(y_spline)
    
    # Mean
    mu <- X%*%beta0 + S%*%etas

    # Latent temporal component
    w <- GPSimul(tPoints, sigma20, phi0)
    
    # Error term
    eps <-  (sqrt(tau20)*zrnorm(TJ[[ind]], 1))[,1]
    
    # Observed process
    Z <- mu + w +  eps
    
    # Append datasets
    Dat = rbind(Dat, data.frame(IndIdx = rep(ind, times = TJ[ind]), t = tPoints, 
                     Z = Z, X=X, S=S, w=w, lon=sx, lat=sy))
  }
  return(list(dat=as_tibble(Dat), pars=list(beta=beta0, 
                                            eta=etas, lambda=lambda20, 
                                            theta=c(sigma20, phi0, tau20))))
}


