# Packages ----------------------------------------------------------------

require(Rcpp)
require(RcppArmadillo)
require(RcppEigen)
require(spBayes)

# Rcpp
sourceCpp("Simulation/Experiment1/Rcpp/tNngpColl_Final_AdaRW_Preds.cpp")

# Data ---------------------------------------------------------------------

datAll <- read_csv("Simulation/Experiment1/Data/sim_1_1e+05_3_20200723_040942_sim.txt")

ncov <- 1

# Priors
alphaS <- 2
betaS <- 2 # Gamma on sigma2
alphaPhi <- 1
betaPhi <- 1 # Uniform on phi
alphaT <- 2
betaT <- 2 # Gamma on tau2
muB <- rep(0, ncov)
vB <- diag(10^6, ncov)

# Initial
logtheta0 <- rep(log(1), 3)
beta0 <- matrix(c(1), nrow = ncov, 1)

# Random walk
rwTheta <- (0.1)^2/length(logtheta0)
rwSigma <- diag(rep(rwTheta, length(logtheta0)))
rwSigmam <- rwSigma

# Chains
M <- 1
burnInP <- 0
burnIn <- M*burnInP+1
verb <- 1

# Wrapper for spLM
priors1 <- list("beta.flat", "phi.Unif"=c(0.001, 30), "sigma.sq.IG"=c(2,2), "tau.sq.IG"=c(2, 2))
starting1 <- list("phi"=0.01, "sigma.sq"=3, "tau.sq"=0.1)
tuning1 <- list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01)
Geostat.full<-function(x, y, data, n.samples)
{
  sp <- spBayes::spLM(data~1, coords=cbind(x, y), tuning=tuning1, starting=starting1, 
                      n.samples=n.samples, priors=priors1, cov.model="exponential", verbose = F)
  return(list(sp))
}


# Test ------------------------------------------------------------------

# Increasing sample size
nsimul <- c(100, 1000, 5000, 10000)
# Number of iterations to estimate the runtime
niter <- 100
# Neighbors
neigh <- 30

# Output object
timeslistAll <- list()

# Start 
tsimul <- Sys.time()
for(i in 1:length(nsimul)){
  
  n <- nsimul[i]
  dat <- datAll[1:n, ]
  
  pte <- 0
  ntr <-  round((1-pte)*n)
  trIdx <- sort(sample(1:n, ntr, replace = F))
  
  indstr <- dat$IndIdx[trIdx]
  Ztr <- dat$Z[trIdx]
  Xtr <- as.matrix(dat[trIdx, 1])
  ttr <- dat$t[trIdx]
  
  times_miter <- matrix(NA, ncol = 2, nrow = niter)
  
  # Function ----------------------------------------------------------------
  for(k in 1:niter){
    
    t8 <- Sys.time()
    outPut <- tNngpCollapsed_NAda(t=ttr, Z=Ztr, 
                                  X=as.matrix(Xtr), 
                                  indLab=as.character(indstr),
                                  M=M, burnIn=burnInP, neigh=neigh,  
                                  beta0=beta0, logtheta0=logtheta0,
                                  rwSigma=rwSigma, rwSigmam=rwSigmam, gamma=0.05, madapt=2/3*M,
                                  alphaS=alphaS, betaS=betaS, alphaPhi=alphaPhi, betaPhi=betaPhi, 
                                  alphaT=alphaT, betaT=betaT,
                                  muB=muB, vB=vB, 
                                  verbose = verb, fileName = "",
                                  n_threads = 1)
    t9 <- Sys.time() - t8
    
    tinit <- Sys.time()
    full <- Geostat.full(x=ttr, y=rep(0,ntr), data=Ztr, n.samples=M)
    tfinalfull <- (Sys.time() - tinit)
    
    times_miter[k,] <- c(t9, tfinalfull)
      
  }
  print("NSiMUL FINISHED")
  print(nsimul[i])
  timeslistAll[[i]] <- times_miter
}
tsimulfinal <- Sys.time() - tsimul
save.image(file = "Simulation/Experiment1/WS/CollapsedVSFullComparison.RData")
