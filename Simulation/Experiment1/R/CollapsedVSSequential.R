# Packages ----------------------------------------------------------------
require(spNNGP)
require(Rcpp)
require(RcppArmadillo)
require(RcppEigen)

# Rcpp
sourceCpp("Simulation/Experiment1/Rcpp/tNngpColl_Final_AdaRW_Preds.cpp")

# Data ---------------------------------------------------------------------

datAll <- read.csv("Simulation/Experiment1/Data/sim_1_1e+05_3_20200723_040942_sim.txt")

ncov <- 1

# Neighbors
neigh <- 10

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

# Prior setting
priors1 <- list("beta.flat",
                "phi.Unif"=c(0.005, 30), "sigma.sq.IG"=c(alphaS,betaS),
                "tau.sq.IG"=c(alphaT, betaT))
starting1 <- list("phi"=1, "sigma.sq"=0.3, "tau.sq"=0.1)
tuning1 <- list("phi"=0.05, "sigma.sq"=0.05, "tau.sq"=0.05)

# Wrapper for the sequential NNGP
Geostat.NNSeq <- function(x, y, data, n.samples, nn, method)
{
  sp <- spNNGP::spNNGP(data~1, coords=cbind(x, y), method=method, 
                       family="gaussian", n.neighbors=nn, 
                       starting=starting1, n.omp.threads = 1,
                       tuning=tuning1, 
                       priors=priors1, n.samples=n.samples, 
                       cov.model="exponential", verbose = F)
  
  return(sp)
}

# Test ------------------------------------------------------------------

# Increasing sample size
ns <- c(2^(0:6)*1000, 100000)
# Iterations to estimate runtime
niter <- 100

# Output
timeslistAll <- list()

tsimul <- Sys.time()
for(i in 1:length(ns)){
  t6 <- matrix(NA, nrow = niter, ncol = 2)
  
  n <- ns[i]
  
  dat <- datAll[1:n, ]
  
  pte <- 0
  ntr <-  round((1-pte)*n)
  trIdx <- sort(sample(1:n, ntr, replace = F))
  
  indstr <- dat$IndIdx[trIdx]
  Ztr <- dat$Z[trIdx]
  Xtr <- as.matrix(as.numeric(dat[trIdx, 1]))
  ttr <- dat$t[trIdx]
  
  Ytr <- rep(0,ntr)
  # Function ----------------------------------------------------------------
  for(k in 1:niter){
    t0 <- Sys.time()
    output <- Geostat.NNSeq(x=ttr, y=Ytr, data=Ztr, 
                            n.samples=M, nn=neigh,
                            method = "latent")
    t1 <- Sys.time() - t0
    # t2 <- Sys.time()
    # c(fit, betas, sigmas, taus, phis) %<-% Geostat.NNSeq(x=ttr, y=Ytr, data=Ztr, 
    #                                                      n.samples=M, nn=neigh,
    #                                                      method = "response")
    # t3 <- Sys.time() - t2
    t4 <- Sys.time()
    output <- tNngpCollapsed_NAda(t=ttr, Z=Ztr, 
                                  X=Xtr, 
                                  indLab=as.character(indstr),
                                  M=M, burnIn=burnInP, neigh=neigh,  
                                  beta0=beta0, logtheta0=logtheta0,
                                  rwSigma=rwSigma, rwSigmam=rwSigmam, gamma=0.05, madapt=2/3*M,
                                  alphaS=alphaS, betaS=betaS, alphaPhi=alphaPhi, betaPhi=betaPhi, 
                                  alphaT=alphaT, betaT=betaT,
                                  muB=muB, vB=vB, 
                                  verbose = verb, fileName = "",
                                  n_threads = 1)
    t5 <- Sys.time() - t4
    t6[k,] <- c(t1, t5)
  }
  timeslistAll[[i]] <- t6
}
tfinalsimulcollseq <- Sys.time() - tsimul

