
# Packages ----------------------------------------------------------------

require(Rcpp)
require(RcppEigen)
require(RcppArmadillo)
require(zeallot)
require(bayesplot)
require(magrittr)
require(tidyverse)
require(lubridate)
require(doParallel)
require(geoR)
require(rgeos)
require(grid)
require(matrixStats)
require(coda)
require(bigreadr)

# Data ---------------------------------------------------------------------

dat <- read_csv("Simulation/Experiment2/Data/sim_5_20000_3_2020-07-14_data.csv")

# Data prep 
n <- nrow(dat)

# Individuals idx
inds <- dat$IndIdx
# Response
Z <- dat$Z
# Times
t <- dat$t
# Design
X <- as.matrix(dat[, grep("X", colnames(dat))])

# Train and test
set.seed(130494)
pte <- 0.3
ntr <-  round((1-pte)*n)
trIdx <- sort(sample(1:n, ntr, replace = F))

indstr <- inds[trIdx]
Ztr <- Z[trIdx]
Xtr <- X[trIdx, ]
ttr <- t[trIdx]

indste <- inds[-trIdx]
Zte <- Z[-trIdx]
Xte <- X[-trIdx, ]
tte <- t[-trIdx] 

# Variogram ---------------------------------------------------------------
lmMod <- lm(Ztr~-1+Xtr)
wEsttr <- lmMod$residuals

uinds <- unique(indstr)
nInds <- length(uinds)
pars <- matrix(NA, 3, nInds)

for(i in 1:ncol(pars))
{
  idxs <- which(indstr==uinds[i])
  curve <- data.frame((wEsttr[idxs]))
  vario <- variog(coords = cbind(ttr[idxs], rep(0, length(idxs))), data = curve,
                  max.dist = 20, message=F)
  
  # grid values of theta vector
  initials <- seq(0.1, 2, length.out=10)
  fit <- variofit(vario, ini.cov.pars = expand.grid(initials, initials), 
                  cov.model = "exponential", message=F)
  sigma <- fit$cov.pars[1]
  phi <- fit$cov.pars[2]
  tau <- fit$nugget
  pars[,i] <- c(sigma, phi, tau)
}

theta0 <- rowMeans(pars)
theta0[theta0==0] <- 0.01

# Params ------------------------------------------------------------------

# Neighbors
neigh <- 10

# Priors
alphaS <- 2
betaS <- 2 # InvGamma on sigma2
alphaPhi <- 1
betaPhi <- 1 # Gamma on phi
alphaT <- 2
betaT <- 2 # InvGamma on tau2
ncov <- ncol(X)
muB <- rep(0, ncov)
vB <- diag(10^6, ncov)

# Chains ------------------------------------------------------------------

# Initial
logtheta0 <- log(theta0)
beta0 <- matrix(lmMod$coefficients, nrow = ncov, 1)

# Random walk
rwTheta <- (0.1)^2/length(logtheta0)
rwSigma <- diag(rep(rwTheta, length(logtheta0)))
rwSigmam <- rwSigma

# Function ----------------------------------------------------------------

Rcpp::sourceCpp("Simulation/Experiment2/Rcpp/tNngpColl_Final_AdaRW_Preds.cpp")

set.seed(7777)

M <- 10000
burnInP <- 0.5
burnIn <- M*burnInP+1
verb <- M/2

today <- Sys.time() %>% str_replace_all(., " ", "_") %>% 
  str_remove_all(., "-|:")
filenm <- paste("Simulation/Experiment2/WS/fit", M, length(uinds), n, ncov, pte, today, sep="_")

estnm <- paste(filenm, "est.RData", sep="_")
algonm <- paste(filenm, "_", sep="")
ncores <- parallel::detectCores()

print(paste("Test prop is", pte))
t0 <- Sys.time()
outPut <- tNngpCollapsed_NAda(t=ttr, Z=Ztr, X=Xtr, indLab=as.character(indstr),
                              M=M, burnIn=burnInP, neigh=neigh, 
                              beta0=beta0, logtheta0=logtheta0, 
                              rwSigma=rwSigma, rwSigmam=rwSigmam, gamma=0.05, madapt=2/3*M,
                              alphaS=alphaS, betaS=betaS, alphaPhi=alphaPhi, betaPhi=betaPhi, 
                              alphaT=alphaT, betaT=betaT,
                              muB=muB, vB=vB,
                              verbose = verb, fileName = algonm,
                              n_threads = max(c(ncores, length(uinds))))
testiimation <- Sys.time()-t0
print(testiimation)

c(thetas, betas, acc, wPredsInLast, loglik) %<-% outPut
rm(outPut)

# Chains burned
betasp <- as.matrix(betas[, (burnIn):(M)])
thetasp <- thetas[, (burnIn):(M)]

# Acceptance
print("Acceptance")
(accRate <- sum(acc)/M)
print("Acceptance Burned")
(accRateburn <- sum(acc[(burnIn+1):M])/(M/2))

# DIC
print("Average Deviance")
(Dbar <- mean(-2*loglik))
print("Variance of Deviance")
(Dvarbar <- var(-2*loglik))
print("DIC Gelman")
(DIC <- Dbar + 1/2*Dvarbar)

# Parameters
chains <- rbind(betasp, thetasp)
rm(betas)
rm(thetas)
row.names(chains) <- c(paste("beta[", 0:(ncov-1), "]", sep = ""), "sigma^2", "phi", "tau^2")

save(dat, trIdx, burnIn, chains, acc, wPredsInLast, loglik, testiimation, verb, algonm, neigh, file = estnm)
rm(dat)
rm(X)
rm(Z)
rm(acc)

print("Stats")
# Stats
a <- t(chains) %>%
  as.data.frame() %>%
  gather(Param, Value) %>%
  group_by(Param) %>%
  summarise(Min = min(Value),
            Q025 = quantile(Value, 0.025),
            Q05 = quantile(Value, 0.05),
            Median = median(Value),
            Mean = mean(Value),
            Q95 = quantile(Value, 0.95),
            Q975 = quantile(Value, 0.975),
            Max = max(Value))
write.table(a, file = paste(algonm, "parsTab.txt"))

# Emptying memory
rm(chains)

print("In Sample Preds")
# In sample predictions
yPredsIn <- eigenYPred(Xtr, betasp, wPredsInLast, thetasp[3, ])
rm(Xtr)

quantPreds <- cbind(rowQuantiles(yPredsIn, probs=c(0.025, 0.975)), rowMeans(yPredsIn))
rm(yPredsIn)
colnames(quantPreds) <- c("Q1", "Q2", "Mean")
predBounds <- quantPreds %>% as.data.frame() %>% select(Q1, Mean, Q2)
rm(quantPreds)

(PredCoverage <- mean(Ztr >= predBounds$Q1 & Ztr <= predBounds$Q2))
(MSE_InSamp <- mean((Ztr-predBounds$Mean)^2))
(rMSE_InSamp <- sqrt(MSE_InSamp))
(RelMSE_InSamp <- MSE_InSamp/var(Ztr))
(PIW_InSamp <- mean(predBounds$Q2 - predBounds$Q1))

print("Saving In Sample Preds")
save(predBounds, PredCoverage, MSE_InSamp, rMSE_InSamp, RelMSE_InSamp, PIW_InSamp, file = paste(filenm, "predsIn.RData", sep="_"))

rm(predBounds)

print("OOS Preds")
# OOS Preds

yPredsOutList <- list()

idx_filt <- unique(inds)

for(i in 1:nInds)
{
  trToPredIdx <- indstr==idx_filt[i]
  teToPredIdx <- indste==idx_filt[i]
  outPredsList <- tNngpCollapsed_Preds(t = ttr[trToPredIdx], 
                                       tPred = tte[teToPredIdx], Xpred = Xte[teToPredIdx,],
                                       neigh = neigh, index=i,
                                       betas = betasp, thetas = thetasp,
                                       wPredsIn = wPredsInLast[trToPredIdx, ],
                                       verbose = (ncol(betasp)/10), 
                                       n_threads = max(c(ncores, nInds)), fileName = algonm)
  yPredsOutList[[i]] <- outPredsList[[2]]
}

rm(wPredsInLast)
rm(outPredsList)

yPredsOutMat <- yPredsOutList %>% reduce(.f=rbind)
rm(yPredsOutList)

outQuantPreds <- cbind(rowQuantiles(yPredsOutMat, probs=c(0.025, 0.975)), rowMeans(yPredsOutMat))
rm(yPredsOutMat)
colnames(outQuantPreds) <- c("Q1", "Q2", "Mean")
outPredBounds <- outQuantPreds %>% as.data.frame() %>% select(Q1, Mean, Q2)
rm(outQuantPreds)

# Validation
(OutPredCoverage <- mean(Zte > outPredBounds$Q1 & Zte < outPredBounds$Q2))
(MSE_OutSamp <- mean((Zte-outPredBounds$Mean)^2))
(rMSE_OutSamp <- sqrt(MSE_OutSamp))
(RelMSE_OutSamp <- MSE_OutSamp/mean((Zte-mean(Ztr))^2))
(PIW_OutSamp <- mean(outPredBounds$Q2 - outPredBounds$Q1))

print("Saving OOS Preds")
save(outPredBounds, OutPredCoverage, MSE_OutSamp, rMSE_OutSamp, RelMSE_OutSamp, PIW_OutSamp, file = paste(filenm, "predsOut.RData", sep="_"))