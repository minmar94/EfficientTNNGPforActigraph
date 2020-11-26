
# Packages ----------------------------------------------------------------

require(Rcpp)
require(RcppEigen)
require(RcppArmadillo)
require(zeallot)
require(bayesplot)
require(magrittr)
require(tidyverse)
require(lubridate)
require(geoR)
require(rgeos)
require(grid)
require(matrixStats)

# Data ---------------------------------------------------------------------

dat <- read_csv("Simulation/Splines/Data/SimSplineGRW_5_20000_3_9_2_20201029_175847_sim.txt")
load("Simulation/Splines/Data/SimSplineGRW_5_20000_3_9_2_20201029_175847_pars.RData")

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
S <- as.matrix(dat[, grep("S", colnames(dat))])
ncov <- ncol(X)
nspl <- ncol(S)
XS <- cbind(X, S)
rm(X)
rm(S)

# Train and test
set.seed(130494)
pte <- 0.3
ntr <-  round((1-pte)*n)
trIdx <- sort(sample(1:n, ntr, replace = F))

indstr <- inds[trIdx]
Ztr <- Z[trIdx]
XStr <- XS[trIdx, ]
ttr <- t[trIdx]

indste <- inds[-trIdx]
Zte <- Z[-trIdx]
XSte <- XS[-trIdx, ]
tte <- t[-trIdx] 

# Variogram ---------------------------------------------------------------
lmMod <- lm(Ztr~-1+XStr)
wEsttr <- lmMod$residuals

uinds <- unique(indstr)
nInds <- length(uinds)
pars <- matrix(NA, 3, length(uinds))

for(i in 1:ncol(pars))
{
  idxs <- which(indstr==uinds[i])
  curve <- data.frame((wEsttr[idxs]))
  vario <- variog(coords = cbind(ttr[idxs], rep(0, length(idxs))), data = curve,
                  max.dist = 20, messages=F)
  
  # grid values of theta vector
  initials <- seq(0.1, 2, length.out=5)
  fit <- variofit(vario, ini.cov.pars = expand.grid(initials, initials), 
                  cov.model = "exponential", messages=F)
  sigma <- fit$cov.pars[1]
  phi <- fit$cov.pars[2]
  tau <- fit$nugget
  pars[,i] <- c(sigma, phi, tau)
}

rm(curve, vario, initials, fit)

theta0 <- apply(pars, 1, median)
theta0[theta0==0] <- 0.01

rm(pars)

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
muB <- rep(0, ncov)
vB <- diag(10^6, ncov)
muE <- rep(0, nspl)

KLon <- diag(2, sqrt(nspl))
KLon[1, 1] <- 1
KLon[sqrt(nspl), sqrt(nspl)] <- 1
KLon[abs(row(KLon)-col(KLon))==1] <- -1
KLat <- diag(2, sqrt(nspl))
KLat[1, 1] <- 1
KLat[sqrt(nspl), sqrt(nspl)] <- 1
KLat[abs(row(KLat)-col(KLat))==1] <- -1
T1 <- kronecker(diag(1, sqrt(nspl)), KLat)
T2 <- kronecker(KLon, diag(1, sqrt(nspl)))

KE <- T1+T2


# Chains ------------------------------------------------------------------

# Initial
logtheta0 <- log(theta0)
beta0 <- matrix(lmMod$coefficients[1:ncov], nrow = ncov, 1)
eta0 <- matrix(lmMod$coefficients[(ncov+1):(ncov+nspl)], nrow = nspl, 1)

# Random walk
rwTheta <- (0.1)^2/length(logtheta0)
rwSigma <- diag(rep(rwTheta, length(logtheta0)))
rwSigmam <- rwSigma

# Function ----------------------------------------------------------------

Rcpp::sourceCpp("Simulation/Splines/Rcpp/tNngpColl_PSpline.cpp")

set.seed(7777)

M <- 10000
burnInP <- 0.5
burnIn <- M*burnInP+1
verb <- floor(M/2)

today <- Sys.time() %>% str_replace_all(., " ", "_") %>% 
  str_remove_all(., "-|:")
filenm <- paste("Simulation/Splines/WS/fitPSplineGRW", M, length(uinds), n, ncov, nspl, pte, today, sep="_")

estnm <- paste(filenm, "est.RData", sep="_")
algonm <- paste(filenm, "_", sep="")
plotnm <- paste(filenm, "plots.pdf", sep="_")
ncores <- parallel::detectCores()

print(paste("Test prop is", pte))
t0 <- Sys.time()
outPut <- tNngpCollapsed_NAda(t=ttr, Z=Ztr, 
                              XS=as.matrix(XStr), ncov=ncov, nspl=nspl, 
                              indLab=as.character(indstr),
                              M=M, burnIn=burnInP, neigh=neigh,  
                              beta0=beta0, logtheta0=logtheta0, eta0=eta0, lambda20=1,
                              rwSigma=rwSigma, rwSigmam=rwSigmam, gamma=0.05, madapt=2/3*M,
                              alphaS=alphaS, betaS=betaS, alphaPhi=alphaPhi, betaPhi=betaPhi, 
                              alphaT=alphaT, betaT=betaT,
                              muB=muB, vB=vB, muE, KE=KE,
                              alphal0 = 1/2, betal0 = 1/2,
                              verbose = verb, fileName = algonm,
                              n_threads = max(c(ncores, nInds)))
testiimation <- Sys.time()-t0
print(testiimation)

# Save output
save(dat, trIdx, burnIn, outPut, testiimation, verb, algonm, neigh, file = estnm)
c(thetas, betas, etas, lambdas, acc, wPredsInLast, loglik) %<-% outPut

rm(outPut)

# Acceptance
print("Acceptance")
(accRate <- sum(acc)/M)
print("Acceptance Burned")
(accRateburn <- sum(acc[(burnIn):M])/(M-burnIn+1))

# DIC
print("Average Deviance")
(Dbar <- mean(-2*loglik))
print("Variance of Deviance")
(Dvarbar <- var(-2*loglik))
print("DIC Gelman")
(DIC <- Dbar + 1/2*Dvarbar)

# Parameters
chains <- rbind(betas[, (burnIn):M], etas[, (burnIn):M], thetas[, (burnIn):M], lambdas[(burnIn):M])

row.names(chains) <- c(paste("beta[", 0:(ncov-1), "]", sep = ""), 
                       paste("eta[", 0:(nspl-1), "]", sep = ""),
                       "sigma^2", "phi", "tau^2", "lambda")

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
write.table(a, paste(algonm, "parstab.txt", sep=""))



# Predictions in Sample ---------------------------------------------------
betasp <- as.matrix(betas[, (M-verb+1):M])
etasp <- as.matrix(etas[, (M-verb+1):M])
thetasp <- as.matrix(thetas[, (M-verb+1):(M)])

yPredsInLast <- eigenYPred(as.matrix(XStr), rbind(betasp, etasp), 
                           wPredsInLast, thetasp[3, ])
yPredsIn <- yPredsInLast
rm(yPredsInLast)
quantPreds <- cbind(rowQuantiles(yPredsIn, probs=c(0.025, 0.975)), Rfast::rowmeans(yPredsIn))
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
save(predBounds, PredCoverage, MSE_InSamp, rMSE_InSamp, RelMSE_InSamp, PIW_InSamp, file = paste(filenm, "predsInLast.RData", sep="_"))

rm(predBounds)

# Predictions OOS-------------------------------------------------------------

yPredsOutList <- list()

for(i in 1:nInds)
{
  trToPredIdx <- indstr==uinds[i]
  teToPredIdx <- indste==uinds[i]
  outPredsList <- tNngpCollapsed_Preds(t = ttr[trToPredIdx], 
                                       tPred = tte[teToPredIdx], XSpred = XSte[teToPredIdx,],
                                       neigh = neigh, index=i,
                                       betas = betasp, etas = etasp, thetas = thetasp,
                                       wPredsIn = wPredsInLast[trToPredIdx, ],
                                       verbose = ncol(betasp)/10, fileName = algonm,
                                       n_threads = max(c(ncores, nInds)))

  yPredsOutList[[i]] <- outPredsList[[2]]
}

yPredsOutMat <- yPredsOutList %>% reduce(.f=rbind)
rm(yPredsOutList)
outQuantPreds <- cbind(rowQuantiles(yPredsOutMat, probs=c(0.025, 0.975)), Rfast::rowmeans(yPredsOutMat))
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
save(outPredBounds, OutPredCoverage, MSE_OutSamp, rMSE_OutSamp, RelMSE_OutSamp, PIW_OutSamp, file = paste(filenm, "predsOutLast.RData", sep="_"))
