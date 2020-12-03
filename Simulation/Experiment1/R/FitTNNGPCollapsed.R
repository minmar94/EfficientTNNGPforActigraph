# Packages ----------------------------------------------------------------
require(Rcpp)
require(RcppEigen)
require(RcppArmadillo)
require(zeallot)
require(tidyverse)
require(Rfast2)


# Data ---------------------------------------------------------------------

datAll <- read_csv("Simulation/Experiment1/Data/sim_1_1e+05_3_20200723_040942_sim.txt")

dat <- datAll
n <- nrow(dat)

# Individuals idx
inds <- dat$IndIdx

set.seed(130494)
pte <- 0.3
ntr <-  round((1-pte)*n)
trIdx <- sort(sample(1:n, ntr, replace = F))

indstr <- inds[trIdx]
Ztr <- dat$Z[trIdx]
Xtr <- as.matrix(dat[trIdx, grepl("^X", colnames(dat))])
ttr <- dat$t[trIdx]

indste <- inds[-trIdx]
Zte <- dat$Z[-trIdx]
Xte <- as.matrix(dat[-trIdx, grepl("^X", colnames(dat))])
tte <- dat$t[-trIdx]

ncov <- ncol(Xtr)
uinds <- unique(indstr)
# Neighbors
neigh <- 10

# Priors
alphaS <- 2
betaS <- 2 # Gamma on sigma2
alphaPhi <- 1
betaPhi <- 1 # Gamma on phi
alphaT <- 2
betaT <- 2 # Gamma on tau2
muB <- rep(0, ncov)
vB <- diag(1, ncov)

# Chains ------------------------------------------------------------------

# Initial
logtheta0 <- c(2, 3, 3)
beta0 <- matrix(c(4), ncov, 1)

# Random walk
rwTheta <- (0.1)^2/length(logtheta0)
rwSigma <- diag(rep(rwTheta, length(logtheta0)))
rwSigmam <- rwSigma

# Function ----------------------------------------------------------------

sourceCpp("Simulation/Experiment1/Rcpp/tNngpColl_Final_AdaRW_Preds.cpp")

set.seed(7777)

M <- 10000
burnInP <- 0.5
burnIn <- M*burnInP+1
verb <- floor(M/2)

today <- Sys.time() %>% str_replace_all(., " ", "_") %>% 
  str_remove_all(., "-|:")
filenm <- paste("Simulation/Experiment1/WS/FitExp1", M, length(uinds), n, ncov, pte, today, sep="_")

estnm <- paste(filenm, "est.RData", sep="_")
algonm <- paste(filenm, "_", sep="")
plotnm <- paste(filenm, "plots.pdf", sep="_")
ncores <- parallel::detectCores()-2

print(paste("Test prop is", pte))
t0 <- Sys.time()
outPut <- tNngpCollapsed_NAda(t=ttr, Z=Ztr, 
                              X=as.matrix(Xtr), 
                              indLab=as.character(indstr),
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


# Results -----------------------------------------------------------------

# Save output
save(dat, trIdx, burnIn, outPut, testiimation, verb, algonm, neigh, file = estnm)
c(thetas, betas, acc, wPredsInLast, loglik) %<-% outPut

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
chains <- rbind(betas[, (burnIn):M], thetas[, (burnIn):M])

row.names(chains) <- c(paste("beta[", 0:(ncov-1), "]", sep = ""), 
                       "sigma^2", "phi", "tau^2")

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
thetasp <- as.matrix(thetas[, (M-verb+1):(M)])

yPredsInLast <- eigenYPred(as.matrix(Xtr), betasp, 
                           wPredsInLast, thetasp[3, ])
yPredsIn <- yPredsInLast
rm(yPredsInLast)
quantPreds <- cbind(rowQuantile(yPredsIn, probs=c(0.025, 0.975)), Rfast::rowmeans(yPredsIn))
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
save(predBounds, PredCoverage, MSE_InSamp, rMSE_InSamp, RelMSE_InSamp, PIW_InSamp, 
     file = paste(filenm, "predsInLast.RData", sep="_"))

rm(predBounds)

# Predictions OOS-------------------------------------------------------------

outPredsList <- tNngpCollapsed_Preds(t = ttr, 
                                     tPred = tte, Xpred = Xte,
                                     neigh = neigh, index=1,
                                     betas = betasp, thetas = thetasp,
                                     wPredsIn = wPredsInLast,
                                     verbose = ncol(betasp)/10, fileName = algonm,
                                     n_threads = max(c(ncores, length(uinds))))

yPredsOutMat <- outPredsList[[2]]
rm(outPredsList)
outQuantPreds <- cbind(rowQuantile(yPredsOutMat, probs=c(0.025, 0.975)), Rfast::rowmeans(yPredsOutMat))
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
save(outPredBounds, OutPredCoverage, MSE_OutSamp, rMSE_OutSamp, RelMSE_OutSamp, PIW_OutSamp,
     file = paste(filenm, "predsOutLast.RData", sep="_"))

