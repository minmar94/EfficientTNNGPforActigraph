require(zeallot)
require(tidyverse)
require(Rfast2)
require(spNNGP)

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

Ytr <- rep(0, ntr)
Yte <- rep(0, n - ntr)

# Neighbors
neigh <- 10

# Priors
alphaS <- 2
betaS <- 2 # Gamma on sigma2
lPhi <- .5
uPhi <- 30 # Uniform on phi
alphaT <- 2
betaT <- 2 # Gamma on tau2

ncores <- parallel::detectCores()-2

# Prior setting
priors1 <- list("beta.flat",
                "phi.Unif"=c(lPhi, uPhi), "sigma.sq.IG"=c(alphaS,betaS),
                "tau.sq.IG"=c(alphaT, betaT))
starting1 <- list("phi"= 5, "sigma.sq"= 0.5, "tau.sq"= 0.5,
                  "beta"=list(4,4,4,4))
tuning1 <- list("phi"=0.035, "sigma.sq"=0.035, "tau.sq"=0.035)
Geostat.NNSeq <- function(x, y, data, n.samples, nn, method)
{
  sp <- spNNGP::spNNGP(Z~-1+., coords=cbind(x, y), method=method, data = data,
                       family="gaussian", n.neighbors=nn,
                       starting=starting1,
                       tuning=tuning1, fit.rep = T, n.omp.threads = ncores,
                       priors=priors1, n.samples=n.samples, n.report = n.samples,
                       cov.model="exponential", verbose = T)
  beta0 <- sp$p.beta.samples
  sigma<-(sp$p.theta.samples[,1])
  tau<-(sp$p.theta.samples[,2])
  phi<-(sp$p.theta.samples[,3])
  
  return(list(fit=sp, beta0=beta0, sigma=sigma, tau=tau,phi=phi))
}

# Function ----------------------------------------------------------------

set.seed(7777)

M <- 100
burnInP <- 0.5
burnIn <- M*burnInP+1
verb <- floor(M/2)

today <- Sys.time() %>% str_replace_all(., " ", "_") %>% 
  str_remove_all(., "-|:")
filenm <- paste("Simulation/Experiment1/WS/fitSequential", M, 1, n, ncov, pte, today, sep="_")

estnm <- paste(filenm, "est.RData", sep="_")
algonm <- paste(filenm, "_", sep="")
plotnm <- paste(filenm, "plots.pdf", sep="_")

print(paste("Test prop is", pte))
t0 <- Sys.time()
outPut_sequential <- Geostat.NNSeq(x=ttr, y=Ytr, data=as.data.frame(cbind(Z=Ztr, Xtr)), n.samples=M, 
                                   nn=neigh, method = "latent")
test_sequential <- Sys.time() - t0
print(test_sequential)


fit_sequential <- outPut_sequential[[1]]
betas_sequential <- outPut_sequential[[2]]
sigmas_sequential <- outPut_sequential[[3]]
taus_sequential <- outPut_sequential[[4]]
phis_sequential <- outPut_sequential[[5]]


# In sample

pred_sequential_InSample <- fit_sequential$y.rep.samples # already burned
outQuantPredsIn <- rowQuantile(pred_sequential_InSample, probs=c(0.025, 0.975))
outQuantPredsIn <- cbind(outQuantPredsIn, rowMeans(pred_sequential_InSample))
(cv_sequential_InSample <- sum(Ztr > outQuantPredsIn[,1] & Ztr < outQuantPredsIn[,2])/(ntr))
(MSE_InSamp_sequential <- sum((Ztr-outQuantPredsIn[,3])^2)/(ntr))
(rMSE_InSamp_sequential <- sqrt(MSE_InSamp_sequential))
(RelMSE_InSamp_sequential <- MSE_InSamp_sequential/var(Ztr))
(PIW_sequential_InSample <- mean(outQuantPredsIn[,2] - outQuantPredsIn[,1]))


# Parameters
chains_sequential <- rbind(t(betas_sequential)[, (burnIn+1):M], t(fit_sequential$p.theta.samples)[, (burnIn+1):M])
#chains <- rbind(betas[, 1:M], thetas[, 1:M])
row.names(chains_sequential) <- c(paste("beta[", 0:(ncov-1), "]", sep = ""), "sigma^2", "phi", "tau^2")

print("Stats")
# Stats
a <- t(chains_sequential) %>%
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
            Max = max(Value)) %>% 
  filter(Q05*Q95>0)
write.table(a, paste(algonm, "parstab.txt", sep=""))

print("Saving In Sample Preds")
save(outQuantPredsIn, cv_sequential_InSample, MSE_InSamp_sequential, RelMSE_InSamp_sequential, rMSE_InSamp_sequential,
     PIW_sequential_InSample, test_sequential,
     file = paste(filenm, "predsInLast.RData", sep="_"))

rm(outQuantPredsIn)


# Predictions OOS-------------------------------------------------------------

t0pred <- Sys.time()
pred_sequential <- predict(object = fit_sequential, X.0 = Xte,  coords.0 = cbind(tte, Yte),
                           verbose = F, n.omp.threads = ncores)
tpred_final_sequential <- Sys.time()-t0pred
outQuantPreds <- rowQuantile(pred_sequential$p.y.0[,(burnIn):M], probs=c(0.025, 0.975))
outQuantPreds <- cbind(outQuantPreds, rowMeans(pred_sequential$p.y.0[,(burnIn):M]))
(cv_OutSamp_sequential <- sum(Zte > outQuantPreds[,1] & Zte < outQuantPreds[,2])/(n - ntr))
(MSE_OutSamp_sequential <- sum((Zte-outQuantPreds[3,])^2)/(n - ntr))
(rMSE_OutSamp_sequential <- sqrt(MSE_OutSamp_sequential))
(RelMSE_OutSamp_sequential <- MSE_OutSamp_sequential/mean((Zte-mean(Ztr))^2))
(PIW_sequential <- mean(outQuantPreds[,2] - outQuantPreds[,1]))


save(outQuantPreds, cv_OutSamp_sequential, MSE_OutSamp_sequential, rMSE_OutSamp_sequential, RelMSE_OutSamp_sequential,
     tpred_final_sequential, PIW_sequential,
     file = paste(filenm, "predsOutLast.RData", sep="_"))
