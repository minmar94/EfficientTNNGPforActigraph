# Source 
source("Simulation/Experiment1/R/NIndOneTraj.R")


# Generation
set.seed(7777)
NIND <- 1
NPOINTS <- 10^5
ncov <- 3
betaind <- rnorm(NIND)
betacov <- rnorm(ncov)
theta <- c(1, 1, 1)

today <- Sys.time() %>% str_replace_all(., " ", "_") %>% 
  str_remove_all(., "-|:")
filenm <- paste("Simulation/Experiment1/Data/sim", NIND, NPOINTS, ncov, today, sep="_")

simnm <- paste(filenm, "sim.txt", sep="_")
parnm <- paste(filenm, "pars.RData", sep="_")

DataSim <- SimData(J = NIND, TJ = rep(NPOINTS, NIND), beta0 = c(betaind, betacov),
                   sigma20 = theta[1], phi0 = theta[2], tau20 = theta[3])

write_csv(DataSim, path=simnm)
save(betaind, betacov, file = parnm)