# Source 
source("Simulation/Splines/R/NIndOneTrajSpline.R")

# Generation
set.seed(7777)

# Set values
NIND <- 1
NPOINTS <- 1000
ncovs <- 3
betaind <- rnorm(NIND, sd=2)
betacov <- rnorm(ncovs, sd=2)
smin = 0
smax = 10
kns = 9
deg = 2
lambda20 = .5
theta <- c(1, 1, 1)

# File output names
today <- Sys.time() %>% str_replace_all(., " ", "_") %>% 
  str_remove_all(., "-|:")
filenm <- paste("Simulation/Splines/Data/SimSplineGRW", NIND, NPOINTS, ncovs, kns, deg, today, sep="_")

simnm <- paste(filenm, "sim.txt", sep="_")
parnm <- paste(filenm, "pars.RData", sep="_")

# Run simulation
simOut <- SimData(J = NIND, TJ = rep(NPOINTS, NIND), 
                  beta0 = c(betaind, betacov),
                  smin = smin, smax = smax, kns = kns, deg = deg, 
                  lambda20 = lambda20,
                  sigma20 = theta[1], phi0 = theta[2], tau20 = theta[3])


# Save simulated dataset and parameters
datSim <- simOut$dat
parSim <- simOut$pars

write_csv(datSim, path=simnm)
save(parSim, file = parnm)

