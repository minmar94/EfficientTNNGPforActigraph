# Efficient estimation of temporal Gaussian processes models for physical activity level in a large scale population study
Optimizing the NNGP methodology for the uni-dimensional/temporal setting.
Application to health study conducted by the Fielding School of Public Health of the University of California Los Angeles.

## Experiment 1

  - <b>Data</b>: contains a simulated dataset for one individual with 10^5 data points and 3 covariates. True parameters used for the simulation are in "\_pars.RData"
  - <b>R</b>:
    - "CollapsedVSFull.R" is the code for comparison between Collapsed NNGP and Full GP for one sampler iteration as the sample size increases
    - "CollapsedVSSequential.R" is the code for comparison between Collapsed NNGP and Sequential NNGP for one sampler iteration as the sample size increases
  - Rcpp
  - WS
    
