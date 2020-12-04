# Efficient estimation of temporal Gaussian processes models for physical activity level in a large scale population study
Optimizing the NNGP methodology for the uni-dimensional/temporal setting.
Application to health study conducted by the Fielding School of Public Health of the University of California Los Angeles.

## Experiment 1

  - <b>Data</b>: contains a simulated dataset for one individual with 10^5 data points and 3 covariates. True parameters used for the simulation are in "\_pars.RData"
  - <b>R</b>:
    - "NIndOneTraj.R" function to simulate a temporal Gaussian process for N individuals
    - "SimulateData.R" simulate the dataset 
    - "CollapsedVSFull.R" comparison between Collapsed NNGP and Full GP for one sampler iteration as the sample size increases
    - "CollapsedVSSequential.R" \rightarrow comparison between Collapsed NNGP and Sequential NNGP for one sampler iteration as the sample size increases
    - "FitTNNGPCollapsed.R" fit of the Collapsed NNGP on the simulated dataset
    - "FitTNNGPSequential.R" fit of the Sequential NNGP on the simulated dataset
  - <b>Rcpp</b>:
    - "MVGen.cpp" auxiliary functions for simulation of a temporal Gaussian process
    - "tNngpColl_Final_AdaRW_Preds.cpp" main functions for estimation and prediction of the Collapsed NNGP 
  - <b>WS</b>: Workspace folder to store the results
    
