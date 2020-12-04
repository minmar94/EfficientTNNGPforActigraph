# Efficient estimation of temporal Gaussian processes models for physical activity level in a large scale population study
Optimizing the NNGP methodology for the uni-dimensional/temporal setting.
Application to health study conducted by the Fielding School of Public Health of the University of California Los Angeles.

## Simulation
### Experiment 1

  - <b>Data</b>: contains a simulated dataset for one individual with 10^5 data points and 3 covariates. True parameters used for the simulation are in "\_pars.RData"
  - <b>R</b>:
    - "NIndOneTraj.R" function to simulate a temporal Gaussian process for N individuals
    - "SimulateData.R" simulate the dataset 
    - "CollapsedVSFull.R" comparison between Collapsed NNGP and Full GP for one sampler iteration as the sample size increases
    - "CollapsedVSSequential.R" \rightarrow comparison between Collapsed NNGP and Sequential NNGP for one sampler iteration as the sample size increases
    - "FitTNNGPCollapsed.R" fit the Collapsed NNGP on the simulated dataset
    - "FitTNNGPSequential.R" fit the Sequential NNGP on the simulated dataset
  - <b>Rcpp</b>:
    - "MVGen.cpp" auxiliary functions for simulation of a temporal Gaussian process
    - "tNngpColl_Final_AdaRW_Preds.cpp" main functions for estimation and prediction of the Collapsed NNGP 
  - <b>WS</b>: Workspace folder to store the results
  
 ### Splines

  - <b>Data</b>: contains a simulated dataset for 5 individuals with 2\cdot 10^5 data points and 3 covariates, with a spline effect on the mean. 
  True parameters used for the simulation are in "\_pars.RData"
  - <b>R</b>:
    - "NIndOneTrajSpline.R" function to simulate a temporal Gaussian process for N individuals with a spline effect on the mean
    - "SplineSimulation.R" simulate the dataset 
    - "SimFitPS_Last.R" fit the Collapsed NNGP on the simulated dataset with penalized splines
    - "SimFitPS_Last.R" fit the Sequential NNGP on the simulated dataset with shrinking splines
  - <b>Rcpp</b>:
    - "MVGenSpline.cpp" auxiliary functions for simulation of a temporal Gaussian process with a spline effect on the mean.
    - "tNngpColl_PSpline.cpp" main functions for estimation and prediction of the Collapsed NNGP with penalized splines
    - "tNngpColl_Shrink.cpp" main functions for estimation and prediction of the Collapsed NNGP with shrinking splines
  - <b>WS</b>: Workspace folder to store the results
    
