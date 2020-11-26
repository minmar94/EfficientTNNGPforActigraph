#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Ziggurat.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppZiggurat)]]
// [[Rcpp::depends(RcppEigen)]]

static Ziggurat::Ziggurat::Ziggurat zigg;

using namespace Eigen; 
using namespace Rcpp;

// Defining useful types
typedef Eigen::Triplet<double> T;
typedef Eigen::VectorXd ColVecd;
typedef Eigen::Map<ColVecd> MapVecd;
typedef Eigen::MatrixXd Matd;
typedef Eigen::Map<Matd> MapMatd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseLU<SpMat> LU;
typedef Eigen::SimplicialLDLT<SpMat, Lower, NaturalOrdering<int>> CholLDL;
typedef Eigen::SimplicialLLT<SpMat, Lower, NaturalOrdering<int>> CholLL;


// Wrapper for random generation of n x j Standard Normal variables
// [[Rcpp::export]] 
MatrixXd zrnorm(int n, int j) {
  MatrixXd x(n, j);
  for (int i=0; i<n; i++) {
    for(int k=0; k<j; k++){
      x(i, k) = zigg.norm();
    }
  }
  return x;
}

// Wrapper for random generation of a ntpoints x 1 vector from a multivariate Normal with Sigma2 as cross-covariance matrix
MatrixXd wrandomgen(MatrixXd Sigma2, int ntpoints){
  return((Sigma2.llt().matrixL())*zrnorm(ntpoints, 1));
}


// Simulation of the realization of a temporal Gaussian process with mean 0 and exponential covariance function of parameters sigma2 and phi over the set t
// [[Rcpp::export]]
ColVecd GPSimul(const ArrayXd& t, const float& sigma2, const float& phi) {
  
  int nts = t.size();
  ColVecd w(nts);
  MatrixXd Sigma2(nts, nts);

  for(int i=0; i<nts; i++){
    Sigma2(i,i) = sigma2;
    for(int j=i+1; j<nts; j++){
      Sigma2(j, i) = sigma2*exp(-phi*(t(j) - t(i))); 
      Sigma2(i, j) = Sigma2(j, i);
    }
  }
  
  w = wrandomgen(Sigma2, nts);
  
  return w;
}