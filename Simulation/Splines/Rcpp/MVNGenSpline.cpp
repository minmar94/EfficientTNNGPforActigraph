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


// [[Rcpp::export]]
arma::mat covFun(const arma::vec& tt, const float& sigma2, const float& phi) {
  
  int ndist = tt.size();
  arma::mat covout(ndist, ndist);
  covout(0,0) = sigma2;
  
  for(int i=0; i<ndist; i++){
    covout(i,i) = sigma2;
    for(int j=i+1; j<ndist; j++){
     covout(j, i) = sigma2*exp(-phi*(tt(j) - tt(i))); 
     covout(i, j) = covout(j, i);
    }
  }
  
  return covout;
}

// [[Rcpp::export]]
MatrixXd wrandomgen(MatrixXd Sigma2, int ntpoints){
  return((Sigma2.llt().matrixL())*zrnorm(ntpoints, 1));
}

// [[Rcpp::export]]
MatrixXd CholEigen(MatrixXd Sigma2){
  return(Sigma2.llt().matrixL());
}

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
  
  w = (Sigma2.llt().matrixL())*zrnorm(nts, 1);
  
  return w;
}