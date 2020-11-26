#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Ziggurat.h>
#include <iostream>
#include <fstream>
#include <vector> 
#include <string> 
#include <algorithm> 
#include <sstream> 
#include <iterator> 
#include <math.h>
#include <chrono>
#include <limits>
#include <random>
// [[Rcpp::depends(RcppZiggurat)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]


// Namespaces
using namespace Eigen; 
using namespace std::chrono;
using namespace std;

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

// Static Objects
static Ziggurat::Ziggurat::Ziggurat zigg;

// Auxiliary functions

// Math

// [[Rcpp::export]]
Matd eigenMapMatMult(const MapMatd A, const MapMatd B){
  Matd C = A * B;
  
  return C;
}

// [[Rcpp::export]]
Matd eigenMapMatSum(const MapMatd A, const MapMatd B){
  Matd C = A + B;
  
  return C;
}

// [[Rcpp::export]]
Matd eigenYPred(const MapMatd X, const MapMatd Beta, const MapMatd w, const MapVecd Tau){
  const int n = w.rows();
  const int M = w.cols();
  
  Matd C(n, M);
  for (int i=0; i<n; i++){
    for (int j=0; j<M; j++){
      C(i, j) = X.row(i) * Beta.col(j) + w(i,j) + sqrt(Tau(j))*zigg.norm();
    }
  }
  
  return C;
}


// Priors
double LogInvGamma(double x, double alpha, double beta) {
  
  double out = alpha*log(beta)-lgamma(alpha)-(alpha+1)*log(x)-beta/x;
  
  return out;
}
double LogGamma(double x, double alpha, double beta) {
  
  double out = alpha*log(beta)+(alpha-1)*log(x)-beta*x-lgamma(alpha);
  
  return out;
}
double LogUnif(double x, double Low, double Up) {
  
  double out = -log(Up-Low);
  
  return out;
}


// Finding neighbors of Out-Of-Sample data
Array2Xi NewNeighbors(const ArrayXd& x, const ArrayXd& y, const int& neigh) {
  
  // Data
  const int& nx = x.size();
  const int& ny = y.size();
  
  // Auxiliary counters
  int basex = 0;
  
  // Neighborhood limits of adjacent observations (transitivity over time)
  Array2Xi neighLims(2, ny);
  
  // Looping on new observations
  for (int i = 0; i<ny; i++){
    
    // Auxiliary condition
    int furthest=ny;
    int ncur=0;
    int jmax=0;
    double dmax = abs(y(i)-x(0));
    bool cond {true};
    
    // Looping on old observations
    for (int j = basex; j<nx && cond; j++){
      
      if (x(j)<y(i))
      {
        jmax = std::max(0, j-neigh+1);
        furthest = j;
        dmax = abs(y(i)-x(jmax));
        ncur = j-jmax+1;
      }else
      {
        double d = abs(y(i)-x(j));
        
        if (ncur<neigh)
        {
          furthest = j;
          ncur++;
          if (d>dmax)
          {
            dmax = d;
            furthest = j+(neigh-ncur);
            cond = false;
          }
        }else
        {
          if (d<dmax)
          {
            furthest = j;
            double dmax2 = abs(y(i)-x(jmax+1));
            
            if (d<dmax2)
            {
              jmax++;
              dmax = dmax2;
            }else
            {
              cond = false;
            }
          }else
          {
            cond = false;
          }
        }
      }
    }
    neighLims.col(i) << furthest-neigh+1, furthest;
    basex = max(furthest-neigh, 0);
  }
  
  // Returning output
  return(neighLims);
}

// Random generations
// [[Rcpp::export]]
Matd zrnorm(int n, int j) {
  
  Matd x(n, j);
  
  for (int i=0; i<n; i++) {
    for(int k=0; k<j; k++){
      x(i, k) = zigg.norm();
    }
  }
  return x;
}

// NNGP Collapsed algorithm N Inds
// [[Rcpp::export]]
Rcpp::List tNngpCollapsed_NAda(const Eigen::Map<ArrayXd> t, const MapVecd Z, 
                               const MapMatd& XS, const int ncov, const int nspl, 
                               const vector<string> indLab, 
                               const int& M, const double& burnIn, const int neigh, 
                               ColVecd beta0, ColVecd logtheta0, ColVecd eta0, double lambda20,
                               const MapMatd rwSigma, Matd rwSigmam, const double& gamma, const double& madapt,
                               const double& alphaS, const double& betaS, const double& alphaPhi, const double& betaPhi, const double& alphaT, const double& betaT,
                               const MapVecd muB, const MapMatd vB, const MapVecd muE, const MapMatd& KE,
                               const double& alphal0, const double& betal0,
                               const int verbose, const string& fileName, 
                               const int n_threads) {
  IOFormat ChainFmt(8, 0, ",", "", "", "\n", "", "");
  
  ofstream monitor;
  string mname = fileName+"monitor.txt";
  ofstream chainsFile;
  ofstream predsFile;
  ofstream rwSigmaFile;
  string sname = fileName+"rwSigma.txt";
  
  // Data infos
  const int& n = Z.size();
  const int& nall = ncov+nspl;
  
  const Matd& tXS = XS.adjoint();
  const Matd& tXSXS = (tXS*XS);
  const ColVecd& tXSZ = (tXS*Z);
  
  // obs per individual, computing lowi.
  // All okay as long as observations are ordered first by individual, second by time
  std::vector<string> inds;
  std::unordered_map<string, int> nks;
  std::unordered_map<string, std::vector<int>> lowj;
  
  // Counting
  int nInd = 0;
  for (int k=0; k<n; k++){
    // lowi
    VectorXi temp(2);
    temp << 0, nks[indLab[k]]-(neigh-1);
    lowj[indLab[k]].push_back(temp.maxCoeff());
    
    // Counts for individual
    nks[indLab[k]]++;
    
    // Etichette uniche
    if (nks[indLab[k]]==1){
      inds.push_back(indLab[k]);
      nInd++;
    }
  }
  
  // Auxiliary
  // Random Walk
  const int& d = rwSigma.cols();
  ArrayXd acc = ArrayXd::Zero(M);
  LLT<Matrix<double, 3, 3>> lltRW(rwSigma);
  Matd Lrw = lltRW.matrixL();
  arma::vec mhUnifs = arma::randu(M-1);
  // Adaptive
  double accCur=0;
  LLT<Matrix<double, 3, 3>> lltRWm(rwSigmam);
  Matd Lrwm = lltRWm.matrixL();
  arma::vec adaUnifs = arma::randu(M-1);
  // Gibbs
  const Matd& FmatB = vB.inverse();
  const ColVecd& fmatB = vB.ldlt().solve(muB);
  Matd FmatE = lambda20*KE;
  ColVecd KEmuE = KE*muE;
  ColVecd fmatE = lambda20*KEmuE;
  Matd Fmat = Matd::Zero(nall, nall);
  Fmat.topLeftCorner(ncov, ncov) = FmatB;
  ColVecd fmat(nall);
  fmat.head(ncov) = fmatB;
  
  const double& alphal = alphal0 + 1/2;
  
  const int mBurned = ceil(M*burnIn);
  
  //Starting points
  double sigma20 = exp(logtheta0(0));
  double phi0 = exp(logtheta0(1));
  double tau20 = exp(logtheta0(2));
  
  // Chain definition
  Array<double, 3, Dynamic> lthetas(3, M);
  ArrayXXd betas(ncov, M);
  ArrayXXd etas(nspl, M);
  ColVecd lambda2s(M);
  Matrix<double, Dynamic, Dynamic> wPreds(n, verbose);
  lthetas.col(0) << logtheta0(0), logtheta0(1), logtheta0(2);
  betas.col(0) = beta0;
  etas.col(0) = eta0;
  ColVecd betaeta(nall);
  betaeta << beta0, eta0;
  ArrayXd loglik(M-mBurned);
  
  // First loop on idividuals
  std::unordered_map<string, SpMat> Omegas;
  std::unordered_map<string, ArrayXd> allDist;
  double logd0 = 0;
  Matd v(n, nall);
  int cselv=0;
  std::vector<int> csels(nInd);
  
  for (int k=0; k<nInd; k++){
    
    csels[k] = cselv;
    // Distances
    const int& nk = nks[inds[k]];
    
    // Covariance matrix and distances
    // Defining covariance object
    ArrayXd distk(nk*neigh-neigh*(neigh+1)/2);
    
    std::vector<T> tripletListC;
    tripletListC.reserve(nk*neigh);
    
    // Initializing auxiliary counter and first value of covariance
    int cntc = 0;
    tripletListC.push_back(T(0, 0, sigma20));
    
    // Looping over rows
    for (int i = 1; i<nk; i++){
      // Assigning values on the diagonal
      tripletListC.push_back(T(i, i, sigma20));
      
      // Looping over cols
      for (int j=lowj[inds[k]][i-1]; j<i; j++){
        // Distances computing
        distk(cntc) = t(cselv+i)-t(cselv+j);
        // Covariance computing
        const double cij=sigma20 * exp(-distk(cntc)*phi0);
        tripletListC.push_back(T(i, j, cij));
        tripletListC.push_back(T(j, i, cij));
        // Updating counter
        cntc++;
      }
    }
    allDist[inds[k]] = distk;
    // Building sparse covariance matrix
    SpMat C(nk, nk);
    
    C.setFromTriplets(tripletListC.begin(), tripletListC.end());
    
    // AD Function
    std::vector<T> tripletListIAt;
    tripletListIAt.reserve(nk*neigh);
    std::vector<T> tripletListSDIA;
    tripletListSDIA.reserve(nk*neigh);
    
    ColVecd Dvec(nk);
    Dvec(0) = sigma20;
    double logD = log(sigma20);
    
    tripletListIAt.push_back(T(0, 0, 1));
    tripletListSDIA.push_back(T(0, 0, 1/Dvec(0)));
    for(int i = 0; i < (nk-1); i++)
    {
      // Auxiliary objects to store necessary matrices blocks
      Matd Cblock1 = C.block(lowj[inds[k]][i], lowj[inds[k]][i], 
                             i-lowj[inds[k]][i]+1, i-lowj[inds[k]][i]+1);
      Matd Cblock2 = C.block(lowj[inds[k]][i], i+1, i-lowj[inds[k]][i]+1, 1);
      ColVecd temp = Cblock1.ldlt().solve(Cblock2);
      
      Dvec(i+1) = (sigma20 - (C.block(i+1, lowj[inds[k]][i], 
                              1, i-lowj[inds[k]][i]+1)*temp)(0));
      logD += log(Dvec(i+1));
      
      for (int j = lowj[inds[k]][i]; j < i+1; j++){
        tripletListIAt.push_back(T(j, i+1, -temp(j-lowj[inds[k]][i])));
        tripletListSDIA.push_back(T(i+1, j, -temp(j-lowj[inds[k]][i])/Dvec(i+1)));
      }
      tripletListIAt.push_back(T(i+1, i+1, 1));
      tripletListSDIA.push_back(T(i+1, i+1, 1/Dvec(i+1)));
    }
    SpMat IAt(nk, nk);
    SpMat SolveDIA(nk, nk);
    IAt.setFromTriplets(tripletListIAt.begin(), tripletListIAt.end());
    SolveDIA.setFromTriplets(tripletListSDIA.begin(), tripletListSDIA.end());
    
    // Computing Omega
    SpMat Omega = (IAt.triangularView<UnitUpper>()*(SolveDIA));
    for (int i = 0; i < nk; i++)
    {
      Omega.coeffRef(i, i) += 1/tau20;
    }
    Omegas[inds[k]] = Omega;
    
    // Solving Omega
    CholLL solverO(Omega);
    
    SpMat L = solverO.matrixL();
    
    // First piece of log-likelihood d0
    logd0 += nk*log(tau20) + logD + 2*L.diagonal().array().log().matrix().sum();
    
    // v for Gibbs on all covs
    v.block(cselv, 0,  nk, nall) = solverO.solve(XS.block(cselv, 0, nk, nall));
    
    cselv += nk;
    
  }
  
  Matd tXSv = (tXS*v);
  ColVecd vtZ = (v.transpose()*Z);

  omp_set_num_threads(n_threads);
  
  monitor.open(mname);
  monitor << "Start looping!" << std::endl;
  monitor.close();
  
  // Chain loops
  for(int m = 1; m < M; m++){
    if (m%verbose==0)
    {
      Rcpp::Rcout << "Iteration" << std::endl << m << std::endl;

      monitor.open(mname.c_str(), std::ios_base::app);
      monitor << "Iteration: " << std::endl << m << std::endl;
      monitor << "Acceptance rate: " << std::endl << accCur/verbose << std::endl;
      monitor.close();
      
      rwSigmaFile.open(sname.c_str());
      rwSigmaFile << rwSigmam;
      rwSigmaFile.close();
      rwSigmaFile.clear();
      
      string cname = fileName + std::to_string(m/verbose) + "_chains.txt";
      chainsFile.open(cname.c_str());
      MatrixXd curChains(lthetas.rows()+betas.rows()+etas.rows()+1, verbose);
      curChains << lthetas.leftCols(m).rightCols(verbose), betas.leftCols(m).rightCols(verbose), 
                   etas.leftCols(m).rightCols(verbose), lambda2s.head(m).tail(verbose).transpose();
      chainsFile << curChains.format(ChainFmt) << "\n";
      chainsFile.close();
      
      accCur=0;
    }
    
    
    // Computing residuals
    const ColVecd& mm = XS*betaeta;
    const ColVecd& r = Z - mm;
    
    // Proposal for m-h on theta: adaptation
    Matd LrwCur = Lrw;
    if (adaUnifs(m-1)>gamma && m>2*d){
      if (m<=madapt)
      {
        Matd lthetasCur = lthetas.leftCols(m-1).rightCols(m*2/3);
        Matd centered = lthetasCur.colwise() - lthetasCur.rowwise().mean();
        rwSigmam = ((centered * centered.adjoint()) / double(lthetasCur.cols() - 1))*pow(2.38, 2)/d;
        lltRWm.compute(rwSigmam);
        Lrwm = lltRWm.matrixL();
      }
      LrwCur = Lrwm;
    }
    
    // Proposal for m-h on theta: proposal
    ColVecd logtheta1 = logtheta0 + LrwCur*zrnorm(3, 1);
    double sigma21 = exp(logtheta1(0));
    double phi1    = exp(logtheta1(1));
    double tau21   = exp(logtheta1(2));
    
    // Initializing objects to store individuals values
    std::unordered_map<string, SpMat> Omega1s;
    double logd1 = 0;
    ColVecd v0(n);
    ColVecd v1(n);
    Matd vNew(n, nall);
    
    // Individual objects
    CholLL solverO;
    string indsk;
    int nk;
    int cselsk;
    std::vector<int> lowjk;
    ArrayXd distk;
    
    // Loop over individuals
    #pragma omp parallel for private(solverO, indsk, nk, cselsk, lowjk, distk) reduction(+:logd1) schedule(auto)
    for (int k=0; k<nInd; k++){
      
      indsk = inds[k];
      nk = nks[indsk];
      lowjk = lowj[indsk];
      cselsk = csels[k];
      distk = allDist[inds[k]];
        
      solverO.compute(Omegas[indsk]);

      // v0 piece computation and predictions
      const ColVecd& tempV0 = solverO.solve(r.segment(cselsk, nk));
      if (m>mBurned)
      {
        SpMat L = solverO.matrixL();
        #pragma omp critical
        {
          wPreds.block(cselsk, (m-1)%verbose, nk, 1) = tempV0/tau20 + L.triangularView<Lower>().solve(zrnorm(nk, 1));
        }

      }

      // Covariance matrix
      std::vector<T> tripletListC1;
      tripletListC1.reserve(nk*neigh);
      
      int cntc1 = 0;
      tripletListC1.push_back(T(0, 0, sigma21));
      for (int i = 1; i<nk; i++){
        tripletListC1.push_back(T(i, i, sigma21));
        for (int j=lowjk[i-1]; j<i; j++){
          const double c1ij = sigma21 * exp(-distk(cntc1)*phi1);
          tripletListC1.push_back(T(i, j, c1ij));
          tripletListC1.push_back(T(j, i, c1ij));
          cntc1++;
        }
      }
      SpMat C1(nk, nk);
      C1.setFromTriplets(tripletListC1.begin(), tripletListC1.end());

      // AD Function
      std::vector<T> tripletListIA1t;
      tripletListIA1t.reserve(nk*neigh);
      std::vector<T> tripletListSDIA1;
      tripletListSDIA1.reserve(nk*neigh);
      
      ColVecd Dvec1(nk);
      Dvec1(0) = sigma21;
      double logD1 = log(sigma21);
      
      tripletListIA1t.push_back(T(0, 0, 1));
      tripletListSDIA1.push_back(T(0, 0, 1/Dvec1(0)));
      for(int i = 0; i < (nk-1); i++)
      {
        Matd C1block1 = C1.block(lowjk[i], lowjk[i], 
                                 i-lowjk[i]+1, i-lowjk[i]+1);
        Matd C1block2 = C1.block(lowjk[i], i+1, i-lowjk[i]+1, 1);
        ColVecd temp = C1block1.ldlt().solve(C1block2);
        
        Dvec1(i+1) = (sigma21 - (C1.block(i+1, lowjk[i], 
                                 1, i-lowjk[i]+1)*temp)(0));
        logD1 += log(Dvec1(i+1));
        
        for (int j = lowjk[i]; j < i+1; j++){
          tripletListIA1t.push_back(T(j, i+1, -temp(j-lowjk[i])));
          tripletListSDIA1.push_back(T(i+1, j, -temp(j-lowjk[i])/Dvec1(i+1)));
        }
        tripletListIA1t.push_back(T(i+1, i+1, 1));
        tripletListSDIA1.push_back(T(i+1, i+1, 1/Dvec1(i+1)));
      }
      SpMat IA1t(nk, nk);
      SpMat SolveDIA1(nk, nk);
      IA1t.setFromTriplets(tripletListIA1t.begin(), tripletListIA1t.end());
      SolveDIA1.setFromTriplets(tripletListSDIA1.begin(), tripletListSDIA1.end());
      
      // Computing Omega
      SpMat Omega1 = (IA1t.triangularView<UnitUpper>()*(SolveDIA1));

      for (int i = 0; i < nk; i++)
      {
        Omega1.coeffRef(i, i) += 1/tau21;
      }

      // Solving Omega
      CholLL solverO1(Omega1);
      SpMat L1 = solverO1.matrixL();
      
      // First piece of log-likelihood 
      logd1 += nk*log(tau21) + logD1 + 2*L1.diagonal().array().log().matrix().sum();

      // v and v0 update
      #pragma omp critical
      {
        Omega1s[indsk] = Omega1;
        v0.segment(cselsk,  nk) = tempV0;
        v1.segment(cselsk,  nk) = solverO1.solve(r.segment(cselsk, nk));
        vNew.block(cselsk, 0,  nk, nall) = solverO1.solve(XS.block(cselsk, 0, nk, nall));
        
      }
    }
    
    // Second piece of likelihood current
    const double&  q0 = r.dot(r)/tau20 - r.dot(v0)/(tau20*tau20);
    // Second piece of likelihood new
    const double&  q1 = r.dot(r)/tau21 - r.dot(v1)/(tau21*tau21);
    
    // Printing out predictions
    if ((m>mBurned))
    {
      loglik(m-mBurned-1) = -n/2*log(2*PI) - q0*0.5 - 0.5*logd0;
      
      if ((m%verbose==0)){
        string pname = fileName + std::to_string(m/verbose)+"_preds.txt";
        predsFile.open(pname.c_str());
        predsFile << wPreds.format(ChainFmt) << "\n";
        predsFile.close();
      }
    }
    
    // MH Step on theta
    const double& Num = - q1*0.5 - 0.5*logd1 -
      log(sigma20) - log(phi0) - log(tau20) +
      LogInvGamma(sigma21, alphaS, betaS) +
      LogGamma(phi1, alphaPhi, betaPhi) +
      LogInvGamma(tau21, alphaT, betaT);
    const double&  Den = -q0*0.5 - 0.5*logd0 -
      log(sigma21) - log(phi1) - log(tau21) +
      LogInvGamma(sigma20, alphaS, betaS) +
      LogGamma(phi0, alphaPhi, betaPhi) +
      LogInvGamma(tau20, alphaT, betaT);
    
    // Acceptance probs
    const double& AccProb = exp(Num-Den);
    
    if (mhUnifs(m-1) < AccProb)
    {
      // Storing accepted values
      logtheta0 = logtheta1;
      sigma20 = sigma21;
      phi0 = phi1;
      tau20 = tau21;
      
      acc(m) = 1;
      accCur += 1;
      
      tXSv = (tXS*vNew);
      vtZ = (vNew.transpose()*Z);
      
      Omegas = Omega1s;
      logd0 = logd1;
    }
    
    // Updating the chain
    lthetas.col(m) << logtheta0(0), logtheta0(1), logtheta0(2);
    
    // Gibbs sampler for Beta and eta
    Fmat.bottomRightCorner(nspl, nspl) = FmatE;
    fmat.tail(nspl) = fmatE;
    
    // Computing be and BE
    const ColVecd& be = (tXSZ)/tau20 - vtZ/(tau20*tau20) + fmat;
    const Matd& BE = tXSXS/tau20 - (tXSv)/(tau20*tau20) + Fmat;
    
    // Simulation
    LLT<Matd> cholBE(BE);
    Matd LBE = cholBE.matrixL();
    betaeta = cholBE.solve(be) + 
      LBE.triangularView<Lower>().solve(zrnorm(nall, 1));
    
    beta0 = betaeta.head(ncov);
    betas.col(m) = beta0;
    eta0 = betaeta.tail(nspl);
    etas.col(m) = eta0;
    
    // Gibbs sampler for lambda
    
    // Computing f.c. parameters (alphal already computed on top)
    const double& invbetal = 1/(betal0 + eta0.dot(KE*eta0)/2);
    
    // Simulation
    std::default_random_engine generator;
    std::gamma_distribution<double> distribution(alphal, invbetal);
    lambda20 = distribution(generator);
    lambda2s(m) = lambda20;
    
    FmatE = lambda20*KE;
    fmatE = lambda20*KEmuE;
    
  }
  // Monitoring
  monitor.open(mname.c_str(), std::ios_base::app);
  monitor << "Iteration: " << std::endl << M << std::endl;
  monitor << "Acceptance rate: " << std::endl << accCur/verbose << std::endl;
  monitor.close();
  
  // Last parameter chains
  string cname = fileName + std::to_string(M/verbose) + "_chains.txt";
  chainsFile.open(cname.c_str());
  MatrixXd curChains(lthetas.rows()+betas.rows()+etas.rows()+1, (M-1)%verbose+1);
  curChains << lthetas.rightCols((M-1)%verbose+1), betas.rightCols((M-1)%verbose+1), 
               etas.rightCols((M-1)%verbose+1), lambda2s.tail((M-1)%verbose+1).transpose();
  chainsFile << curChains.format(ChainFmt) << "\n";
  chainsFile.close();
  
  // Last predictions
  const ColVecd& mm = XS*betaeta;
  const ColVecd& r = Z - mm;
  ColVecd v0(n);
  CholLL solverO;
  string indsk;
  int nk;
  int cselsk;

  // Loop over individuals
  #pragma omp parallel for private(solverO, indsk, nk, cselsk) schedule(auto)
  for (int k=0; k<nInd; k++){
    indsk = inds[k];
    nk = nks[indsk];
    cselsk = csels[k];
    solverO.compute(Omegas[indsk]);

    // v0 piece computation and predictions
    ColVecd tempV0 = solverO.solve(r.segment(cselsk, nk));
    SpMat L = solverO.matrixL();
  
    #pragma omp critical
    {
      wPreds.block(cselsk, (M-1)%verbose, nk, 1) = tempV0/tau20 + L.triangularView<Lower>().solve(zrnorm(nk, 1));
      v0.segment(cselsk,  nk) = tempV0;
    }
  }
  const double&  q0 = r.dot(r)/tau20 - r.dot(v0)/(tau20*tau20);
  loglik(M-mBurned-1) = -n/2*log(2*PI) - q0*0.5 - 0.5*logd0;
  
  string pname = fileName+std::to_string(M/verbose)+"_preds.txt";
  predsFile.open(pname.c_str());
  predsFile << wPreds.format(ChainFmt) << "\n";
  predsFile.close();
  
  // Computing covariance parameters from logs
  ArrayXXd thetas = lthetas.array().exp();
  
  // Monitoring
  monitor.open(mname.c_str(), std::ios_base::app);
  monitor << "Finished" << std::endl;
  monitor.close();
  
  // Returning the output
  return Rcpp::List::create(thetas, betas, etas, lambda2s, acc, wPreds, loglik);
}



// [[Rcpp::export]]
Rcpp::List tNngpCollapsed_Preds(const Eigen::Map<ArrayXd> t, 
                                const Eigen::Map<ArrayXd> tPred, const MapMatd XSpred,
                                const int& neigh, const int& index,
                                const Eigen::Map<ArrayXXd> betas, const Eigen::Map<ArrayXXd> etas, const Eigen::Map<ArrayXXd> thetas, 
                                const MapMatd wPredsIn,
                                const int verbose, const int n_threads, const string& fileName) {
  ofstream monitor;
  string mname = fileName+"PredsMonitor.txt";
  
  omp_set_num_threads(n_threads);
  
  // Monitoring
  monitor.open(mname.c_str(), std::ios_base::app);
  monitor << "Individual OOS: " << std::endl << index << std::endl;
  monitor.close();
  
  // Data
  const int& nPred = tPred.size();
  const int& M = thetas.cols();
  const int& nall = XSpred.cols();
  
  // Auxiliary
  Array2Xi NNeigh(2, nPred);
  NNeigh = NewNeighbors(t, tPred, neigh);
  std::unordered_map<int, ColVecd> dPreds;
  std::unordered_map<int, ColVecd> dNeighs;
  
  // Computing distances for covariances
  for (int i=0; i<nPred; i++){
    int ll = NNeigh(0, i);
    int top = NNeigh(1, i);
    
    // Distances s0 to N0
    ColVecd tmps0(neigh);
    
    // Distances N0 to N0
    ColVecd tmpN0((neigh-1)*neigh/2);
    int counter = 0;
    
    // Computations
    // Looping over neighbor set (rows)
    for(int k = ll; k <= top; k++){
      // Vector of distances of neigbors from new obs
      tmps0(k-ll) = abs(tPred(i) - t(k));
      
      // Looping over neighbor set (cols)
      for (int j = k+1; j<=top; j++){
        // Vector of distances between neighbors
        tmpN0(counter) = abs(t(j)-t(k));
        counter++;
      }
    }
    
    // Storing
    dPreds[i] = tmps0;
    dNeighs[i] = tmpN0;
  }
  
  // Predictions initialization
  ArrayXXd wPredsOut(nPred, M);
  ArrayXXd yPredsOut(nPred, M);
  
  // Loop on parameter chains
  for(int m = 0; m < M; m++){
    if (m%verbose==0)
    {
      Rcpp::Rcout << "Iteration" << std::endl << m << std::endl;
    }
    
    // Current parameters
    double sigma2 = thetas(0, m);
    double phi = thetas(1, m);
    double tau2 = thetas(2, m);
    ColVecd betaetam(nall);
    betaetam << betas.col(m), etas.col(m);
    ColVecd wm = wPredsIn.col(m);
    
    // Loop on value to predict
    #pragma omp parallel for schedule(auto)
    for(int i = 0; i<nPred; i++){
      
      // Auxiliary
      int ll = NNeigh(0, i);
      int top = NNeigh(1, i);
      const ColVecd dists0i = dPreds[i];
      
      // InitializingC(s0, N0)
      ColVecd covs0i(neigh);
      // Initializing C(N0, N0)
      Matd CN0i(neigh, neigh);
      
      // Computations
      // First value on diagonal
      CN0i(0, 0) = sigma2;
      // Auxiliary counter
      int skippedEl = 0;
      // Looping over rows
      for(int k = 0; k < neigh; k++){
        // C(s0, N0)
        covs0i(k) = sigma2 * exp(-dists0i(k)*phi);
        
        // Values on diagonal
        CN0i(k, k) = sigma2;
        // Updating counter
        skippedEl += k+1;
        // Looping over columns
        for (int j = k+1; j < neigh; j++){
          CN0i(k, j) = sigma2 * exp(-dNeighs[i](j+k*neigh-skippedEl)*phi);
          CN0i(j, k) = CN0i(k, j);
        }
      }
      
      // Computing mean and variance for latent process
      double mi = covs0i.dot(CN0i.llt().solve(wm.segment(ll, top-ll+1)));
      double vi = sigma2 - covs0i.dot(CN0i.llt().solve(covs0i));
      
      // Latent gaussian process
      wPredsOut(i, m) = mi + sqrt(vi)*zigg.norm();
      
      // Process
      yPredsOut(i, m) = XSpred.row(i)*betaetam + wPredsOut(i, m) + sqrt(tau2)*zigg.norm();
      
    }
  }
  
  // Returning output
  return Rcpp::List::create(wPredsOut, yPredsOut);
  
}

