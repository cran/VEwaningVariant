#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "VEwaningVariant.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

List estTheta(const arma::vec& E, const arma::vec& U, const arma::vec& R, 
              const arma::vec& lag, const arma::ivec& A, 
              const arma::ivec& Gam, const arma::ivec& Psi, 
              const arma::ivec& delta, const int wgt, 
              const arma::mat& fR, const arma::vec& fEX, 
              const arma::vec& Psiprob, 
              const double minWgt, const double maxWgt,
              const arma::mat& censor_expXbetaG0, 
              const arma::mat& censor_expXbetaG1, 
              const arma::mat& censor_LambdaTB,
              const arma::mat& censor_LambdaTU,
              const arma::mat& censor_LambdaR, 
              const arma::mat& unBlind_expXbeta,
              const arma::mat& unBlind_LambdaT,
              const arma::vec& timesB,
              const arma::vec& timesU, 
              const int gFunc, const arma::ivec& group, 
              const arma::vec& v, const arma::vec& theta,
              const int type) {

  List eb(3), eu(3);
  NumericVector scoreb, scoreu;
  NumericVector gradb, gradu;
  NumericVector meatb, meatu;

  if (type == 2 || type == 3) {
    // type =2 blinded only, =3 both phases
    eb = estimateb(E, U, R, lag, A, Gam, Psi, delta,
                   wgt, fEX, Psiprob, minWgt, maxWgt,
                   unBlind_expXbeta,
                   censor_expXbetaG0,
                   censor_expXbetaG1,
                   unBlind_LambdaT,
                   censor_LambdaTB,
                   censor_LambdaR,
                   timesB,  
                   gFunc,
                   group,
                   v,
                   theta,
                   type);
  } else {
    vec sb(theta.n_elem);
    mat gb(theta.n_elem, theta.n_elem);
    mat mb(theta.n_elem, theta.n_elem);
    eb(0) = sb.zeros();
    eb(1) = gb.zeros();
    eb(2) = mb.zeros();
  }

  if (type == 1 || type == 3) {
    // type =1 unblinded only, =3 both phases
    eu = estimateu(E, U, R, lag, A, Gam, Psi, delta,
                   wgt, fR, fEX, Psiprob, minWgt, maxWgt,
                   censor_expXbetaG0,
                   censor_expXbetaG1,
                   censor_LambdaTU,
                   censor_LambdaR,
                   timesU,  
                   gFunc,
                   group,
                   v,
                   theta,
                   type);
  } else {
    vec su(theta.n_elem);
    mat gu(theta.n_elem, theta.n_elem);
    mat mu(theta.n_elem, theta.n_elem);
    eu(0) = su.zeros();
    eu(1) = gu.zeros();
    eu(2) = mu.zeros();
  }

  List toReturn(3);

  scoreb = eb(0);
  scoreu = eu(0);
  toReturn(0) = scoreb + scoreu;

  gradb = eb(1);
  gradu = eu(1);
  toReturn(1) = gradb + gradu;

  meatb = eb(2);
  meatu = eu(2);
  toReturn(2) = meatb + meatu;

  return toReturn;
}


