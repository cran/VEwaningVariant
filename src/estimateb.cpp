#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "VEwaningVariant.h"
#include "gFunction.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

List estimateb(const arma::vec& E, const arma::vec& U, 
               const arma::vec& R, const arma::vec& lag, 
               const arma::ivec& A, const arma::ivec& Gam, 
               const arma::ivec& Psi, const arma::ivec& delta,
               const int wgt, 
               const arma::vec& fEX, const arma::vec& Psiprob, 
               const double minWgt, const double maxWgt,
               const arma::mat& unBlind_expXbeta,
               const arma::mat& censor_expXbetaG0, 
               const arma::mat& censor_expXbetaG1, 
               const arma::mat& unBlind_LambdaT,
               const arma::mat& censor_LambdaT,
               const arma::mat& censor_LambdaR, 
               const arma::vec& times, const int gFunc, 
               const arma::ivec& group, const arma::vec& v, 
               const arma::vec& theta, const int type) {

  int n = E.n_elem;
  int m = theta.n_elem;
  int nTimes = times.n_elem;

  vec estb;
  mat inflmat, gradb;

  vec KR0, KR1;
  mat KU01, KR01, prebT;
  mat KU0, KU1;

  mat ZmZbar;
  mat meatb;
  double YbSum;
  rowvec Zbar;
  vec dNmYdLambdaHat;
  vec dNbt, Ybt, expgE, tEl;
  vec expgR, tRl;
  List gFuncE, gFuncR;

  estb.zeros(m);
  inflmat.zeros(n,m);
  gradb.zeros(m,m);

  for (int it = 0; it < nTimes; ++it) {

    // wgt == 2 when censoring and others
    // wgt == 1 when only censoring
    // wgt == 0 when no model dependent weighting
    if (wgt != 0) {
      KU0 = KU(censor_expXbetaG0.col(0), censor_expXbetaG1.col(0), 
               censor_LambdaT(it,0), censor_LambdaR.col(0), R, times(it), n);

      KU1 = KU(censor_expXbetaG0.col(1), censor_expXbetaG1.col(1), 
               censor_LambdaT(it,1), censor_LambdaR.col(1), R, times(it), n+1);

      KU01 = join_rows(KU0, KU1);
    } else {
      KU01.zeros(1,1);
    }

    if (wgt == 2) {
      KR0 = KR(unBlind_expXbeta.col(0), unBlind_LambdaT(it,0), n);

      KR1 = KR(unBlind_expXbeta.col(1), unBlind_LambdaT(it,1), n+1);

      KR01 = join_rows(KR0, KR1);
    } else {
      KR01.zeros(1,1);
    }

    prebT = preb(E, U, R, lag, A, delta,
                 wgt, times(it), KR01, fEX, KU01,
                 minWgt, maxWgt, group);

    //  exponential term for tilYb(t)
    // (t_m - E_i - lag)
    // {n}
    tEl = -((E + lag) - times(it));

    gFuncE = gFunction(gFunc, tEl, theta, v);
    vec gfe0 = gFuncE(0);

    expgE = exp(theta(0) + gfe0);

    Ybt = prebT.col(1) + prebT.col(2) % expgE;

    // {n}
    dNbt = prebT.col(0);

    //   Get Zb, ZbminusZbar, and blinded part of estimating equation

    // {1}
    YbSum = sum(Ybt);

    // {n} dN - dLambda Y
    dNmYdLambdaHat = dNbt - Ybt * (sum(dNbt) / YbSum);

    mat gfe = gFuncE(1);
    if (type == 1) gfe.col(0).zeros();

    uvec indA0 = find(A == 0);

    gfe.rows(indA0).zeros();

    // {1 x d}
    Zbar = (Ybt.t() * gfe) / YbSum;

    // {n x d}
    ZmZbar = gfe;
    ZmZbar.each_row() -= Zbar;

    estb += (dNbt.t() * ZmZbar).t();

    inflmat += ZmZbar.each_col() % dNmYdLambdaHat;

    ZmZbar.each_col() %= sqrt(Ybt) / sqrt(YbSum);

    mat tmat = ZmZbar.t() * ZmZbar;
    gradb += tmat * sum(dNbt);

  }

  meatb = inflmat.t() * inflmat;

  List toReturn(3);
  toReturn(0) = estb;
  toReturn(1) = gradb;
  toReturn(2) = meatb;

  return toReturn;

}

