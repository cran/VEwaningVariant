#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "VEwaningVariant.h"
#include "gFunction.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

List estimateu(const arma::vec& E, const arma::vec& U, 
               const arma::vec& R, const arma::vec& lag, 
               const arma::ivec& A, const arma::ivec& Gam, 
               const arma::ivec& Psi, const arma::ivec& delta,
               const int wgt, const arma::mat& fR, 
               const arma::vec& fEX, const arma::vec& Psiprob, 
               const double minWgt, const double maxWgt,
               const arma::mat& censor_expXbetaG0, 
               const arma::mat& censor_expXbetaG1, 
               const arma::mat& censor_LambdaT,
               const arma::mat& censor_LambdaR, const arma::vec& times, 
               const int gFunc, const arma::ivec& group, 
               const arma::vec& v, const arma::vec& theta,
               const int type) {

  int n = E.n_elem;
  int nTimes = times.n_elem;
  int m = theta.n_elem;

  double YuSum;

  vec dNmYdLambdaHat, dNut, estu, KU0, KU1;
  vec tEl, tRl, Yut;

  rowvec Zbar;

  mat gradu, inflmat, KU01, meatu, preuT;
  mat Z, ZmZbar;

  List gFuncE, gFuncR;

  estu.zeros(m);
  inflmat.zeros(n,m);
  gradu.zeros(m,m);

  uvec indA0 = find(A == 0);
  uvec indA1 = find(A == 1);

  for (int it = 0; it < nTimes; ++it) {

    // wgt == 2 when censoring and others
    // wgt == 1 when only censoring
    // wgt == 0 when no model dependent weighting
    if (wgt != 0) {
      KU0 = KU(censor_expXbetaG0.col(0), censor_expXbetaG1.col(0), censor_LambdaT(it,0), 
               censor_LambdaR.col(0), R, times(it), n);

      KU1 = KU(censor_expXbetaG0.col(1), censor_expXbetaG1.col(1), censor_LambdaT(it,1), 
               censor_LambdaR.col(1), R, times(it), n+1);

      KU01 = join_rows(KU0, KU1);
    } else {
      KU01.zeros(1,1);
    }

    // preu returns an {n x 3} matrix
    // first column is dN
    // second column is w0
    // third column is w1
    preuT = preu(E, U, R, lag, A, Gam, Psi, delta,
                 wgt, times(it), fR, fEX, Psiprob, KU01,  
                 minWgt, maxWgt, group);

    // {t_m - R_i - lag}
    // {n}
    tRl = -(R + lag - times(it));

    // {t_m - E_i - lag}
    // {n}
    tEl = -(E + lag - times(it));

    // gFunction() returns a list of two elements
    // first element is g(t) for each individual {n}
    // second element is g'(t) for each individual {n x nTheta}
    gFuncR = gFunction(gFunc, tRl, theta, v);
    vec gvR = gFuncR(0);

    gFuncE = gFunction(gFunc, tEl, theta, v);
    vec gvE = gFuncE(0);

    // w0 exp(g(t-R-lag)) + w1 exp(g(t-E-lag))
    // {n}
    Yut = preuT.col(1) % exp(gvR) + preuT.col(2) % exp(gvE);

    // {n}

    dNut = preuT.col(0);

    YuSum = sum(Yut);
 
    ////   Get Zu, ZuminusZbar, and unblinded part of estimating equation

    // {n} dN - dLambda Y
    dNmYdLambdaHat = dNut - Yut * (sum(dNut) / YuSum);

    ZmZbar.zeros(n,m);

    mat ge = gFuncE(1);
    mat gr = gFuncR(1);
    if (type == 1) {
      ge.col(0).zeros();
      gr.col(0).zeros();
    }

    ge.rows(indA0).zeros();
    gr.rows(indA1).zeros();

    // {n x d}
    Z = ge + gr;

    // {1 x d}
    Zbar = (Yut.t() * Z) / YuSum;

    // {n x d}
    ZmZbar = Z;
    ZmZbar.each_row() -= Zbar;

    estu += (dNut.t() * ZmZbar).t();

    inflmat += ZmZbar.each_col() % dNmYdLambdaHat;

    ZmZbar.each_col() %= sqrt(Yut) / sqrt(YuSum);

    mat tmat = ZmZbar.t() * ZmZbar;
    gradu += tmat * sum(dNut);

  }

  meatu = inflmat.t() * inflmat;

  List toReturn(3);
  toReturn(0) = estu;
  toReturn(1) = gradu;
  toReturn(2) = meatu;

  return toReturn;
}
