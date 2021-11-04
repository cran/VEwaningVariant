#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// expXbetaG0: vector of exp(x beta) from dropout fit with G = 0 {n+2}
//             last two elements of the vector are the means for A = 0 and A = 1
// expXbetaG1: vector of exp(x beta) from dropout fit with G = 1 {n+2}
//             last two elements of the vector are the means for A = 0 and A = 1
// Lambdat: Lambda at current t
// LambdaR: Lambda evaluated at each R {n}
// R: R {n}
// time: current time point
// iMean: the index of the mean value to use for stabilization

vec KU(const vec& expXbetaG0, 
       const vec& expXbetaG1,
       const double Lambdat, 
       const vec& LambdaR, 
       const vec& R, 
       const double time,
       const int iMean) {

  int n = expXbetaG0.n_elem;

  vec s3, lep;
  vec KU, KU_mean;
  double tmp;

  s3.zeros(n);
  KU.zeros(n-2);
  lep = LambdaR % expXbetaG0 + (Lambdat - LambdaR) % expXbetaG1;

  for (int i = 0; i < n-2; ++i) {
    if (time < R(i)) {
      KU(i) = exp(-(Lambdat * expXbetaG0(i)));
    } else {
      KU(i) = exp(-lep(i));
      tmp = LambdaR(i) * expXbetaG0(iMean) + 
            (Lambdat - LambdaR(i)) * expXbetaG1(iMean);
      s3(i) = exp(-tmp);
    }

  }

  KU_mean = exp(-(Lambdat * expXbetaG0(iMean))) + s3.head(n-2);

  KU = KU_mean / KU;

  return KU;

}
