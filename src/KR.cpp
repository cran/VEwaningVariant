#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// expXbeta: vector of exp(x beta) from unblinded fit {n+2}
//           last two elements of the vector are the means for A = 0 and A = 1
// Lambdat: Lambda at current t
// iMean: integer indicating which of the last two elements to use as the
//  average value
vec KR(const vec& expXbeta, const double LambdaT, const int iMean) {

  int n = expXbeta.n_elem;
  vec krVec, res(n-2);

  // survival probability at each time point if each participant randomized
  // to placebo
  // S(t,X_i)
  // {n}
  krVec = exp(- expXbeta * LambdaT);

  //  return stabilized weights

  res = krVec(iMean) / krVec.head(n-2);

  return res;

}
