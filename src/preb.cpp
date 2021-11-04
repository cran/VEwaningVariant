#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

mat preb(const vec& E, 
         const vec& U, 
         const vec& R, 
         const vec& lag, 
         const ivec& A, 
         const ivec& delta,
         const int wgt,
         const double time, 
         const mat& KR, 
         const vec& fEX,  
         const mat& KU,  
         const double minWgt,
         const double maxWgt,
         const ivec& group) {

  int n = E.n_elem;

  vec Y, dNt, wt0, wt1, irt; 
  mat w01;

  // I(E_i < t <= U_i) i = 1:n
  // {n}
  Y.zeros(n);
  Y.elem(find( (E < time) && (time < (U+1e-8)) )).ones();

  //  Construct numerator indicators for dtilNb and tilYb

  // {n} I(R_i >= t)
  irt.zeros(n);
  irt.elem(find(R > (time-1e-8))).ones();

  //  stabilized weights

  if (wgt == 0) {
    // no model based weighting
    w01.ones(n,2);
    w01.each_col() %= irt;
  } else if (wgt == 1) {
    // weighting only depends on censoring
    w01 = KU;
    w01.each_col() %= irt;
    mat w01m = clamp(w01, minWgt, maxWgt);
    w01 = w01m;
  } else {
    // weighting depends on entry, psi, unblinding, and censoring models
    vec tmp = (fEX % irt);
    w01 = KR % KU;
    w01.each_col() %= tmp;
    mat w01m = clamp(w01, minWgt, maxWgt);
    w01 = w01m;
  }

  // (1-A_i)I(R_i >= t) fEX KRfR0
  wt0 = w01.col(0);
  wt0.elem(find(A == 1)).zeros();

  // A_i I(E_i + lag <= t) I(R_i >= t) fEx KRfR1
  wt1 = w01.col(1);
  wt1.elem(find((A == 0) || ((E + lag) > time))).zeros();

  // I(U_i == t; Delta = nu)
  // {n}
  uvec ind = find( (U > (time-1e-8)) && 
                   (U < (time+1e-8)) && 
                   (group == 1) );

  dNt.zeros(n);
  dNt.elem(ind) = wt0.elem(ind) + wt1.elem(ind);

  wt0 = wt0 % Y;
  wt1 = wt1 % Y;

  return join_rows(dNt, wt0, wt1);
}
