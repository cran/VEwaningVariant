#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

mat preu(const vec& E, 
         const vec& U, 
         const vec& R, 
         const vec& lag, 
         const ivec& A, 
         const ivec& Gam, 
         const ivec& Psi, 
         const ivec& delta,
         const int wgt,
         const double time, 
         const mat& fR, 
         const vec& fEX,  
         const vec& Psiprob, 
         const mat& KU,  
         const double minWgt,
         const double maxWgt,
         const ivec& group) {

  int n = E.n_elem;

  vec Y, dNt, wt0, wt1, tmp, irt;

  mat w01;

  // Observed at-risk indicator at time t
  // Y = I(E_i < t <= U_i)
  Y.zeros(n);
  Y(find((E < time) && (time < (U+1e-8)))).ones();

  irt.zeros(n);
  irt.elem(find(Gam == 1)).ones();

  if (wgt == 0) {
    // no model based weighting
    w01.ones(n,2);
    w01.each_col() %= irt;
  } else if(wgt == 1) {
    // weighting only depends on censoring
    w01 = KU;
    w01.each_col() %= irt;
    mat w01m = clamp(w01, minWgt, maxWgt);
    w01 = w01m;
  } else {
    // weighting depends on entry, psi, unblinding, and censoring models
    tmp = fEX % Psiprob % irt;
    w01 = fR % KU;
    w01.each_col() %= tmp;
    mat w01m = clamp(w01, minWgt, maxWgt);
    w01 = w01m;

  }

  // (1-A_i) I(t - R >= lag) {I(Gamma = 1, Psi = 1) w0/h01}
  wt0 = w01.col(0);
  wt0.elem(find((time-R) < lag || A == 1 || Psi == 0)).zeros();

  // A_i I(t > R) {{I(Gamma = 1) w1/h11}
  wt1 = w01.col(1);
  wt1.elem(find(time < (R+1e-8) || A == 0)).zeros();
    
  //  dtilNu(t) and Ytilu(t)
    
  // Indicator that a participant is observed to be infected at time t
  // dNt = I(U_i == t, Delta == 1)
  // {n}
  uvec ind = find( (U > (time-1e-8)) && 
                   (U < (time+1e-8)) && 
                   (group == 1) );
  dNt.zeros(n);
  dNt.elem(ind) = wt0.elem(ind) + wt1.elem(ind);

  wt0 = wt0 % Y;
  wt1 = wt1 % Y;

  return join_rows(dNt, wt0,  wt1);
}
