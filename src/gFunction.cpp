#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

List gPiecewise(const vec& u, const vec& theta, const vec& v) {

  int n = u.n_elem;
  int m = theta.n_elem;
  int nv = v.n_elem;

  vec gu, temp;
  mat gutheta;

  gu.zeros(n);
  gutheta.zeros(n,m);

  gutheta.col(0).ones();

  for (int i = 2; i < nv; ++i) {

    temp.zeros(n);
    temp.elem(find( (u > v(i-1)) && 
                    (u < (v(i)+1e-8)) )).ones();
    gutheta.col(i-1) = temp;
    gu = gu + theta(i-1) * gutheta.col(i-1);
  }

  List toReturn(2);
  toReturn(0) = gu;
  toReturn(1) = gutheta;

  return toReturn;

}   

List gLinear(const vec& u, const vec& theta) {

  int n = u.n_elem;
  int m = theta.n_elem;

  vec gu;
  mat gutheta(n,m);
  gu = theta(1)*u;

  gutheta.col(0).ones();
  gutheta.col(1) = u;

  List toReturn(2);
  toReturn(0) = gu;
  toReturn(1) = gutheta;

  return toReturn;

}   

List gSplineLinear(const vec& u, const vec& theta, const vec& knots) {

  int n = u.n_elem;
  int m = theta.n_elem;
  int nv = knots.n_elem;

  vec gu, temp;
  mat gutheta;

  gutheta.zeros(n,m);

  gutheta.col(0).ones();

  gu = theta(1)*u;
  gutheta.col(1) = u;

  for (int i = 0; i < nv; ++i) {
    temp = u - knots(i);
    temp.elem(find(temp < 0.0)).zeros();
    gu = gu + theta(i+2)*temp;
    gutheta.col(i+2) = temp;
  }

  List toReturn(2);
  toReturn(0) = gu;
  toReturn(1) = gutheta;

  return toReturn;

}

List gSplineCubic(const vec& u, const vec& theta, const vec& knots) {

  int n = u.n_elem;
  int m = theta.n_elem;
  int nv = knots.n_elem;

  vec gu, temp;
  mat gutheta;

  gu.zeros(n);
  gutheta.zeros(n,m);

  gutheta.col(0).ones();

  gutheta.col(1) = u;
  gu = gu + theta(1)*gutheta.col(1);

  gutheta.col(2) = pow(u, 2);
  gu = gu + theta(2)*gutheta.col(2);

  gutheta.col(3) = pow(u, 3);
  gu = gu + theta(3)*gutheta.col(3);

  for (int k = 0; k < nv; ++k) {
    temp = pow(u - knots(k),3);
    temp.elem(find(temp < 0.0)).zeros();
    gu = gu + theta(k+4)*temp;
    gutheta.col(k+4) = temp;
  }

  List toReturn(2);
  toReturn(0) = gu;
  toReturn(1) = gutheta;

  return toReturn;

}


// 1 - piecewise; 2 - linear; 3 - linear spline; 4 - cubic spline
// [[Rcpp::export]]
List gFunction(const int gFunc, const arma::vec& u, const arma::vec& theta, 
               const arma::vec& knots) {

  if (gFunc == 1) {
    return gPiecewise(u, theta, knots);
  } else if (gFunc == 2) {
    return gLinear(u, theta);
  } else if (gFunc == 3) {
    return gSplineLinear(u, theta, knots);
  } else {
    return gSplineCubic(u, theta, knots);
  }
}
