// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// estTheta
List estTheta(const arma::vec& E, const arma::vec& U, const arma::vec& R, const arma::vec& lag, const arma::ivec& A, const arma::ivec& Gam, const arma::ivec& Psi, const arma::ivec& delta, const int wgt, const arma::mat& fR, const arma::vec& fEX, const arma::vec& Psiprob, const double minWgt, const double maxWgt, const arma::mat& censor_expXbetaG0, const arma::mat& censor_expXbetaG1, const arma::mat& censor_LambdaTB, const arma::mat& censor_LambdaTU, const arma::mat& censor_LambdaR, const arma::mat& unBlind_expXbeta, const arma::mat& unBlind_LambdaT, const arma::vec& timesB, const arma::vec& timesU, const int gFunc, const arma::ivec& group, const arma::vec& v, const arma::vec& theta, const int type);
RcppExport SEXP _VEwaningVariant_estTheta(SEXP ESEXP, SEXP USEXP, SEXP RSEXP, SEXP lagSEXP, SEXP ASEXP, SEXP GamSEXP, SEXP PsiSEXP, SEXP deltaSEXP, SEXP wgtSEXP, SEXP fRSEXP, SEXP fEXSEXP, SEXP PsiprobSEXP, SEXP minWgtSEXP, SEXP maxWgtSEXP, SEXP censor_expXbetaG0SEXP, SEXP censor_expXbetaG1SEXP, SEXP censor_LambdaTBSEXP, SEXP censor_LambdaTUSEXP, SEXP censor_LambdaRSEXP, SEXP unBlind_expXbetaSEXP, SEXP unBlind_LambdaTSEXP, SEXP timesBSEXP, SEXP timesUSEXP, SEXP gFuncSEXP, SEXP groupSEXP, SEXP vSEXP, SEXP thetaSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type E(ESEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type Gam(GamSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const int >::type wgt(wgtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fR(fRSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type fEX(fEXSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Psiprob(PsiprobSEXP);
    Rcpp::traits::input_parameter< const double >::type minWgt(minWgtSEXP);
    Rcpp::traits::input_parameter< const double >::type maxWgt(maxWgtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type censor_expXbetaG0(censor_expXbetaG0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type censor_expXbetaG1(censor_expXbetaG1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type censor_LambdaTB(censor_LambdaTBSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type censor_LambdaTU(censor_LambdaTUSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type censor_LambdaR(censor_LambdaRSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type unBlind_expXbeta(unBlind_expXbetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type unBlind_LambdaT(unBlind_LambdaTSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type timesB(timesBSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type timesU(timesUSEXP);
    Rcpp::traits::input_parameter< const int >::type gFunc(gFuncSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(estTheta(E, U, R, lag, A, Gam, Psi, delta, wgt, fR, fEX, Psiprob, minWgt, maxWgt, censor_expXbetaG0, censor_expXbetaG1, censor_LambdaTB, censor_LambdaTU, censor_LambdaR, unBlind_expXbeta, unBlind_LambdaT, timesB, timesU, gFunc, group, v, theta, type));
    return rcpp_result_gen;
END_RCPP
}
// gFunction
List gFunction(const int gFunc, const arma::vec& u, const arma::vec& theta, const arma::vec& knots);
RcppExport SEXP _VEwaningVariant_gFunction(SEXP gFuncSEXP, SEXP uSEXP, SEXP thetaSEXP, SEXP knotsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type gFunc(gFuncSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type knots(knotsSEXP);
    rcpp_result_gen = Rcpp::wrap(gFunction(gFunc, u, theta, knots));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VEwaningVariant_estTheta", (DL_FUNC) &_VEwaningVariant_estTheta, 28},
    {"_VEwaningVariant_gFunction", (DL_FUNC) &_VEwaningVariant_gFunction, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_VEwaningVariant(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}