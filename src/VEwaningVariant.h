arma::mat preu(const arma::vec& E, const arma::vec& U, const arma::vec& R, 
               const arma::vec& lag, const arma::ivec& A, const arma::ivec& Gam, 
               const arma::ivec& Psi, const arma::ivec& delta, const int wgt, 
               const double time, const arma::mat& fR, const arma::vec& fEX,  
               const arma::vec& Psiprob, const arma::mat& KU, const double minWgt,  
               const double maxWgt, const arma::ivec& group);

arma::mat preb(const arma::vec& E, const arma::vec& U, const arma::vec& R, 
               const arma::vec& lag, const arma::ivec& A, const arma::ivec& delta,
               const int wgt, const double time, const arma::mat& KR, 
               const arma::vec& fEX,  const arma::mat& KU, const double minWgt,
               const double maxWgt, const arma::ivec& group);

arma::vec KU(const arma::vec& expXbetaG0, const arma::vec& expXbetaG1,
              const double Lambdat, const arma::vec& LambdaR, 
              const arma::vec& R, const double time, const int iMean);

arma::vec KR(const arma::vec& expXbeta, const double Lambdat, const int iMean);

Rcpp::List estimateb(const arma::vec& E, const arma::vec& U, 
               const arma::vec& R, const arma::vec& lag, 
               const arma::ivec& A, const arma::ivec& Gam, 
               const arma::ivec& Psi, const arma::ivec& delta,
               const int wgt, 
               const arma::vec& fEX, const arma::vec& Psiprob, 
               const double minWgt, const double maxWgt,
               const arma::mat& unBlind_expXbeta,
               const arma::mat& dropout_expXbetaG0, 
               const arma::mat& dropout_expXbetaG1, 
               const arma::mat& unBlind_LambdaT,
               const arma::mat& dropout_LambdaT,
               const arma::mat& dropout_LambdaR, 
               const arma::vec& times, const int gFunc, 
               const arma::ivec& group, const arma::vec& v, 
               const arma::vec& theta, const int type);

Rcpp::List estimateu(const arma::vec& E, const arma::vec& U, const arma::vec& R, const arma::vec& lag, 
               const arma::ivec& A, const arma::ivec& Gam, const arma::ivec& Psi, const arma::ivec& delta,
               const int wgt, const arma::mat& fR, const arma::vec& fEX, const arma::vec& Psiprob, 
               const double minWgt, const double maxWgt,
               const arma::mat& dropout_expXbetaG0, 
               const arma::mat& dropout_expXbetaG1, 
               const arma::mat& dropout_LambdaT,
               const arma::mat& dropout_LambdaR, const arma::vec& times, 
               const int gFunc, const arma::ivec& group, const arma::vec& v, 
               const arma::vec& theta, const int type);

