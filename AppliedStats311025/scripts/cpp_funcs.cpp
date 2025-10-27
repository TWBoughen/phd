#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;


NumericVector subsetNumVec(NumericVector x, IntegerVector index) {  
  // Length of the index vector
  int n = index.size();
  // Initialize output vector
  NumericVector out(n);
  
  // Subtract 1 from index as C++ starts to count at 0
  // Loop through index vector and extract values of x at the given positions
  for (int i = 0; i < n; i++) {
    out[i] = x[index[i]];
  }
  
  // Return output
  return out;
}


double rho(double lambda, NumericVector pars){
  double a = pars[0];
  double eps = pars[1];
  double k0 = pars[2];
  double b = pars[3];
  
  IntegerVector k0_vec_int = seq_len(k0);
  NumericVector k0_vec = as<NumericVector>(k0_vec_int);
  NumericVector logfracs = log(pow(k0_vec-1,a) + eps) -
    log(lambda + pow(k0_vec-1, a)+eps);
  NumericVector plogsums = cumsum(logfracs);
  
  double out = sum(exp(plogsums)) + (pow(k0,a)+eps)*exp(plogsums(k0-1))/(lambda-b);
  return(out);
}


double rho_adj(double lambda, NumericVector pars){
  return 1.0-rho(lambda, pars);
}

double find_lambda(NumericVector pars){
  
  double upper = 20;
  double tol = 1e-7;
  double lower = pars[3]+tol;
  int max_iter=100;
  double mid;
  
  for(int i;i<=max_iter;i++){
    mid= (lower+upper)/2.0;
    double f_mid = rho_adj(mid, pars);
    if(std::abs(f_mid)<tol){
      break;
    }else if(f_mid>0){
      upper=mid;
    }else{
      lower=mid;
    }
  }
  
  return mid;
}






// [[Rcpp::export]]
Rcpp::NumericMatrix rmvnorm_rcpp(int n, Rcpp::NumericVector mean, Rcpp::NumericMatrix sigma) {
  int p = mean.size();
  Rcpp::NumericMatrix result(n, p);
  
  // Convert Rcpp::NumericVector and Rcpp::NumericMatrix to Armadillo types
  arma::vec mean_arma = Rcpp::as<arma::vec>(mean);
  arma::mat sigma_arma = Rcpp::as<arma::mat>(sigma);
  
  // Cholesky decomposition of the covariance matrix
  arma::mat L = chol(sigma_arma, "lower");
  
  // Generate random samples from standard normal distribution
  arma::mat z = arma::randn(n, p);
  
  // Transform the standard normal samples to multivariate normal
  arma::mat result_arma = z * L.t() + repmat(mean_arma.t(), n, 1);
  
  // Convert result back to Rcpp::NumericMatrix
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      result(i, j) = result_arma(i, j);
    }
  }
  
  return result;
}






// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


bool check_positive(NumericVector pars,bool strict=true){
  if(strict){
    return std::all_of(pars.begin(),pars.end(),[](double x) {return x>0;});
  }else{
    return std::all_of(pars.begin(),pars.end(),[](double x) {return x>=0;});
  }
  
}


NumericVector llh_x(NumericVector x, NumericVector pars, double lambda, double minx=0){
  double a = pars[0];
  double eps = pars[1];
  double k0 = pars[2];
  double b = pars[3];
  IntegerVector k0vec_int = seq_len(k0-minx);
  NumericVector k0vec = as<NumericVector>(k0vec_int) +minx-1;
  NumericVector result(x.size());
  for(int i=0;i<x.size();i++){
    if(x[i]==0){
      result[i] = log(lambda) - log(lambda+eps);
    }else if(x[i]<k0){
      IntegerVector xvec_int = seq(minx,x[i]-1);
      NumericVector xvec = as<NumericVector>(xvec_int);
      result[i] = log(lambda) - log(lambda + eps + pow(x[i],a)) +  
        sum(log(1 - lambda/(lambda + eps + pow(xvec,a))));
    }else if(x[i]>=k0){
      result[i] = sum(log(1 - lambda/(lambda + eps + pow(k0vec,a))))+
        Rf_lbeta(x[i]-k0 + (pow(k0,a) + eps)/b,1+lambda/b) - Rf_lbeta((pow(k0,a) + eps)/b,lambda/b);
    }
  }
  return result;
}


double llh(NumericMatrix dat, NumericVector pars, double lambda, double minx=0){
  if(!check_positive(pars)){
    return -INFINITY;
  }
  return(sum(dat(_,1) * llh_x(dat(_,0),pars,lambda,minx)));
}


double prior(NumericVector pars){
  return dgamma(pars,1.0,1.0,true)[0] +
    dgamma(pars, 1.0,1.0,true)[1]+
    // dunif(pars, 4.99,10000.0,true)[2]+
    dgamma(pars, 1.0,1.0,true)[3];
}

double posterior(NumericMatrix dat, NumericVector pars, double lambda,double minx=0, bool auto_cond = false){
  if(auto_cond){
    minx = min(dat(_,0))-1;
  }
  return(llh(dat, pars, lambda, minx) + prior(pars));
}
//[[Rcpp::export]]
NumericVector posterior_k0(NumericMatrix dat, NumericVector k0, NumericVector pars, double lambda,double minx=0, bool auto_cond=false){
  NumericVector result(k0.size());
  for(int i=0;i<k0.size();i++){
    if(k0[i]<=min(dat(_,0)) | k0[i]<=5.0){
      result[i] = -INFINITY;
    }else if(k0[i]>=max(dat(_,0))){
      result[i] = -INFINITY;
    }else{
      pars[2] = k0[i];
      double lambda2 = find_lambda(pars);
      result[i] = posterior(dat, pars, lambda, minx, auto_cond);
    }
  }
  return result;
}



//[[Rcpp::export]]
List polylin_mcmc(double N, NumericMatrix dat, NumericVector init_pars,
                  NumericMatrix sig, double k0_jump=2, double burn_in=5e3,
                  double lookback = 200, double cov_scale=0.25,double minx=0, bool auto_cond=false,bool gibbs=true){
  NumericMatrix par_mat(N+1, 4);
  NumericVector lambda_vec(N+1);
  NumericVector post_vec(N+1);
  NumericVector old_pars = init_pars;
  NumericVector accepted(N+1);
  accepted(0) = 0.0;
  par_mat(0,_) = init_pars;
  lambda_vec(0) = find_lambda(init_pars);
  post_vec(0) = posterior(dat, init_pars, lambda_vec(0), minx, auto_cond);
  double st = 1;
  double old_post = post_vec(0);
  double old_lambda = lambda_vec(0);
  double acc_rate=0;
  double min_cov_idv = 0.00000001;
  NumericVector min_cov_vec = {min_cov_idv,0,0,0,
                               0,min_cov_idv,0,0,
                               0,0,min_cov_idv,0,
                               0,0,0,min_cov_idv};
  NumericMatrix min_cov(4,4,min_cov_vec.begin());
  
  NumericVector k0_probs;
  for(int i=1;i<=N;i++){
    if(i % static_cast<int>(lookback) ==0){
      Rcout<< i <<" : (" << old_pars[0]<<", "<<
        old_pars[1]<<", "<<
          old_pars[2]<<", "<<
            old_pars[3]<< ") - "<<100*acc_rate<<
              '\n';
    }
    
    accepted[i] = 0.0;
    //propose new values
    if(i>=burn_in & i % static_cast<int>(lookback)){
      IntegerVector prev = i-seq_len(lookback);
      NumericMatrix prev_mat(prev.size(), par_mat.ncol());
      double tot_accepted = 0.0;
      for(int j=0;j<prev.size();j++){
        prev_mat(j,_) = par_mat(prev[j],_);
        tot_accepted = tot_accepted + accepted[prev[j]];
      }
      acc_rate = tot_accepted/lookback;
      st = st * exp(-cov_scale*(0.234-acc_rate));
      arma::uvec cols = {0,1,2,3};
      NumericMatrix cov_mat = Rcpp::wrap(st*arma::cov(Rcpp::as<arma::mat>(prev_mat).cols(cols))+
        Rcpp::as<arma::mat>(min_cov));
      sig = cov_mat;
    }
    
    double cur_k0 = par_mat(i-1,2);
    double cur_k0_idx;
    IntegerVector poss_k0_int = seq_len(2*k0_jump +1);
    // NumericVector poss_k0 = as<NumericVector>(poss_k0_int) - k0_jump-1 + cur_k0;
    
    
    for(double idx=0;idx<dat(_,0).size();idx++){
      // Rcout<<dat(idx,0)<<std::endl;
      if(dat(idx,0)==cur_k0){
        cur_k0_idx = idx;
      }
    }
    
    IntegerVector poss_k0_idx = (poss_k0_int) - k0_jump-1 + cur_k0_idx;
    NumericVector poss_k0(poss_k0_idx.size());
    for(int k=0;k<poss_k0_idx.size();k++){
      if(poss_k0_idx[k]<=2){
        poss_k0[k] = -1;
      }else if(poss_k0_idx[k]>=dat(_,0).size()){
        poss_k0[k] = -2;
      }else{
        poss_k0[k] = dat(poss_k0_idx[k],0);
      }
    }
    NumericVector y = posterior_k0(dat, poss_k0,par_mat(i-1,_), old_lambda, minx,auto_cond);
    k0_probs = exp(y-max(y));
    
    
    if(sum(k0_probs>0)){
      
      double new_k0 = sample(poss_k0, 1, false, k0_probs)[0];
      
      NumericVector new_pars0 = rmvnorm_rcpp(1, {old_pars[0], old_pars[1],new_k0, old_pars[3]},sig);
  
      NumericVector new_pars = NumericVector::create(new_pars0(0), new_pars0(1),
                                                     new_k0, new_pars0(3));
      double new_lambda = find_lambda(new_pars);
      double new_post = posterior(dat, new_pars, new_lambda, minx, auto_cond);
      
      
      
      double logA = new_post-old_post;
      
      
      double U = log(runif(1))[0];
      if(U<logA){//accepting
        old_post = new_post;
        old_lambda = new_lambda;
        old_pars = new_pars;
        accepted[i] = 1;
      }
    }
    post_vec[i] = old_post;
    lambda_vec[i] = old_lambda;
    par_mat(i,_) = old_pars;
    
  }
  return List::create(
    Named("pars") = par_mat,
    Named("lambdas") = lambda_vec,
    Named("post") = post_vec,
    Named("sig") = sig
  );
}

// igp code


const double lnan(const double l) {
  return (l != l) ? -INFINITY : l;
}

const double lr1() {
  return log(runif(1)[0]);
}

template <class T>
const NumericVector tv(T x) {
  return NumericVector::create(x);
}

const IntegerVector ti(const int x) {
  return IntegerVector::create(x);
}

const LogicalVector tl(const bool x) {
  return LogicalVector::create(x);
}

const double ldnorm(const double x, const double mean, const double sd) {
  return dnorm(tv(x), mean, sd, true)[0];
}

const double ldgamma(const double x, const double shape, const double rate) {
  return dgamma(tv(x), shape, 1.0 / rate, true)[0];
}

const NumericVector lpgamma(const NumericVector x, const double shape, const double rate) {
  return pgamma(x, shape, 1.0 / rate, true, true);
}

const double ldbeta(const double x, const double a, const double b) {
  return dbeta(tv(x), a, b, true)[0];
}

const double ldunif(const double x, const double a, const double b) {
  return dunif(tv(x), a, b, true)[0];
}

template <class T>
void update(T & par_curr,
            T par_prop,
            double & lpost_curr,
            const double lpost_prop,
            double & llik_curr,
            const double llik_prop,
            double & s,
            const int i,
            const int burnin,
            const double invt = 1.0,
            const double factor = 10.0) {
  // Metropolis-Hastings step
  // invt is for Metropolis-coupled MCMC / parallel tempering, not power posterior
  const double lalpha = invt * (lpost_prop - lpost_curr);
  const bool accept_reject = lalpha > lr1();
  par_curr = accept_reject ? par_prop : par_curr;
  lpost_curr = accept_reject ? lpost_prop : lpost_curr;
  llik_curr = accept_reject ? llik_prop : llik_curr;
  if (i < burnin) {
    s = sqrt(s * s + (accept_reject ? 3.0 : (-1.0)) * s * s / factor / sqrt(i + 1.0));
  }
}

const double sd_curr(const vec x, const int i) {
  return sqrt(as_scalar(cov(x.head(i))));
}

const double cor_curr(const vec x, const vec y, const int i) {
  return as_scalar(cor(x.head(i), y.head(i)));
}

const int sample_1(const IntegerVector seq) {
  // sample 1 value from seq uniformly
  return sample(seq, 1, false)[0];
  // true or false shouldn't matter
}

const int sample_w(const IntegerVector seq, const NumericVector weights) {
  // sample 1 value from seq with weights
  return sample(seq, 1, true, weights)[0];
}

const double lse(const NumericVector x) {
  return log(sum(exp(x)));
}

DataFrame df_scalars(const int iter,
                     const int thin,
                     const int burn,
                     const int freq,
                     const bool mc3_or_marg) {
  return
  DataFrame::create(Named("iter") = ti(iter),
                    Named("thin") = ti(thin),
                    Named("burn") = ti(burn),
                    Named("freq") = ti(freq),
                    Named("mc3_or_marg") = tl(mc3_or_marg));
}

const IntegerVector for_Smix(const int u, const int x_min) {
  const uvec x0 = regspace<uvec>(1, 9);
  uvec x(x0);
  x = join_cols(x, x0 * pow(10, 1));
  x = join_cols(x, x0 * pow(10, 2));
  x = join_cols(x, x0 * pow(10, 3));
  x = join_cols(x, x0 * pow(10, 4));
  x = join_cols(x, x0 * pow(10, 5));
  IntegerVector x1 = wrap(x);
  x1 = x1[(x1 <= u) & (x1 >= x_min)];
  x1.insert(x1.end(), u);
  return x1;
}

const bool ispm1(const double x, const double precision = 1.0e-10) {
  // TRUE if x=+/-1.0 up to floating-point prec
  return abs(abs(x) - 1.0) < precision;
}

const double sqrt1mx2(const double x) {
  // sqrt(1 - x^2)
  return ispm1(x) ? 0.0 : sqrt(1.0 - x * x);
}

const double intdiv(const int a, const int b) {
  return (double) a / (double) b;
}

const double odds(const double p) {
  return p / (1.0 - p);
}

const NumericVector pm(const NumericVector v, const bool plus_or_minus) {
  const int n = v.size();
  const double factor = (double) plus_or_minus * 2.0 - 1.0;
  const NumericVector u = head(v, n - 1) + factor * tail(v, n - 1);
  return u;
}




const double llik_igpd(const NumericVector par,
                       const IntegerVector x,
                       const IntegerVector count,
                       const int u,
                       const double phiu) {
  if (x.size() != count.size()) {
    stop("llik_igp: lengths of x & count have to be equal.");
  }
  if (is_true(any(x <= 0))) {
    stop("llik_igpd: all of x has to be +ve integers.");
  }
  if (par.size() != 2) {
    stop("llik_igpd: length of par has to be 2.");
  }
  const double shape = par[0],
                          sigma = par[1],
                                     sigmau = sigma + shape * u;
  const LogicalVector above = x > u;
  const NumericVector
  xu(x[above]), cu(count[above]);
  const double nu(sum(cu));
  double l;
  if (u <= 1 || u <= min(x) || u >= max(x) ||
      sigma <= 0.0 || sigmau <= 0.0) {
    l = -INFINITY;
  }
  else {
    NumericVector yu, zu;
    if (shape != 0.0) {
      yu = 1.0 + shape / sigmau * (xu - 1.0 - u);
      zu = 1.0 + shape / (sigmau + shape * (xu - 1.0 - u));
      if (is_true(any(yu <= 0.0))) {
        l = -INFINITY;
      }
      else {
        l = sum(cu * log(1.0 - pow(zu, -1.0 / shape))) - 1.0 / shape * sum(cu * log(yu));
      }
    }
    else {
      yu = (xu - 1.0 - u) / sigmau;
      l = nu * log(1.0 - exp(-1.0 / sigmau)) - sum(cu * yu);
    }
    l += nu * log(phiu);
  }
  return lnan(l);
}

// [[Rcpp::export]]
const double lpost_igpd(const NumericVector par,
                        const IntegerVector x,
                        const IntegerVector count,
                        const int u,
                        const double m_shape,
                        const double s_shape,
                        const double a_sigma,
                        const double b_sigma,
                        const double phiu) {
  if (x.size() != count.size()) {
    stop("lpost_igpd: lengths of x & count have to be equal.");
  }
  const double shape = par[0], sigma = par[1];
  double l;
  if (u <= 1 || u <= min(x) || u >= max(x) ||
      sigma <= 0.0) {
    l = -INFINITY;
  }
  else {
    l =
      llik_igpd(par, x, count, u, phiu) +
      ldnorm(shape, m_shape, s_shape) +
      ldgamma(sigma, a_sigma, b_sigma);
  }
  return lnan(l);
}




// [[Rcpp::export]]
void igp_mcmc(int N, IntegerVector x, IntegerVector count, int u){
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
































