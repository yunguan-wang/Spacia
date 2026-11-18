//--------------------------------------------------------------
// Header (header)
//--------------------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;


//--------------------------------------------------------------
// Functions (Functions_cpp)
//--------------------------------------------------------------
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return repmat(mu, 1, n).t() + Y * chol(sigma);
}

double log_sum_exp(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  double maxv = x(maxi);
  if (!(maxv > -datum::inf)) {
    return -datum::inf;
  }
  double cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) & (x(i) > -datum::inf)) {
      cumsum += exp(x(i) - maxv);
    }
  }
  return maxv + log1p(cumsum);
}


double c_rtexp_onesided(double a)
{
  bool stop = false;
  double lambda = 0.5 * (a + sqrt(pow(a, 2) + 4));
  double z;
  
  while(!stop)
  {
    z = a - log(R::runif(0,1))/lambda;
    stop = (-2* log(R::runif(0,1)) > pow(z - lambda, 2));
  }
  return z;
}

// [[Rcpp::export]]
arma::vec ctruncnorm(arma::vec mean, 
                     arma::vec sd, 
                     bool upper){
  int m = mean.size();
  vec p = normcdf(- mean / sd);
  vec res = zeros(m);
  if(upper){
    // right truncated
    p += (1.0 - p) % randu(m);
    for(int j = 0; j < m; j++) {
      if(- mean(j) / sd(j) > 3.48672170399){
        res(j) = sd(j) * c_rtexp_onesided(- mean(j) / sd(j)) + mean(j);
      } else{
        res(j) = R::qnorm(p(j), mean(j), sd(j), true, false);
      }
    }
  } else{
    // left truncated
    p = p % randu(m);
    for(int j = 0; j < m; j++) {
      if(- mean(j) / sd(j) < -3.48672170399){
        res(j) = - sd(j) * c_rtexp_onesided(mean(j) / sd(j)) + mean(j);
      } else{
        res(j) = R::qnorm(p(j), mean(j), sd(j), true, false);
      }
    }
  }
  return res;
}


// Allocated once per chain. Fixed designs and label indices do not change.
struct GibbsWorkspace {
  mat X_delta, X1_delta, V_beta, Xb_design;
  vec mu_z, m_beta, mu_u, probit_prob, m_b, ones_bags, ones_instances;
  uvec positive_labels, negative_labels;
  GibbsWorkspace(const mat& Xb, const mat& Xbeta, const vec& y)
      : X_delta(y.n_elem, Xbeta.n_cols, fill::zeros),
        X1_delta(y.n_elem, Xbeta.n_cols + 1, fill::ones),
        V_beta(Xbeta.n_cols + 1, Xbeta.n_cols + 1),
        Xb_design(join_rows(ones<vec>(Xb.n_rows), Xb)),
        mu_z(y.n_elem), m_beta(Xbeta.n_cols + 1), mu_u(Xb.n_rows),
        probit_prob(Xb.n_rows), m_b(Xb.n_cols + 1),
        ones_bags(y.n_elem, fill::ones), ones_instances(Xb.n_rows, fill::ones),
        positive_labels(find(y == 1)), negative_labels(find(y == 0)) {}
};

void gibbs_step(
    const arma::mat& Xb,
    const arma::mat& Xbeta,
    const arma::vec& y,
    const arma::vec& ninst,
    
    const arma::vec& hp_mu_beta,
    const arma::vec& hp_mu_b,
    const arma::mat& hp_Sig_beta,
    const arma::mat& hp_Sig_b,
    
    arma::vec& beta,
    arma::vec& b,
    arma::vec& delta,
    arma::vec& u,
    arma::vec& z,
    
    const arma::mat& hp_Sig_beta_inv,
    const arma::mat& hp_Sig_b_inv,
    const arma::mat& V_b, const arma::mat& chol_V_b, GibbsWorkspace& ws){
  
  int n = y.size(), db = Xb.n_cols,dbeta = Xbeta.n_cols;
  
  // set X_delta, mu_z
  mat& X_delta = ws.X_delta;
  X_delta.zeros();
  vec& mu_z = ws.mu_z;
  int pos = 0;
  for(int i = 0; i < n; i++){
    if(accu(delta.subvec(pos, pos + ninst[i] - 1)) != 0){
      X_delta.row(i) = sum(Xbeta.rows(pos + find(delta.subvec(pos, pos + ninst[i] - 1) == 1)),
                  0); // column sum
    }
    pos += ninst[i];
  }
  mu_z = beta(0) + X_delta * beta.tail(dbeta);
  
  // Rcout << "mu_z=\n" << mu_z << "\n";
  
  // update z
  const vec& one_colvec = ws.ones_bags;
  //// y == 1
  const uvec& idx1 = ws.positive_labels;
  z(idx1) = ctruncnorm(mu_z(idx1), one_colvec(idx1), true);
  //// y == 0
  const uvec& idx0 = ws.negative_labels;
  z(idx0) = ctruncnorm(mu_z(idx0), one_colvec(idx0), false);
  // z[y==1]<-rtruncnorm(sum(y==1),mean=mu_z[y==1],sd=1,a=0,b=Inf)
  //   z[y==0]<-rtruncnorm(sum(y==0),mean=mu_z[y==0],sd=1,a=-Inf,b=0)
  
  
  // update beta
  mat& X1_delta = ws.X1_delta;
  X1_delta.cols(1, dbeta) = X_delta;
  // X1_delta<-cbind(1,X_delta)
  mat& V_beta = ws.V_beta;
  V_beta = inv_sympd(hp_Sig_beta_inv + X1_delta.t() * X1_delta);
  vec& m_beta = ws.m_beta;
  m_beta = V_beta * (hp_Sig_beta_inv * hp_mu_beta + X1_delta.t() * z);
  beta = m_beta + (randn(1, dbeta + 1) * chol(V_beta)).t();
  
  // update delta
  vec& mu_u = ws.mu_u;
  mu_u = b(0) + Xb * b.tail(db);
  vec& probit_prob = ws.probit_prob;
  probit_prob = normcdf(mu_u);
  
  // Rcout << "mu_u=\n" << mu_u << "\n";
  pos = 0;
  double A, B;
  for(int i = 0; i < n; i++){

    // Preserve upstream's frozen within-bag delta state. This total is
    // invariant inside the sender loop; do not change the transition kernel.
    const double selected_total = accu(
        Xbeta.rows(pos + find(delta.subvec(pos, pos + ninst[i] - 1) == 1)) * beta.tail(dbeta));
    for(int j = 0; j < ninst[i]; j++){
      // This sender has not yet been updated, so its current delta equals
      // the frozen value. Use the original total throughout this bag.
      const double contribution = dot(Xbeta.row(pos + j), beta.tail(dbeta));
      double tmp = z(i) - beta(0) - selected_total + delta(pos + j) * contribution;
      
      A = exp(-0.5 * pow(tmp - contribution, 2.0));
      B = exp(-0.5 * pow(tmp, 2.0));
      
      double prim_prob = (A * probit_prob(pos + j)) / 
        (A * probit_prob(pos + j) + B * (1 - probit_prob(pos + j)));
      
      if(std::isnan(prim_prob)){
        prim_prob = 0;
      }
      // Rcout << "prim_prob=\n" << prim_prob << "\n";
      if(conv_to<double>::from(randu(1)) < prim_prob){
        delta(pos + j) = 1;
      } else{
        delta(pos + j) = 0;
      }  
    }
    pos += ninst[i];
  }
  
  // update u
  //// delta == 1

  const vec& one_colvec2 = ws.ones_instances;
  uvec idx_delta1 = find(delta == 1);
  u(idx_delta1) = ctruncnorm(mu_u(idx_delta1), one_colvec2(idx_delta1), true);
  //// delta == 0
  uvec idx_delta2 = find(delta == 0);
  u(idx_delta2) = ctruncnorm(mu_u(idx_delta2), one_colvec2(idx_delta2), false);
  // u [delta==1]<-rtruncnorm(sum(delta==1),mean=mu_u[delta==1],sd=1,a=0,b=Inf);
  // u[delta==0]<-rtruncnorm(sum(delta==0),mean=mu_u[delta==0],sd=1,a=-Inf,b=0)

  // // update u
  // //// delta == 1
  // one_colvec = ones(X.n_rows);
  // uvec idx_delta1 = find(delta == 1);
  // u(idx) = ctruncnorm(mu_u(idx), one_colvec(idx), true);
  // //// delta == 0
  // idx = find(delta == 0);
  // u(idx) = ctruncnorm(mu_u(idx), one_colvec(idx), false);
  // // u [delta==1]<-rtruncnorm(sum(delta==1),mean=mu_u[delta==1],sd=1,a=0,b=Inf);
  // // u[delta==0]<-rtruncnorm(sum(delta==0),mean=mu_u[delta==0],sd=1,a=-Inf,b=0)
  
  // update b
  // one_colvec = ones(X.n_rows);
  vec& m_b = ws.m_b;
  m_b = V_b * (hp_Sig_b_inv * hp_mu_b + ws.Xb_design.t() * u);
  b = m_b + (randn(1, db + 1) * chol_V_b).t();
  // Rcout << "here3\n";
  
}

// [[Rcpp::export]]
Rcpp::List MICProB_1Gibbs_cpp(
    arma::mat Xb,
    arma::mat Xbeta,
    arma::vec y,
    arma::vec ninst,

    arma::vec hp_mu_beta,
    arma::vec hp_mu_b,
    arma::mat hp_Sig_beta,
    arma::mat hp_Sig_b,

    arma::vec beta,
    arma::vec b,
    arma::vec delta,
    arma::vec u,
    arma::vec z,

    arma::mat hp_Sig_beta_inv,
    arma::mat hp_Sig_b_inv,
    arma::mat V_b){

  GibbsWorkspace ws(Xb, Xbeta, y);
  gibbs_step(Xb, Xbeta, y, ninst, hp_mu_beta, hp_mu_b, hp_Sig_beta,
             hp_Sig_b, beta, b, delta, u, z, hp_Sig_beta_inv,
             hp_Sig_b_inv, V_b, chol(V_b), ws);
  return List::create(Named("beta")=beta, Named("b")=b,
                      Named("delta")=delta, Named("u")=u, Named("z")=z);
}

// Run a whole chain with one R/C++ transition. Keep upstream RNG draw order,
// initial output row, thinning schedule, and frozen within-bag delta updates.
// [[Rcpp::export]]
Rcpp::List MICProB_chain_cpp(
    const arma::mat& Xb, const arma::mat& Xbeta,
    const arma::vec& y, const arma::vec& ninst,
    const arma::vec& hp_mu_beta, const arma::vec& hp_mu_b,
    const arma::mat& hp_Sig_beta, const arma::mat& hp_Sig_b,
    arma::vec beta, arma::vec b, arma::vec delta,
    const arma::mat& hp_Sig_beta_inv, const arma::mat& hp_Sig_b_inv,
    const arma::mat& V_b, int ntotal, int nwarm, int nthin,
    bool return_delta) {
  if (ntotal <= nwarm || nwarm < 0 || nthin <= 0 ||
      (ntotal - nwarm) % nthin != 0)
    Rcpp::stop("Require ntotal > nwarm >= 0, nthin > 0, and post-warmup iterations divisible by nthin.");
  const int nsave = (ntotal - nwarm) / nthin;
  arma::mat beta_post(nsave+1, beta.n_elem), b_post(nsave+1, b.n_elem);
  arma::mat delta_post;
  if (return_delta) {
    delta_post.set_size(nsave+1, delta.n_elem);
    delta_post.row(0) = delta.t();
  }
  beta_post.row(0) = beta.t();
  b_post.row(0) = b.t();
  arma::vec u(delta.n_elem, arma::fill::zeros), z(y.n_elem, arma::fill::zeros);
  arma::vec pip(delta.n_elem, arma::fill::zeros);
  const arma::mat chol_V_b = chol(V_b);
  GibbsWorkspace ws(Xb, Xbeta, y);
  for (int iter=1; iter<=ntotal; ++iter) {
    if (iter % 100 == 0) Rcpp::checkUserInterrupt();
    gibbs_step(Xb, Xbeta, y, ninst, hp_mu_beta, hp_mu_b,
               hp_Sig_beta, hp_Sig_b, beta, b, delta, u, z,
               hp_Sig_beta_inv, hp_Sig_b_inv, V_b, chol_V_b, ws);
    const int post_iter = iter - nwarm;
    if (post_iter > 0 && post_iter % nthin == 0) {
      const int row = post_iter / nthin;
      beta_post.row(row) = beta.t();
      b_post.row(row) = b.t();
      pip += delta;
      if (return_delta) delta_post.row(row) = delta.t();
    }
  }
  Rcpp::List out = Rcpp::List::create(
      Named("beta")=beta_post, Named("b")=b_post, Named("pip")=pip/nsave);
  if (return_delta) out["delta"] = delta_post;
  return out;
}
