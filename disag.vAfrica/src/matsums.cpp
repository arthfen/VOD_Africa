#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec Arma_colSums(const arma::mat& x) {
  return arma::sum(x, 0);
}

// [[Rcpp::export]]
double fitCpp(arma::vec rho, const arma::colvec& Qty, const arma::mat& R, const arma::mat& RtR, const double& r, const double& N, const arma::vec& b0, const arma::mat& S1, const arma::mat& S2, const arma::mat& S3, const arma::mat& S4, const arma::mat& S5, const arma::mat& S6, const arma::mat& S7, const arma::mat& S8, const arma::mat& S9, const arma::mat& S10, const arma::mat& S11, const arma::mat& S12, const arma::mat& S13, const arma::mat& S14, const arma::mat& S15, const arma::mat& S16, const double& sldiagA, const int& it){
	rho = arma::exp(rho);
	if(min(rho) <= 1e-10) return arma::datum::inf;
	if(max(rho) > 1e10) return arma::datum::inf;

	arma::mat Sl = rho(1) * S1 + rho(2) * S2 + rho(3) * S3 + rho(4) * S4 + rho(5) * S5 + rho(6) * S6 + rho(7) * S7 + rho(8) * S8 + rho(9) * S9 + rho(10) * S10 + rho(11) * S11 + rho(12) * S12 + rho(13) * S13 + rho(14) * S14 + rho(15) * S15 + rho(16) * S16;
	double ldetSl;
	try {
		ldetSl = arma::log_det_sympd(Sl.rows(3, Sl.n_rows - 1).cols(3, Sl.n_cols - 1));
	} catch(...) {
		return arma::datum::inf;
	}

	arma::mat V = RtR/rho(0) + Sl;
	double ldetV;
	try {
		ldetV = arma::log_det_sympd(V);
	} catch(...) {
		return arma::datum::inf;
	}

	arma::mat g0 = 2*R.t()*Qty/rho(0) - 2*Sl*b0; // gradient
	arma::mat H0 = -2*V; // approximate Hessian
	arma::mat hat_beta;
	try {
		hat_beta = b0 - solve(H0, g0, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);
	} catch(...) {
		return arma::datum::inf;
	}

	arma::mat logLik = - Qty.t()*Qty/rho(0) - r/rho(0);
	logLik += - b0.t() * Sl * b0;
	logLik += 2*Qty.t()*R/rho(0)*(hat_beta - b0) - 2*b0.t()*Sl*(hat_beta - b0);
	logLik += - (hat_beta - b0).t() * V * (hat_beta - b0);
	logLik += - sldiagA - N * std::log(rho(0));
	logLik += ldetSl;
	logLik += - N * std::log(2 * M_PI) - ldetV;
	
	logLik += - 100/pow(2, (it-1)) * (hat_beta - b0).t() * (hat_beta - b0); // attempt to penalize solutions that fall very far
	
	return -0.5 * arma::as_scalar(logLik);
}

// [[Rcpp::export]]
arma::vec hbCpp(arma::vec rho, const arma::colvec& Qty, const arma::mat& R, const arma::mat& RtR, const double& r, const double& N, const arma::vec& b0, const arma::mat& S1, const arma::mat& S2, const arma::mat& S3, const arma::mat& S4, const arma::mat& S5, const arma::mat& S6, const arma::mat& S7, const arma::mat& S8, const arma::mat& S9, const arma::mat& S10, const arma::mat& S11, const arma::mat& S12, const arma::mat& S13, const arma::mat& S14, const arma::mat& S15, const arma::mat& S16){
	rho = arma::exp(rho);

	arma::mat Sl = rho(1) * S1 + rho(2) * S2 + rho(3) * S3 + rho(4) * S4 + rho(5) * S5 + rho(6) * S6 + rho(7) * S7 + rho(8) * S8 + rho(9) * S9 + rho(10) * S10 + rho(11) * S11 + rho(12) * S12 + rho(13) * S13 + rho(14) * S14 + rho(15) * S15 + rho(16) * S16;
	arma::mat V = RtR/rho(0) + Sl;
	arma::mat g0 = 2*R.t()*Qty/rho(0) - 2*Sl*b0; // gradient
	arma::mat H0 = -2*V; // approximate Hessian
	arma::mat hat_beta = b0 - solve(H0, g0, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);

	return hat_beta;
}

// [[Rcpp::export]]
double ll_at_min(arma::vec rho, const arma::colvec& Qty, const arma::mat& R, const arma::mat& RtR, const double& r, const double& N, const arma::vec& b0, const arma::mat& S1, const arma::mat& S2, const arma::mat& S3, const arma::mat& S4, const arma::mat& S5, const arma::mat& S6, const arma::mat& S7, const arma::mat& S8, const arma::mat& S9, const arma::mat& S10, const arma::mat& S11, const arma::mat& S12, const arma::mat& S13, const arma::mat& S14, const arma::mat& S15, const arma::mat& S16, const double& sldiagA){
	rho = arma::exp(rho);
	if(min(rho) <= 1e-10) return arma::datum::inf;
	if(max(rho) > 1e10) return arma::datum::inf;

	arma::mat Sl = rho(1) * S1 + rho(2) * S2 + rho(3) * S3 + rho(4) * S4 + rho(5) * S5 + rho(6) * S6 + rho(7) * S7 + rho(8) * S8 + rho(9) * S9 + rho(10) * S10 + rho(11) * S11 + rho(12) * S12 + rho(13) * S13 + rho(14) * S14 + rho(15) * S15 + rho(16) * S16;
	double ldetSl;
	try {
		ldetSl = arma::log_det_sympd(Sl.rows(3, Sl.n_rows - 1).cols(3, Sl.n_cols - 1));
	} catch(...) {
		return arma::datum::inf;
	}

	arma::mat V = RtR/rho(0) + Sl;
	double ldetV;
	try {
		ldetV = arma::log_det_sympd(V);
	} catch(...) {
		return arma::datum::inf;
	}

	arma::mat logLik = - Qty.t()*Qty/rho(0) - r/rho(0);
	logLik += - b0.t() * Sl * b0;
	logLik += - sldiagA - N * std::log(rho(0));
	logLik += ldetSl;
	logLik += - N * std::log(2 * M_PI) - ldetV;
	
	return -0.5 * arma::as_scalar(logLik);
}

// [[Rcpp::export]]
arma::vec getIndex(const arma::vec& v, const int b){
	arma::vec res(v.size());

	int j = 0;
	for (int i = 0; i != v.size(); i++) {
		if (v[i] == b){
			res[j] = i;
			j++;
		}
	}
	res = res.head(j);

	return res;
}

// [[Rcpp::export]]
arma::uvec getRowIndex_NA(const arma::mat& m){
	int n = m.n_rows;
	arma::uvec res = find_nonfinite(m);
	res -= floor(res/n)*n;

	return unique(res) + 1;
}

arma::vec g(const arma::vec& x){
	arma::vec v = arma::exp(x);
	return v;
}

arma::mat gmat(arma::mat& x){
	arma::mat v = arma::exp(x);
	return v;
}

arma::vec dg(const arma::vec& x){
	arma::vec v = arma::exp(x);
	return v;
}

// [[Rcpp::export]]
Rcpp::List tpmm_agg_noblock_fast(const arma::mat& Xs, const arma::vec& groups, const arma::vec& blocks, const arma::vec& b) {
	int m = Xs.n_cols;
	int n = groups.size();

	// declare the output vectors
	arma::vec gzd(n);
	arma::mat AGX(n, m);
	arma::vec MMt(n);
	arma::mat X_tmp;

	// declare the temporary vectors
	arma::uvec idx;
	arma::vec zd_tmp = Xs * b;

	for(int i = 0; i != n; i++) {
		idx = arma::find(blocks == groups[i]);

		X_tmp = Xs.rows(idx);
		X_tmp.each_col() %= dg(zd_tmp.rows(idx));
		AGX.row(i) = arma::sum(X_tmp, 0);
		gzd.row(i) = arma::sum(g(zd_tmp.rows(idx)));
		MMt.row(i) = idx.size();
	}

	return Rcpp::List::create(Rcpp::Named("gz.") = gzd, Rcpp::Named("AGX") = AGX, Rcpp::Named("MMt") = MMt);
}

// [[Rcpp::export]]
arma::vec tpmm_vec_q_noblock(const arma::mat& Xs, const arma::mat& b, int type) {
	// declare the temporary vectors
	arma::mat gzd = Xs * b;
//	gzd = gmat(gzd);

	arma::vec out;
	if(type == 1){
		arma::vec P_050 = {0.50};
		out = arma::quantile(gzd, P_050, 1); // type == 1 : median predictions
	} else if(type == 2){
		arma::vec P_005 = {0.05};
		out = arma::quantile(gzd, P_005, 1); // type == 2 : 5% percentile
	} else if(type == 3){
		arma::vec P_095 = {0.95};
		out = arma::quantile(gzd, P_095, 1); // type == 3 : 95% percentile
	} else if(type == 4){
		out = arma::stddev(gzd, 0, 1);  // type == 4 : standard deviation
	}	

	return out;
}

// [[Rcpp::export]]
arma::vec tpmm_vec_noblock(const arma::mat& Xs, const arma::mat& b) {
//	return g(Xs * b);
	return Xs * b;
}

// [[Rcpp::export]]
arma::mat tpmm_mat_noblock(const arma::mat& Xs, const arma::mat& b) {
//	return g(Xs * b);
	return Xs * b;
}

// [[Rcpp::export]]
arma::vec tpmv_vec_noblock(const arma::mat& Xs, const arma::vec& b) {
//	return g(Xs * b);
	return Xs * b;
}
