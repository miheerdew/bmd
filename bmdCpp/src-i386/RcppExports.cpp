// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// initializeC
Rcpp::IntegerVector initializeC(int n, const arma::vec& cors, float alpha, bool conserv);
RcppExport SEXP bmdCpp_initializeC(SEXP nSEXP, SEXP corsSEXP, SEXP alphaSEXP, SEXP conservSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cors(corsSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type conserv(conservSEXP);
    rcpp_result_gen = Rcpp::wrap(initializeC(n, cors, alpha, conserv));
    return rcpp_result_gen;
END_RCPP
}
// bh_rejectC
Rcpp::IntegerVector bh_rejectC(const arma::vec& pvals, double alpha, bool conserv);
RcppExport SEXP bmdCpp_bh_rejectC(SEXP pvalsSEXP, SEXP alphaSEXP, SEXP conservSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pvals(pvalsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type conserv(conservSEXP);
    rcpp_result_gen = Rcpp::wrap(bh_rejectC(pvals, alpha, conserv));
    return rcpp_result_gen;
END_RCPP
}
// updateC
Rcpp::IntegerVector updateC(const arma::mat& first, const arma::mat& second, const arma::rowvec& f4ColSums, const arma::mat& f2, const arma::mat& f3, const arma::mat& cross, float alpha, bool conserv);
RcppExport SEXP bmdCpp_updateC(SEXP firstSEXP, SEXP secondSEXP, SEXP f4ColSumsSEXP, SEXP f2SEXP, SEXP f3SEXP, SEXP crossSEXP, SEXP alphaSEXP, SEXP conservSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type first(firstSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type second(secondSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type f4ColSums(f4ColSumsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f3(f3SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cross(crossSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type conserv(conservSEXP);
    rcpp_result_gen = Rcpp::wrap(updateC(first, second, f4ColSums, f2, f3, cross, alpha, conserv));
    return rcpp_result_gen;
END_RCPP
}
// pvalsC
Rcpp::NumericVector pvalsC(const arma::mat& first, const arma::mat& second, const arma::rowvec& f4ColSums, const arma::mat& f2, const arma::mat& f3, const arma::mat& cross);
RcppExport SEXP bmdCpp_pvalsC(SEXP firstSEXP, SEXP secondSEXP, SEXP f4ColSumsSEXP, SEXP f2SEXP, SEXP f3SEXP, SEXP crossSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type first(firstSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type second(secondSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type f4ColSums(f4ColSumsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f3(f3SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cross(crossSEXP);
    rcpp_result_gen = Rcpp::wrap(pvalsC(first, second, f4ColSums, f2, f3, cross));
    return rcpp_result_gen;
END_RCPP
}