// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


arma::uvec bh_reject(const arma::vec &pvals, double alpha, bool conserv){
  int m = pvals.n_elem;
  
  if(conserv) {
    double mults = 0;
    for(int j = 1; j <= m; j++) {
      mults += 1.0/j ;
    }
    alpha = alpha/mults;
  }
  
  arma::vec sorted = sort(pvals);
  
  int j;
  for(j = m; j > 0; j--){
    if(sorted(j-1) <= alpha*j/m) break;
  }
  
  if(j==0){
    return arma::uvec();
  } else {
    return(find(pvals <= alpha*j/m));
  }
}

arma::vec pvals(const arma::mat& first,
          const arma::mat& second,
          const arma::rowvec& f4sum,
          const arma::mat& f2,
          const arma::mat& f3,
          const arma::mat& cross){
  arma::vec corsums = arma::sum(cross, 1);
  int n = first.n_rows;
  arma::vec sRowSum = arma::sum(second, 1); // n x 1
  arma::mat sRowSum2 = cross * arma::pow(second,2).t(); // p x n
  arma::vec star1 = f2.t() * arma::pow(sRowSum,2); // p x 1
  arma::vec star2 = f4sum.t() % arma::pow(corsums,2); //p x 1
  arma::vec star3 = 2 * corsums % arma::sum(f2 % sRowSum2.t()).t(); // p x 1
  arma::vec star4 = arma::sum(arma::pow(sRowSum2,2),1); // p x 1
  arma::vec dagger1 = corsums % (f3.t() * sRowSum); // p x 1
  
  arma::mat aux = first % sRowSum2.t(); // n x p
  arma::rowvec dagger2 = arma::sum(aux.each_col() % sRowSum); // 1 x p
  
  arma::vec vars = (star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2.t()) / (n - 1);
  
  arma::vec pvals = sqrt(n) * corsums / arma::sqrt(vars);
  pvals.transform([](double zstat) {return R::pnorm(zstat,0,1,false,false);});
  
  return pvals;
}

arma::uvec initialize(int n, const arma::vec& cors, float alpha, bool conserv){
  arma::vec pvalues = atanh(cors)*sqrt(n-3); //fischer transform
  pvalues.transform([](double stat) {return R::pnorm(stat,0,1,false,false);});
  return bh_reject(pvalues, alpha, conserv);
}

//' @title
//' inititalizeC
//' 
//' @param n is the sample size
//' 
//' @param cors is the vector of correlations to the initital node.
//'
//' @param alpha is the significance level for the BH test
//'
//' @param conserv (=TRUE or FALSE). If conserv=TRUE, BHY test is used at alpha 
//' otherwise BH test is used at alpha.
//' @details
//' \code{initializeC} takes a vector of correlations, firscher transforms it 
//'  (to get a z-statistic) and computes the pvalues. The BH (or BHY) procedure is 
//'  used on these pvales to reject a set of indices which are then returned
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector initializeC(int n, const arma::vec& cors, float alpha, bool conserv){
  arma::uvec res = initialize(n, cors, alpha, conserv);
  res += 1;
  return Rcpp::IntegerVector(res.begin(), res.end());
}

//' @title
//' bh_rejectC
//' @description
//' The rcpp version of bhreject 
//' 
//' @param pvals is the vector of pvalues
//'
//' @param alpha is the significance level for the BH test
//' 
//' @param conserv (=TRUE or FALSE). If conserv=TRUE, BHY test is used at level alpha 
//' otherwise BH test is used at level alpha.
//' @details
//' \code{bh_rejectC} takes a vector of p-values and applies the Benjamini Hochberg (Yekutelli)
//' multiple testing procedure at significane level alpha, and returns the indices 
//' of the pvalues which were rejected.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector bh_rejectC(const arma::vec &pvals, double alpha, bool conserv){
  arma::uvec res = bh_reject(pvals, alpha, conserv);
  res += 1;
  return Rcpp::IntegerVector(res.begin(), res.end());
}

//' @title
//' updateC
//' @description
//' The rcpp version of update5. update5 is first step for a one sided update.
//' One starts with a subset A of X or Y indices, and wishes to get a set B=update5(A)
//' at the other side. Call the otherside as the first side (or the side being tested)
//' and side corresponding to A as the second side.
//' 
//' @param first (either X or Y) is the matrix corresponding to the first side.
//' 
//' @param second (e.g X[,A] or Y[,A]) is the sub-matrix (of the second side)
//' corresponding to the columns of the A.
//' 
//' @param f4ColSums is the rowvector columnSums(first^4)
//' 
//' @param f2 is the matrix first^2
//' 
//' @param f3 is the matrix first^3
//'
//' @param cross is the matrix cor(first, second)
//'
//' @param alpha is the significance level for the BH test
//'
//' @param conserv (=TRUE or FALSE). If conserv=TRUE, BHY test is used at alpha 
//' otherwise BH test is used at alpha.
//' @details
//' \code{updateC} computes the z-statistic corresponding to the sum of 
//' correlations to the second matrix for each column of the first matrix
//' (using the variance formula from SteigerHakstian1982). It converts the 
//' z-statistic to pvalues and rejects the columns of first matrix which 
//' were significant. The indicies corresponding to the significant columns are
//' returned.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector updateC(const arma::mat& first,
                            const arma::mat& second,
                            const arma::rowvec& f4ColSums,
                            const arma::mat& f2,
                            const arma::mat& f3,
                            const arma::mat& cross,
                            float alpha,
                            bool conserv){
  arma::uvec res = bh_reject(pvals(first, second, f4ColSums, f2, f3, cross), alpha, conserv);
  res += 1;
  return Rcpp::IntegerVector(res.begin(), res.end());
}


//' @title
//' pvalsC
//' 
//' @description
//' The p-values in the one-sided updated step. first corresponds to the side 
//' (X or Y) the p-values will be generated for. This is side opposite to 
//' that of the index set A that is being used.
//' 
//' @param first (either X or Y) is the matrix corresponding to the first side.
//' 
//' @param second (e.g X[,A] or Y[,A]) is the sub-matrix (of the second side)
//' corresponding to the columns of the A.
//' 
//' @param f4ColSums is the rowvector columnSums(first^4)
//' 
//' @param f2 is the matrix first^2
//' 
//' @param f3 is the matrix first^3
//'
//' @param cross is the matrix cor(first, second)
//'
//' @details
//' \code{pvalsC} computes the z-statistic corresponding to the sum of 
//' correlations to the second matrix for each column of the first matrix
//' (using the variance formula from SteigerHakstian1982). 
//' The pvalue vector corresponding to the columns of first matrix is returned.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pvalsC(const arma::mat& first,
                           const arma::mat& second,
                           const arma::rowvec& f4ColSums,
                           const arma::mat& f2,
                           const arma::mat& f3,
                           const arma::mat& cross){
  arma::vec res = pvals(first, second, f4ColSums, f2, f3, cross);
  return Rcpp::NumericVector(res.begin(), res.end());
}