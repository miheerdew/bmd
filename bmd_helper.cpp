#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]


uvec bh_reject0(const vec &pvals, double alpha, bool conserv){
  int m = pvals.n_elem;

  if(conserv) {
    double mults = 0;
    for(int j = 1; j <= m; j++) {
      mults += 1.0/j ;
    }
    alpha = alpha/mults;
  }

  vec sorted = sort(pvals);

  int j;
  for(j = m; j > 0; j--){
    if(sorted(j-1) <= alpha*j/m) break;
  }

  if(j==0){
    return uvec();
  } else {
    return(find(pvals <= alpha*j/m));
  }
}

vec pvals(const mat& first,
          const mat& second,
          const rowvec& f4sum,
          const mat& f2,
          const mat& f3,
          const mat& cross){
  vec corsums = sum(cross, 1);
  int n = first.n_rows;
  vec sRowSum = sum(second, 1); // n x 1
  mat sRowSum2 = cross * pow(second,2).t(); // p x n
  vec star1 = f2.t() * pow(sRowSum,2); // p x 1
  vec star2 = f4sum.t() % pow(corsums,2); //p x 1
  vec star3 = 2 * corsums % sum(f2 % sRowSum2.t()).t(); // p x 1
  vec star4 = sum(pow(sRowSum2,2),1); // p x 1
  vec dagger1 = corsums % (f3.t() * sRowSum); // p x 1

  mat aux = first % sRowSum2.t(); // n x p
  rowvec dagger2 = sum(aux.each_col() % sRowSum); // 1 x p

  vec vars = (star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2.t()) / (n - 1);

  vec pvals = sqrt(n) * corsums / sqrt(vars);
  pvals.transform([](double zstat) {return R::pnorm(zstat,0,1,false,false);});

  return pvals;
}

uvec initialize(int n, const vec& cors, float alpha, bool conserv){
  vec pvalues = atanh(cors)*sqrt(n-3); //fischer transform
  pvalues.transform([](double stat) {return R::pnorm(stat,0,1,false,false);});
  return bh_reject0(pvalues, alpha, conserv);
}


// [[Rcpp::export]]
Rcpp::IntegerVector initializeC(int n, const vec& cors, float alpha, bool conserv){
  uvec res = initialize(n, cors, alpha, conserv);
  res += 1;
  return Rcpp::IntegerVector(res.begin(), res.end());
}

// [[Rcpp::export]]
Rcpp::IntegerVector bh_rejectC(const vec &pvals, double alpha, bool conserv){
  uvec res = bh_reject0(pvals, alpha, conserv);
  res += 1;
  return Rcpp::IntegerVector(res.begin(), res.end());
}

// [[Rcpp::export]]
Rcpp::IntegerVector updateC(const mat& first,
                            const mat& second,
                            const rowvec& f4ColSums,
                            const mat& f2,
                            const mat& f3,
                            const mat& cross,
                            float alpha,
                            bool conserv){
  uvec res = bh_reject0(pvals(first, second, f4ColSums, f2, f3, cross), alpha, conserv);
  res += 1;
  return Rcpp::IntegerVector(res.begin(), res.end());
  }

// [[Rcpp::export]]
Rcpp::NumericVector pvalsC(const mat& first,
                           const mat& second,
                           const rowvec& f4ColSums,
                           const mat& f2,
                           const mat& f3,
                           const mat& cross){
  vec res = pvals(first, second, f4ColSums, f2, f3, cross);
  return Rcpp::NumericVector(res.begin(), res.end());
}
