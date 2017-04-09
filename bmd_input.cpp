#include <RcppArmadillo.h>

using namespace arma;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

uvec bh_reject(vec pvals, double alpha, bool conserv){
  uint m = pvals.n_elem;

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

// [[Rcpp::export]]
Rcpp::IntegerVector bh_rejectC(vec pvals, double alpha, bool conserv = true){
  return Rcpp::wrap(bh_reject(pvals, alpha, conserv)+1);
}

class BmdInput {

public:
  const mat X, Y;
  mat cormat;
  unsigned dx, dy, n;
  const double alpha;

private:
  const mat X2, X3, Y2, Y3; //ORDER DEPENDECY
  const rowvec X4sum, Y4sum;

public:

  BmdInput(const mat& X_scaled, const mat& Y_scaled, double alpha_=0.05) :
    X(X_scaled), Y(Y_scaled),
    X2(pow(X,2)), X3(X2%X), X4sum(sum(pow(X2,2))),
    Y2(pow(Y,2)), Y3(Y2%Y), Y4sum(sum(pow(Y2,2))),
    dx(X_scaled.n_cols),
    dy(Y_scaled.n_cols),
    n(X_scaled.n_rows),
    alpha(alpha_)
    {
      Rprintf("Computing cross-correlation matrix\n");
      cormat = cor(X_scaled, Y_scaled);
    };

  Rcpp::NumericMatrix cross_cors(uvec A, uvec B){
    mat ret = cormat.submat(A-1,B-1);
    return Rcpp::wrap(ret);
  }

  Rcpp::NumericVector vars_wrap(uvec A, bool test_x){
    set_sides(A-1, test_x);
    return Rcpp::wrap(vars());
  }

  Rcpp::NumericMatrix get_X(){
    return Rcpp::wrap(X);
  }

  Rcpp::NumericMatrix get_Y(){
    return Rcpp::wrap(Y);
  }

  Rcpp::NumericVector pvals_wrap(uvec A, bool test_x){
    set_sides(A-1, test_x);
    return Rcpp::wrap(pvals());
  }

  Rcpp::IntegerVector update_wrap(uvec A, bool test_x, bool conserv){
    return Rcpp::wrap(update(A-1,test_x, conserv)+1);
  }

  Rcpp::IntegerVector init_wrap(uint u, bool test_x, bool conserv){
    return Rcpp::wrap(init(u-1, test_x, conserv)+1);
  }

  uvec update(uvec A, bool test_x, bool conserv){
    set_sides(A, test_x);
    return bh_reject(pvals(), alpha, conserv);
  }

  uvec init(uint u, bool test_x, bool conserv){
    vec cors;
    if (test_x) {
      cors = cormat.col(u);
    } else {
      cors = cormat.row(u).t();
    }
    vec pvalues = atanh(cors)*sqrt(n-3); //fischer transform
    pvalues.transform([](double stat) {return R::pnorm(stat,0,1,false,false);});
    return bh_reject(pvalues, alpha, conserv);
  }

protected :
  const mat *first, *f2, *f3; const rowvec *f4sum;
  mat cross, second;
  vec corsums;

  vec pvals(){
    vec pvals = sqrt(n) * corsums / sqrt(vars());
    pvals.transform([](double zstat) {return R::pnorm(zstat,0,1,false,false);});
    return pvals;
  }

  void set_sides(uvec B, bool test_x){
    if (test_x){
      cross = cormat.cols(B);
      first = &X; f2 = &X2; f3=&X3; f4sum = &X4sum;
      second = Y.cols(B);
    } else {
      cross = cormat.rows(B).t();
      first = &Y; f2 = &Y2; f3=&Y3; f4sum = &Y4sum;
      second = X.cols(B);
    }
    corsums = sum(cross, 1);
  }

  vec vars(){
    vec sRowSum = sum(second, 1); // n x 1
    mat sRowSum2 = cross * pow(second,2).t(); // p x n
    vec star1 = f2->t() * pow(sRowSum,2); // p x 1
    vec star2 = f4sum->t() % pow(corsums,2); //p x 1
    vec star3 = 2 * corsums % sum((*f2) % sRowSum2.t()).t(); // p x 1
    vec star4 = sum(pow(sRowSum2,2),1); // p x 1
    vec dagger1 = corsums % (f3->t() * sRowSum); // p x 1

    mat aux = (*first) % sRowSum2.t(); // n x p
    rowvec dagger2 = sum(aux.each_col() % sRowSum); // 1 x p

    return (star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2.t()) / (n - 1);
  }

};

RCPP_MODULE(yada) {
 Rcpp::class_<BmdInput>( "BmdInput" )
    .constructor<mat,mat>()
    .method("cross_cors", &BmdInput::cross_cors, "Usage: cross_cors(A, B)")
    .method("vars", &BmdInput::vars_wrap, "Usage: vars(A, test_x)")
    .method("pvals", &BmdInput::pvals_wrap, "Usage : pvals(A, test_x)")
    .method("init", &BmdInput::init_wrap, "Usage : init(comm_id, test_x)")
    .method("update", &BmdInput::update_wrap)
    .field_readonly( "n", &BmdInput::n)
    .field_readonly( "dx", &BmdInput::dx)
    .field_readonly( "dy", &BmdInput::dy)
    .property( "X", &BmdInput::get_X)
    .property( "Y", &BmdInput::get_Y)
   ;
 }
