#include <RcppArmadillo.h>

using namespace arma;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]


class BmdInput {

  public:

  BmdInput(mat X_scaled, mat Y_scaled) :
    X(X_scaled), Y(Y_scaled),
    dx(X_scaled.n_cols),
    dy(Y_scaled.n_cols),
    n(X_scaled.n_rows),
    cormat(cor(X_scaled, Y_scaled))
    {};

  Rcpp::NumericMatrix cross_cors(uvec A, uvec B){
    mat ret = cormat.submat(A-1,B-1);
    return Rcpp::wrap(ret);
  }

  Rcpp::NumericVector vars_wrapper(uvec A, bool test_x){
    rowvec ret = vars(A-1,test_x);
    return Rcpp::wrap(ret);
  }

  Rcpp::NumericMatrix get_X(){
    return Rcpp::wrap(X);
  }

  Rcpp::NumericMatrix get_Y(){
    return Rcpp::wrap(Y);
  }

protected :
  inline rowvec ri4(){
    return sum(pow(*first,4));
  }

  inline rowvec rjjhh(uword j, uword h){
    double val = as_scalar(sum(square(second->col(j) % second->col(h))));
    rowvec ret(first->n_cols);
    ret.fill(val);
    return ret;
  }

  inline rowvec riijh(uword j, uword h){
    mat sq = square(*first);
    return sum(sq.each_col() % (second->col(j) % second->col(h)) );
  }

  inline rowvec ri3j(uword j){
    mat f3 = pow(*first,3);
    return sum(f3.each_col() % second->col(j));
  }

  inline rowvec rijjh(uword j, uword h){
    return sum( first->each_col() %
                ( square(second->col(j)) % second->col(h) )
              );
  }

  inline rowvec r(uword j){
    rowvec res;
    if(test_x){
      res = cormat.col(j).t();
    } else {
      res = cormat.row(j);
    }
    return res;
  }

  inline rowvec cov_rij_rih(uword j, uword h){
    //See Steiger and Hakstian 1982 equation (5.1) (with k=i)
    return
      (riijh(j, h) +
      0.25*r(j)%r(h)%( ri4() + riijh(j,j) + riijh(h,h) + rjjhh(j,h) )
      - 0.5*r(j)%( ri3j(h) + rijjh(j, h))
      - 0.5*r(h)%( ri3j(j) + rijjh(h, j)))/(n-1);
  }

  rowvec vars(uvec B, bool test_x){
    set_sides(test_x);
    rowvec res((*first).n_cols, fill::zeros);
    for(uword j : B){
      for(uword h: B){
        res += cov_rij_rih(j,h);
      }
    }
    return res;
  }

public:
  const mat X, Y, cormat;
  unsigned dx, dy, n;

private:
  const mat* first;
  const mat* second;
  bool test_x;

  void set_sides(bool test_x_){
    test_x  = test_x_;
    first = (test_x) ? &X : &Y;
    second = (test_x) ? &Y : &X;
  }
};

RCPP_MODULE(yada) {
 Rcpp::class_<BmdInput>( "BmdInput" )
    .constructor<mat,mat>()
    .method("cross_cors", &BmdInput::cross_cors)
    .method("vars", &BmdInput::vars_wrapper)
    .field_readonly( "n", &BmdInput::n)
    .field_readonly( "dx", &BmdInput::dx)
    .field_readonly( "dy", &BmdInput::dy)
    .property( "X", &BmdInput::get_X)
    .property( "Y", &BmdInput::get_Y)
   ;
 }
