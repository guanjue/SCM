

#include <Rcpp.h>
#include <math.h> 
using namespace Rcpp;

// [[Rcpp::export]]
double get_nloglik_s(double y, int x, double a0, double a1, double b0, double b1, double lam) {
  double tmp1, tmp2, nloglik_s, x1, x2;

  x1 = a0+a1*x;
  x2 = b0+b1*x;

  if (y!=0) 
  {
    tmp1 = - lgamma(y+1) + lgamma(x1+y) - lgamma(x1) + y*log(x2);
    tmp2 = - (x1+y)*log(x2+1) + log(1 - pow(((x2+1)/((x2)*(1+lam)+1)),(x1+y)));
  }
  else 
  {
  tmp1 = 0;
  {
  if(lam>1)
  tmp2 = - (x1)*log(x2+1) + log(1 - pow(((x2+1)/((x2)*(1+lam)+1)),(x1)) + pow(((x2+1)/((x2)*lam+1)),(x1)));
  else
  tmp2 = - (x1)*log((x2)*lam+1) + log(pow((((x2)*lam+1)/(x2+1)),(x1)) - pow((((x2)*lam+1)/((x2)*(1+lam)+1)),(x1)) + 1);
  }
  }
  nloglik_s = -(tmp1+tmp2);

  return nloglik_s;
}

// [[Rcpp::export]]
double get_nloglik_c(double lam, NumericMatrix parg_all, NumericVector y_c, NumericVector matx_c) {
  double nloglik_c = 0;
  int n = y_c.size();

  for (int i = 0; i < n; i++) {
    nloglik_c += get_nloglik_s(y_c[i], matx_c[i], parg_all(i,0), parg_all(i,1), parg_all(i,2), parg_all(i,3), lam);
  }

  return nloglik_c;
}


// [[Rcpp::export]]
double get_nloglik_g(NumericVector parg, NumericVector parc_all, NumericVector y_g, NumericVector matx_g) {
  double nloglik_g = 0;
  int n = y_g.size();

  for (int i = 0; i < n; i++) {
    nloglik_g += get_nloglik_s(y_g[i], matx_g[i], parg[0], parg[1], parg[2], parg[3], parc_all[i]);
  }
  
  return nloglik_g;
}


// [[Rcpp::export]]
DP() {
  
}







