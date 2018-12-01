

#include <Rcpp.h>
#include <math.h> 
using namespace Rcpp;

// [[Rcpp::export]]
double get_nloglik_s(double y, int x, double a0, double a1, double b0, double b1, double lam) {
  double tmp1, tmp2, nloglik_s;

  if (y!=0) 
  	tmp1 = - lgamma(y+1) + lgamma(a0+a1*x+y) - lgamma(a0+a1*x) + y*log(b0+b1*x);
  else 
  	tmp1 = 0;
  
  if(y==0) {if(lam>1)
	          tmp2 = - (a0+a1*x)*log(b0+b1*x+1) + log(1 - pow(((b0+b1*x+1)/((b0+b1*x)*(1+lam)+1)),(a0+a1*x)) + pow(((b0+b1*x+1)/((b0+b1*x)*lam+1)),(a0+a1*x)));
            else
	          tmp2 = - (a0+a1*x)*log((b0+b1*x)*lam+1) + log(pow((((b0+b1*x)*lam+1)/(b0+b1*x+1)),(a0+a1*x)) - pow((((b0+b1*x)*lam+1)/((b0+b1*x)*(1+lam)+1)),(a0+a1*x)) + 1);}
  else 
	  tmp2 = - (a0+a1*x+y)*log(b0+b1*x+1) + log(1 - pow(((b0+b1*x+1)/((b0+b1*x)*(1+lam)+1)),(a0+a1*x+y)));
  
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









