// slphyx@Shift-Enter

#include <Rmath.h>
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// modulus function
int mod(int a, const int n){
  // division reminder
  // a - dividend
  // n - divisor
  return(a - floor(a/n)*n);
}

// designed to be used with deSolve
// j - array index 0-4
// t - month
// [[Rcpp::export]]
double foi3(arma::vec Y, double t, int j, List params){
  //  inf_state_names = c("UV0I0", "UV0IR", "LV0I0", "LV0IR","UV1I0", "UV1IR", "LV1I0", "LV1IR","UV2I0", "UV2IR",
  //      "LV2I0", "LV2IR","UV3I0", "UV3IR", "LV3I0", "LV3IR")

  //position inf_states = {3, 7, 4, 8, 12, 16, 13, 17, 20, 24, 21, 25, 28, 32, 29, 33} - 1


  arma::mat contact = as<arma::mat>(params["contact"]);

  arma::vec contact_adj = params["contact_adj"];

  arma::vec age_grps = params["age_grps"];
  arma::uvec inf_state_pos = {2, 6, 3, 7, 11, 15, 12, 16, 19, 23, 20, 24, 27, 31, 28, 32};
  int n_age = age_grps.size();


  arma::mat pop(n_age,34);
  std::copy(Y.begin(), Y.end(), pop.begin());

  arma::uvec row_indx = arma::linspace<arma::uvec>(0,(n_age-1),n_age);
  arma::mat inf = pop.submat(row_indx,inf_state_pos);
  arma::vec N = sum(pop,1);

  // double shift = 6;
  double c_season;
  double xi = params["xi"];
  double phi = params["phi"];
  double ns = params["ns"];

  // convert t into month in a year (1-12)
  int tt = mod(t, 12);
  if(tt==0) tt = 12;

  //amp (2^(ns - 1)))*(Cos@(2*Pi*(tt - phi)/12) + 1)^ns + 1 - amp
  // amp -> xi
  // phi -> phi
  c_season = (xi*pow(2,(ns-1)))*pow((cos(2*M_PI*(tt-phi)/12)+1),ns) + 1 - xi;

  // if((tt+shift)>12) {
  //   c_season = xi*cos(0.018181818*M_PI*(tt + shift - 12 - 1));
  // } else {
  //   c_season = xi*cos(0.018181818*M_PI*(tt + shift - 1));
  // }

  arma::vec valpha = as<arma::vec>(wrap(params["alpha"]));

  arma::vec c_infect = as<arma::vec>(wrap(inf*valpha));

  double res = 0.0;


  for(int i = 0; i <= (n_age-1); i++){
    res = res +((N(i)!=0)?contact_adj(j-1)*contact(j-1 , i)/N(i):0.0) * c_season * c_infect(i);
  }

  return(res);

}
