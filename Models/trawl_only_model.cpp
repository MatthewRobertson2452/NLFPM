
#include <TMB.hpp> 
#include <iostream>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace density;
  Type nll = 0.0;
  Type zero = 0.0;
  Type one = 1.0;
  Type two = 2.0;
  
  //input data
  DATA_INTEGER(n);
  DATA_INTEGER(nyrs);
  DATA_IVECTOR(iyear);
  DATA_VECTOR(pa);
  
  PARAMETER_VECTOR(iye);
  PARAMETER(logrw_var);
  
  Type rw_var = exp(logrw_var);
  vector<Type> mu(nyrs);
  
  //Observation Model
  int iy;
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    //id = idex(i);
    
   mu(iy) = one-exp(-iye(iy));
    
   nll -= dbinom(pa(i), one, mu(iy), true);
    
  }
  
  //RW on these likelihoods results in better estimates
  vector<Type> del_iye=log(iye);
  nll -= dnorm(del_iye(0),Type(10.0),rw_var, true);
  for(int i = 1;i < nyrs;++i){
    nll -= dnorm(del_iye(i),del_iye(i-1),rw_var, true);
  }
  
  
  REPORT(mu);
  REPORT(iye);
  
  ADREPORT(iye);
  
  return nll;
}

