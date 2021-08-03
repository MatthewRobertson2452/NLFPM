
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
  DATA_INTEGER(ndex);
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(idex);
  DATA_SCALAR(k);
  DATA_VECTOR(pa);
  
  PARAMETER_VECTOR(iye);
  PARAMETER(lbeta);
  PARAMETER(lchi);
  PARAMETER(logrw_var);
  
  Type rw_var = exp(logrw_var);
  Type beta=exp(lbeta);
  Type chi=exp(lchi);
  vector<Type> mu(nyrs);
  matrix<Type> new_mu(nyrs,ndex);
  
  //Observation Model
  int id, iy;
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    id = idex(i);
    
    mu(iy) = one-exp(-iye(iy));
    
    if(id==0){new_mu(iy,0)=mu(iy);}
    if(id>0){new_mu(iy,id)=(k*pow(mu(iy),beta))/(pow(chi,beta)+pow(mu(iy),beta));}
    
    if(isNA(pa(i))==false){
      if(id==0){nll -= dbinom(pa(i), one, new_mu(iy,0), true);}
      if(id>0){nll -= dbinom(pa(i), one, new_mu(iy,id), true);}
    }
    
  }
  
  //RW on these likelihoods results in better estimates
  vector<Type> del_iye=log(iye);
  nll -= dnorm(del_iye(0),Type(10.0),rw_var, true);
  for(int i = 1;i < nyrs;++i){
    nll -= dnorm(del_iye(i),del_iye(i-1),rw_var, true);
  }
  
  
  REPORT(mu);
  REPORT(iye);
  REPORT(beta);
  REPORT(chi);
  REPORT(new_mu);
  
  ADREPORT(iye);
  
  return nll;
}

