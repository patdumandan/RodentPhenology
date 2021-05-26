//BB spline model for Portal rodents 

data { 
  int num_data; //rows of observations 
  int num_basis1; //no. of basis (order-1) 
  int num_basis2; //no. of basis (order-1)
  
  int Y[num_data]; //response variable (# of reproductive)
  int n[num_data]; //total no.of indivs in plot
   
  vector[num_data] X3; //lag_NDVI
  vector[num_data] X4; //lag_temp
  vector[num_data] X5; //lag warm precip
  vector[num_data] X6;// lag cool precip
  
  matrix[num_basis1, num_data] B1; //matrix of coefficients of splines(rows), length of X1 (columns)
  matrix[num_basis2, num_data] B2; //matrix of coefficients of splines(rows), length of X2 (columns)

} 
 
parameters { 
 row_vector[num_basis1] a_raw; // smooth terms for month
  row_vector[num_basis2] b_raw; // smooth terms for year
  real a0; //intercept
  real ndvi_eff;
  real cool_ppt_eff;
  real temp_eff;
  real warm_ppt_eff;
  
  real<lower=0> sigma; //error term for shape params
  real<lower=0> tau; // for noncentered parameterization of spline coefficients (month)
 // real<lower=0> phi; // for noncentered parameterization of spline coefficients (year)
  
 vector <lower=0, upper=1> [num_data] pred_repro;//breeding odds as a parameter
} 
 
transformed parameters { 
  row_vector[num_basis1] a; //noncentered parameters of splines
  row_vector[num_basis2] b; //noncentered parameters of splines
  
  //beta dist of probability as a deterministic function
  vector <lower=0, upper=1> [num_data] Y_hat; // mean of response variable
 
 //beta shape params
  vector <lower=0> [num_data] a1;
  vector <lower=0> [num_data] b1;
  
  a = a_raw*tau;  
  b = b_raw*tau; 
  
  Y_hat = inv_logit(a0 + ndvi_eff*X3 +temp_eff*X4+ warm_ppt_eff*X5+ cool_ppt_eff*X6 + to_vector(a*B1)+ to_vector(b*B2)); 
  

a1=Y_hat*sigma;
b1=(1-Y_hat)*sigma;
}

model { 
  a0~ normal(0,1);
  a_raw ~ normal(0, 1); 
  b_raw ~ normal(0, 1);
  tau ~ normal(0, 1); 
//  phi ~ normal(0, 1);
  sigma ~ normal(0, 1); 
  ndvi_eff~normal(0,1);
  warm_ppt_eff~normal(0,1);
  cool_ppt_eff~normal(0,1);
  temp_eff~normal(0,1);
  pred_repro ~ beta(a1, b1);
  Y~ binomial(n, pred_repro);
}