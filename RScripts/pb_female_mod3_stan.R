#pb female exclosure####

require(bayesplot)
library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(splines)
library(rstan)
library(mgcv)
library(dotwhisker)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

portal_reprod=read.csv("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/ReproData/reproductive_full_data.csv")

pb_plot=portal_reprod%>%filter(species=="PB")
PB_male_con=pb_plot%>%filter(treatment=="control", sex=="male")
PB_male_ex=pb_plot%>%filter(treatment=="exclosure", sex=="male")
PB_female_con=pb_plot%>%filter(treatment=="control", sex=="female")
PB_female_ex=pb_plot%>%filter(treatment=="exclosure", sex=="female")

# PB female exclosure####
X1= PB_female_ex$month
X2= PB_female_ex$year
X3= PB_female_ex$ndvis_lag
X4= PB_female_ex$temps_lag_mean
X5= PB_female_ex$ppts_lag_warm
X6= PB_female_ex$ppts_lag_cool
X7= PB_female_ex$ppts_cool
X8= PB_female_ex$ppts_warm
X9= PB_female_ex$temps_mean
X10= PB_female_ex$ndvis
X11= PB_female_ex$pps
X12= PB_female_ex$pbs
X13= PB_female_ex$dipos
num_data = length(X1)
B1 = t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = FALSE)) 
B2 = t(bs(X2, df=NULL, knots=NULL, degree=3, intercept = FALSE))
num_basis1 = nrow(B1)
num_basis2 = nrow(B2)
Y_pbfex = PB_female_ex$reproductive
n_pbfex=PB_female_ex$abundance

dat_list1=list(X1= PB_female_ex$month,
               X2= PB_female_ex$year,
               X3= PB_female_ex$ndvis_lag,
               X4= PB_female_ex$temps_lag_mean,
               X5= PB_female_ex$ppts_lag_warm,
               X6= PB_female_ex$ppts_lag_cool,
               X7= PB_female_ex$ppts_cool,
               X8= PB_female_ex$ppts_warm,
               X9= PB_female_ex$temps_mean,
               X10= PB_female_ex$ndvis,
               X11= PB_female_ex$pps,
               X12= PB_female_ex$pbs,
               X13= PB_female_ex$dipos,
               num_data = length(X1),
               B1 = t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = FALSE)), 
               B2 = t(bs(X2, df=NULL, knots=NULL, degree=3, intercept = FALSE)),
               num_basis1 = nrow(B1),
               num_basis2 = nrow(B2),
               Y_pbfex = PB_female_ex$reproductive,
               n_pbfex=PB_female_ex$abundance)

pbf_ex1<-stan(model_code="
data { 
  int num_data; //rows of observations 
  int num_basis1; //no. of basis (order-1) 
  int num_basis2; //no. of basis (order-1)
  
  int Y_pbfex[num_data]; //response variable (# of reproductive)
  int n_pbfex[num_data]; //total no.of indivs in plot
   
  vector[num_data] X3; //lag_NDVI
  vector[num_data] X4; //lag_temp
  vector[num_data] X5; //lag warm precip
  vector[num_data] X6;// lag cool precip
  vector[num_data] X7;// cool precip
  vector[num_data] X8;// warm  precip
  vector[num_data] X9;// mean  precip
  vector[num_data] X10;// ndvi precip
  vector[num_data] X11;// PP biomass
  vector[num_data] X12;// PB biomass
  vector[num_data] X13;// DIPO biomass
  
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
  real ndvi_lag_eff;
  real cool_ppt_lag_eff;
  real temp_lag_eff;
  real warm_ppt_lag_eff;
  real pp_eff;
  real pb_eff;
  real dipo_eff;
  
  
  real<lower=0> sigma; //error term for shape params
  real<lower=0> tau; // for noncentered parameterization of spline coefficients (month)
  real<lower=0> phi; // for noncentered parameterization of spline coefficients (year)
 vector <lower=0, upper=1> [num_data] pred_repro;//breeding odds as a parameter
} 
 
transformed parameters { 
  row_vector[num_basis1] a; //noncentered parameters of splines (month)
  row_vector[num_basis2] b; //noncentered parameters of splines (year)
  
  //beta dist of probability as a deterministic function
  vector <lower=0, upper=1> [num_data] Y_hat; // mean of response variable
 
 //beta shape params
  vector <lower=0> [num_data] a1;
  vector <lower=0> [num_data] b1;
  
  a = a_raw*tau;  
  b = b_raw*phi; 
  
  Y_hat = inv_logit(a0 + ndvi_eff*X10 +temp_eff*X9+ warm_ppt_eff*X8+ cool_ppt_eff*X7 + 
  cool_ppt_lag_eff*X6 +warm_ppt_lag_eff*X5+ temp_lag_eff*X4+ ndvi_lag_eff*X3 + 
  pp_eff*X11+ pb_eff*X12+ dipo_eff*X13+
  to_vector(a*B1)+ to_vector(b*B2)); 
  

a1=Y_hat*sigma;
b1=(1-Y_hat)*sigma;
}

model { 
  a0~ normal(0,1);
  a_raw ~ normal(0, 1); 
  b_raw ~ normal(0, 1);
  tau ~ normal (0,1); 
  phi~ normal (0,1);
  sigma ~ normal (0,1); 
  ndvi_eff~normal(0,0.1);
  warm_ppt_eff~normal(0,0.1);
  cool_ppt_eff~normal(0,0.1);
  temp_eff~normal(0,0.1);
  ndvi_lag_eff~normal(0,0.1);
  warm_ppt_lag_eff~normal(0,0.1);
  cool_ppt_lag_eff~normal(0,0.1);
  temp_lag_eff~normal(0,0.1);
  pp_eff~normal(0,0.1);
  pb_eff~normal(0,0.1);
  dipo_eff~normal(0,0.1);
  
  pred_repro ~ beta(a1, b1);
  Y_pbfex~ binomial(n_pbfex, pred_repro);
}

generated quantities{

vector [num_data]log_lik;

for (i in 1:num_data){

log_lik[i] = beta_binomial_lpmf( Y_pbfex[i]| n_pbfex[i], a1[i], b1[i]);
}

}
",
iter=5000,
data =dat_list1, control=list(adapt_delta=0.99, max_treedepth=12))

saveRDS(pbf_ex1, "pbf_ex1_output.RDS")
