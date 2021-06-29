#FINAL MODEL FOR REPRODUCTIVE PHENOLOGY####

#model specs:
#abiotic predictors: year, NDVI, mean temp, warm and cool precip (lag: 0 and 1)
#horizon for cumulative values of precip: 90 (precip in the last 3 months)
#biotic predictors: biomass of intra- and interspecific competitors (lag: 0) 

require(splines)
require(loo)
require(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#MODEL 1: abiotic variables only####
dat_list1=list(X1= PP_female_ex$month,
               X2= PP_female_ex$year,
               X3= PP_female_ex$ndvis_lag,
               X4= PP_female_ex$temps_lag_mean,
               X5= PP_female_ex$ppts_lag_warm,
               X6= PP_female_ex$ppts_lag_cool,
               X7= PP_female_ex$ppts_cool,
               X8= PP_female_ex$ppts_warm,
               X9= PP_female_ex$temps_mean,
               X10= PP_female_ex$ndvis,
               num_data = length(X1),
               B1 = t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = FALSE)), 
               B2 = t(bs(X2, df=NULL, knots=NULL, degree=3, intercept = FALSE)),
               num_basis1 = nrow(B1),
               num_basis2 = nrow(B2),
               Y = PP_female_ex$reproductive,
               n=PP_female_ex$abundance)

mod1<-stan(model_code="
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
  vector[num_data] X7;// cool precip
  vector[num_data] X8;// warm  precip
  vector[num_data] X9;// mean  precip
  vector[num_data] X10;// ndvi precip
  
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
  to_vector(a*B1)+ to_vector(b*B2)); 
  

a1=Y_hat*sigma;
b1=(1-Y_hat)*sigma;
}

model { 
  a0~ normal(0,1);
  a_raw ~ normal(0, 1); 
  b_raw ~ normal(0, 1);
  tau ~ normal(0, 1); 
  phi~ normal(0,1);
  sigma ~ normal(0, 1); 
  ndvi_eff~normal(0,1);
  warm_ppt_eff~normal(0,1);
  cool_ppt_eff~normal(0,1);
  temp_eff~normal(0,1);
  ndvi_lag_eff~normal(0,1);
  warm_ppt_lag_eff~normal(0,1);
  cool_ppt_lag_eff~normal(0,1);
  temp_lag_eff~normal(0,1);
  pred_repro ~ beta(a1, b1);
  Y~ binomial(n, pred_repro);
}

generated quantities{

vector [num_data]log_lik;

for (i in 1:num_data){

log_lik[i] = beta_binomial_lpmf( Y[i]| n[i], a1[i], b1[i]);
}

}
",
iter=2000,
data =dat_list1, control=list(adapt_delta=0.99, max_treedepth=12))


#MODEL 2: all abiotic variables + intraspecific competition####

X1= PP_female_ex$month
X2= PP_female_ex$year
X3= PP_female_ex$ndvis_lag
X4= PP_female_ex$temps_lag_mean
X5= PP_female_ex$ppts_lag_warm
X6= PP_female_ex$ppts_lag_cool
X7= PP_female_ex$ppts_cool
X8= PP_female_ex$ppts_warm
X9= PP_female_ex$temps_mean
X10= PP_female_ex$ndvis
num_data= length(X1)
B1= t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = FALSE)) 
B2= t(bs(X2, df=NULL, knots=NULL, degree=3, intercept = FALSE))
num_basis1= nrow(B1)
num_basis2= nrow(B2)
Y= PP_female_ex$reproductive
n= PP_female_ex$abundance

dat_list2=list(X1= PP_female_ex$month,
               X2= PP_female_ex$year,
               X3= PP_female_ex$ndvis_lag,
               X4= PP_female_ex$temps_lag_mean,
               X5= PP_female_ex$ppts_lag_warm,
               X6= PP_female_ex$ppts_lag_cool,
               X7= PP_female_ex$ppts_cool,
               X8= PP_female_ex$ppts_warm,
               X9= PP_female_ex$temps_mean,
               X10= PP_female_ex$ndvis,
               X11= PP_female_ex$pps,
               num_data = length(X1),
               B1 = t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = FALSE)), 
               B2 = t(bs(X2, df=NULL, knots=NULL, degree=3, intercept = FALSE)),
               num_basis1 = nrow(B1),
               num_basis2 = nrow(B2),
               Y = PP_female_ex$reproductive,
               n=PP_female_ex$abundance)

mod2<-stan(model_code="
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
  vector[num_data] X7;// cool precip
  vector[num_data] X8;// warm  precip
  vector[num_data] X9;// mean  precip
  vector[num_data] X10;// ndvi precip
  vector[num_data] X11;// PP biomass
  
  
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
  cool_ppt_lag_eff*X6 +warm_ppt_lag_eff*X5+ temp_lag_eff*X4+ ndvi_lag_eff*X3 + pp_eff*X11+
  to_vector(a*B1)+ to_vector(b*B2)); 
  

a1=Y_hat*sigma;
b1=(1-Y_hat)*sigma;
}

model { 
  a0~ normal(0,1);
  a_raw ~ normal(0, 1); 
  b_raw ~ normal(0, 1);
  tau ~ normal(0, 1); 
  phi~ normal(0,1);
  sigma ~ normal(0, 1); 
  ndvi_eff~normal(0,1);
  warm_ppt_eff~normal(0,1);
  cool_ppt_eff~normal(0,1);
  temp_eff~normal(0,1);
  ndvi_lag_eff~normal(0,1);
  warm_ppt_lag_eff~normal(0,1);
  cool_ppt_lag_eff~normal(0,1);
  temp_lag_eff~normal(0,1);
  pp_eff~normal(0,1);
  pred_repro ~ beta(a1, b1);
  Y~ binomial(n, pred_repro);
}

generated quantities{

vector [num_data]log_lik;

for (i in 1:num_data){

log_lik[i] = beta_binomial_lpmf( Y[i]| n[i], a1[i], b1[i]);
}

}
",
iter=2000,
data =dat_list2, control=list(adapt_delta=0.99, max_treedepth=12))

print(mod2, pars=c ("a0", "ndvi_eff", "warm_ppt_eff", "cool_ppt_eff", "temp_eff",
                    "ndvi_lag_eff", "warm_ppt_lag_eff", "cool_ppt_lag_eff", "temp_lag_eff", "pp_eff"))

ppex_f=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+pps+temps_mean+ppts_warm+ppts_cool+
                 ndvis_lag+ temps_lag_mean+ppts_lag_warm+ppts_lag_cool,
                 data=PP_female_ex, weights=abundance, family=binomial(link="logit"))
summary(ppex_f)


#MODEL 3: all abiotic variables + intra and interspecific competition####
dat_list3=list(X1= PP_female_ex$month,
               X2= PP_female_ex$year,
               X3= PP_female_ex$ndvis_lag,
               X4= PP_female_ex$temps_lag_mean,
               X5= PP_female_ex$ppts_lag_warm,
               X6= PP_female_ex$ppts_lag_cool,
               X7= PP_female_ex$ppts_cool,
               X8= PP_female_ex$ppts_warm,
               X9= PP_female_ex$temps_mean,
               X10= PP_female_ex$ndvis,
               X11= PP_female_ex$pps,
               X12= PP_female_ex$pbs,
               X13= PP_female_ex$dipos,
               num_data = length(X1),
               B1 = t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = FALSE)), 
               B2 = t(bs(X2, df=NULL, knots=NULL, degree=3, intercept = FALSE)),
               num_basis1 = nrow(B1),
               num_basis2 = nrow(B2),
               Y = PP_female_ex$reproductive,
               n=PP_female_ex$abundance)

mod3<-stan(model_code="
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
  tau ~ normal(0, 1); 
  phi~ normal(0,1);
  sigma ~ normal(0, 1); 
  ndvi_eff~normal(0,1);
  warm_ppt_eff~normal(0,1);
  cool_ppt_eff~normal(0,1);
  temp_eff~normal(0,1);
  ndvi_lag_eff~normal(0,1);
  warm_ppt_lag_eff~normal(0,1);
  cool_ppt_lag_eff~normal(0,1);
  temp_lag_eff~normal(0,1);
  pp_eff~normal(0,1);
  pB_eff~normal(0,1);
  dipo_eff~normal(0,1);
  
  pred_repro ~ beta(a1, b1);
  Y~ binomial(n, pred_repro);
}

generated quantities{

vector [num_data]log_lik;

for (i in 1:num_data){

log_lik[i] = beta_binomial_lpmf( Y[i]| n[i], a1[i], b1[i]);
}

}
",
iter=2000,
data =dat_list3, control=list(adapt_delta=0.99, max_treedepth=12))

#MODEL VISUALS AND EVALS (PPC and LOOIC)####

ff<-extract(mod2)
Y_hat_med <- array(NA, length(Y)) #median estimate
Y_hat_ub <- array(NA, length(Y)) #upper boundary
Y_hat_lb <- array(NA, length(Y)) #lower boundary

for (i in 1:length(Y)) {
  Y_hat_med[i] <- median(ff$pred_repro[,i]);
  Y_hat_lb[i] <- quantile(ff$pred_repro[,i],probs = 0.025)
  Y_hat_ub[i] <- quantile(ff$pred_repro[,i],probs = 0.975)
}

prop=PP_female_ex$proportion
plot(X1,prop, xaxt="n", ylab="P(breeding)", xlab="month", main="PP female (exclosure)") #plot raw data
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12))
lines(smooth.spline(X1, Y_hat_med), col="blue")
lines(smooth.spline(X1, Y_hat_ub), lty=2, col="red")
lines(smooth.spline(X1, Y_hat_lb), lty=2, col="red")

yrep=extract(mod2)$pred_repro
Y_hat_med <- array(NA, length(Y)) #median estimate
for (i in 1:length(Y)) {
  Y_hat_med[i] <- mean(ff$pred_repro[,i])}

#plot posterior draws####
post2=rstan::extract(mod2)$pred_repro
post2=as.data.frame(post2)
post2=t(post2)
t2=cbind(PP_female_ex$month, post2)
t2=as.data.frame(t2)
t2=t2%>%
  rename("month"="V1")
t2=reshape2::melt(t2, id=c("month"))

#plot posterior draws####
plot(t2$value~t2$month, type="l", col="grey", ylim=c(0,1), xaxt="n", ylab="P(breeding)", xlab="month",main="PP female(exclosure)")
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12))
points(PP_female_ex$proportion~PP_female_ex$month, col="blue", pch=16)
lines(smooth.spline(X1, Y_hat_med), col="red")

#looic calculations####
loglik2=extract_log_lik(mod2)
loo(loglik2)
