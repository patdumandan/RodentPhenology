portal_reprod=read.csv("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/ReproData/reproductive_full_data.csv")

pb_plot=portal_reprod%>%filter(species=="PB")
PB_male_con=pb_plot%>%filter(treatment=="control", sex=="male")
PB_male_ex=pb_plot%>%filter(treatment=="exclosure", sex=="male")
PB_female_con=pb_plot%>%filter(treatment=="control", sex=="female")
PB_female_ex=pb_plot%>%filter(treatment=="exclosure", sex=="female")

#PB male exclosure####
XY1= PB_male_ex$month
XY2= PB_male_ex$year
XY3= PB_male_ex$ndvis_lag
XY4= PB_male_ex$temps_lag_mean
XY5= PB_male_ex$ppts_lag_warm
XY6= PB_male_ex$ppts_lag_cool
XY7= PB_male_ex$ppts_cool
XY8= PB_male_ex$ppts_warm
XY9= PB_male_ex$temps_mean
XY10= PB_male_ex$ndvis
XY11= PB_male_ex$pps
XY12= PB_male_ex$pbs
XY13= PB_male_ex$dipos
num_data_XY = length(XY1)
B1_XY = t(bs(XY1, df=NULL, knots=NULL, degree=3, intercept = FALSE)) 
B2_XY = t(bs(XY2, df=NULL, knots=NULL, degree=3, intercept = FALSE))
num_basis1_XY = nrow(B1_XY)
num_basis2_XY = nrow(B2_XY)
Y_pbmex = PB_male_ex$reproductive
n_pbmex=PB_male_ex$abundance

dat_list3=list(XY1= PB_male_ex$month,
               XY2= PB_male_ex$year,
               XY3= PB_male_ex$ndvis_lag,
               XY4= PB_male_ex$temps_lag_mean,
               XY5= PB_male_ex$ppts_lag_warm,
               XY6= PB_male_ex$ppts_lag_cool,
               XY7= PB_male_ex$ppts_cool,
               XY8= PB_male_ex$ppts_warm,
               XY9= PB_male_ex$temps_mean,
               XY10= PB_male_ex$ndvis,
               XY11= PB_male_ex$pps,
               XY12= PB_male_ex$pbs,
               XY13= PB_male_ex$dipos,
               num_data_XY = length(XY1),
               B1_XY = t(bs(XY1, df=NULL, knots=NULL, degree=3, intercept = FALSE)), 
               B2_XY = t(bs(XY2, df=NULL, knots=NULL, degree=3, intercept = FALSE)),
               num_basis1_XY = nrow(B1_XY),
               num_basis2_XY = nrow(B2_XY),
               Y_pbmex = PB_male_ex$reproductive,
               n_pbmex=PB_male_ex$abundance)

pbm_ex1<-stan(model_code="
data { 
  int num_data_XY; //rows of observations 
  int num_basis1_XY; //no. of basis (order-1) 
  int num_basis2_XY; //no. of basis (order-1)
  
  int Y_pbmex[num_data_XY]; //response variable (# of reproductive)
  int n_pbmex[num_data_XY]; //total no.of indivs in plot
   
  vector[num_data_XY] XY3; //lag_NDVI
  vector[num_data_XY] XY4; //lag_temp
  vector[num_data_XY] XY5; //lag warm precip
  vector[num_data_XY] XY6;// lag cool precip
  vector[num_data_XY] XY7;// cool precip
  vector[num_data_XY] XY8;// warm  precip
  vector[num_data_XY] XY9;// mean  precip
  vector[num_data_XY] XY10;// ndvi precip
  vector[num_data_XY] XY11;// PP biomass
  vector[num_data_XY] XY12;// PB biomass
  vector[num_data_XY] XY13;// DIPO biomass
  
  matrix[num_basis1_XY, num_data_XY] B1_XY; //matrix of coefficients of splines(rows), length of X1 (columns)
  matrix[num_basis2_XY, num_data_XY] B2_XY; //matrix of coefficients of splines(rows), length of X2 (columns)

} 
 
parameters { 
 row_vector[num_basis1_XY] a_raw; // smooth terms for month
  row_vector[num_basis2_XY] b_raw; // smooth terms for year
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
 vector <lower=0, upper=1> [num_data_XY] pred_repro;//breeding odds as a parameter
} 
 
transformed parameters { 
  row_vector[num_basis1_XY] a; //noncentered parameters of splines (month)
  row_vector[num_basis2_XY] b; //noncentered parameters of splines (year)
  
  //beta dist of probability as a deterministic function
  vector <lower=0, upper=1> [num_data_XY] Y_hat; // mean of response variable
 
 //beta shape params
  vector <lower=0> [num_data_XY] a1;
  vector <lower=0> [num_data_XY] b1;
  
  a = a_raw*tau;  
  b = b_raw*phi; 
  
  Y_hat = inv_logit(a0 + ndvi_eff*XY10 +temp_eff*XY9+ warm_ppt_eff*XY8+ cool_ppt_eff*XY7 + 
  cool_ppt_lag_eff*XY6 +warm_ppt_lag_eff*XY5+ temp_lag_eff*XY4+ ndvi_lag_eff*XY3 + 
  pp_eff*XY11+ pb_eff*XY12+ dipo_eff*XY13+
  to_vector(a*B1_XY)+ to_vector(b*B2_XY)); 
  

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
  pb_eff~normal(0,1);
  dipo_eff~normal(0,1);
  
  pred_repro ~ beta(a1, b1);
  Y_pbmex~ binomial(n_pbmex, pred_repro);
}

generated quantities{

vector [num_data_XY]log_lik;

for (i in 1:num_data_XY){

log_lik[i] = beta_binomial_lpmf( Y_pbmex[i]| n_pbmex[i], a1[i], b1[i]);
}

}
",
iter=3000,
data =dat_list3) #, control=list(adapt_delta=0.99, max_treedepth=12))

#PB MALE CONTROL####


XYY1= PB_male_con$month
XYY2= PB_male_con$year
XYY3= PB_male_con$ndvis_lag
XYY4= PB_male_con$temps_lag_mean
XYY5= PB_male_con$ppts_lag_warm
XYY6= PB_male_con$ppts_lag_cool
XYY7= PB_male_con$ppts_cool
XYY8= PB_male_con$ppts_warm
XYY9= PB_male_con$temps_mean
XYY10= PB_male_con$ndvis
XYY11= PB_male_con$pps
XYY12= PB_male_con$pbs
XYY13= PB_male_con$dipos
num_data_XYY = length(XYY1)
B1_XYY = t(bs(XYY1, df=NULL, knots=NULL, degree=3, intercept = FALSE)) 
B2_XYY = t(bs(XYY2, df=NULL, knots=NULL, degree=3, intercept = FALSE))
num_basis1_XYY = nrow(B1_XYY)
num_basis2_XYY = nrow(B2_XYY)
Y_pbmcon = PB_male_con$reproductive
n_pbmcon=PB_male_con$abundance

dat_list4=list(XYY1= PB_male_con$month,
               XYY2= PB_male_con$year,
               XYY3= PB_male_con$ndvis_lag,
               XYY4= PB_male_con$temps_lag_mean,
               XYY5= PB_male_con$ppts_lag_warm,
               XYY6= PB_male_con$ppts_lag_cool,
               XYY7= PB_male_con$ppts_cool,
               XYY8= PB_male_con$ppts_warm,
               XYY9= PB_male_con$temps_mean,
               XYY10= PB_male_con$ndvis,
               XYY11= PB_male_con$pps,
               XYY12= PB_male_con$pbs,
               XYY13= PB_male_con$dipos,
               num_data_XYY = length(XYY1),
               B1_XYY = t(bs(XYY1, df=NULL, knots=NULL, degree=3, intercept = FALSE)), 
               B2_XYY = t(bs(XYY2, df=NULL, knots=NULL, degree=3, intercept = FALSE)),
               num_basis1_XYY = nrow(B1_XYY),
               num_basis2_XYY = nrow(B2_XYY),
               Y_pbmcon = PB_male_con$reproductive,
               n_pbmcon=PB_male_con$abundance)

pbm_con1<-stan(model_code="
data { 
  int num_data_XYY; //rows of observations 
  int num_basis1_XYY; //no. of basis (order-1) 
  int num_basis2_XYY; //no. of basis (order-1)
  
  int Y_pbmcon[num_data_XYY]; //response variable (# of reproductive)
  int n_pbmcon[num_data_XYY]; //total no.of indivs in plot
   
  vector[num_data_XYY] XYY3; //lag_NDVI
  vector[num_data_XYY] XYY4; //lag_temp
  vector[num_data_XYY] XYY5; //lag warm precip
  vector[num_data_XYY] XYY6;// lag cool precip
  vector[num_data_XYY] XYY7;// cool precip
  vector[num_data_XYY] XYY8;// warm  precip
  vector[num_data_XYY] XYY9;// mean  precip
  vector[num_data_XYY] XYY10;// ndvi precip
  vector[num_data_XYY] XYY11;// PP biomass
  vector[num_data_XYY] XYY12;// PB biomass
  vector[num_data_XYY] XYY13;// DIPO biomass
  
  matrix[num_basis1_XYY, num_data_XYY] B1_XYY; //matrix of coefficients of splines(rows), length of X1 (columns)
  matrix[num_basis2_XYY, num_data_XYY] B2_XYY; //matrix of coefficients of splines(rows), length of X2 (columns)

} 
 
parameters { 
 row_vector[num_basis1_XYY] a_raw; // smooth terms for month
  row_vector[num_basis2_XYY] b_raw; // smooth terms for year
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
 vector <lower=0, upper=1> [num_data_XYY] pred_repro;//breeding odds as a parameter
} 
 
transformed parameters { 
  row_vector[num_basis1_XYY] a; //noncentered parameters of splines (month)
  row_vector[num_basis2_XYY] b; //noncentered parameters of splines (year)
  
  //beta dist of probability as a deterministic function
  vector <lower=0, upper=1> [num_data_XYY] Y_hat; // mean of response variable
 
 //beta shape params
  vector <lower=0> [num_data_XYY] a1;
  vector <lower=0> [num_data_XYY] b1;
  
  a = a_raw*tau;  
  b = b_raw*phi; 
  
  Y_hat = inv_logit(a0 + ndvi_eff*XYY10 +temp_eff*XYY9+ warm_ppt_eff*XYY8+ cool_ppt_eff*XYY7 + 
  cool_ppt_lag_eff*XYY6 +warm_ppt_lag_eff*XYY5+ temp_lag_eff*XYY4+ ndvi_lag_eff*XYY3 + 
  pp_eff*XYY11+ pb_eff*XYY12+ dipo_eff*XYY13+
  to_vector(a*B1_XYY)+ to_vector(b*B2_XYY)); 
  

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
  pb_eff~normal(0,1);
  dipo_eff~normal(0,1);
  
  pred_repro ~ beta(a1, b1);
  Y_pbmcon~ binomial(n_pbmcon, pred_repro);
}

generated quantities{

vector [num_data_XYY]log_lik;

for (i in 1:num_data_XYY){

log_lik[i] = beta_binomial_lpmf( Y_pbmcon[i]| n_pbmcon[i], a1[i], b1[i]);
}

}
",
iter=3000,
data =dat_list4) #, control=list(adapt_delta=0.99, max_treedepth=12)) 

# Note: 1 divergent transition, low ESS for both male models

post3=rstan::extract(pbm_con1)
post4=rstan::extract(pbm_ex1)
post3=as.data.frame(post3)
post4=as.data.frame(post4)

p3=post3%>%select(ndvi_eff,ndvi_lag_eff, temp_eff, temp_lag_eff,warm_ppt_eff, warm_ppt_lag_eff, cool_ppt_eff,
                  cool_ppt_lag_eff, pp_eff, pb_eff, dipo_eff)

p4=post4%>%select(ndvi_eff,ndvi_lag_eff, temp_eff, temp_lag_eff,warm_ppt_eff, warm_ppt_lag_eff, cool_ppt_eff,
                  cool_ppt_lag_eff, pp_eff, pb_eff, dipo_eff)

combined_pbm1 <- rbind(mcmc_intervals_data(p3),mcmc_intervals_data(p4))
combined_pbm1$model<- rep(c("Control", "Exclosure"), each = ncol(p3))
theme_set(bayesplot::theme_default())
pos <- position_nudge(y = ifelse(combined_pbm1$model == "Exclosure", 0, 0.1))

ggplot(combined_pbm1, aes(x = m, y = parameter, color = model)) + xlab("coefficient estimate")+
  geom_point(position = pos) +ggtitle("PB males (abiotic+intra+interspecific competition)")#+ geom_vline(xintercept=0, linetype="dotted")+
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos)

