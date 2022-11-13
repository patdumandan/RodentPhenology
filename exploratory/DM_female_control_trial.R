#this code is the spline model with a beta-binomial distribution####
#proportion is a parameter estimated by the model, makes use of the no.of reproductive and total indiv.data

library("splines")
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dm_plot=var_lag2%>%filter(species=="DM")

DM_male_con=dm_plot%>%filter(treatment=="control", sex=="male")
DM_male_ex=dm_plot%>%filter(treatment=="exclosure", sex=="male")
DM_female_con=dm_plot%>%filter(treatment=="control", sex=="female")
DM_female_ex=dm_plot%>%filter(treatment=="exclosure", sex=="female")

X1 =DM_female_con$month
X2= DM_female_con$year
X3= DM_female_con$lag_ndvi
X4= DM_female_con$lag_ppt
num_data = length(X1)
B1 = t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = FALSE)) # creating the B-splines, degree=3 for cubic spline
B2 = t(bs(X2, df=NULL, knots=NULL, degree=3, intercept = FALSE))
num_basis1 = nrow(B1)
num_basis2 = nrow(B2)
Y = DM_female_con$reproductive
n=DM_female_con$abundance

dat_list2=list(X1 =DM_female_con$month,
               X2= DM_female_con$year,
               X3= DM_female_con$lag_ndvi,
               X4= DM_female_con$lag_ppt,
               num_data = length(X1),
               B1 = t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = FALSE)), # creating the B-splines, degree=3 for cubic spline
               B2 = t(bs(X2, df=NULL, knots=NULL, degree=3, intercept = FALSE)),
               num_basis1 = nrow(B1),
               num_basis2 = nrow(B2),
               Y = DM_female_con$reproductive,
               n=DM_female_con$abundance)

mod_dm<-stan(model_code="
data { 
  int num_data; //rows of observations 
  int num_basis1; //no. of basis (order-1) 
  int num_basis2; //no. of basis (order-1)
  
  int Y[num_data]; //response variable (# of reproductive)
  int n[num_data]; //total no.of indivs in plot
   
  vector[num_data] X3; //lag_NDVI
  vector[num_data] X4; //lag_ppt
  
  matrix[num_basis1, num_data] B1; //matrix of coefficients of splines(rows), length of X1 (columns)
  matrix[num_basis2, num_data] B2; //matrix of coefficients of splines(rows), length of X1 (columns)

} 
 
parameters { 
 row_vector[num_basis1] a_raw; // smooth terms for month
  row_vector[num_basis2] b_raw; // smooth terms for year
  real a0; //intercept
  real ndvi_eff;
  real ppt_eff;
  
  real<lower=0> sigma; //error term for shape params
  real<lower=0> tau; // for noncentered parameterization of spline coefficients
  real<lower=0> phi; // for noncentered parameterization of spline coefficients (year)
  
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
  b = b_raw*phi; 
  
  Y_hat = inv_logit(a0 +ndvi_eff*X3+ ppt_eff*X4+ to_vector(a*B1)+ to_vector(b*B2)); 
  

a1=Y_hat*sigma;
b1=(1-Y_hat)*sigma;
}

model { 
  a_raw ~ normal(0, 1); 
  b_raw ~ normal(0, 1);
  tau ~ normal(0, 1); 
  phi ~ normal(0, 1);
  sigma ~ normal(0, 1); 
  ndvi_eff~normal(0,1);
  pred_repro ~ beta(a1, b1);
  Y~ binomial(n, pred_repro);
}",
iter=200,
data =dat_list2)

#plotting regression lines over raw data####
ff<-extract(mod_dm)
Y_hat_med <- array(NA, length(Y)) #median estimate
Y_hat_ub <- array(NA, length(Y)) #upper boundary
Y_hat_lb <- array(NA, length(Y)) #lower boundary

for (i in 1:length(Y)) {
  Y_hat_med[i] <- median(ff$pred_repro[,i]);
  Y_hat_lb[i] <- quantile(ff$pred_repro[,i],probs = 0.025)
  Y_hat_ub[i] <- quantile(ff$pred_repro[,i],probs = 0.975)
}

prop=DM_female_con$proportion
plot(X1,prop, xaxt="n", pch=16, ylab="P(breeding)", xlab="month") #plot raw data
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12))
lines(smooth.spline(X1, Y_hat_med), col="blue")
lines(smooth.spline(X1, Y_hat_ub), lty=2, col="red")
lines(smooth.spline(X1, Y_hat_lb), lty=2, col="red")

yrep=extract(mod_dm)$pred_repro
Y_hat_med <- array(NA, length(Y)) #median estimate
for (i in 1:length(Y)) {
  Y_hat_med[i] <- mean(ff$Y_hat[,i])}

#plot posterior draws####
post1=rstan::extract(mod_dm)$pred_repro
post1=as.data.frame(post1)
post1=t(post1)
t3=cbind(DM_female_con$month, post1)
t3=as.data.frame(t3)
t3=t3%>%
  rename("month"="V1")
t3=reshape2::melt(t3, id=c("month"))

#plot posterior draws####
plot(t3$value~t3$month, type="l", col="grey", ylim=c(0,1), xaxt="n")
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12))
points(DM_female_con$proportion~DM_female_con$month, col="blue", pch=16)
lines(smooth.spline(X1, Y_hat_med), col="red")

#model output summary####
print(mod_dm, pars=c("a0","ndvi_eff", "ppt_eff"))

#sort of the same results when using mgcv package but high uncertainty with the one in stan
dmmod=mgcv::gam(proportion~s(month, bs="cc")+s(year)+lag_ndvi+lag_ppt, 
                weights=abundance, data=DM_female_con, family=binomial(link="logit"))
summary(dmmod)
