library("splines")
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

PB_male_con=pb_plot%>%filter(treatment=="control", sex=="male")
PB_male_ex=pb_plot%>%filter(treatment=="exclosure", sex=="male")
PB_female_con=pb_plot%>%filter(treatment=="control", sex=="female")
PB_female_ex=pb_plot%>%filter(treatment=="exclosure", sex=="female")

dat_list=list(X1= PB_female_con$month,
              X2= PB_female_con$lag_ndvi,
              X3= PB_female_con$years,
              B=  t(bs(X1, df = NULL, knots = NULL, degree=3, intercept = TRUE)), # creating B-splines using splines package
              B3= t(bs(X3, df = NULL, knots = NULL, degree=3, intercept = TRUE)),
              num_data = length(X1),
              num_basis = nrow(B),
              num_basis3 = nrow(B3),
              Y = PB_female_con$reproductive,
              n=PB_female_con$abundance)

smbb<-stan(model_code="
 data { 
   int num_data; //rows of observations 
   int num_basis; //no. of basis (order-1) 
   int num_basis3; //no. of basis (order-1) 
   int Y[num_data]; //response variable (e.g., no.of breeding obs.)
   int n[num_data]; //total no.of indivs in plot
   vector[num_data] X1; //month
   vector[num_data] X2; //ndvi 
   vector[num_data] X3; //standardized year
   matrix[num_basis, num_data] B; //matrix of coefficients of splines, length of X1 
   matrix[num_basis3, num_data] B3;//matrix of coefficients of splines, length of X3 
 } 
  
 parameters { 
   real a0; //intercept
   
   row_vector[num_basis] a_raw; // smooth terms for month
   row_vector[num_basis] b_raw; // smooth terms for year

   real<lower=0> tau; // for noncentered parameterization of spline coefficients a
   real<lower=0> phi; // for noncentered parameterization of spline coefficients b
   
   real ndvi_eff;
   real<lower=0> sigma; //variance for shape params 
   
   vector <lower=0, upper=1> [num_data] pred_repro;
 } 
  
 transformed parameters { 
   
   row_vector[num_basis] a; //coefficients for month spline
   row_vector[num_basis3] b; // coefficients for year spline
   
   vector <lower=0, upper=1> [num_data] Y_hat; // mean of response variable
   vector <lower=0> [num_data] a1;
   vector <lower=0> [num_data] b1;
   
   a = a_raw*tau;
   b = b_raw*phi;
   
   Y_hat = inv_logit(a0+ ndvi_eff*X2+ to_vector(a*B)+ to_vector(b*B3)); 
   
 
 a1=Y_hat*sigma;
 b1=(1-Y_hat)*sigma;
 }
 
 model { 
   a_raw ~ normal(0, 1); 
   b_raw~normal(0,1);
   tau ~ normal(0, 1); 
   sigma ~ normal(0, 1); 
   ndvi_eff~normal(0,1);
   pred_repro ~ beta(a1, b1);
   Y~ binomial(n, pred_repro);
 }",
iter=2000, control=list(max_treedepth=12), 
data =dat_list)

Y=PB_female_con$proportion
#plotting regression lines over raw data####
ff<-extract(smbb)
Y_hat_med <- array(NA, length(Y)) #median estimate
Y_hat_ub <- array(NA, length(Y)) #upper boundary
Y_hat_lb <- array(NA, length(Y)) #lower boundary

for (i in 1:length(Y)) {
   Y_hat_med[i] <- median(ff$Y_hat[,i]);
   Y_hat_lb[i] <- quantile(ff$Y_hat[,i],probs = 0.025)
   Y_hat_ub[i] <- quantile(ff$Y_hat[,i],probs = 0.975)
}

prop=PB_female_con$proportion
plot(X1,prop, xaxt="n") #plot raw data
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12))
lines(smooth.spline(X1, Y_hat_med), col="blue")
lines(smooth.spline(X1, Y_hat_ub), lty=2, col="red")
lines(smooth.spline(X1, Y_hat_lb), lty=2, col="red")

yrep=extract(smbb)$Y_hat
Y_hat_med <- array(NA, length(Y)) #median estimate
for (i in 1:length(Y)) {
   Y_hat_med[i] <- mean(ff$Y_hat[,i])}

#plot posterior draws####
post1=rstan::extract(smbb)$pred_repro
post1=as.data.frame(post1)
post1=t(post1)
t3=cbind(PB_female_con$month, post1)
t3=as.data.frame(t3)
t3=t3%>%
   rename("month"="V1")
t3=reshape2::melt(t3, id=c("month"))

#plot posterior draws####
plot(t3$value~t3$month, type="l", col="grey", ylim=c(0,1), xaxt="n")
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12))
points(PB_female_ex$proportion~PB_female_ex$month, col="blue", pch=16)
lines(smooth.spline(X1, Y_hat_med), col="red")

#model output summary####
print(smbb, pars=c("a0","ndvi_eff"))


