#this code is for model with splines assuming normal distribution####
#model uses the proportion data

library("splines")
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

PB_male_con=pb_plot%>%filter(treatment=="control", sex=="male")
PB_male_ex=pb_plot%>%filter(treatment=="exclosure", sex=="male")
PB_female_con=pb_plot%>%filter(treatment=="control", sex=="female")
PB_female_ex=pb_plot%>%filter(treatment=="exclosure", sex=="female")

X1 <- PB_female_con$month
X2= PB_female_con$lag_ndvi
X3= PB_female_con$lag_ppt
B <- t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = TRUE)) # creating the B-splines, degree=3 for cubic spline
num_data <- length(X1)
num_basis <- nrow(B)
Y = PB_female_con$proportion

mod2<-stan(model_code="
data { 
  int num_data; //rows of observations 
  int num_basis; //no. of basis (order-1) 
  real <lower=0,upper=1>Y[num_data]; //response variable (e.g., no.of breeding obs.)
  vector[num_data] X1; //month
  vector[num_data] X2; //ndvi 
 // vector[num_data] X3; //precip 
  matrix[num_basis, num_data] B; //matrix of coefficients of splines(rows), length of X1 (columns)
} 
 
parameters { 
  row_vector[num_basis] a_raw; // smooth terms for month
  real a0; //intercept
  real ndvi_eff;
  real<lower=0> sigma; //error term 
  real<lower=0> tau; // for noncentered parameterization of spline coefficients
} 
 
transformed parameters { 
  row_vector[num_basis] a; //noncentered parameters of splines
  vector[num_data] Y_hat; 
  
  a = a_raw*tau;  

 Y_hat=a0*X1 + ndvi_eff*X2 + to_vector(a*B);
 }

model { 
  a_raw ~ normal(0, 1); 
  tau ~ normal(0, 1); 
  sigma ~ normal(0, 1); 
  ndvi_eff~normal(0,1);
  Y~ normal(Y_hat, sigma);
  }",
iter=300, control=list(adapt_delta=0.95), 
data =list(X1 =X1, # generating inputs
           X2=X2,
           B =B, # creating the B-splines
           num_data=num_data,
           num_basis=num_basis,
           Y=Y))

#plotting splines####
ff<-extract(sm)
Y_hat_med <- array(NA, length(Y)) #median estimate
Y_hat_ub <- array(NA, length(Y)) #upper boundary
Y_hat_lb <- array(NA, length(Y)) #lower boundary

for (i in 1:length(Y)) {
  Y_hat_med[i] <- median(ff$Y_hat[,i]);
  Y_hat_lb[i] <- quantile(ff$Y_hat[,i],probs = 0.025)
  Y_hat_ub[i] <- quantile(ff$Y_hat[,i],probs = 0.975)
}

plot(X1,Y, xaxt="n") #plot raw data
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12))
lines(smooth.spline(X1, Y_hat_med), col="blue")
lines(smooth.spline(X1, Y_hat_lb), lty=2, col="red") #0.025
lines(smooth.spline(X1, Y_hat_ub), lty=2, col="red") #0.975

yrep=extract(sm)$Y_hat
Y_hat_med <- array(NA, length(Y)) #median estimate
for (i in 1:length(Y)) {
  Y_hat_med[i] <- mean(ff$Y_hat[,i])}

#extract posterior draws####
post1=rstan::extract(sm)$Y_hat
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
points(PB_female_con$proportion~PB_female_con$month, col="blue", pch=16)
lines(smooth.spline(X1, Y_hat_med), col="red")

#model output summary####
print(sm, pars=c("a0", "ndvi_eff"))
