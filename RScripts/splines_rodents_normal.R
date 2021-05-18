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
X2= PB_female_con$years
B <- t(bs(X1, knots=seq(-5,5,1), degree=3, intercept = TRUE)) # creating the B-splines
num_data <- length(X1)
num_basis <- nrow(B)
Y = PB_female_con$proportion

sm<-stan(model_code="
data { 
  int num_data; //rows of observations 
  int num_basis; //no. of basis (order-1) 
  real <lower=0,upper=1>Y[num_data]; //response variable (e.g., no.of breeding obs.)
  vector[num_data] X1; //month
 // vector[num_data] X2; //ndvi 
  matrix[num_basis, num_data] B; //matrix of coefficients of splines(rows), length of X1 (columns)
} 
 
parameters { 
  row_vector[num_basis] a_raw; // smooth terms for month
  real a0; //intercept
//  real ndvi_eff;
  real<lower=0> sigma; //error term 
  real<lower=0> tau; // for noncentered parameterization of spline coefficients
} 
 
transformed parameters { 
  row_vector[num_basis] a; //noncentered parameters of splines
  vector <lower=0, upper=1> [num_data] repro_mu; //estimated P(breeding) 
  vector <lower=0> [num_data] a1;
  vector <lower=0> [num_data] b1;
  a = a_raw*tau;  

 for (i in 1:num_data){
 
 repro_mu=inv_logit(a0*X1 + to_vector(a*B));
 }
 
 a1=repro_mu*tau;
 b1=(1-a1)*tau;
}

model { 
  a_raw ~ normal(0, 1); 
  tau ~ cauchy(0, 1); 
  sigma ~ cauchy(0, 1); 
 // ndvi_eff~normal(0,1);


  for (i in 1:num_data){
  
  Y[i]~ beta(a1, b1);
  }

}",
iter=300, control=list(adapt_delta=0.95), 
data =list(X1 =X1, # generating inputs
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
  Y_hat_lb[i] <- quantile(ff$Y_hat[,i],probs = 0.25)
  Y_hat_ub[i] <- quantile(ff$Y_hat[,i],probs = 0.75)
}

plot(X1,Y, xaxt="n") #plot raw data
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12))
lines(smooth.spline(X, Y_hat_med), col="blue")

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
lines(smooth.spline(X, Y_hat_med), col="red")

nd=lag(all_bmass3$ndvi, k=1)
