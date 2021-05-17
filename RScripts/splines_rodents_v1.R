library("splines")
library("rstan")

X <- PB_female_con$month
B <- t(bs(X, knots=seq(-5,5,1), degree=3, intercept = TRUE)) # creating the B-splines
num_data <- length(X)
num_basis <- nrow(B)
Y = PB_female_con$proportion

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

sm<-stan(model_code="
data { 
  int num_data; //rows of observations 
  int num_basis; //no. of basis (order-1) 
  vector[num_data] Y; //response variable (e.g., proportion)
  vector[num_data] X; //predictor (e.g., month)
  matrix[num_basis, num_data] B; //matrix of coefficients of splines(rows), length of X (columns)
} 
 
parameters { 
  row_vector[num_basis] a_raw; // spline coefficients
  real a0; //intercept
  real<lower=0> sigma; //error term 
  real<lower=0> tau; // for noncentered parameterization of spline coefficients
} 
 
transformed parameters { 
  row_vector[num_basis] a; //noncentered parameters of splines
  vector[num_data] Y_hat; // mean of response variable
  
  a = a_raw*tau;  
  Y_hat = a0*X + to_vector(a*B); 
} 
 
model { 
  a_raw ~ normal(0, 1); 
  tau ~ cauchy(0, 1); 
  sigma ~ cauchy(0, 1); 
  Y ~ normal(Y_hat, sigma); 
}",
         iter=300, control=list(adapt_delta=0.95), 
         data =list(X =X, # generating inputs
                    B =B, # creating the B-splines
                    num_data=num_data,
                    num_knots=num_knots,
                    spline_degree=spline_degree,
                    num_basis=num_basis,
                    a0 =a0, # intercept
                    a =a, # coefficients of B-splines
                    Y_true =Y_true, # generating the output
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

plot(X,Y) #plot raw data
polygon(c(rev(X), X), c(rev(Y_hat_lb), Y_hat_ub), col = 'grey80', border = NA) #plot uncertainty in estimates
lines(X[order(X)], Y_hat_med[order(Y_hat_med)], col="Red", lw=2) #fitted curve
lines(X, Y, col="blue",lw=2) #true curve

yrep=extract(sm)$Y_hat

post1=rstan::extract(sm)$Y_hat
post1=as.data.frame(post1)
post1=t(post1)
t3=cbind(PB_female_con$month, post1)
t3=as.data.frame(t3)
t3=t3%>%
  rename(year=V2, month=V1)
t3=reshape2::melt(t3, id=c("month"))
plot(t3$value~t3$month, type="l", col="grey", ylim=c(0,1))
points(PB_female_con$proportion~PB_female_con$month, col="blue", pch=16)
mean_eff=apply(t3,1,mean)
lines()