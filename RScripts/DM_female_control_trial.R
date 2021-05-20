dmbb<-stan(model_code="
data { 
  int num_data; //rows of observations 
  int num_basis; //no. of basis (order-1) 
  int Y[num_data]; //response variable (e.g., no.of breeding obs.)
  int n[num_data]; //total no.of indivs in plot
  vector[num_data] X1; //month
  vector[num_data] X2; //ndvi 
  matrix[num_basis, num_data] B; //matrix of coefficients of splines(rows), length of X1 (columns)
} 

parameters{

 row_vector[num_basis] a_raw; // smooth terms for month
 real a0; //intercept
 real ndvi_eff;
 real<lower=0> sigma; //variance of likelihood
 real<lower=0> tau; // for noncentered parameterization of spline coefficients??
 real <lower=0, upper=1> reprod_pred; //breeding odds
}

transformed parameters{
 row_vector[num_basis] a; //coefficients of splines
 vector <lower=0, upper=1> [num_data] reprod_mu;
 vector <lower=0> [num_data] a1;
 vector <lower=0> [num_data] b1;
 
 a = a_raw*tau;  // NCP of spline coefs??
 
 for (i in 1:num_data){
 
 reprod_mu[i] =inv_logit(a0*X1 +ndvi_eff*X2[i]+ to_vector(a*B));
 }
 
 a1=reprod_mu * sigma;
 b1=(1-reprod_mu) *sigma;
}

model { 
  a_raw ~ normal(0, 1); 
  tau ~ normal(0, 1); 
  sigma ~ normal(0, 1); 
  ndvi_eff~normal(0,1);
  Y_hat ~ beta(a1, b1);
  Y~ binomial(n, Y_hat);
}",
iter=200, control=list(adapt_delta=0.95), 
data =list(X1=X1, # generating inputs
           X2=X2,
           n=DM_female_con$abundance,
           B=B, # creating the B-splines
           num_data=num_data,
           num_basis=num_basis,
           Y=Y))
saveRDS(dmbb, "trial_DM.RDS")

Y=DM_female_con$proportion

#plotting regression lines over raw data####
ff<-extract(dmbb)
Y_hat_med <- array(NA, length(Y)) #median estimate
Y_hat_ub <- array(NA, length(Y)) #upper boundary
Y_hat_lb <- array(NA, length(Y)) #lower boundary

for (i in 1:length(Y)) {
  Y_hat_med[i] <- median(ff$Y_hat[,i]);
  Y_hat_lb[i] <- quantile(ff$Y_hat[,i],probs = 0.025)
  Y_hat_ub[i] <- quantile(ff$Y_hat[,i],probs = 0.975)
}

prop=DM_female_con$proportion
plot(X1,prop, xaxt="n") #plot raw data
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12))
lines(smooth.spline(X1, Y_hat_med), col="blue")
lines(smooth.spline(X1, Y_hat_ub), lty=2, col="red")
lines(smooth.spline(X1, Y_hat_lb), lty=2, col="red")

yrep=extract(dmbb)$Y_hat
Y_hat_med <- array(NA, length(Y)) #median estimate
for (i in 1:length(Y)) {
  Y_hat_med[i] <- mean(ff$Y_hat[,i])}

#plot posterior draws####
post1=rstan::extract(dmbb)$Y_hat
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
print(dmbb, pars=c("a0", "ndvi_eff"))

ggplot(DM_female_con, aes(y=proportion, x=month, col=treatment)) +
  geom_point() + 
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("DM females (control")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
