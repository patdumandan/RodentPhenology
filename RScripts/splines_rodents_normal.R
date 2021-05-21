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

X1 =PB_female_con$month
X2= PB_female_con$year
X3= PB_female_con$lag_ndvi
X4= PB_female_con$lag_ppt
num_data = length(X1)
B1 = t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = FALSE)) # creating the B-splines, degree=3 for cubic spline
B2 = t(bs(X2, df=NULL, knots=NULL, degree=3, intercept = FALSE))
num_basis1 = nrow(B1)
num_basis2 = nrow(B2)
Y = PB_female_con$proportion

dat_list=list(X1 =PB_female_con$month,
              X2= PB_female_con$year,
              X3= PB_female_con$lag_ndvi,
              X4= PB_female_con$lag_ppt,
              num_data = length(X1),
              B1 = t(bs(X1, df=NULL, knots=NULL, degree=3, intercept = FALSE)), # creating the B-splines, degree=3 for cubic spline
              B2 = t(bs(X2, df=NULL, knots=NULL, degree=3, intercept = FALSE)),
              num_basis1 = nrow(B1),
              num_basis2 = nrow(B2),
              Y = PB_female_con$proportion)

mod2<-stan(model_code="
data { 
  int num_data; //rows of observations 
  int num_basis1; //no. of basis (order-1) 
  int num_basis2; //no. of basis (order-1)
  real <lower=0,upper=1>Y[num_data]; //response variable (e.g., no.of breeding obs.)
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
  real ppt_eff;
  real<lower=0> sigma; //error term for shape params 
  real<lower=0> tau; // for noncentered parameterization of spline coefficients (month)
  real<lower=0> phi; // for noncentered parameterization of spline coefficients (year)
  
} 
 
transformed parameters { 
  row_vector[num_basis1] a; //noncentered parameters of splines
  row_vector[num_basis2] b; //noncentered parameters of splines
  
  vector[num_data] Y_hat; 
  
  a = a_raw*tau;  
  b = b_raw*phi;
  
 Y_hat=a0 + ndvi_eff*X3 +ppt_eff*X4 + to_vector(a*B1)+ to_vector(b*B2);
 }

model { 
  a_raw ~ normal(0, 1); 
  b_raw ~ normal(0, 1); 
  phi~normal(0,1);
  tau ~ normal(0, 1); 
  sigma ~ normal(0, 1); 
  ndvi_eff~normal(0,1);
  ppt_eff~normal(0,1);
  
  Y~ normal(Y_hat, sigma);
  }",iter=300,
           data =dat_list, chains=2)

#plotting splines####
ff<-extract(mod2)
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

yrep=extract(mod2)$Y
Y_hat_med <- array(NA, length(Y)) #median estimate
for (i in 1:length(Y)) {
  Y_hat_med[i] <- mean(ff$Y_hat[,i])}

#extract posterior draws####
post1=rstan::extract(mod2)$Y
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

ggplot(PB_female_con, aes(y=proportion, x=month, col=treatment)) +
  geom_point() + 
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PB females in control")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#model output summary####
print(mod2, pars=c("a0", "ndvi_eff"))
