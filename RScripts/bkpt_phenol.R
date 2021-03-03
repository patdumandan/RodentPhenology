require(rstan)
options(mc.cores = parallel::detectCores())

bkpt_con=stan(model_code="
data {
  int<lower=1> N;                      // sample size
  vector<lower=0,upper=1>[N] proportion;        // response 
  vector [N] month;                       // predictor 1
  vector [N] year; // predictor 2
}

parameters{
  real alpha;
  real beta1;
  real beta2;
  real <lower=1, upper=12>bkpt;
  real <lower=0> sigma;
  real year_eff;
}

transformed parameters{
  real prop_mu [N];

for ( i in 1:N) {
  if (month[i] < bkpt)
    
    prop_mu[i]=alpha+beta1* (month[i]*bkpt)+ year_eff*year[i];
  
  else
    prop_mu[i]= alpha+ beta2 * (month[i]*bkpt+ year_eff*year[i]);
  
}
  }

model {
  
  alpha~normal(0,10);
  beta1~ normal(0,10);
  beta2~normal (0,10);
  bkpt~ uniform(1,12);
  sigma~cauchy(0,2.5);
  year_eff~normal(0,10);
  proportion~normal(prop_mu, sigma); }

generated quantities {
real prop_mu_pred [N];

for (j in 1:N)
prop_mu_pred[j]=normal_rng(prop_mu[j],sigma);

}


", data=datlist, chains=4, iter=1000)

datlist=list(N=length(PB_gam_con$X),reprodn=PB_gam_con$reproductive, 
             year=PB_gam_con$years,
             month=PB_gam_con$month, proportion=PB_gam_con$proportion)

saveRDS(bkpt_con, "bkpt_con_v1.RDS")
bkpt_con=readRDS("bkpt_con_v1.RDS")
print(bkpt_con, pars=c("alpha", "beta1", "beta2", "bkpt", "year_eff"))

post=rstan::extract(bkpt_con)$prop_mu_pred
post=as.data.frame(post)
post=t(post)
t1=cbind(PB_gam_con$month,PB_gam_con$year, post)
t1=as.data.frame(t1)
t1=t1%>%
  rename(year=V2, month=V1)
t2=reshape2::melt(t1, id=c("year", "month"))
plot(t2$value~t2$month, type="l", col="grey", ylim=c(0,1))
points(PB_gam_con$proportion~PB_gam_con$month, col="blue", pch=16)

abline(v=2.85, col="black", lty="dashed")
rect(1,0,5,1, col = rgb(0.5,0.5,0.5,1/4))

#exclosure####
bkpt_ex=stan(model_code="
data {
  int<lower=1> N;                      // sample size
  vector<lower=0,upper=1>[N] proportion;        // response 
  vector [N] month;                       // predictor
  vector [N] year;
}

parameters{
  real alpha;
  real beta1;
  real beta2;
  real year_eff;
  real <lower=1, upper=12>bkpt;
  real <lower=0> sigma;
  
}

transformed parameters{
  real prop_mu [N];

for ( i in 1:N) {
  if (month[i] < bkpt){
    
    prop_mu[i]=alpha+beta1* (month[i]*bkpt)+year_eff*year[i];
  }
  else{
    prop_mu[i]= alpha+ beta2 * (month[i]*bkpt)+year_eff*year[i];
  }
}
  }

model {
  
  alpha~normal(0,10);
  beta1~ normal(0,10);
  beta2~normal (0,10);
  bkpt~ uniform(1,12);
  sigma~cauchy(0,2.5);
  year_eff~normal(0,10);
  
  proportion~normal(prop_mu, sigma); }
generated quantities {
real prop_mu_pred [N];

for (j in 1:N)
prop_mu_pred[j]=normal_rng(prop_mu[j],sigma);

}


", data=datlist, chains=4, iter=1000)

datlist=list(N=length(PB_gam_ex$X),reprodn=PB_gam_ex$reproductive, month=PB_gam_ex$month, 
             proportion=PB_gam_ex$proportion, year=PB_gam_ex$years)

saveRDS(bkpt_ex, "bkpt_ex_v1.RDS")
re
print(bkpt_ex, pars=c("alpha", "beta1", "beta2", "bkpt", "year_eff"))

post=rstan::extract(bkpt_ex)$prop_mu_pred
post=as.data.frame(post)
post=t(post)
t1=cbind(PB_gam_ex$month,PB_gam_ex$year, post)
t1=as.data.frame(t1)
t1=t1%>%
  rename(year=V2, month=V1)
t2=reshape2::melt(t1, id=c("year", "month"))
plot(t2$value~t2$month, type="l", col="grey", ylim=c(0,1))
points(PB_gam_ex$proportion~PB_gam_ex$month, col="blue", pch=16)

abline(v=1.4, col="black", lty="dashed")
rect(1,0,2,1, col = rgb(0.5,0.5,0.5,1/4))
