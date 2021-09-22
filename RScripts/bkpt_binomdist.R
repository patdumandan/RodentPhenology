library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(dotwhisker)
library(ggpubr)
library(rstan)

options(mc.cores = parallel::detectCores())

portal_reprod=read.csv("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/ReproData/reproductive_full_data.csv")

pb_plot=portal_reprod%>%filter(species=="PB")
PB_male_con=pb_plot%>%filter(treatment=="control", sex=="male")
PB_male_ex=pb_plot%>%filter(treatment=="exclosure", sex=="male")
PB_female_con=pb_plot%>%filter(treatment=="control", sex=="female")
PB_female_ex=pb_plot%>%filter(treatment=="exclosure", sex=="female")


pp_plot=portal_reprod%>%filter(species=="PP")
PP_male_con=pp_plot%>%filter(treatment=="control", sex=="male")
PP_male_ex=pp_plot%>%filter(treatment=="exclosure", sex=="male")
PP_female_con=pp_plot%>%filter(treatment=="control", sex=="female")
PP_female_ex=pp_plot%>%filter(treatment=="exclosure", sex=="female")


datlist=list(N=length(PB_female_con$X),
             abundance=PB_female_con$abundance,
             reproductive=PB_female_con$reproductive,
             month=PB_female_con$month,
             years=PB_female_con$year)

bkpt_con=stan(model_code="
data {
  int<lower=1> N;                      // sample size
  int <lower=0>  abundance [N]; 
  int <lower=0> reproductive [N];
//  vector<lower=0,upper=1>[N] proportion;        // response 
  vector [N] month;                       // predictor month
  vector [N] years; // predictor year
}

parameters{
  real alpha; // intercept
  real beta1;//slope 1
  real beta2;//slope 2
  real <lower=1, upper=12>bkpt;
  real <lower=0> sigma;
  
 // real<lower=0,upper=1> prop_mu [N]; // probability of success in range (0,1)
  real year_eff; //year effect
}

transformed parameters{
  real <lower=0,upper=1> prop_mu [N];

for ( i in 1:N) {
  if (month[i] < bkpt)
    
    prop_mu[i]=inv_logit(alpha+beta1* (month[i]*bkpt)+ year_eff*years[i]);
  
  else
    prop_mu[i]= inv_logit(alpha+ beta2 * (month[i]*bkpt)+ year_eff*years[i]);
  
}
  }

model {
  
  alpha~normal(0,1);
  beta1~ normal(0,1);
  beta2~normal (0,1);
  bkpt~ normal(0,1);
  sigma~cauchy(0,2.5);
  year_eff~normal(0,1);
  
 }

", data=datlist, chains=4, iter=1000)

#visualization####
post=rstan::extract(bkpt_con)$prop_mu
post=as.data.frame(post)
post=t(post)
t1=cbind(PB_female_con$month,PB_female_con$year, post)
t1=as.data.frame(t1)
t1=t1%>%
  rename(year=V2, month=V1)
t2=reshape2::melt(t1, id=c("year", "month"))
plot(t2$value~t2$month, type="l", col="grey", ylim=c(0,1))
points(PB_female_con$proportion~PB_female_con$month, col="blue", pch=16)

print(bkpt_con, pars=c("alpha", "beta1", "beta2", "bkpt"))
abline(v=4, type="l", lty=2)
rect(xleft=1, ybottom=0, ytop=1, xright=6, col = rgb(0.5,0.5,0.5,1/4))
