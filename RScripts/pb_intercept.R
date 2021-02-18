require(portalr)
require(dplyr)
require(ggplot2)
require(rstan)
require(brms)
require(rstanarm)
require(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#Data Curation####
####load cleaned individual-level data
Portal_data_indiv=summarize_individual_rodents(
  path = get_default_data_path(),
  clean = TRUE,
  type = "Rodents",
  length = "all",
  unknowns = FALSE,
  time = "date",
  fillweight = FALSE,
  min_plots = 1,
  min_traps = 1,
  download_if_missing = TRUE,
  quiet = FALSE
)%>%filter(!(treatment=="removal")& !is.na(treatment)& !is.na(sex))%>%
  mutate(month=as.character(month))

# add columns for reproductive traits
Portal_data=load_rodent_data(clean = TRUE)
Portal_rodent=Portal_data[["rodent_data"]]%>%mutate(month=factor(month))%>%
  filter(!is.na(sex))

#create full dataset without removal plots
portal1=left_join(Portal_data_indiv, Portal_rodent)%>%
  filter(!(treatment=="removal")& !is.na(treatment)& !is.na(sex))%>%
  mutate(month=as.integer(month), Month=recode(month, "1"= "Jan", "2"="Feb", "3"="Mar", "4"="Apr",
                                                 "5"="May","6"="Jun", "7"="Jul", "8"="Aug", "9"="Sept",
                                                 "10"="Oct","11"="Nov", "12"="Dec"))%>%
  select(period, month, Month, day, year, plot, stake,
         treatment, species, sex, reprod, vagina, nipples,lactation, pregnant, testes,hfl,wgt, tag)

str(portal1)


### Data Manipulation ####


portal_male=portal1%>%filter(sex=="M") #49% of individuals are males
head(portal_male)

####determine threshold for breeding adult male individuals

target_repro=c("S", "M", "R")
repro_male=portal_male%>%
  filter(testes==c("S", "M", "R"))
head(repro_male)

#BA=repro_male%>%
#  filter(species=="BA")%>%
#  arrange(wgt)

#size thresholds
BA=repro_male%>%
  filter(species=="BA", wgt >=6)

DM=repro_male%>%
  filter(species=="DM", wgt >=15)

DO=repro_male%>%
  filter(species=="DO", wgt >=29)

DS=repro_male%>%
  filter(species=="DS", wgt >=12)

NEA=repro_male%>%
  filter(species=="NA", wgt >=121)

OL=repro_male%>%
  filter(species=="OL", wgt >=19)

OT=repro_male%>%
  filter(species=="OT", wgt >=10)

PH=repro_male%>%
  filter(species=="PH", wgt >=18)

PL=repro_male%>%
  filter(species=="PL", wgt >=20)

PB=repro_male%>%
  filter(species=="PB", wgt >=16)

PP=repro_male%>%
  filter(species=="PP", wgt >=10)

PE=repro_male%>%
  filter(species=="PE", wgt >=7)

PI=repro_male%>%
  filter(species=="PI", wgt >=15)

PF=repro_male%>%
  filter(species=="PF", wgt >=4)

PM=repro_male%>%
  filter(species=="PM", wgt >=11)

RM=repro_male%>%
  filter(species=="RM", wgt >=4)

RF=repro_male%>%
  filter(species=="RF", wgt >=11)

RO=repro_male%>%
  filter(species=="RO", wgt >=6)

SF=repro_male%>%
  filter(species=="SF", wgt >=39)

SH=repro_male%>%
  filter(species=="SH", wgt >=51)

SO=repro_male%>%
  filter(species=="SO", wgt >=68)

#create full dataset with all species

full_repro_male_dat=rbind(SO,SH,SF,RO,RM,RF,PP,PM,PL,PH,PF,PE,PB,OT,OL,NEA,DS,DO,DM,BA, PI)
full_repro_male_dat=as.data.frame(full_repro_male_dat)
head(full_repro_male_dat)

####calculate proportion of reproductive individuals

#get count of reproductive males for each species per month per year per trt
repro_dat=full_repro_male_dat%>%
  group_by(month, year, treatment, species)%>%
  summarise(reproductive=n())

#get total observed abundance for each species per month per year per trt
total_rodents=portal_male%>%
  group_by(month,year, treatment, species)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion=right_join(repro_dat, total_rodents)%>%
  mutate(proportion=reproductive/abundance, month=as.integer(month))%>%
  arrange(proportion)
head(total_proportion)

#length(unique(total_proportion$species)) #21 spp
#max(total_proportion$proportion, na.rm=T) #1

#summary stats####
#total indivs across all years: 58345
#no. of males: 28755
#% of males: 49%
#no of reproductive males: 4971
# % of males in repro.mode across all years: 18%

####Testing for Autocorrelation####

y=PB_dat_M$proportion
mod_lm=lm(proportion~year, data=PB_dat_M)
yrep=mod_lm$residuals
require(lmtest)
dwtest(yrep~PB_dat_M$year) #autocorrelation exists (DW=0.06, p < 0.05)

### Data Visualization ####

#sample using PB data only
PB_dat_M=total_proportion%>%filter(species=="PB", !(treatment=="spectabs"))
plot(PB_dat_M$proportion~PB_dat_M$year)
####plot PB proportion across all years per trt type
PB_dat_M%>%
  ggplot(mapping = aes(x=year, y=proportion, colour=treatment))+
  geom_point()+geom_smooth(se = FALSE, method = 'loess')

####plot PB proportion across all years per trt type per month
PB_dat_M%>%
  ggplot(mapping = aes(x=year, y=proportion, colour=treatment))+
  geom_point()+geom_smooth(se = FALSE, method = 'lm')+
  facet_wrap(~month)
PB_dat_M=as.data.frame(PB_dat_M)
####plot PB proportion across all years per month per trt type
#control
PB_dat_M_con=total_proportion%>%filter(species=="PB", treatment=="control")
PB_dat_M_con%>%
  ggplot(mapping = aes(x=year, y=proportion))+
  geom_point()+geom_smooth(se = FALSE, method = 'lm')+
  facet_wrap(~month)+ggtitle("PB male control")

#exclosure
PB_dat_M_ex=total_proportion%>%filter(species=="PB", treatment=="exclosure")
PB_dat_M_ex%>%
  ggplot(mapping = aes(x=year, y=proportion))+
  geom_point()+geom_smooth(se = FALSE, method = 'lm')+
  facet_wrap(~month)+ggtitle("PB male exclosure")

PB_dat_M[is.na(PB_dat_M)] <- 0 #set non-detects to 0
PB_dat_M$trt<-ifelse(PB_dat_M$treatment=="control", 0, 1) 
PB_dat_M$years=(PB_dat_M$year-mean(PB_dat_M$year))/(2*sd(PB_dat_M$year)) #standardize year
PBprop=PB_dat_M$proportion
PBrep=PB_dat_M$reproductive
PB_dat_M$mon_cos= cos(2*pi*(PB_dat_M$month/12))
PB_dat_M$mon_sin= sin(2*pi*(PB_dat_M$month/12))


dat_list=list(
  N=length(PB_dat_M$month),
  y=PB_dat_M$reproductive,
  n=PB_dat_M$abundance,
  year=PB_dat_M$years,
  treatment=PB_dat_M$trt,
  mon_cos=PB_dat_M$mon_cos, 
  mon_sin=PB_dat_M$mon_sin,
  Nmon=length(unique(PB_dat_M$month)),
  Nsp=length(unique(PB_dat_M$species)))


#Beta-binomial model variations####
#intercept only model####
mod_int5=stan(model_code="
   data{
  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // reproductive indivs
  int <lower=0>  n[N];       // total males
 // vector [N] year;// year
//  vector[N] treatment;// treatment
 // int month[N]; //ID of each month
// int Nmon; //no.of months
 // int species[N]; //species ID
// int Nsp; //no.of species
 }
                
 parameters {
  real alpha;// intercept
  //real year_eff; //slope year
 // real trt_eff; //slope treatment effect
  //real<lower=0> sigma_mon[Nmon];//error for random intercept (month)
  //real <lower=0> mon_non;//non-centered error term for species
  //real<lower=0> sigma_sp[Nsp];//error for random intercept (species)
  //real <lower=0> sp_non;//non-centered error term for month
  real <lower=0> phi;
  real <lower=0, upper=1> pred_repro[N] ;//proportion of reproductive event 
              }
   
  transformed parameters{
  vector <lower=0, upper=1> [N] repro_mu; //so we can add statement describing proportion (not able to do in parameters block)
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  //vector [Nmon] alpha_mon; //random intercept per species
   //vector [Nsp] alpha_sp; //random intercept per species
  //vector [Nmon] yr_mon; //random slope per month for year effect
  //vector [Nmon] trt_mon;//random slope per month for treatment effect

  
  //for (j in 1:Nmon) {
  
  //alpha_mon[j]= mon_non*sigma_mon[j];}
  
 //  for (k in 1:Nsp) {
  
  //alpha_sp[k]= sp_non*sigma_sp[k];}
  
  //model:
  
  for (i in 1:N){
  
  repro_mu[i]= inv_logit(alpha);
  }
  
  A = repro_mu * phi;
  B = (1 - repro_mu)* phi;
  
  }
 model {
  //priors
  alpha~normal(0,1);
 // year_eff~ normal (0,1);
  //trt_eff~ normal (0,1);
  //mon_non~ normal(0,1);
  //sigma_mon~ normal(0,1);
  phi ~normal(0,1);
  //sp_non~ normal(0,1);
  //sigma_sp~ normal(0,1);
  
  //model likelihood:
  
  pred_repro ~ beta(A, B); // survival estimate, beta dist.
  y~binomial(n, pred_repro); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate
 
 }
  
  generated quantities {
  
  real pred_y [N];//predictions on proportions
  real log_lik [N];// for looic calculations
  
    pred_y = beta_rng(A, B);
    
    for (x in 1:N){
    log_lik[x]= beta_lpdf(pred_repro[x]| A[x], B[x]);}
   
  }    ", data=dat_list, chains=4, iter=3000)
saveRDS(mod_int5, "PB_intercept.RDS")
m1_loo=extract(mod_int5)$log_lik
loo(m1_loo) #-559.9 +/- 26.2
print(mod_int6, pars=c("alpha", "phi"))

#year only####
mod_int9=stan(model_code="
   data{
  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // reproductive indivs
  int <lower=0>  n[N];       // total males
  vector [N] year;// year
  //vector[N] treatment;// treatment
 // int month[N]; //ID of each month
// int Nmon; //no.of months
 // int species[N]; //species ID
// int Nsp; //no.of species
 }
                
 parameters {
  real alpha;// intercept
  real year_eff; //slope year
  real trt_eff; //slope treatment effect
  //real<lower=0> sigma_mon[Nmon];//error for random intercept (month)
  //real <lower=0> mon_non;//non-centered error term for species
  //real<lower=0> sigma_sp[Nsp];//error for random intercept (species)
  //real <lower=0> sp_non;//non-centered error term for month
  real <lower=0> phi;
  real <lower=0, upper=1> pred_repro[N] ;//proportion of reproductive event 
              }
   
  transformed parameters{
  vector <lower=0, upper=1> [N] repro_mu; //so we can add statement describing proportion (not able to do in parameters block)
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  //vector [Nmon] alpha_mon; //random intercept per species
   //vector [Nsp] alpha_sp; //random intercept per species
  //vector [Nmon] yr_mon; //random slope per month for year effect
  //vector [Nmon] trt_mon;//random slope per month for treatment effect

  
  //for (j in 1:Nmon) {
  
  //alpha_mon[j]= mon_non*sigma_mon[j];}
  
 //  for (k in 1:Nsp) {
  
  //alpha_sp[k]= sp_non*sigma_sp[k];}
  
  //model:
  
  for (i in 1:N){
  
  repro_mu[i]= inv_logit(alpha+ year_eff*year[i]);
  }
  
  A = repro_mu * phi;
  B = (1 - repro_mu)* phi;
  
  }
 model {
  //priors
  alpha~normal(0,1);
  year_eff~ normal (0,1);
  trt_eff~ normal (0,1);
  //mon_non~ normal(0,1);
  //sigma_mon~ normal(0,1);
  phi ~normal(0,1);
  //sp_non~ normal(0,1);
  //sigma_sp~ normal(0,1);
  
  //model likelihood:
  
  pred_repro ~ beta(A, B); // survival estimate, beta dist.
  y~binomial(n, pred_repro); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate
 
 }
  
  generated quantities {
  
  real pred_y [N];//predictions on proportions
  real log_lik [N];// for looic calculations
  
    pred_y = beta_rng(A, B);
    
    for (x in 1:N){
    log_lik[x]= beta_lpdf(pred_repro[x]| A[x], B[x]);}
   
  }    ", data=dat_list, chains=2, iter=300)

#year and treatment effect model####
mod_int6=stan(model_code="
   data{
  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // reproductive indivs
  int <lower=0>  n[N];       // total males
  vector [N] year;// year
  vector[N] treatment;// treatment
 // int month[N]; //ID of each month
// int Nmon; //no.of months
 // int species[N]; //species ID
// int Nsp; //no.of species
 }
                
 parameters {
  real alpha;// intercept
  real year_eff; //slope year
  real trt_eff; //slope treatment effect
  //real<lower=0> sigma_mon[Nmon];//error for random intercept (month)
  //real <lower=0> mon_non;//non-centered error term for species
  //real<lower=0> sigma_sp[Nsp];//error for random intercept (species)
  //real <lower=0> sp_non;//non-centered error term for month
  real <lower=0> phi;
  real <lower=0, upper=1> pred_repro[N] ;//proportion of reproductive event 
              }
   
  transformed parameters{
  vector <lower=0, upper=1> [N] repro_mu; //so we can add statement describing proportion (not able to do in parameters block)
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  //vector [Nmon] alpha_mon; //random intercept per species
   //vector [Nsp] alpha_sp; //random intercept per species
  //vector [Nmon] yr_mon; //random slope per month for year effect
  //vector [Nmon] trt_mon;//random slope per month for treatment effect

  
  //for (j in 1:Nmon) {
  
  //alpha_mon[j]= mon_non*sigma_mon[j];}
  
 //  for (k in 1:Nsp) {
  
  //alpha_sp[k]= sp_non*sigma_sp[k];}
  
  //model:
  
  for (i in 1:N){
  
  repro_mu[i]= inv_logit(alpha+ year_eff*year[i]+trt_eff*treatment[i]);
  }
  
  A = repro_mu * phi;
  B = (1 - repro_mu)* phi;
  
  }
 model {
  //priors
  alpha~normal(0,1);
  year_eff~ normal (0,1);
  trt_eff~ normal (0,1);
  //mon_non~ normal(0,1);
  //sigma_mon~ normal(0,1);
  phi ~normal(0,1);
  //sp_non~ normal(0,1);
  //sigma_sp~ normal(0,1);
  
  //model likelihood:
  
  pred_repro ~ beta(A, B); // survival estimate, beta dist.
  y~binomial(n, pred_repro); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate
 
 }
  
  generated quantities {
  
  real pred_y [N];//predictions on proportions
  real log_lik [N];// for looic calculations
  
    pred_y = beta_rng(A, B);
    
    for (x in 1:N){
    log_lik[x]= beta_lpdf(pred_repro[x]| A[x], B[x]);}
   
  }    ", data=dat_list, chains=4, iter=3000)
saveRDS(mod_int6, "PB_intyrtrt.RDS")
m2_loo=extract(mod_int6)$log_lik
mod_int6=readRDS("./model_output/PB_intyrtrt.RDS")
yrep2=rstan::extract(mod_int6)$pred_y
bayesplot::ppc_dens_overlay(y, yrep2)
y=PB_dat_M$proportion
con_pb=yrep2[,which(PB_dat_M$treatment=="control"& PB_dat_M$month==3)]
con_pbmat=as.matrix(con_pb)
con_pbs=con_pbmat[1:300,]
matplot(t(con_pbs), type="l", col="grey", main="PB control (March)")
mean_con_pb=apply(con_pb, 2, mean)
con_pb_obs=PB_dat_M%>%filter(treatment=="control"& month==3)
lines(mean_con_pb~c(1:length(mean_con_pb)), col="white")
points(con_pb_obs$proportion, col="black", cex=1 )



loo(m2_loo) #-558.4+26
print(mod_int6, pars=c("alpha", "phi", "year_eff", "trt_eff"))

#year, treatment effect and trigonometric functions for a circular variable(month)####
mod_int7=stan(model_code="
   data{
  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // reproductive indivs
  int <lower=0>  n[N];       // total males
  vector [N] year;// year
  vector [N]mon_cos;//cosine month
  vector [N]mon_sin;//sine of month
  vector[N] treatment;// treatment
 // int month[N]; //ID of each month
// int Nmon; //no.of months
 // int species[N]; //species ID
// int Nsp; //no.of species
 }
                
 parameters {
  real alpha;// intercept
  real year_eff; //slope year
  real trt_eff; //slope treatment effect
  real monc_eff;
  real mons_eff;
    real monc_eff2;
  real mons_eff2;
  //real<lower=0> sigma_mon[Nmon];//error for random intercept (month)
  //real <lower=0> mon_non;//non-centered error term for species
  //real<lower=0> sigma_sp[Nsp];//error for random intercept (species)
  //real <lower=0> sp_non;//non-centered error term for month
  real <lower=0> phi;
  real <lower=0, upper=1> pred_repro[N] ;//proportion of reproductive event 
              }
   
  transformed parameters{
  vector <lower=0, upper=1> [N] repro_mu; //so we can add statement describing proportion (not able to do in parameters block)
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  //vector [Nmon] alpha_mon; //random intercept per species
   //vector [Nsp] alpha_sp; //random intercept per species
  //vector [Nmon] yr_mon; //random slope per month for year effect
  //vector [Nmon] trt_mon;//random slope per month for treatment effect

  
  //for (j in 1:Nmon) {
  
  //alpha_mon[j]= mon_non*sigma_mon[j];}
  
 //  for (k in 1:Nsp) {
  
  //alpha_sp[k]= sp_non*sigma_sp[k];}
  
  //model:
  
  for (i in 1:N){
  
  repro_mu[i]= inv_logit(alpha+ year_eff*year[i]+trt_eff*treatment[i]+monc_eff*mon_cos[i]+mons_eff*mon_sin[i]+
  +monc_eff2*mon_cos[i]+mons_eff2*mon_sin[i]);
  }
  
  A = repro_mu * phi;
  B = (1 - repro_mu)* phi;
  
  }
 model {
  //priors
  alpha~normal(0,1);
  year_eff~ normal (0,1);
  trt_eff~ normal (0,1);
  monc_eff~normal(0,1);
  mons_eff~normal(0,1);
   monc_eff2~normal(0,1);
  mons_eff2~normal(0,1);
  //mon_non~ normal(0,1);
  //sigma_mon~ normal(0,1);
  phi ~normal(0,1);
  //sp_non~ normal(0,1);
  //sigma_sp~ normal(0,1);
  
  //model likelihood:
  
  pred_repro ~ beta(A, B); // survival estimate, beta dist.
  y~binomial(n, pred_repro); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate
 
 }
  
  generated quantities {
  
  real pred_y [N];//predictions on proportions
  real log_lik [N];// for looic calculations
  
    pred_y = beta_rng(A, B);
    
    for (x in 1:N){
    log_lik[x]= beta_lpdf(pred_repro[x]| A[x], B[x]);}
   
  }    ", data=dat_list, chains=2, iter=300)
saveRDS(mod_int7, "PB_intyrtrtmon.RDS")
PB_intyrtrtmon=readRDS("./model_output/PB_intyrtrtmon.RDS")
m3_loo=extract(PB_intyrtrtmon)$log_lik
loo(m3_loo) #-593.7+/- 27.7
y=PB_dat_M$proportion
yrep=extract(mod_int7)$repro_mu
yrep=as.data.frame(yrep)
ppc_dens_overlay(y, yrep[1:500,])
print(mod_int7, pars=c("alpha", "phi", "year_eff", "trt_eff", "monc_eff", "mons_eff"))

#PB control
con_pb=yrep[which(PB_dat_M$treatment=="control"& PB_dat_M$month==3)]
con_pbmat=as.matrix(con_pb)
con_pbs=con_pbmat[1:1000,]
matplot(t(con_pbs), type="l", col="grey", main="PB control (March)", ylim=c(0,1))
mean_con_pb=apply(con_pb, 2, mean)
con_pb_obs=PB_dat_M%>%filter(treatment=="control"& month==3)
lines(mean_con_pb~c(1:length(mean_con_pb)), col="white")
points(con_pb_obs$proportion, col="black", cex=2 )

#PB exclosure
ex_pb=yrep[which(PB_dat_M$treatment=="exclosure"& PB_dat_M$month==3)]
matplot(t(ex_pb), type="l", col="grey", main="PB exclosure (March)")
mean_ex_pb=apply(ex_pb, 2, mean)
ex_pb_obs=PB_dat_M%>%filter(treatment=="exclosure"& month==3)
lines(mean_ex_pb~c(1:length(mean_ex_pb)), col="white")
points(ex_pb_obs$proportion, col="black", cex=0.5)



#beta-binomial with autocorrelation structure (brms)####
#define custom response distribution (Buerkner)

total_prop=read.csv("./reconfigured_data/raw_cleaned/reprod_propn_male.csv")

beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)

stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"
stanvars <- stanvar(scode = stan_funs, block = "functions")

total_prop$trt<-ifelse(total_prop$treatment=="control", 0, 1) 
total_prop$years=(total_prop$year-mean(total_prop$year))/(2*sd(total_prop$year)) #standardize year
total_prop[is.na(total_prop)] <- 0 #set non-detects to 0


#subset for PB only
PB_dat_M=total_prop%>%filter(species=="PB", !(treatment=="spectabs"))

PB_dat_M[is.na(PB_dat_M)] <- 0 #set non-detects to 0
PB_dat_M$trt<-ifelse(PB_dat_M$treatment=="control", 0, 1) 
PB_dat_M$years=(PB_dat_M$year-mean(PB_dat_M$year))/(2*sd(PB_dat_M$year)) #standardize year
PBprop=PB_dat_M$proportion
PBrep=PB_dat_M$reproductive

pb_mod=brm(reproductive | vint (abundance) ~ ar(p=1, cov=TRUE)+trt, data= PB_dat_M,
           family=beta_binomial2, stanvar=stanvars, chains=2, iter=100)
saveRDS(pb_mod, "autocor_brms.RDS")

stancode(pb_mod)
expose_functions(pb_mod, vectorize = TRUE)
posterior_predict_beta_binomial2 <- function(i, prep,...) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  trials <- prep$data$vint1[i]
  beta_binomial2_rng(mu, phi, trials)
}

pp_check(pb_mod)

log_lik_beta_binomial2 <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, trials)
}

loo(pb_mod)
