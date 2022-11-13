total_proportion=read.csv("./reconfigured_data/raw_cleaned/reprod_propn_male.csv", stringsAsFactors = FALSE)
total_proportion[is.na(total_proportion)] <- 0 #set non-detects to 0
total_proportion$trt<-ifelse(total_proportion$treatment=="control", 0, 1) 
total_proportion$years=(total_proportion$year-mean(total_proportion$year))/(2*sd(total_proportion$year)) #standardize year


PB_dat_M_con=total_proportion%>%filter(species=="PB", treatment=="control")
PB_dat_M=total_proportion%>%filter(species=="PB", !(treatment=="spectabs"))
PBprop=total_proportion$proportion
PBrep=total_proportion$reproductive
total_proportion$mon_cos= cos( 0.5235988*total_proportion$month)
total_proportion$mon_sin= sin(0.5235988*total_proportion$month)
total_proportion$spcode=as.integer(as.factor(total_proportion$species))
str(total_proportion)
PB_dat_M[is.na(PB_dat_M)] <- 0 #set non-detects to 0
PB_dat_M$trt<-ifelse(PB_dat_M$treatment=="control", 0, 1) 
PB_dat_M$years=(PB_dat_M$year-mean(PB_dat_M$year))/(2*sd(PB_dat_M$year)) #standardize year
PBprop=PB_dat_M$proportion
PBrep=PB_dat_M$reproductive
PB_dat_M$mon_cos= cos(2*pi*(PB_dat_M$month/12))
PB_dat_M$mon_sin= sin(2*pi*(PB_dat_M$month/12))

dat_list=list(
  N=length(total_proportion$month),
  y=total_proportion$reproductive,
  n=total_proportion$abundance,
  year=total_proportion$years,
  treatment=total_proportion$trt,
  mon_cos=total_proportion$mon_cos, 
  species=total_proportion$spcode,
  mon_sin=total_proportion$mon_sin,
  Nmon=length(unique(total_proportion$month)),
  Nsp=length(unique(total_proportion$species)))

mod_all_M=stan(model_code="
   data{
  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // reproductive indivs
  int <lower=0>  n[N];       // total males
  vector [N] year;// year
  vector [N]mon_cos;//cosine month
  vector [N]mon_sin;//sine of month
  vector[N] treatment;// treatment
  int species[N]; //species ID
  int Nsp; //no.of species
 }
                
 parameters {
  real alpha;// intercept
  real year_eff; //slope year
  real trt_eff; //slope treatment effect
  real monc_eff;
  real mons_eff;
 
  real<lower=0> sigma_sp[Nsp];//error for random intercept (species)
  real <lower=0> sp_non;//non-centered error term for month
  real <lower=0> phi;
  real <lower=0, upper=1> pred_repro[N] ;//proportion of reproductive event 
              }
   
  transformed parameters{
  vector <lower=0, upper=1> [N] repro_mu; //so we can add statement describing proportion (not able to do in parameters block)
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
   vector [Nsp] alpha_sp; //random intercept per species
 
  for (k in 1:Nsp) {
  
  alpha_sp[k]= sp_non*sigma_sp[k];}
  
  //model:
  
  for (i in 1:N){
  
  repro_mu[i]= inv_logit(alpha+ alpha_sp[species[i]]+ year_eff*year[i]+trt_eff*treatment[i]+monc_eff*mon_cos[i]+mons_eff*mon_sin[i]);
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
  phi ~normal(0,1);
  sp_non~ normal(0,1);
  sigma_sp~ normal(0,1);
  
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

y=PB_dat_M$proportion

PB_intyrtrtmon=readRDS("./model_output/PB_intyrtrtmon.RDS")
yrep=rstan::extract(PB_intyrtrtmon)$pred_y
yrep=as.data.frame(yrep)
ppc_dens_overlay(y, yrep)
print(PB_intyrtrtmon, pars=c("alpha", "phi", "year_eff", "trt_eff", "monc_eff", "mons_eff"))


#PB control
con_pb=yrep[which(PB_dat_M$treatment=="control"& PB_dat_M$month==3)]
con_pbmat=as.matrix(con_pb)
con_pbs=con_pbmat[1:1000,]
matplot(t(con_pbs), type="l", col="grey", main="PB control (March)", ylim=c(0,1))
mean_con_pb=apply(con_pb, 2, mean)
con_pb_obs=PB_dat_M%>%filter(treatment=="control"& month==3)
lines(mean_con_pb~c(1:length(mean_con_pb)), col="white")
points(con_pb_obs$proportion, col="black", cex=2, pch=21)
