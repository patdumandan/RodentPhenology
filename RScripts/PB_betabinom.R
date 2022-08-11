#data manipulation####
PB_props$trt<-ifelse(PB_props$variable=="repro_con", 0, 1)
PB_props$years=(PB_props$year-mean(PB_props$year))/(2*sd(PB_props$year))
PB_props[is.na(PB_props)] <- 0 #set non-detects to 0

dat_list=list(
  N=length(PB_props$month),
  y=PB_props$value,
  n=PB_props$count,
  year=PB_props$years,
  treatment=PB_props$trt,
  month=PB_props$month,
  Nmon=12
)

#data analysis####
pb_prop_male_mod=stan(model_code = "

data{
  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // reproductive indivs
  int <lower=0>  n[N];       // total males
  vector [N] year;// year
   vector [N]mon_cos;//cosine month
  vector [N]mon_sin;//sine of month
  vector[N] treatment;// treatment
  int month[N]; //ID of each month
  int Nmon; //no.of months
 }
                
 parameters {
  real alpha;// intercept
  real year_eff; //slope year
 // real trt_eff; //slope treatment effect
 real<lower=0> sigma_mon[Nmon];//error for random intercept (month)
  real <lower=0> mon_non;//non-centered error term for month
  real <lower=0> phi;
  real monc_eff;
  real mons_eff;
  real <lower=0, upper=1> pred_repro[N] ;//proportion of reproductive event 
              }
   
  transformed parameters{
  vector <lower=0, upper=1> [N] repro_mu; //mean estimated proportion of reproductive event
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  vector [Nmon] alpha_mon; //random intercept per species
  //vector [Nmon] yr_mon; //random slope per family for year effect
  //vector [Nmon] trt_mon;//random slope per family for treatment effect

  
  for (j in 1:Nmon) {
  
  alpha_mon[j]= mon_non*sigma_mon[j];
  }
  
  //model:
  
  for (i in 1:N){
  
  repro_mu[i]= inv_logit(alpha+alpha_mon[month[i]]*month[i]+ year_eff*year[i]);
  }
  
  A = repro_mu * phi;
  B = (1 - repro_mu)* phi;
  
  }
 model {
  //priors
  alpha~normal(0,1);
  year_eff~ normal (0,1);
//  trt_eff~ normal (0,1);
  mon_non~ normal(0,1);
  sigma_mon~ normal(0,1);
  phi ~normal(0,1);
  monc_eff~normal(0,1);
  mons_eff~normal(0,1);

  
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
   
  }", data=dat_list, chains=2, iter=300)

#model output####
print(all_prop_male_mod, pars=c("alpha","alpha_mon", "alpha_sp","trt_eff", "year_eff"))
saveRDS(all_prop_male_mod, "all_sp_betabinom_mod.RDS")

y=PB_dat_M$proportion
yrep2=rstan::extract(pb_prop_male_mod)$repro_mu
con_pb=yrep2[,which(PB_dat_M$treatment=="control"& PB_dat_M$month==3)]
con_pbmat=as.matrix(con_pb)
con_pbs=con_pbmat[1:300,]

matplot(t(con_pbs), type="l", col="grey", main="PB control (March)")
mean_con_pb=apply(con_pb, 2, mean)
con_pb_obs=PB_dat_M%>%filter(treatment=="control"& month==3)
lines(mean_con_pb~c(1:length(mean_con_pb)), col="blue")
points(con_pb_obs$proportion, col="black", cex=1 )
ppc_dens_overlay(y, yrep2)
