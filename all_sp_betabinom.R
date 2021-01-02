#DATA MANIPULATION####
all_props$trt<-ifelse(all_props$trt=="control", 0, 1)
all_props$years=(all_props$year-mean(all_props$year))/(2*sd(all_props$year))
all_props[is.na(all_props)] <- 0 #set non-detects to 0
all_props$spcode=as.integer(all_props$species)

dat_list=list(
  N=length(all_props$month),
  y=all_props$repro,
  n=all_props$count,
  year=all_props$years,
  treatment=all_props$trt,
  species=all_props$spcode,
  month=all_props$month, 
  Nmon=length(unique(all_props$month)),
  Nsp=length(unique(all_props$species))
)

all_prop_male_mod=stan(model_code = "

data{
  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // reproductive indivs
  int <lower=0>  n[N];       // total males
  vector [N] year;// year
  vector[N] treatment;// treatment
  int month[N]; //ID of each month
  int Nmon; //no.of months
  int species[N]; //species ID
  int Nsp; //no.of species
 }
                
 parameters {
  real alpha;// intercept
  real year_eff; //slope year
  real trt_eff; //slope treatment effect
  real<lower=0> sigma_mon[Nmon];//error for random intercept (month)
  real <lower=0> mon_non;//non-centered error term for species
  real<lower=0> sigma_sp[Nsp];//error for random intercept (species)
  real <lower=0> sp_non;//non-centered error term for month
  real <lower=0> phi;
  real <lower=0, upper=1> pred_repro[N] ;//proportion of reproductive event 
              }
   
  transformed parameters{
  vector <lower=0, upper=1> [N] repro_mu; //mean estimated proportion of reproductive event
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  vector [Nmon] alpha_mon; //random intercept per species
   vector [Nsp] alpha_sp; //random intercept per species
  //vector [Nmon] yr_mon; //random slope per family for year effect
  //vector [Nmon] trt_mon;//random slope per family for treatment effect

  
  for (j in 1:Nmon) {
  
  alpha_mon[j]= mon_non*sigma_mon[j];
  }
  
   for (k in 1:Nsp) {
  
  alpha_sp[k]= sp_non*sigma_sp[k];
  }
  
  //model:
  
  for (i in 1:N){
  
  repro_mu[i]= inv_logit(alpha+alpha_mon[month[i]]+alpha_sp[species[i]]+ year_eff*year[i]+trt_eff*treatment[i]);
  }
  
  A = repro_mu * phi;
  B = (1 - repro_mu)* phi;
  
  }
 model {
  //priors
  year_eff~ normal (0,1);
  trt_eff~ normal (0,1);
  mon_non~ normal(0,1);
  sigma_mon~ normal(0,1);
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
   
  }", data=dat_list, chains=4, iter=3000)

#MODEL OUTPUT####
print(all_prop_male_mod, pars=c("alpha","alpha_mon", "alpha_sp","trt_eff", "year_eff"))
saveRDS(all_prop_male_mod, "all_sp_betabinom_mod.RDS")
