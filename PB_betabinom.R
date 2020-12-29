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

pb_prop_male_mod=stan(model_code = "

data{
  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // reproductive indivs
  int <lower=0>  n[N];       // total males
  vector [N] year;// year
  vector[N] treatment;// treatment
  int month[N]; //ID of each month
  int Nmon; //no.of months
 }
                
 parameters {
  real alpha;// intercept
  real year_eff; //slope year
  real trt_eff; //slope treatment effect
  real<lower=0> sigma_mon[Nmon];//error for random intercept (month)
  real <lower=0> mon_non;//non-centered error term for month
  real <lower=0> phi;
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
  
  repro_mu[i]= inv_logit(alpha+alpha_mon[month[i]]+ year_eff*year[i]+trt_eff*treatment[i]);
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
