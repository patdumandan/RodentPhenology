PB_props$treat<-ifelse(PB_props$trt=="control", 0, 1)
PB_props$years=(PB_props$year-mean(PB_props$year))/(2*sd(PB_props$year))
PB_props[is.na(PB_props)] <- 0 #set non-detects to 0

dat_list_beta=list(
  N=length(PB_props$month),
  y=PB_props$proportion,
  year=PB_props$years,
  treatment=PB_props$treat,
  month=PB_props$month,
  Nmon=12
)

PB_beta=betareg(proportion~years+treat, link="logit",data=all_sp)

PB_beta=stan(model_code="
data {
  int <lower=0> N;
  real <lower=0,upper=1> y[N];//proportion
  vector [N] year;// year
  vector[N] treatment;// treatment
  int month[N]; //ID of each month
  int Nmon; //no.of months
}
parameters {
  real alpha;
  real year_eff;
  real trt_eff;
  real<lower=0> sigma_mon[Nmon];//error for random intercept (month)
  real <lower=0> mon_non;//non-centered error term for month
  real <lower=0> phi;
}

transformed parameters{
  vector [Nmon] alpha_mon; //random intercept per month
  vector <lower=0, upper=1> [N] prop_mu; //estimated survival 
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  
  for (j in 1:Nmon) {
  
  alpha_mon[j]= mon_non*sigma_mon[j];
  }

  for (i in 1:N){
  prop_mu[i] = inv_logit(alpha + alpha_mon[month[i]]+year_eff * year[i]+ trt_eff* treatment[i]);
}

 A = prop_mu * phi;
 B = (1 - prop_mu )* phi;

}

model {

  year_eff~ normal (0,1);
  trt_eff~ normal (0,1);
  mon_non~ normal(0,1);
  sigma_mon~ normal(0,1);
  phi~ normal(0,1);
   
  for (f in 1:N){
  
  y[f]~ beta(A, B);
  }
}
generated quantities {
  
  real log_lik [N];//predictions
  
    log_lik = beta_rng(A, B);
   
  }

", data=dat_list_beta, chains=1, iter=300)
