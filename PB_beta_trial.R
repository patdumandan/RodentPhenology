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


PB_beta=stan(model_code="
data {
  int<lower=0> N;
  int <lower=0,upper=1> y[N];//proportion
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
}

transformed parameters{
  vector [Nmon] alpha_mon; //random intercept per month

  for (j in 1:Nmon) {
  
  alpha_mon[j]= mon_non*sigma_mon[j];
  }

}
model {
for (i in 1:N){
  y[i] ~ bernoulli_logit(alpha + alpha_mon[month[i]]+year_eff * year[i]+ trt_eff* treatment[i]);
}
  year_eff~ normal (0,1);
  trt_eff~ normal (0,1);
  mon_non~ normal(0,1);
  sigma_mon~ normal(0,1);

}
", data=dat_list_log, chains=2, iter=300)
