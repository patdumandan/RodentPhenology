bkpt_con=stan(model_code="
data {
  int<lower=0> N;                      // sample size
  int <lower=0> abundance [N];
  int <lower=0> reproductive [N];
  vector [N] month;                       // estimate within-year breakpoint
}
parameters{
  real alpha; //intercept
  real beta1; //slope 1
  real beta2; //slope1
  
  real<lower=0,upper=1> proportion[N];        // response 
  
  real <lower=1, upper=12>bkpt;
  real <lower=0> sigma;

}
transformed parameters{
  vector <lower=0, upper=1> [N] prop_mu;
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  
  
for ( i in 1:N) {
  if (month[i] < bkpt)
    
    prop_mu[i]=inv_logit(alpha+beta1* (month[i]*bkpt));
  
  else
    prop_mu[i]= inv_logit(alpha+ beta2 * (month[i]*bkpt));
  
}

A= prop_mu*sigma;
B= (1-prop_mu)*sigma;

  }
model {
  
  alpha~normal(0,10);
  beta1~ normal(0,10);
  beta2~normal (0,10);
  bkpt~ normal (0,10);
  sigma~cauchy(0,2.5);
  
  proportion~ beta(A,B);
  
  reproductive~binomial(abundance, prop_mu); }
  

", data=datlist, chains=4, iter=1000)

datlist=list(N=length(PB_female_con$X),
             abundance=PB_female_con$abundance,
             reproductive=PB_female_con$reproductive,
             proportion=PB_female_con$proportion,
             month=PB_female_con$month)
