#DATA MANIPULATION####
all_female_sp$trt<-ifelse(all_female_sp$trt=="control", 0, 1) 
all_female_sp$years=(all_female_sp$year-mean(all_female_sp$year))/(2*sd(all_female_sp$year)) #standardize year
all_female_sp[is.na(all_female_sp)] <- 0 #set non-detects to 0
all_female_sp$spcode=as.integer(all_female_sp$species) 

dat_fem_list=list(
  N=length(all_female_sp$month),
  y=all_female_sp$repro,
  n=all_female_sp$count,
  year=all_female_sp$years,
  treatment=all_female_sp$trt,
  species=all_female_sp$spcode,
  month=all_female_sp$month, 
  Nmon=length(unique(all_female_sp$month)),
  Nsp=length(unique(all_female_sp$species))
)

#MODEL INFRASTRUCTURE ####
all_prop_female_mod=stan(model_code = "

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
  //vector [Nmon] yr_mon; //random slope per month for year effect
  //vector [Nmon] trt_mon;//random slope per month for treatment effect

  
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
   
  }", data=dat_fem_list, chains=4, iter=3000, control=list(adapt_delta=0.99))

#MODEL OUTPUT####
print(all_prop_female_mod, pars=c("alpha","alpha_mon", "alpha_sp","trt_eff", "year_eff"))
#3 DIVERGENT TRANSITIONS
saveRDS(all_prop_female_mod, "all_female_sp_betabinom_mod.RDS")


#MODEL EVALUATION###

post_female=rstan::extract(all_prop_female_mod)$log_lik
loo(post_female) #-843.6 +/- 58.9


#DATA VISUALIZATION####
all_female_sp%>%
  ggplot(mapping=aes(x=reorder(Month,month), y=repro, fill=trt))+
  geom_col(position=position_dodge())+
  xlab("Month")+ylab("No.of reproductive individuals")+
  facet_wrap(~species)+theme(axis.text.x = element_text(angle = 45, hjust=1),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(fill="set-up")+
  scale_fill_grey()
plot_all_raw

all_female_sp%>%
  ggplot(mapping=aes(x=year, y=repro))+
  geom_col(position = position_dodge())+labs(fill="set-up")+
  xlab("Year")+ylab("No.ofreproductive individuals")+
  facet_wrap(~trt)+
  ggtitle("Female rodents")

