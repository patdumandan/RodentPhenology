allmt$years=(allmt$year-mean(allmt$year))/(2*sd(allmt$year))
allmt$spcode=as.integer(allmt$species) 

require(rstan)

male_trt=stan(model_code="
          data {
          
          int <lower=0> N; //no.of rows
          real <lower=0, upper=1> proportion[N];//proportion of reproductive indivs
          vector[N] years;
          
          }    
              
          parameters {
          
          real alpha;
          real beta;
          }
              
          model {
          
          proportion~bernoulli_logit(alpha+beta*years);
          
          
          alpha~normal(0,1);
          beta~normal(0,1);
          }", data= list(N=length(allmt$proportion), years=allmt$years, proportion=allmt$proportion),
              chains=2, iter=100)
              
              
              
              ")