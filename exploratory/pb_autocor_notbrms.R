total_proportion=read.csv("./reconfigured_data/raw_cleaned/reprod_propn_male.csv")
PB_dat_M=total_proportion%>%filter(species=="PB", !(treatment=="spectabs"))
PB_dat_M[is.na(PB_dat_M)] <- 0 #set non-detects to 0
PB_dat_M$trt<-ifelse(PB_dat_M$treatment=="control", 0, 1) 
PB_dat_M$years=(PB_dat_M$year-mean(PB_dat_M$year))/(2*sd(PB_dat_M$year)) #standardize year
PBprop=PB_dat_M$proportion
PBrep=PB_dat_M$reproductive

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


pbyr_autocor=stan(model_code=
"functions { 
  matrix cov_matrix_ar1(real ar, real sigma, int nrows) { 
    matrix[nrows, nrows] mat; 
    vector[nrows - 1] gamma; 
    mat = diag_matrix(rep_vector(1, nrows)); 
    for (i in 2:nrows) { 
      gamma[i - 1] = pow(ar, i - 1); 
      for (j in 1:(i - 1)) { 
        mat[i, j] = gamma[i - j]; 
        mat[j, i] = gamma[i - j]; 
      } 
    } 
    return sigma^2 / (1 - ar^2) * mat; 
  }
} 

data { 
  int<lower=1> N;  // total number of observations 
//  vector [N] year;// year
  int <lower=0> y[N];       // reproductive indivs
  int <lower=0>  n[N];       // total males
}

transformed data {
  vector[N] se2 = rep_vector(0, N); 
} 
parameters { 
  real alpha;
//  real year_eff;
  real<lower=0> sigma;  // residual SD 
  real <lower=-1,upper=1> phi;  // autoregressive effects 
  real <lower=0, upper=1> pred_repro[N] ;//proportion of reproductive indivs
  real <lower=0>psi;//overdispersion param
} 

transformed parameters{
  vector [N] repro_mu; //so we can add statement describing proportion (not able to do in parameters block)
  vector[N] A;
  vector [N] B;
 
//model:
  
  for (i in 1:N){
  
  repro_mu[i]= inv_logit(alpha);
  }
  
  A = repro_mu * psi;
  B = (1 - repro_mu)* psi;
}

model {

  matrix[N, N] res_cov_matrix;
  matrix[N, N] Sigma; 
  res_cov_matrix = cov_matrix_ar1(phi, sigma, N);
  Sigma = res_cov_matrix + diag_matrix(se2);
  Sigma = cholesky_decompose(Sigma); 

//likelihood:

  alpha~normal(0,1);
//  year_eff~ normal (0,1);
  psi~normal(0,1);
  sigma ~ cauchy(0,5);
  
  pred_repro ~beta (A,B);
  y~ binomial(n, pred_repro);
}

generated quantities {
  
  real pred_y [N];//predictions on proportions
  real log_lik [N];// for looic calculations
  
  pred_y = beta_rng(A, B);
  
  for (x in 1:N){
    log_lik[x]= beta_lpdf(pred_repro[x]| A[x], B[x]);}
  
}", data=dat_list, chains=2, iter=300)
launch_shinystan(pbyr_autocor)
saveRDS(pbyr_autocor, "pbyr_autocor1.RDS")

pbyr_autocor1=readRDS("pbyr_autocor1.RDS")
yrep=rstan::extract(pbyr_autocor1)$pred_y

con_pb=yrep[,which(PB_dat_M$treatment=="control"& PB_dat_M$month==3)]
con_pbmat=as.matrix(con_pb)
con_pbs=con_pbmat[1:300,]
matplot(t(con_pbs), type="l", col="grey", main="PB control (March)")
mean_con_pb=apply(con_pb, 2, mean)
con_pb_obs=PB_dat_M%>%filter(treatment=="control"& month==3)
lines(mean_con_pb~c(1:length(mean_con_pb)), col="white")
points(con_pb_obs$proportion, col="black", cex=1 )


ppc_dens_overlay(y, yrep)


require(lmtest)
res=rstan::extract(pbyr_autocor)$pred_y
res1=apply(res,2,mean)
resq=cbind(res1, PB_dat_M$proportion)
str(resq)
resq=as.data.frame(resq)
resq$diff=resq$V2-resq$res1
dwtest(resq$diff~PB_dat_M$year) #1.99; P=0.42
