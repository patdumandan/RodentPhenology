dat_list=list(
  N=length(PB_dat_M$month),
  y=PB_dat_M$reproductive,
  n=PB_dat_M$abundance,
  Kar=1,
  Kma=1,
  N_tg=1,
  begin_tg=1,
  end_tg=413,
  nobs_tg=length(PB_dat_M$year),
  year=PB_dat_M$years,
  treatment=PB_dat_M$trt,
  mon_cos=PB_dat_M$mon_cos, 
  mon_sin=PB_dat_M$mon_sin,
  Nmon=length(unique(PB_dat_M$month)),
  Nsp=length(unique(PB_dat_M$species)))



pb_autocor=stan(model_code = "
functions {
  /* multi-normal log-PDF for time-series covariance structures 
  * assuming homogoneous variances
  * Args: 
    *   y: response vector 
  *   mu: mean parameter vector
  *   sigma: residual standard deviation
  *   chol_cor: cholesky factor of the correlation matrix
  *   se2: square of user defined standard errors 
  *     should be set to zero if none are defined 
  *   nobs: number of observations in each group 
  *   begin: the first observation in each group 
  *   end: the last observation in each group 
  * Returns: 
    *   sum of the log-PDF values of all observations 
  */ 
    real normal_time_hom_lpdf(vector y, vector mu, real sigma, matrix chol_cor, 
                              vector se2, int[] nobs, int[] begin, int[] end) {
      int I = size(nobs);
      int has_se = max(se2) > 0;
      vector[I] lp; 
      for (i in 1:I) { 
        matrix[nobs[i], nobs[i]] L;
        L = sigma * chol_cor[1:nobs[i], 1:nobs[i]];
        if (has_se) {
          // need to add 'se' to the correlation matrix itself
          L = multiply_lower_tri_self_transpose(L);
          L += diag_matrix(se2[begin[i]:end[i]]);
          L = cholesky_decompose(L);
        }
        lp[i] = multi_normal_cholesky_lpdf(
          y[begin[i]:end[i]] | mu[begin[i]:end[i]], L
        );
      }                        
      return sum(lp); 
    }
  /* multi-normal log-PDF for time-series covariance structures 
  * assuming heterogenous variances
  * Args: 
    *   y: response vector 
  *   mu: mean parameter vector
  *   sigma: residual standard deviation vector
  *   chol_cor: cholesky factor of the correlation matrix
  *   se2: square of user defined standard errors 
  *     should be set to zero if none are defined 
  *   nobs: number of observations in each group 
  *   begin: the first observation in each group 
  *   end: the last observation in each group 
  * Returns: 
    *   sum of the log-PDF values of all observations 
  */ 
    real normal_time_het_lpdf(vector y, vector mu, vector sigma, matrix chol_cor, 
                              vector se2, int[] nobs, int[] begin, int[] end) {
      int I = size(nobs);
      int has_se = max(se2) > 0;
      vector[I] lp; 
      for (i in 1:I) { 
        matrix[nobs[i], nobs[i]] L;
        L = diag_pre_multiply(sigma[begin[i]:end[i]], 
                              chol_cor[1:nobs[i], 1:nobs[i]]);
        if (has_se) {
          // need to add 'se' to the correlation matrix itself
          L = multiply_lower_tri_self_transpose(L);
          L += diag_matrix(se2[begin[i]:end[i]]);
          L = cholesky_decompose(L);
        }
        lp[i] = multi_normal_cholesky_lpdf(
          y[begin[i]:end[i]] | mu[begin[i]:end[i]], L
        );
      }                        
      return sum(lp); 
    }
  /* multi-student-t log-PDF for time-series covariance structures 
  * assuming homogoneous variances
  * Args: 
    *   y: response vector 
  *   nu: degrees of freedom parameter 
  *   mu: mean parameter vector
  *   sigma: scale parameter
  *   chol_cor: cholesky factor of the correlation matrix
  *   se2: square of user defined standard errors 
  *     should be set to zero if none are defined 
  *   nobs: number of observations in each group 
  *   begin: the first observation in each group 
  *   end: the last observation in each group 
  * Returns: 
    *   sum of the log-PDF values of all observations 
  */ 
    real student_t_time_hom_lpdf(vector y, real nu, vector mu, real sigma, 
                                 matrix chol_cor, vector se2, int[] nobs, 
                                 int[] begin, int[] end) { 
      int I = size(nobs);
      int has_se = max(se2) > 0;
      vector[I] lp; 
      for (i in 1:I) { 
        matrix[nobs[i], nobs[i]] Cov; 
        Cov = sigma * chol_cor[1:nobs[i], 1:nobs[i]];
        Cov = multiply_lower_tri_self_transpose(Cov);
        if (has_se) {
          Cov += diag_matrix(se2[begin[i]:end[i]]);
        }
        lp[i] = multi_student_t_lpdf(
          y[begin[i]:end[i]] | nu, mu[begin[i]:end[i]], Cov
        );
      }                        
      return sum(lp); 
    }
  /* multi-student-t log-PDF for time-series covariance structures 
  * assuming heterogenous variances
  * Args: 
    *   y: response vector 
  *   nu: degrees of freedom parameter 
  *   mu: mean parameter vector
  *   sigma: scale parameter vector
  *   chol_cor: cholesky factor of the correlation matrix
  *   se2: square of user defined standard errors 
  *     should be set to zero if none are defined 
  *   nobs: number of observations in each group 
  *   begin: the first observation in each group 
  *   end: the last observation in each group 
  * Returns: 
    *   sum of the log-PDF values of all observations 
  */ 
    real student_t_time_het_lpdf(vector y, real nu, vector mu, vector sigma, 
                                 matrix chol_cor, vector se2, int[] nobs, 
                                 int[] begin, int[] end) { 
      int I = size(nobs);
      int has_se = max(se2) > 0;
      vector[I] lp; 
      for (i in 1:I) { 
        matrix[nobs[i], nobs[i]] Cov; 
        Cov = diag_pre_multiply(sigma[begin[i]:end[i]], 
                                chol_cor[1:nobs[i], 1:nobs[i]]);
        Cov = multiply_lower_tri_self_transpose(Cov);
        if (has_se) {
          Cov += diag_matrix(se2[begin[i]:end[i]]);
        }
        lp[i] = multi_student_t_lpdf(
          y[begin[i]:end[i]] | nu, mu[begin[i]:end[i]], Cov
        );
      }                        
      return sum(lp); 
    }
  /* scale and correlate time-series residuals 
  * Args: 
    *   zerr: standardized and independent residuals
  *   sderr: standard deviation of the residuals
  *   chol_cor: cholesky factor of the correlation matrix
  *   nobs: number of observations in each group 
  *   begin: the first observation in each group 
  *   end: the last observation in each group 
  * Returns: 
    *   vector of scaled and correlated residuals
  */ 
    vector scale_time_err(vector zerr, real sderr, matrix chol_cor, 
                          int[] nobs, int[] begin, int[] end) { 
      vector[rows(zerr)] err; 
      for (i in 1:size(nobs)) { 
        err[begin[i]:end[i]] = 
          sderr * chol_cor[1:nobs[i], 1:nobs[i]] * zerr[begin[i]:end[i]];
      }                        
      return err; 
    }
  /* compute the cholesky factor of an AR1 correlation matrix
  * Args: 
    *   ar: AR1 autocorrelation 
  *   nrows: number of rows of the covariance matrix 
  * Returns: 
    *   A nrows x nrows matrix 
  */ 
    matrix cholesky_cor_ar1(real ar, int nrows) { 
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
      return cholesky_decompose(1 / (1 - ar^2) * mat); 
    }
  /* compute the cholesky factor of a MA1 correlation matrix
  * Args: 
    *   ma: MA1 autocorrelation 
  *   nrows: number of rows of the covariance matrix 
  * Returns: 
    *   A nrows x nrows MA1 covariance matrix 
  */ 
    matrix cholesky_cor_ma1(real ma, int nrows) { 
      matrix[nrows, nrows] mat; 
      mat = diag_matrix(rep_vector(1 + ma^2, nrows)); 
      if (nrows > 1) { 
        mat[1, 2] = ma; 
        for (i in 2:(nrows - 1)) { 
          mat[i, i - 1] = ma; 
          mat[i, i + 1] = ma; 
        } 
        mat[nrows, nrows - 1] = ma; 
      } 
      return cholesky_decompose(mat); 
    }
  /* compute the cholesky factor of an ARMA1 correlation matrix
  * Args: 
    *   ar: AR1 autocorrelation 
  *   ma: MA1 autocorrelation 
  *   nrows: number of rows of the covariance matrix 
  * Returns: 
    *   A nrows x nrows matrix 
  */ 
    matrix cholesky_cor_arma1(real ar, real ma, int nrows) { 
      matrix[nrows, nrows] mat; 
      vector[nrows] gamma; 
      mat = diag_matrix(rep_vector(1 + ma^2 + 2 * ar * ma, nrows)); 
      gamma[1] = (1 + ar * ma) * (ar + ma); 
      for (i in 2:nrows) { 
        gamma[i] = gamma[1] * pow(ar, i - 1); 
        for (j in 1:(i - 1)) { 
          mat[i, j] = gamma[i - j]; 
          mat[j, i] = gamma[i - j]; 
        } 
      } 
      return cholesky_decompose(1 / (1 - ar^2) * mat); 
    }
  /* compute the cholesky factor of a compound symmetry correlation matrix
  * Args: 
    *   cosy: compound symmetry correlation
  *   nrows: number of rows of the covariance matrix 
  * Returns: 
    *   A nrows x nrows covariance matrix 
  */ 
    matrix cholesky_cor_cosy(real cosy, int nrows) { 
      matrix[nrows, nrows] mat; 
      mat = diag_matrix(rep_vector(1, nrows)); 
      for (i in 2:nrows) { 
        for (j in 1:(i - 1)) { 
          mat[i, j] = cosy; 
          mat[j, i] = mat[i, j];
        } 
      } 
      return cholesky_decompose(mat); 
    }
  
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
  
}

 data{
  int<lower=1> N; // no.of obs
  int <lower=0> y[N];       // reproductive indivs
  int <lower=0>  n[N];       // total males
   int<lower=0> Kar;  // AR order 
  vector [N] year;// year
  int<lower=0> Kma;  // MA order
  int<lower=1> N_tg;
  //int<lower=1> begin_tg[N_tg];
  //int<lower=1> end_tg[N_tg];
  int<lower=1> nobs_tg[N_tg];
 }

 transformed data {
  int max_lag = max(Kar, Kma);
  int max_nobs_tg = max(nobs_tg);  // maximum dimension of the autocorrelation matrix               
 }          
 
parameters {
  real alpha;// intercept
  real year_eff; //slope year
  vector<lower=-1,upper=1>[Kar] ar;  // autoregressive coefficients
   //vector[N] zerr;  // unscaled residuals
  //real<lower=0> sderr;  // SD of residuals
  real<lower=0> phi;
   real <lower=0, upper=1> pred_repro[N] ;//proportion of reproductive event 
  }
  
transformed parameters {
  vector <lower=0, upper=1> [N] repro_mu; //so we can add statement describing proportion (not able to do in parameters block)
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  
  // cholesky factor of the autocorrelation matrix
  matrix[max_nobs_tg, max_nobs_tg] chol_cor;
  vector[N] err;  // actual residuals
  // compute residual covariance matrix
  chol_cor = cholesky_cor_ar1(ar[1], max_nobs_tg);
  // compute correlated time-series residuals
//  err = scale_time_err(zerr, sderr, chol_cor, nobs_tg, begin_tg, end_tg);
 
 
  
   for (i in 1:N){
  
  repro_mu[i]= inv_logit(alpha+ year_eff*year[i]);
  }
  
  A = repro_mu * phi;
  B = (1 - repro_mu)* phi;
  
  }

model {
  //priors
  alpha~normal(0,1);
  year_eff~ normal (0,1);
 // sderr~normal(0,2.5);
  //zerr~normal(0,2.5);
  phi~ normal(0,1);
  year_eff~normal(0,1);
  
  pred_repro ~ beta(A, B); // survival estimate, beta dist.
  y~binomial(n, pred_repro); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate
  }

generated quantities {
  
  real pred_y [N];//predictions on proportions
  real log_lik [N];// for looic calculations
  
    pred_y = beta_rng(A, B);
    
    for (x in 1:N){
    log_lik[x]= beta_lpdf(pred_repro[x]| A[x], B[x]);}
   
  }", data=dat_list, chains=2, iter=100) 

repro=as.vector(PB_dat_M$reproductive)
repro[413]    
