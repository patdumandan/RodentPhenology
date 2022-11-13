pb_aut=stan(model_code="
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
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  // data for custom integer vectors
  int vint1[N];
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data needed for ARMA correlations
  int<lower=0> Kar;  // AR order
  int<lower=0> Kma;  // MA order
  // see the functions block for details
  int<lower=1> N_tg;
  int<lower=1> begin_tg[N_tg];
  int<lower=1> end_tg[N_tg];
  int<lower=1> nobs_tg[N_tg];
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  int max_lag = max(Kar, Kma);
  int max_nobs_tg = max(nobs_tg);  // maximum dimension of the autocorrelation matrix
  // no known standard errors specified by the user
  vector[N] se2 = rep_vector(0.0, N);
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector<lower=-1,upper=1>[Kar] ar;  // autoregressive coefficients
  vector[N] zerr;  // unscaled residuals
  real<lower=0> sderr;  // SD of residuals
  real<lower=0> phi;
}
transformed parameters {
  // cholesky factor of the autocorrelation matrix
  matrix[max_nobs_tg, max_nobs_tg] chol_cor;
  vector[N] err;  // actual residuals
  // compute residual covariance matrix
  chol_cor = cholesky_cor_ar1(ar[1], max_nobs_tg);
  // compute correlated time-series residuals
  err = scale_time_err(zerr, sderr, cholcor, nobs_tg, begin_tg, end_tg);
}
model {
  // likelihood including all constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b + err;
    for (n in 1:N) {
      // apply the inverse link function
      mu[n] = inv_logit(mu[n]);
    }
    for (n in 1:N) {
      target += beta_binomial2_lpmf(Y[n] | mu[n], phi, vint1[n]);
    }
  }
  // priors including all constants
  target += student_t_lpdf(Intercept | 3, 0, 2.5);
  target += student_t_lpdf(sderr | 3, 0, 2.5);
  target += std_normal_lpdf(zerr);
  target += gamma_lpdf(phi | 0.01, 0.01);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}")