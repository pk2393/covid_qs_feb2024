functions{
  real protection(real a, real frac, real g1, real g2, real time){
    return(a*(frac*exp(-g1*time)+(1-frac)*exp(-g2*time)));
  }
  real hazard_func(real c){
    return (1-exp(-c));
  }
}

data{
  int N;
  int M;
  array[N] int pos_feb2024;
  vector[N] imm_last;
  vector[N] infect_period;
  vector[N] vax_period;
  vector[N] time_since_imm;
  matrix[N, M] X;
  array[N] int sex;

  int K;
  int A;
  array[N] int pref;
  array[N] int age_group;
  array[N] int job_type;
  
  array[2] int ns_sex;
  array[A] int ns_pop_m_age;
  array[A] int ns_pop_f_age;
  array[47, A] int ns_pop_pref_m;
  array[47, A] int ns_pop_pref_f;
  array[K, A] int ns_job_f;
  array[K, A] int ns_job_m;
}

parameters{
  vector[M] beta_risk_background;

  real<lower=0, upper=1> v_infect_xbb;
  real<lower=0, upper=1> v_infect_preX;
  real<lower=0, upper=1> v_infect_wuhan;
  
  real<lower=0, upper=1> v_vax_xbb;
  real<lower=0, upper=1> v_vax_bivalent;
  real<lower=0, upper=1> v_vax_wuhan;
  
  real<lower=0> h1;
  real<lower=0> h2;
  real<lower=0, upper=1> frac_short_infect;
  real<lower=0, upper=1> frac_short_vax;

  simplex[N] alpha;
  real<lower=0, upper=1e-2> weight_empty;
  
  vector[46] RE_pref_log_raw;
  real<lower=0> sig_RE_pref;
  
  real<lower=0> tau;
  
}

transformed parameters{
  real<lower=0> g1 = log(2)/h1;
  real<lower=0> g2 = log(2)/h2;
  vector[N] protection_individual;
  for(n in 1:N){ 
      if(imm_last[n]==1&&infect_period[n]==1){
      }else if(imm_last[n]==1&&infect_period[n]==2){ //infection by xbb sublineages
        protection_individual[n]=protection(v_infect_xbb, frac_short_infect, g1, g2, time_since_imm[n]);
      }else if(imm_last[n]==1&&infect_period[n]==3){ //infection by pre-XBB Omicron
        protection_individual[n]=protection(v_infect_preX, frac_short_infect, g1, g2, time_since_imm[n]);
      }else if(imm_last[n]==1&&infect_period[n]==4){ //infection by pre-Omicron
        protection_individual[n]=protection(v_infect_wuhan, frac_short_infect, g1, g2, time_since_imm[n]);
      }      else if(imm_last[n]==2&&vax_period[n]==1){ // vax by xbb.1.5
        protection_individual[n]=protection(v_vax_xbb, frac_short_vax, g1, g2, time_since_imm[n]);
      }else if(imm_last[n]==2&&vax_period[n]==2){ //vax by bivalent
        protection_individual[n]=protection(v_vax_bivalent, frac_short_vax, g1, g2, time_since_imm[n]);
      }else if(imm_last[n]==2&&vax_period[n]==3){ //vax by monovalent 
        protection_individual[n]=protection(v_vax_wuhan, frac_short_vax, g1, g2, time_since_imm[n]);
      }else{
        protection_individual[n]=0;
      }
  }
 
  vector[47] RE_pref_log;
  RE_pref_log = append_row(RE_pref_log_raw, -sum(RE_pref_log_raw));
  vector[N] RE_pref_individual;
  for(n in 1:N){
    RE_pref_individual[n]= exp(RE_pref_log[pref[n]]);
  }

  vector<lower=0>[N] risk = exp(X*beta_risk_background);
  vector<lower=0>[N] lambda; 
  for(n in 1:N){
    lambda[n]=risk[n]*RE_pref_individual[n];
  }
  vector<lower=0>[N] lambda_effective;
  for(n in 1:N){
    lambda_effective[n]=(1-protection_individual[n])*lambda[n];
  }

  vector<lower=0, upper=1>[N] prob_infection;
  for(n in 1:N){
    prob_infection[n]=hazard_func(lambda_effective[n]);
    }
  
  vector[2] weight_sex_ns = rep_vector(0, 2);
  vector[A] weight_age_ns_f = rep_vector(0, A);
  vector[A] weight_age_ns_m = rep_vector(0, A);
  matrix[47, A] weight_pref_ns_f = rep_matrix(weight_empty, 47, A);
  matrix[47, A] weight_pref_ns_m = rep_matrix(weight_empty, 47, A);
  matrix[K, A] weight_job_ns_f = rep_matrix(weight_empty, K, A);
  matrix[K, A] weight_job_ns_m = rep_matrix(weight_empty, K, A);

  for(n in 1:N){
    int sex_i=sex[n];
    int age_group_i=age_group[n]-1;
    int pref_i=pref[n];
    int job_i=job_type[n];
    if(sex_i==0){
      weight_sex_ns[sex_i+1]+=alpha[n];
      weight_age_ns_m[age_group_i]+=alpha[n];
      weight_pref_ns_m[pref_i, age_group_i]+=alpha[n];
      weight_job_ns_m[job_i, age_group_i]+=alpha[n];
    }else{
      weight_sex_ns[sex_i+1]+=alpha[n];
      weight_age_ns_f[age_group_i]+=alpha[n];
      weight_pref_ns_f[pref_i, age_group_i]+=alpha[n];
      weight_job_ns_f[job_i, age_group_i]+=alpha[n];
    }
  }

  weight_sex_ns /= sum(weight_sex_ns);
  weight_age_ns_f /= sum(weight_age_ns_f);
  weight_age_ns_m /= sum(weight_age_ns_m);
  
  for (a in 1:A) {
     weight_pref_ns_f[, a] /= sum(weight_pref_ns_f[, a]);
     weight_pref_ns_m[, a] /= sum(weight_pref_ns_m[, a]);
     weight_job_ns_f[, a] /= sum(weight_job_ns_f[, a]);
     weight_job_ns_m[, a] /= sum(weight_job_ns_m[, a]);
  }
}

model{
  // priors: immunity
  v_infect_xbb ~ beta(10, 10);
  v_infect_preX ~ beta(10, 10);
  v_infect_wuhan ~ beta(10, 10); 
  
  v_vax_xbb ~ beta(10, 10);
  v_vax_bivalent ~ beta(10, 10);
  v_vax_wuhan ~ beta(10, 10);
  
  h1~normal(60, 10); 
  h2~normal(360, 120);  
  frac_short_infect ~ beta(20, 20);
  frac_short_vax ~ beta(20, 20);
  
  // prior information on protection
  log(1-protection(v_vax_xbb, frac_short_vax, g1, g2, 30)) ~ normal(log(1-0.442), (log(1-0.217)-log(1-0.603))/3.92); // VE XBB-monovalent  @ 0-2 mos = 0.442(0.217 - 0.603) -> JN1: may be slightly lower
  log(1-protection(v_vax_xbb, frac_short_vax, g1, g2, 90)) ~ normal(log(1-0.241), (log(1-0.007)-log(1-0.429))/3.92); // VE XBB-monovalent  @ 2-4 mos = 0.241(-0.007 - 0.429)
  log(1-protection(v_vax_xbb, frac_short_vax, g1, g2, 150)) ~ normal(log(1-0.267), (log(1+0.275)-log(1-0.579))/3.92); // VE XBB-monovalent  @ 4-6 mos = 0.267(-0.275 - 0.579)
  log(1-protection(v_vax_bivalent, frac_short_vax, g1, g2, 30)) ~ normal(log(1-0.022), (log(1+0.357)-log(1-0.295))/3.92); // VE bivalent: @ 0-2mo 0.022(-0.357, 0.295)
  log(1-protection(v_vax_bivalent, frac_short_vax, g1, g2, 90)) ~ normal(log(1-0.151), (log(1+0.554)-log(1-0.536))/3.92); // VE bivalent: @ 2-4mo 0.151(-0.554, 0.536)
  log(1-protection(v_infect_xbb, frac_short_infect, g1, g2, 90)) ~ normal(log(1-0.493), (log(1-0.292)-log(1-0.636))/3.92);  // Protection from infection @ 0-6 mos = 0.493 (0.292 - 0.636)

  log(1-protection(v_vax_xbb, frac_short_vax, g1, g2, 30)) ~ normal(log(1-0.45), (log(1-0.3)-log(1-0.6))/3.92); //  range 30-60: young = 40%, elder = 50%

  log(1-(1-protection(v_vax_xbb, frac_short_vax, g1, g2, 80))/(1-protection(v_vax_bivalent, frac_short_vax, g1, g2, 674))) ~ normal(log(1-0.49),(log(1-0.19)-log(1-0.68))/3.92); // 0.49(0.19-0.68)

  // lambda : intercept
  beta_risk_background[1] ~ normal(0, 10);
  // lambda: others
  tau ~ cauchy(0,1);
  beta_risk_background[2:M] ~ double_exponential(0, tau);
  
  // priors: RE by prefecture 
  sig_RE_pref ~  inv_gamma(2, 2); 
  RE_pref_log ~ normal(0, sig_RE_pref);

  // weighting
  weight_empty~uniform(0, 1e-2); // small weight to fill potentially empty elements
  
  vector[N] alpha_prior = rep_vector(5, N); // pbeta(1/4K) <=.01
  // vector[N] alpha = rep_vector(3.85, N); // pbeta(1/5K) <=.01
  alpha ~ dirichlet(alpha_prior);
  
  // fitting
  // infection event in Feb 2024
  for(n in 1:N){
    target+=bernoulli_lpmf(pos_feb2024[n]|prob_infection[n]);
  }
  
  // Census
  target+=multinomial_lpmf(ns_sex|weight_sex_ns);
  target+=multinomial_lpmf(ns_pop_m_age|weight_age_ns_m);
  target+=multinomial_lpmf(ns_pop_f_age|weight_age_ns_f);
  
  for(a in 1:A){
      target+=multinomial_lpmf(ns_job_m[:, a]|weight_job_ns_m[,a]);
      target+=multinomial_lpmf(ns_job_f[:, a]|weight_job_ns_f[,a]);
      
      target+=multinomial_lpmf(ns_pop_pref_m[:, a]|weight_pref_ns_m[,a]);
      target+=multinomial_lpmf(ns_pop_pref_f[:, a]|weight_pref_ns_f[,a]);
      }
}

generated quantities{
  vector[181] decay_infect_xbb;
  for(n in 1:181){
      decay_infect_xbb[n]=protection(v_infect_xbb, frac_short_infect, g1, g2, n-1);
  }
  vector[181] decay_infect_preX;
  for(n in 1:181){
      decay_infect_preX[n]=protection(v_infect_preX, frac_short_infect, g1, g2, n-1);
  }
  vector[181] decay_infect_wuhan;
  for(n in 1:181){
      decay_infect_wuhan[n]=protection(v_infect_wuhan, frac_short_infect, g1, g2, n-1);
  }
  vector[181] decay_vax_xbb;
  for(n in 1:181){
      decay_vax_xbb[n]=protection(v_vax_xbb, frac_short_vax, g1, g2, n-1);
  }
  vector[181] decay_vax_bivalent;
  for(n in 1:181){
      decay_vax_bivalent[n]=protection(v_vax_bivalent, frac_short_vax, g1, g2, n-1);
  }
  vector[181] decay_vax_wuhan;
  for(n in 1:181){
      decay_vax_wuhan[n]=protection(v_vax_wuhan, frac_short_vax, g1, g2, n-1);
  }
  
}
