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

  array[N] int pref;
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
  
  // fitting
  // infection event in Feb 2024
  for(n in 1:N){
    target+=bernoulli_lpmf(pos_feb2024[n]|prob_infection[n]);
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

  


