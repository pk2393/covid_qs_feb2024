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
  simplex[N] alpha;
  real<lower=0, upper=1e-2> weight_empty;
}

transformed parameters{
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
  // weighting
  weight_empty~uniform(0, 1e-2); // small weight to fill potentially empty elements
  
  vector[N] alpha_prior = rep_vector(5, N); // pbeta(1/5K) <=.01
  // vector[N] alpha = rep_vector(3.85, N); // pbeta(1/5K) <=.01
  alpha ~ dirichlet(alpha_prior);
  
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
