rm(list=ls())
library(dplyr);library(magrittr);library(lubridate)
library(cmdstanr)

readr::read_csv("synthetic_data_census.csv") -> df_synthetic_census

df_pop_sex <- readr::read_csv(file = "./data/census_pop_sex.csv")
df_pop_pref_m <- readr::read_csv(file = "./data/census_pop_pref_age_male.csv")
df_pop_pref_f <- readr::read_csv(file = "./data/census_pop_pref_age_female.csv")
df_job_m <- readr::read_csv(file = "./data/census_pop_job_age_male.csv")
df_job_f <- readr::read_csv(file = "./data/census_pop_job_age_female.csv")

data_stan=list(
  N=nrow(df_synthetic_census),
  
  age_group=df_synthetic_census$age_group,
  pref=df_synthetic_census$pref, 
  job_type=df_synthetic_census$job_type,   
  sex=df_synthetic_census$sex,
  
  K=nrow(df_job_f),
  A=ncol(df_job_f),
  ns_sex=c(sum(df_pop_sex$pop_m), sum(df_pop_sex$pop_f)),
  ns_pop_m_age=df_pop_sex$pop_m, 
  ns_pop_f_age=df_pop_sex$pop_f, 
  ns_pop_pref_m=as.matrix(df_pop_pref_m), 
  ns_pop_pref_f=as.matrix(df_pop_pref_f), 
  ns_job_m=as.matrix(df_job_m), 
  ns_job_f=as.matrix(df_job_f)

)

library(cmdstanr)

stanfile=paste0("./stan_code/test_census_model.stan")

model <- cmdstanr::cmdstan_model(stanfile, force_recompile=T)

niter=1e3
w=.5

fit<-model$sample(
  data=data_stan, 
  iter_warmup=niter*w, 
  iter_sampling=niter, 
  chains=4,
  parallel_chains = 4, 
  refresh=10,
  max_treedepth = 15,
  adapt_delta = .8,
  seed = 1
)

fit$cmdstan_diagnose()
