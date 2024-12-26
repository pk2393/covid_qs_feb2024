rm(list=ls())
library(dplyr);library(magrittr);library(lubridate)
library(cmdstanr)

df_analysis <- readr::read_csv(file = "../data/data_questionnaire.csv")%>%
  filter(infect_period2!=1)

df_pop_sex <- readr::read_csv(file = "../data/census_pop_sex.csv")
df_pop_pref_m <- readr::read_csv(file = "../data/census_pop_pref_age_male.csv")
df_pop_pref_f <- readr::read_csv(file = "../data/census_pop_pref_age_female.csv")
df_job_m <- readr::read_csv(file = "../data/census_pop_job_age_male.csv")
df_job_f <- readr::read_csv(file = "../data/census_pop_job_age_female.csv")

covariates<-df_analysis%>%
  mutate(intercept=1)%>%
  select(intercept, 
         age_30_39, age_40_49, age_50_59, age_60_69, age_70,
         sex, 
         diabetes, neoplasm, immunosuppression, respiratory, cardiovascular, 
         cerebrovascular, liver, obesity, alcohol, smoke, 
         hh_size, infection_history
         )%>%as.matrix()

data_stan=list(
  N=nrow(df_analysis),
  M=dim(covariates)[2],
  
  pos_feb2024=df_analysis$pos_feb2024, 
  imm_last=df_analysis$imm_last, 
  infect_period=df_analysis$infect_period2, 
  vax_period=df_analysis$vax_period2,
  time_since_imm=df_analysis$time_since_imm, 
  age_group=df_analysis$age_group, 
  pref=df_analysis$pref, 
  job_type=df_analysis$job_type,
  sex=df_analysis$sex,
  X=covariates, 
  
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

#####
stanfile=paste0("../stan_code/stan_model.stan")

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
