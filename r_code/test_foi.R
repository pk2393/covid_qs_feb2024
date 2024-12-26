rm(list=ls())
library(dplyr);library(magrittr);library(lubridate)
library(cmdstanr)

readr::read_csv("../data/synthetic_data_covariates.csv") -> df_synthetic_covariates

df_synthetic_covariates%>%
    mutate(
      age_30_39 = case_when(
        age_group ==3 ~ 1,
        T ~ 0
      ),
      age_40_49 = case_when(
        age_group ==4 ~ 1,
        T ~ 0
      ),
      age_50_59 = case_when(
        age_group ==5 ~ 1,
        T ~ 0
      ),
      age_60_69 = case_when(
        age_group ==6 ~ 1,
        T ~ 0
      ),
      age_70 = case_when(
        age_group ==7 ~ 1,
        T ~ 0
      )
      )->df_analysis_synthetic_covariates


covariates<-df_analysis_synthetic_covariates%>%
  mutate(intercept=1)%>%
  select(intercept, 
         age_30_39, age_40_49, age_50_59, age_60_69, age_70,
         sex, 
         diabetes, neoplasm, immunosuppression, respiratory, cardiovascular, 
         cerebrovascular, liver, obesity, alcohol, smoke, 
         hh_size, infection_history
         )%>%as.matrix()

data_stan=list(
  N=nrow(df_analysis_synthetic_covariates),
  M=dim(covariates)[2],
  pos_feb2024=df_analysis_synthetic_covariates$pos_feb2024, 
  imm_last=df_analysis_synthetic_covariates$imm_last, 
  infect_period=df_analysis_synthetic_covariates$infect_period2, 
  vax_period=df_analysis_synthetic_covariates$vax_period2,
  time_since_imm=df_analysis_synthetic_covariates$time_since_imm, 
  age_group=df_analysis_synthetic_covariates$age_group,
  pref=df_analysis_synthetic_covariates$pref,
  X=covariates
)

library(cmdstanr)

stanfile=paste0("./stan_code/test_foi_model.stan")

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
  adapt_delta = .825,
  seed = 3 
)

fit$cmdstan_diagnose()
