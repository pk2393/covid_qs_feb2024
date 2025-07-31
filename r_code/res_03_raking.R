rm(list=ls())
library(dplyr);library(magrittr);library(lubridate); library(posterior)

load(file = "path_to_RData") # load cmdstan fit object

###
### male vs female distribution, national
###
weight_sex_posterior <- subset_draws(sampling_result_2, variable = "weight_sex_ns")
weight_sex_posterior_summary <- summarize_draws(weight_sex_posterior,
                                                rhat,           
                                                ess_bulk,        
                                                ess_tail,        
                                                mean,
                                                median,
                                                ~quantile(.x, probs = c(.5, 0.025, 0.975)))

DescTools::MultinomCI(pop_sex, conf.level = .95, method = "goodman") -> ci_census

weight_sex_posterior_summary%>%
  mutate(census_ci_mean = ci_census[,1], 
         census_ci_lwr95 = ci_census[,2], 
         census_ci_upr95 = ci_census[,3] 
  )%>%
  readr::write_csv(file = "path_csv")

###
### male and female distribution, national, by age group
###
weight_male_age_posterior <- subset_draws(sampling_result_2, variable = "weight_age_ns_m")
weight_male_age_posterior_summary <- summarize_draws(weight_male_age_posterior,
                                                rhat,           
                                                ess_bulk,        
                                                ess_tail,        
                                                mean,
                                                median,
                                                ~quantile(.x, probs = c(.5, 0.025, 0.975)))

DescTools::MultinomCI(pop_m, conf.level = .95, method = "goodman") -> ci_census

weight_male_age_posterior_summary %>%
  mutate(census_ci_mean = ci_census[,1], 
         census_ci_lwr95 = ci_census[,2], 
         census_ci_upr95 = ci_census[,3] 
  )%>%
  readr::write_csv(file = "path_csv")


##
weight_female_age_posterior <- subset_draws(sampling_result_2, variable = "weight_age_ns_f")
weight_female_age_posterior_summary <- summarize_draws(weight_female_age_posterior,
                                                     rhat,           
                                                     ess_bulk,        
                                                     ess_tail,        
                                                     mean,
                                                     median,
                                                     ~quantile(.x, probs = c(.5, 0.025, 0.975)))

DescTools::MultinomCI(pop_f, conf.level = .95, method = "goodman") -> ci_census

weight_female_age_posterior_summary  %>%
  mutate(census_ci_mean = ci_census[,1], 
         census_ci_lwr95 = ci_census[,2], 
         census_ci_upr95 = ci_census[,3] )%>%
  readr::write_csv(file = "path_csv")

###
### male and female distribution, national, by age group
###
weight_male_age_pref_posterior <- subset_draws(sampling_result_2, variable = "weight_pref_ns_m")
weight_male_age_pref_posterior_summary <- summarize_draws(weight_male_age_pref_posterior,
                                                     rhat,           
                                                     ess_bulk,        
                                                     ess_tail,        
                                                     mean,
                                                     median,
                                                     ~quantile(.x, probs = c(.5, 0.025, 0.975)))

weight_male_age_pref_posterior_summary$`50%`%>%matrix(nrow=47, ncol=6)-> weight_male_age_pref_median_df
weight_male_age_pref_posterior_summary$`2.5%`%>%matrix(nrow=47, ncol=6)-> weight_male_age_pref_lwr_df
weight_male_age_pref_posterior_summary$`97.5%`%>%matrix(nrow=47, ncol=6)-> weight_male_age_pref_upr_df

male_age_distr_pref = list()

for(col in 1:6){
  tmp_df <- matrix(NA, nrow = 47, ncol = 4)%>%data.frame()
  colnames(tmp_df) <- c("pref_idx", "median", "lwr_95CrI", "upr_95CrI")
  
  tmp_df[,1] <- 1:47
  tmp_df[,2] <- weight_male_age_pref_median_df[, col]
  tmp_df[,3] <- weight_male_age_pref_lwr_df[, col]
  tmp_df[,4] <- weight_male_age_pref_upr_df[, col]
  
  DescTools::MultinomCI(df_pop_pref_m[,col], conf.level = .95, method = "goodman") -> ci_census
  
  tmp_df%<>%
    mutate(
      census_mean = ci_census[,1],
      census_lwr95 = ci_census[,2],
      census_upr95 = ci_census[,3])%>%
    mutate(age_group = paste0((col+1)*10, "-", (col+1)*10+9), .after = pref_idx)%>%
    mutate(overlap = case_when(
      lwr_95CrI <= census_upr95 & census_lwr95 <= upr_95CrI ~ T, 
      T ~ F
    ))
  
  male_age_distr_pref[[paste0("age", (col+1)*10)]] <- tmp_df
}

male_pref_df <- male_age_distr_pref[[paste0("age", (1+1)*10)]]
for(col in 2:6){
  male_pref_df <- rbind(male_pref_df, male_age_distr_pref[[paste0("age", (col+1)*10)]])
}

male_pref_df%>%readr::write_csv(file = "path_csv")


#######################
#######################

weight_female_age_pref_posterior <- subset_draws(sampling_result_2, variable = "weight_pref_ns_f")
weight_female_age_pref_posterior_summary <- summarize_draws(weight_female_age_pref_posterior,
                                                          rhat,           
                                                          ess_bulk,        
                                                          ess_tail,        
                                                          mean,
                                                          median,
                                                          ~quantile(.x, probs = c(.5, 0.025, 0.975)))

weight_female_age_pref_posterior_summary$`50%`%>%matrix(nrow=47, ncol=6)-> weight_female_age_pref_median_df
weight_female_age_pref_posterior_summary$`2.5%`%>%matrix(nrow=47, ncol=6)-> weight_female_age_pref_lwr_df
weight_female_age_pref_posterior_summary$`97.5%`%>%matrix(nrow=47, ncol=6)-> weight_female_age_pref_upr_df

female_age_distr_pref = list()

for(col in 1:6){
  tmp_df <- matrix(NA, nrow = 47, ncol = 4)%>%data.frame()
  colnames(tmp_df) <- c("pref_idx", "median", "lwr_95CrI", "upr_95CrI")
  
  tmp_df[,1] <- 1:47
  tmp_df[,2] <- weight_female_age_pref_median_df[, col]
  tmp_df[,3] <- weight_female_age_pref_lwr_df[, col]
  tmp_df[,4] <- weight_female_age_pref_upr_df[, col]
  
  DescTools::MultinomCI(df_pop_pref_f[,col], conf.level = .95, method = "goodman") -> ci_census
  
  tmp_df%<>%
    mutate(
      census_mean = ci_census[,1],
      census_lwr95 = ci_census[,2],
      census_upr95 = ci_census[,3])%>%
    mutate(age_group = paste0((col+1)*10, "-", (col+1)*10+9), .after = pref_idx)%>%
    mutate(overlap = case_when(
      lwr_95CrI <= census_upr95 & census_lwr95 <= upr_95CrI ~ T, 
      T ~ F
    ))
  
  female_age_distr_pref[[paste0("age", (col+1)*10)]] <- tmp_df
}

female_pref_df <- female_age_distr_pref[[paste0("age", (1+1)*10)]]
for(col in 2:6){
  female_pref_df <- rbind(female_pref_df, female_age_distr_pref[[paste0("age", (col+1)*10)]])
}

female_pref_df%>%readr::write_csv(file = "path_csv")


############
############ jobs
############
weight_male_job_posterior <- subset_draws(sampling_result_2, variable = "weight_job_ns_m")
weight_male_job_posterior_summary <- summarize_draws(weight_male_job_posterior,
                                                          rhat,           
                                                          ess_bulk,        
                                                          ess_tail,        
                                                          mean,
                                                          median,
                                                          ~quantile(.x, probs = c(.5, 0.025, 0.975)))

weight_male_job_posterior_summary$`50%`%>%matrix(nrow=21, ncol=6)-> weight_male_job_median_df
weight_male_job_posterior_summary$`2.5%`%>%matrix(nrow=21, ncol=6)-> weight_male_job_lwr_df
weight_male_job_posterior_summary$`97.5%`%>%matrix(nrow=21, ncol=6)-> weight_male_job_upr_df


male_job_distr = list()

for(col in 1:6){
  # col=1
  
  tmp_df <- matrix(NA, nrow = 21, ncol = 4)%>%data.frame()
  colnames(tmp_df) <- c("job_category", "median", "lwr_95CrI", "upr_95CrI")
  
  tmp_df[,1] <- 1:21
  tmp_df[,2] <- weight_male_job_median_df[, col]
  tmp_df[,3] <- weight_male_job_lwr_df[, col]
  tmp_df[,4] <- weight_male_job_upr_df[, col]
  
  DescTools::MultinomCI(df_job_m[,col], conf.level = .95, method = "goodman") -> ci_census
  
  tmp_df%<>%
    mutate(
      census_mean = ci_census[,1],
      census_lwr95 = ci_census[,2],
      census_upr95 = ci_census[,3])%>%
    mutate(age_group = paste0((col+1)*10, "-", (col+1)*10+9), .after = job_category)%>%
    mutate(overlap = case_when(
      lwr_95CrI <= census_upr95 & census_lwr95 <= upr_95CrI ~ T, 
      T ~ F
    ))
  
  
  male_job_distr[[paste0("age", (col+1)*10)]] <- tmp_df
}

male_job_distr_df <- male_job_distr[[paste0("age", (1+1)*10)]]
for(col in 2:6){
  male_job_distr_df <- rbind(male_job_distr_df, male_job_distr[[paste0("age", (col+1)*10)]])
}

male_job_distr_df%>%readr::write_csv(file = "path_csv")


##########
weight_female_job_posterior <- subset_draws(sampling_result_2, variable = "weight_job_ns_f")
weight_female_job_posterior_summary <- summarize_draws(weight_female_job_posterior,
                                                     rhat,           
                                                     ess_bulk,        
                                                     ess_tail,        
                                                     mean,
                                                     median,
                                                     ~quantile(.x, probs = c(.5, 0.025, 0.975)))

weight_female_job_posterior_summary$`50%`%>%matrix(nrow=21, ncol=6)-> weight_female_job_median_df
weight_female_job_posterior_summary$`2.5%`%>%matrix(nrow=21, ncol=6)-> weight_female_job_lwr_df
weight_female_job_posterior_summary$`97.5%`%>%matrix(nrow=21, ncol=6)-> weight_female_job_upr_df


female_job_distr = list()

for(col in 1:6){
  tmp_df <- matrix(NA, nrow = 21, ncol = 4)%>%data.frame()
  colnames(tmp_df) <- c("job_category", "median", "lwr_95CrI", "upr_95CrI")
  
  tmp_df[,1] <- 1:21
  tmp_df[,2] <- weight_female_job_median_df[, col]
  tmp_df[,3] <- weight_female_job_lwr_df[, col]
  tmp_df[,4] <- weight_female_job_upr_df[, col]
  
  DescTools::MultinomCI(df_job_f[,col], conf.level = .95, method = "goodman") -> ci_census
  
  tmp_df%<>%
    mutate(
      census_mean = ci_census[,1],
      census_lwr95 = ci_census[,2],
      census_upr95 = ci_census[,3])%>%
    mutate(age_group = paste0((col+1)*10, "-", (col+1)*10+9), .after = job_category)%>%
    mutate(overlap = case_when(
      lwr_95CrI <= census_upr95 & census_lwr95 <= upr_95CrI ~ T, 
      T ~ F
    ))
  

  female_job_distr[[paste0("age", (col+1)*10)]] <- tmp_df
}


female_job_distr_df <- female_job_distr[[paste0("age", (1+1)*10)]]
for(col in 2:6){
  female_job_distr_df <- rbind(female_job_distr_df, female_job_distr[[paste0("age", (col+1)*10)]])
}

female_job_distr_df%>%readr::write_csv(file = "path_csv")
