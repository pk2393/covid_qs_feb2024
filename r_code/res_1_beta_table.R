rm(list=ls())
library(dplyr);library(magrittr);library(lubridate); library(posterior)

load(file = "path_to_RData") # load cmdstan fit object

beta_risk_draws <- subset_draws(sampling_result_2, variable = "beta_risk_background")

beta_risk_summary <- summarize_draws(beta_risk_draws,
                                               mean,
                                               median,
                                               ~quantile(.x, probs = c(.5, 0.025, 0.975)))

beta_risk_summary_exp <- beta_risk_summary %>%
  mutate(
    mean = exp(mean),
    median = exp(median),
    `50%` = exp(`50%`),
    `2.5%` = exp(`2.5%`),
    `97.5%` = exp(`97.5%`)
  )

beta_risk_summary_exp %>% print(n = 50)

beta_risk_summary_formatted <- beta_risk_summary_exp %>%
  mutate(
    formatted = paste0(round(median, 3), " (", round(`2.5%`, 3), "-", round(`97.5%`, 3), ")")
  )

beta_risk_summary_formatted %>%readr::write_csv(file = "path_csv") 

