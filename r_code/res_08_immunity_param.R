rm(list=ls())
library(dplyr);library(magrittr);library(posterior)
load(file = "path_to_RData") # load cmdstan fit object

variable_list = c("v_infect_xbb", "v_infect_preX", "v_infect_wuhan", 
                  "v_vax_xbb", "v_vax_bivalent", "v_vax_wuhan", 
                  "h1", "h2", "frac_short_infect", "frac_short_vax"
)

imm_table <- matrix(NA, nrow = length(variable_list), ncol = 4)%>%data.frame()

colnames(imm_table) <- c("variable", "median", "lwr_95CrI", "upr_95CrI")
imm_table[,1] <- variable_list

for(n in 1:length(variable_list)){
  # n=1
  for(n in 1:6){
    subset_draws(sampling_result_2, variable = variable_list[n])%>%quantile(c(.5, .025, .975))-> imm_table[n, 2:4]
    
    imm_table[n, 2:4] <- imm_table[n, 2:4]*100
  }
  for(n in 7:8){
    subset_draws(sampling_result_2, variable = variable_list[n])%>%quantile(c(.5, .025, .975)) -> imm_table[n, 2:4]
  }
  
  for(n in 9:10){
    subset_draws(sampling_result_2, variable = variable_list[n])%>%quantile(c(.5, .025, .975)) -> imm_table[n, 2:4]
    
    imm_table[n, 2:4] <- imm_table[n, 2:4]*100
  }
}


imm_table$estimate <- apply(imm_table, 1, function(row) {
  sprintf("%.1f (%.1f, %.1f)", 
          round(as.numeric(row["median"]), 1), 
          round(as.numeric(row["lwr_95CrI"]), 1), 
          round(as.numeric(row["upr_95CrI"]), 1))
  }
)

imm_table %>% readr::write_csv(file = "path_csv")



