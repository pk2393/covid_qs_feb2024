rm(list=ls())
library(dplyr);library(magrittr);library(lubridate); library(posterior)

load(file = "path_to_RData") # load cmdstan fit object

###
### Plots
###
type_decay=c("decay_infect_xbb", "decay_infect_preX", "decay_infect_wuhan", "decay_vax_xbb", "decay_vax_bivalent", "decay_vax_wuhan")

title_decay = c("Infection, XBB.1.5", "Infection, Omicron pre-XBB", "Infection, pre-Omicron", 
                "Vaccination, XBB.1.5", "Vaccination, Bivalent", "Vaccination, Wuhan")

decay_list=list()
decay_data_list = list()

for(i in 1:6){
  mcmc_draws <- subset_draws(sampling_result_2, variable = paste0(type_decay[i]))
  variables <- dimnames(mcmc_draws)[[3]]
  
  quantile_vec <- c(0.025, 0.5, 0.975)
  
  quantile_data<-data.frame()
  
  for (j in 1:length(variables)){
    var<-variables[j]
    var_data <- as.vector(mcmc_draws[, , var])
    quantiles <- quantile(var_data, probs = quantile_vec)
    
    quantile_data <- rbind(quantile_data, 
                           data.frame(Time = j-1,
                                      Variable = var, 
                                      Lower = quantiles[1],
                                      Median = quantiles[2],
                                      Upper = quantiles[3]))
  }
  
  decay_data_list[[paste0(type_decay[i])]] <- quantile_data
  
  
  ggplot(quantile_data, aes(x = Time, y = Median*100)) +
    geom_ribbon(aes(ymin = Lower*100, ymax = Upper*100), fill = "red", alpha = 0.2) +
    geom_line(size = 1, color= "red") +
    theme_minimal() +
    labs(x = "Time from immunization (days)", y = "RR reduction (%)", title = paste0(title_decay[i])) +
    scale_y_continuous(limits = c(0, 100)) +  
    theme(
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      axis.text = element_text(size = 12),
      axis.line = element_line(color = "black", size = 0.5),  
      axis.ticks = element_line(color = "black", size = 0.5), 
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16)
    )->decay_list[[paste0(type_decay[i])]]
}

gridExtra::grid.arrange(decay_list[[paste0(type_decay[1])]], decay_list[[paste0(type_decay[2])]], 
                        decay_list[[paste0(type_decay[3])]], decay_list[[paste0(type_decay[4])]], 
                        decay_list[[paste0(type_decay[5])]], decay_list[[paste0(type_decay[6])]], 
                        ncol=3)



tiff("path_fig", width = 12, height = 6, units = "in", res = 1000, compression = "lzw")
gridExtra::grid.arrange(decay_list[[paste0(type_decay[1])]], decay_list[[paste0(type_decay[2])]], 
                        decay_list[[paste0(type_decay[3])]], decay_list[[paste0(type_decay[4])]], 
                        decay_list[[paste0(type_decay[5])]], decay_list[[paste0(type_decay[6])]], 
                        ncol=3)
dev.off()



###
### Table
###
decay_first_mo <- matrix(NA, nrow=6, ncol = 4)%>%data.frame()
colnames(decay_first_mo) <- c("type", "lwr_95CrI", "median", "upr_95CrI")
decay_first_mo[,1] <- c("infect_xbb_1mo", "infect_preX_1mo", "infect_wuhan_1mo", 
                        "vax_xbb_1mo", "vax_bivalent_1mo", "vax_wuhan_1mo")
# type_decay

for(i in 1:6){
  decay_data_list[[paste0(type_decay[i])]]%>%
    filter(Time == 30) ->tmp_range

  tmp_range[3:5]*100 -> decay_first_mo[i, 2:4]
}

decay_first_mo[, 2:4]<- round(decay_first_mo[, 2:4], 1)

decay_first_mo$estimate <- apply(decay_first_mo, 1, function(row) {
  sprintf("%.1f (%.1f, %.1f)", 
          round(as.numeric(row["median"]), 1), 
          round(as.numeric(row["lwr_95CrI"]), 1), 
          round(as.numeric(row["upr_95CrI"]), 1))
}
)

decay_first_mo %>% readr::write_csv(file = "path_csv")


