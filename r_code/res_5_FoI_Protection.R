rm(list=ls())
library(posterior);library(tidyr);library(ggplot2)
library(patchwork);library(pammtools)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE) 
brewerPalette <- brewer.pal(9, "YlOrRd")


df_analysis <- readr::read_csv(file = "./CODES/data_questionnaire.csv")%>%
  filter(infect_period2!=1)
df_analysis%<>%
  mutate(ID = 1:nrow(.))

load(file = "path_to_RData") # load cmdstan fit object

protect_20<-protect_30<-protect_40<-
  protect_50<-protect_60<-protect_70<-data.frame(index=1:20)

lambda_20<-lambda_30<-lambda_40<-
  lambda_50<-lambda_60<-lambda_70<-data.frame(index=1:20)

protection_individual <- subset_draws(sampling_result_2, variable = "protection_individual")
weight_adj <- subset_draws(sampling_result_2, variable = "alpha")
lambda <- subset_draws(sampling_result_2, variable = "lambda")


bins_p <- seq(0, 1, length.out = 21)       # 20 bins for protection
bins_l <- seq(0, 0.2, length.out = 21)     # 20 bins for lambda
n_p    <- length(bins_p) - 1
n_l    <- length(bins_l) - 1

## posterior draws
Z_p   <- as_draws_matrix(protection_individual)
Z_l   <- as_draws_matrix(lambda)
W     <- as_draws_matrix(weight_adj) 
niter <- nrow(Z_p)



#########
quantile_vec <- c(0.5, 0.25, 0.75, 0.1, 0.9, 0.025, 0.975)
ag_cols_template <- data.frame(
  fit = 0, lwr50 = 0, upr50 = 0, lwr80 = 0, upr80 = 0, lwr95 = 0, upr95 = 0
)
# Protection
protect_aggr <- data.frame(bin_lwr = bins_p[1:n_p], bin_upr = bins_p[2:(n_p+1)])
for (age in c(20, 30, 40, 50, 60, 70)) {
  protect_aggr <- cbind(protect_aggr, setNames(ag_cols_template, paste0("ag", age, "_", names(ag_cols_template))))
}

lambda_aggr <- data.frame(bin_lwr = bins_l[1:n_l], bin_upr = bins_l[2:(n_l+1)])
for (age in c(20, 30, 40, 50, 60, 70)) {
  lambda_aggr <- cbind(lambda_aggr, setNames(ag_cols_template, paste0("ag", age, "_", names(ag_cols_template))))
}

array_list_all <- list()
array_list_no_expos <- list()
array_list_last_infect <- list()
array_list_last_vaccine <- list()

for(ag in 1:6){
  # ag=1
  id_all         <- df_analysis$ID[df_analysis$age_group == ag+1]
  id_no_expos    <- df_analysis$ID[df_analysis$age_group == ag+1 & df_analysis$imm_last == 0]
  id_last_infect <- df_analysis$ID[df_analysis$age_group == ag+1 & df_analysis$imm_last == 1]
  id_last_vax    <- df_analysis$ID[df_analysis$age_group == ag+1 & df_analysis$imm_last == 2]
  
  strata_data <- list(
    all = list(Z_p = Z_p[, id_all], Z_l = Z_l[, id_all], W = W[, id_all]),
    no_expos = list(Z_p = Z_p[, id_no_expos], Z_l = Z_l[, id_no_expos], W = W[, id_no_expos]),
    last_infect = list(Z_p = Z_p[, id_last_infect], Z_l = Z_l[, id_last_infect], W = W[, id_last_infect]),
    last_vax = list(Z_p = Z_p[, id_last_vax], Z_l = Z_l[, id_last_vax], W = W[, id_last_vax])
  )
  
  for(stratum_name in names(strata_data)) {
    current_data <- strata_data[[stratum_name]]
    if(ncol(current_data$W) == 0) next 

    bin_with_cut <- function(x, breaks) {
      as.integer(cut(x, breaks = breaks, include.lowest = TRUE, right = TRUE))
    }
    prot_idx   <- t(apply(current_data$Z_p, 1, bin_with_cut, breaks = seq(0, 1, length.out = 21)))
    lambda_idx <- t(apply(current_data$Z_l, 1, bin_with_cut, breaks = c(seq(0, 0.2, by = 0.01), Inf)))
    lambda_idx[lambda_idx > 20] <- 20L
    
    arr_2d <- array(0, dim = c(n_p, n_l, niter))
    for (it in seq_len(niter)) {

      df_it <- data.frame(
        w = current_data$W[it, ]%>%unlist()%>%as.vector(),
        p = prot_idx[it, ]%>%unlist()%>%as.vector(),
        l = lambda_idx[it, ]%>%unlist()%>%as.vector()
      )
      
      tab <- xtabs(w ~ p + l, data = df_it)
      mat <- matrix(0, nrow = n_p, ncol = n_l)
      mat[as.integer(rownames(tab)), as.integer(colnames(tab))] <- as.vector(tab)
      arr_2d[, , it] <- mat
    }
    
    if (stratum_name == "all") array_list_all[[paste0("ag", ag+1)]] <- arr_2d
    if (stratum_name == "no_expos") array_list_no_expos[[paste0("ag", ag+1)]] <- arr_2d
    if (stratum_name == "last_infect") array_list_last_infect[[paste0("ag", ag+1)]] <- arr_2d
    if (stratum_name == "last_vax") array_list_last_vaccine[[paste0("ag", ag+1)]] <- arr_2d
    
    if (stratum_name == "all") {
      p_dist_iter <- apply(arr_2d, c(1, 3), sum)
      l_dist_iter <- apply(arr_2d, c(2, 3), sum) 
      
      p_dist_norm <- sweep(p_dist_iter, 2, colSums(p_dist_iter), "/")
      l_dist_norm <- sweep(l_dist_iter, 2, colSums(l_dist_iter), "/")

      q_protect <- t(apply(p_dist_norm, 1, quantile, probs = quantile_vec, na.rm = TRUE))
      q_lambda  <- t(apply(l_dist_norm, 1, quantile, probs = quantile_vec, na.rm = TRUE))
      
      start_col <- 2 + (ag - 1) * 7 + 1
      protect_aggr[, start_col:(start_col + 6)] <- q_protect
      lambda_aggr[, start_col:(start_col + 6)]  <- q_lambda
    }
  } 
}

#### Protection 1D
protect_long <- protect_aggr %>%
  pivot_longer(cols = starts_with("ag"), names_to = c("AgeGroup", ".value"), 
               names_pattern = "ag(\\d+)_(fit|lwr50|upr50|lwr80|upr80|lwr95|upr95)") %>%
  mutate(Range = paste0("(", round(bin_lwr, 3), "-", round(bin_upr, 3), ")"))

plot_protect_list <- list()
age_groups <- unique(protect_long$AgeGroup)
age_groups_upper = as.character(c(paste0("-", as.numeric(age_groups[1:5])+rep(9, 5)), "+"))

for (i in 1:length(age_groups)) {
  age_group = age_groups[i]
  age_group_upper=age_groups_upper[i]
  
  df <- subset(protect_long, AgeGroup == age_group)
  
  plot <- ggplot(df, aes(x = bin_lwr, y = fit*100)) +
    geom_step(aes(x = bin_lwr, y = fit*100), color = brewerPalette[9], size = 1) +  
    pammtools::geom_stepribbon(aes(ymin = lwr95*100, ymax = upr95*100), fill = brewerPalette[4], alpha = 0.2) + 
    pammtools::geom_stepribbon(aes(ymin = lwr80*100, ymax = upr80*100), fill = brewerPalette[6], alpha = 0.3) + 
    pammtools::geom_stepribbon(aes(ymin = lwr50*100, ymax = upr50*100), fill = brewerPalette[9], alpha = 0.4) +  
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent(seq(0, 1, by = 0.1))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + 
    coord_cartesian(ylim=c(0, 75))+
    theme_minimal() +
    labs(title = paste(age_group, age_group_upper), x = "Relative Risk Reduction (%)", y = "Proportion (%)") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.ticks = element_line(),
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5)
    )
  
  plot_protect_list[[i]] <- plot
}

final_protect_plot <- wrap_plots(plot_protect_list, ncol = 3)+ 
  plot_annotation(title = "Immune Protection by Age Group") & 
  theme(plot.title = element_text(size = 16, face = "bold"))
final_protect_plot

# tiff("path_prot_tiff", width = 12, height = 6, units = "in", res = 1000, compression = "lzw")
# final_protect_plot
# dev.off()

###
### FoI 1D
###
lambda_long <- lambda_aggr %>%
  pivot_longer(cols = starts_with("ag"), names_to = c("AgeGroup", ".value"), 
               names_pattern = "ag(\\d+)_(fit|lwr50|upr50|lwr80|upr80|lwr95|upr95)") %>%
  mutate(Range = paste0("(", round(bin_lwr, 3), "-", round(bin_upr, 3), ")"))

plot_lambda_list <- list()
age_groups <- unique(lambda_long$AgeGroup)
age_groups_upper = as.character(c(paste0("-", as.numeric(age_groups[1:5])+rep(9, 5)), "+"))

for (i in 1:length(age_groups)) {
  age_group = age_groups[i]
  age_group_upper=age_groups_upper[i]
  
  df <- subset(lambda_long, AgeGroup == age_group)
  
  plot <- ggplot(df, aes(x = bin_lwr, y = fit*100)) +
    geom_step(aes(x = bin_lwr, y = fit*100), color = brewerPalette[9], size = 1) +  
    pammtools::geom_stepribbon(aes(ymin = lwr95*100, ymax = upr95*100), fill = brewerPalette[4], alpha = 0.2) + 
    pammtools::geom_stepribbon(aes(ymin = lwr80*100, ymax = upr80*100), fill = brewerPalette[6], alpha = 0.3) + 
    pammtools::geom_stepribbon(aes(ymin = lwr50*100, ymax = upr50*100), fill = brewerPalette[9], alpha = 0.4) +  
    scale_x_continuous(breaks = seq(0, 0.2, by = 0.02)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + 
    coord_cartesian(ylim=c(0, 75))+
    # ylim(0, 60) +
    theme_minimal() +
    labs(title = paste(age_group, age_group_upper), x = "Force of Infection", y = "Proportion (%)") +
    theme(
      legend.position = "none", 
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.ticks = element_line(),
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5)
    )
  
  plot_lambda_list[[i]] <- plot
}

final_lambda_plot <- wrap_plots(plot_lambda_list, ncol = 3)+ 
  plot_annotation(title = "Force of Infection by Age Group") & 
  theme(plot.title = element_text(size = 16, face = "bold"))
final_lambda_plot

# tiff("path_lambda_tiff", width = 12, height = 6, units = "in", res = 1000, compression = "lzw")
# final_lambda_plot
# dev.off()


###
### 2D heatmap
###

library(reshape2); library(viridis); library(gridExtra)

bins_p <- seq(0, 1, length.out = 21)[1:20] * 100    # 0,5,10,…,95
bins_l <- seq(0, 0.2, length.out = 21)[1:20]        # 0.00,0.01,…,0.19

age_vec  <- 1:6
ag_title <- c("20-29","30-39","40-49","50-59","60-69","70+")

get_arr <- function(ag, nm) {
  key <- paste0("ag", ag+1)
  switch(nm,
         all          = array_list_all[[key]],
         no_expos     = array_list_no_expos[[key]],
         last_infect  = array_list_last_infect[[key]],
         last_vaccine = array_list_last_vaccine[[key]]
  )
}
denom_age <- sapply(age_vec, function(ag) sum(get_arr(ag, "all")))


all_density <- unlist(lapply(age_vec, function(ag) {
  mat_all <- apply(get_arr(ag, "all"), c(1,2), sum)
  as.vector(mat_all / denom_age[ag]) 
}))

threshold_val <- quantile(all_density, 0.975, na.rm = TRUE)
global_min    <- min(all_density, na.rm = TRUE)
threshold_val <- 0.012

heatmap_list <- vector("list", length(age_vec))
heatmap_list_imm_last <- vector("list", length(age_vec) * 3)

name_list      <- c("all","no_expos","last_infect","last_vaccine")
title_2nd_list <- c("", "None", "Infection", "Vaccine")

for (i in seq_along(age_vec)) {
  ag <- age_vec[i]
  denom <- denom_age[ag]
  
  for (n_idx in seq_along(name_list)) {
    arr  <- get_arr(ag, name_list[n_idx])
    mat  <- apply(arr, c(1,2), sum)
    dens <- mat / denom  
    df   <- reshape2::melt(dens)
    
    names(df) <- c("ProtBin","LamBin","Density")
    df$Protection     <- bins_p[df$ProtBin]
    df$Lambda         <- bins_l[df$LamBin]
    df$Density_Capped <- pmin(df$Density, threshold_val)
    
    p <- ggplot(df, aes(x = Protection + 2.5, y = Lambda + 0.005, fill = Density_Capped)) +
      geom_tile(width = 5, height = 0.01) +
      scale_fill_viridis_c(option = "inferno", name = "Density",
                           limits = c(global_min, threshold_val)) +
      scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100), expand = c(0,0)) +
      scale_y_continuous(breaks = seq(0, 0.2, 0.05), limits = c(0, 0.22), expand = c(0,0)) +
      labs(title = paste0(ag_title[i], " ", title_2nd_list[n_idx]),
           x = "Relative Risk Reduction (%)", y = "Force of Infection") +
      theme_minimal() +
      theme(plot.title = element_text(size = 24, face = "bold"),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x  = element_text(size = 12, hjust = 1),
            axis.text.y  = element_text(size = 12),
            panel.grid = element_blank(),
            axis.ticks = element_line(),
            axis.ticks.length = unit(0.2, "cm")) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
    
    if (name_list[n_idx] == "all") {
      heatmap_list[[i]] <- p
    } else {
      heatmap_list_imm_last[[3*(i-1) + (n_idx-1)]] <- p
    }
  }
}

# tiff("path_lambda_protection_heatmap_tiff",
#      width = 18, height = 12, units = "in",
#      res = 1000, compression = "lzw")
# grid.arrange(
#   grobs = heatmap_list,
#   ncol  = 3, nrow = 2
# )
# dev.off()

# tiff("path_heatmap_by_last_imm_tiff", 
#      width = 12, height = 18, units = "in", res = 1000, compression = "lzw")
# grid.arrange(
#   grobs = heatmap_list_imm_last,
#   ncol  = 3, nrow = 6
# )
# dev.off()
