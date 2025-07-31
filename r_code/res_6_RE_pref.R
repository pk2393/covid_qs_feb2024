rm(list=ls())
library(posterior);library(tidyr);library(ggplot2)
library(patchwork);library(pammtools)

df_analysis <- readr::read_csv(file = "./CODES/data_questionnaire.csv")%>%
  filter(infect_period2!=1)%>%
  mutate(ID = 1:nrow(.))

load(file = "path_to_RData") # load cmdstan fit object

pref_info <-readr::read_csv("pref_info.csv")

RE_pref_log <- subset_draws(sampling_result_2, variable = "RE_pref_log")

RE_pref_log_long <- as_draws_df(RE_pref_log) %>%
  pivot_longer(cols = everything(), names_to = "prefecture", values_to = "value") %>%
  mutate(prefecture = factor(prefecture, levels = paste0("RE_pref_log[", 1:47, "]"))) %>%
  filter(is.na(prefecture) ==F)%>%
  mutate(
    pref_index = as.numeric(gsub("RE_pref_log\\[|\\]", "", prefecture)), 
    pref_name = stringr::str_to_title(pref_info$pref_list[pref_index])                        
  )

palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


ggplot(RE_pref_log_long, aes(x = prefecture, y = value, fill = prefecture)) +
  geom_violin(fill = palette[2], trim = TRUE, scale = "area", alpha = 0.5, 
              position = position_dodge(width = 2), linewidth = 0.05) +
  stat_summary(
    fun.data = function(x) {
      return(data.frame(
        y = median(x),
        ymin = quantile(x, 0.025),
        ymax = quantile(x, 0.975)
      ))
    },
    geom = "pointrange",
    color = palette[9], 
    size = .5
  ) +
  labs(
    title = "Effect by Prefecture on the Risk of Infection",
    x = "Prefecture",
    y = "Posterior Distribution"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    legend.position = "none"
  )+
  scale_x_discrete(
    labels = RE_pref_log_long$pref_name
  )+
  theme_minimal() + 
  theme(
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = +.4, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 12), 
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none" 
  ) +
  scale_fill_manual(values = rep("lightgray", 47)) + 
  scale_y_continuous(expand = c(0.25, 0.25)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") -> plot_RE_pref

plot_RE_pref

# tiff("path_tiff", width = 18, height = 6, units = "in", res = 1000, compression = "lzw")
# plot_RE_pref
# dev.off()
