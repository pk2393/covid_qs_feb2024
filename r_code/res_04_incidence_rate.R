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

prob_infection <- subset_draws(sampling_result_2, variable = "prob_infection")
weight_adj <- subset_draws(sampling_result_2, variable = "alpha")

W     <- as_draws_matrix(weight_adj)
prob_mat <- as_draws_matrix(prob_infection)

############################
df_analysis%<>%
  mutate(ID = 1:nrow(.))

df_analysis%>%
  pull(pos_feb2024)-> pos_vec

####
#### national, main result
####
df_analysis%>%
  filter(age_group==2)%>%pull(ID) -> ID_20

df_analysis%>%
  filter(age_group==3)%>%pull(ID) -> ID_30

df_analysis%>%
  filter(age_group==4)%>%pull(ID) -> ID_40

df_analysis%>%
  filter(age_group==5)%>%pull(ID) -> ID_50

df_analysis%>%
  filter(age_group==6)%>%pull(ID) -> ID_60

df_analysis%>%
  filter(age_group==7)%>%pull(ID) -> ID_70

W_20 <- W[, ID_20, drop=F]
W_20_normalized <- sweep(W_20, 1, rowSums(W_20), "/")
W_30 <- W[, ID_30, drop=F]
W_30_normalized <- sweep(W_30, 1, rowSums(W_30), "/")
W_40 <- W[, ID_40, drop=F]
W_40_normalized <- sweep(W_40, 1, rowSums(W_40), "/")
W_50 <- W[, ID_50, drop=F]
W_50_normalized <- sweep(W_50, 1, rowSums(W_50), "/")
W_60 <- W[, ID_60, drop=F]
W_60_normalized <- sweep(W_60, 1, rowSums(W_60), "/")
W_70 <- W[, ID_70, drop=F]
W_70_normalized <- sweep(W_70, 1, rowSums(W_70), "/")

df_res_all = data.frame(
  ag = c("20-29", "30-39", "40-49", "50-59", "60-69", "70+"), 
  Median = 0, Lower = 0,   Upper = 0
)

pos_vec[ID_20]%>%as.numeric()-> pos_vec_20
pos_vec[ID_30]%>%as.numeric()-> pos_vec_30
pos_vec[ID_40]%>%as.numeric()-> pos_vec_40
pos_vec[ID_50]%>%as.numeric()-> pos_vec_50
pos_vec[ID_60]%>%as.numeric()-> pos_vec_60
pos_vec[ID_70]%>%as.numeric()-> pos_vec_70

W_20_normalized %*% pos_vec_20 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[1, 2:4]
W_30_normalized %*% pos_vec_30 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[2, 2:4]
W_40_normalized %*% pos_vec_40 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[3, 2:4]
W_50_normalized %*% pos_vec_50 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[4, 2:4]
W_60_normalized %*% pos_vec_60 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[5, 2:4]
W_70_normalized %*% pos_vec_70 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[6, 2:4]


gg <- df_res_all %>%
  ggplot(aes(x = ag, y = Median*100)) +
  geom_point(size = 3, color = "blue") +  
  geom_errorbar(aes(ymin = Lower*100, ymax = Upper*100), width = 0.2, color = "blue") +  
  labs(
    title = "Incidence Rates by Age Group with 95% CrI",
    x = "Age Group",
    y = "Incidence Rate (%)"
  ) +
  scale_y_continuous(
    limits = c(0, 7),  
    breaks = scales::pretty_breaks(n = 10)  
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.ticks.x = element_line(color = "black", size = 1), 
    axis.ticks.y = element_line(color = "black", size = 1), 
    axis.text.x = element_text(size = 16, hjust = 0.5),
    axis.text.y = element_text(size = 16),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_blank(),      
    axis.line.x = element_line(color = "black"),  
    axis.line.y = element_line(color = "black")   
  )+geom_text(
    aes(label = paste0(round(Median*100, 2), " (", round(Lower*100, 2), " - ", round(Upper*100, 2), ")")),
    vjust = -4, size = 4, color = "black"  
  )
gg


tiff("path_tiff", width = 9, height = 6, units = "in", res = 1000, compression = "lzw")
gg
dev.off()
