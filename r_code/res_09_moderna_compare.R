rm(list=ls())
library(posterior);library(tidyr);library(ggplot2)
df_analysis <- readr::read_csv(file = "./CODES/data_questionnaire.csv")%>%
  filter(infect_period2!=1)
df_analysis%<>%
  mutate(ID = 1:nrow(.))

load(file = "path_to_RData") # load cmdstan fit object


prob_infection <- subset_draws(sampling_result_2, variable = "prob_infection")
weight_adj <- subset_draws(sampling_result_2, variable = "alpha")

set.seed(99)
W     <- as_draws_matrix(weight_adj)
prob_mat <- as_draws_matrix(prob_infection)
Urand      <- matrix(runif(length(prob_mat)), nrow = nrow(prob_mat))
Y_rep  <- (Urand < prob_mat) * 1L  

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
  med_20 = 0,   lwr_20 = 0,  upr_20 = 0,
  med_30 = 0,   lwr_30 = 0,  upr_30 = 0,
  med_40 = 0,   lwr_40 = 0,  upr_40 = 0,
  med_50 = 0,   lwr_50 = 0,  upr_50 = 0,
  med_60 = 0,   lwr_60 = 0,  upr_60 = 0,
  med_70 = 0,   lwr_70 = 0,  upr_70 = 0
)

pos_vec[ID_20]%>%as.numeric()-> pos_vec_20
pos_vec[ID_30]%>%as.numeric()-> pos_vec_30
pos_vec[ID_40]%>%as.numeric()-> pos_vec_40
pos_vec[ID_50]%>%as.numeric()-> pos_vec_50
pos_vec[ID_60]%>%as.numeric()-> pos_vec_60
pos_vec[ID_70]%>%as.numeric()-> pos_vec_70

W_20_normalized %*% pos_vec_20 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[1, 1:3]
W_30_normalized %*% pos_vec_30 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[1, 4:6]
W_40_normalized %*% pos_vec_40 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[1, 7:9]
W_50_normalized %*% pos_vec_50 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[1, 10:12]
W_60_normalized %*% pos_vec_60 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[1, 13:15]
W_70_normalized %*% pos_vec_70 %>% quantile(probs = c(.5, .025, .975)) -> df_res_all[1, 16:18]


counts_20_59 <- 
  (W_20_normalized %*% pos_vec_20)*12731*1e3 +
  (W_30_normalized %*% pos_vec_30)*13351*1e3 +
  (W_40_normalized %*% pos_vec_40)*16702*1e3 +
  (W_50_normalized %*% pos_vec_50)*18053*1e3 
# https://www.stat.go.jp/data/jinsui/2024np/index.html

counts_20_59%>% quantile(probs = c(.5, .025, .975)) 

counts_60over <- 
  (W_60_normalized %*% pos_vec_60)*14804*1e3  +
  (W_70_normalized %*% pos_vec_70)*28917*1e3
counts_60over%>% quantile(probs = c(.5, .025, .975)) 
