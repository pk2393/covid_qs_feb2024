rm(list=ls())
library(posterior);library(tidyr);library(ggplot2)

df_analysis <- readr::read_csv(file = "./CODES/data_questionnaire.csv")%>%
  filter(infect_period2!=1)
df_analysis%<>%
  mutate(ID = 1:nrow(.))

load(file = "path_to_RData") # load cmdstan fit object

library(posterior)
df_analysis%<>%mutate(ID_new = 1:nrow(.))

df_analysis%>%
  select(pos_feb2024, infection_history, infection_last, vax_num, vax_last, pref, sex, age_group)%>%
  mutate(
    antiN_Mar2024 = case_when(
      pos_feb2024 == 0 & infection_history == 0 ~ 0, 
      T ~ 1
      ), 
    antiS_Mar2024 = case_when(
      pos_feb2024 == 0 & infection_history == 0 & vax_num ==0 ~ 0, 
      T ~ 1
      ), 
    antiN_Jan2024 = case_when(
      infection_history == 0 ~ 0, 
      T ~ 1
    ),
    antiS_Jan2024 = case_when(
      infection_history == 0 & vax_num == 0 ~ 0,
      infection_history == 0 & vax_num >=1 & vax_last >="2024/1/15" ~ 0,
      T ~ 1
    )
    ) -> df_seropos


####### 
weight_adj_draws <- subset_draws(sampling_result_2, variable = "alpha")
weight_adj_draws

dim(weight_adj_draws)[1] -> niter

df_analysis%>%
  filter(sex == 0 & age_group==2)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_male_20
df_analysis%>%
  filter(sex == 0 & age_group==3)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_male_30
df_analysis%>%
  filter(sex == 0 & age_group==4)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_male_40
df_analysis%>%
  filter(sex == 0 & age_group==5)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_male_50
df_analysis%>%
  filter(sex == 0 & age_group==6)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_male_60
df_analysis%>%
  filter(sex == 0 & age_group==7)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_male_70


male_antiS_Mar2024 = matrix(NA, nrow = 6, ncol = niter*4)
male_antiN_Mar2024 = matrix(NA, nrow = 6, ncol = niter*4)

male_antiS_Jan2024 = matrix(NA, nrow = 6, ncol = niter*4)
male_antiN_Jan2024 = matrix(NA, nrow = 6, ncol = niter*4)

for(a in 1:niter){
  for(b in 1:4){
    # a=1
    # b=1
    weight_adj_draws[a, b, c(idx_male_20)] %>%as.vector() -> weight_20
    # Mar2024
    weight_20%*%as.vector(df_seropos$antiS_Mar2024[c(idx_male_20)]) / sum(weight_20) -> male_antiS_Mar2024[1, 4*(a-1)+b]
    weight_20%*%as.vector(df_seropos$antiN_Mar2024[c(idx_male_20)]) / sum(weight_20) -> male_antiN_Mar2024[1, 4*(a-1)+b]
    # Jan2024
    weight_20%*%as.vector(df_seropos$antiS_Jan2024[c(idx_male_20)]) / sum(weight_20) -> male_antiS_Jan2024[1, 4*(a-1)+b]
    weight_20%*%as.vector(df_seropos$antiN_Jan2024[c(idx_male_20)]) / sum(weight_20) -> male_antiN_Jan2024[1, 4*(a-1)+b]
    
    weight_adj_draws[a, b, c(idx_male_30)] %>%as.vector() -> weight_30
    # Mar2024
    weight_30%*%as.vector(df_seropos$antiS_Mar2024[c(idx_male_30)]) / sum(weight_30) -> male_antiS_Mar2024[2, 4*(a-1)+b]
    weight_30%*%as.vector(df_seropos$antiN_Mar2024[c(idx_male_30)]) / sum(weight_30) -> male_antiN_Mar2024[2, 4*(a-1)+b]
    # Jan2024
    weight_30%*%as.vector(df_seropos$antiS_Jan2024[c(idx_male_30)]) / sum(weight_30) -> male_antiS_Jan2024[2, 4*(a-1)+b]
    weight_30%*%as.vector(df_seropos$antiN_Jan2024[c(idx_male_30)]) / sum(weight_30) -> male_antiN_Jan2024[2, 4*(a-1)+b]
    
    weight_adj_draws[a, b, c(idx_male_40)] %>%as.vector() -> weight_40
    # Mar2024
    weight_40%*%as.vector(df_seropos$antiS_Mar2024[c(idx_male_40)]) / sum(weight_40) -> male_antiS_Mar2024[3, 4*(a-1)+b]
    weight_40%*%as.vector(df_seropos$antiN_Mar2024[c(idx_male_40)]) / sum(weight_40) -> male_antiN_Mar2024[3, 4*(a-1)+b]
    # Jan2024
    weight_40%*%as.vector(df_seropos$antiS_Jan2024[c(idx_male_40)]) / sum(weight_40) -> male_antiS_Jan2024[3, 4*(a-1)+b]
    weight_40%*%as.vector(df_seropos$antiN_Jan2024[c(idx_male_40)]) / sum(weight_40) -> male_antiN_Jan2024[3, 4*(a-1)+b]
    
    weight_adj_draws[a, b, c(idx_male_50)] %>%as.vector() -> weight_50
    # Mar2024
    weight_50%*%as.vector(df_seropos$antiS_Mar2024[c(idx_male_50)]) / sum(weight_50) -> male_antiS_Mar2024[4, 4*(a-1)+b]
    weight_50%*%as.vector(df_seropos$antiN_Mar2024[c(idx_male_50)]) / sum(weight_50) -> male_antiN_Mar2024[4, 4*(a-1)+b]
    # Jan2024
    weight_50%*%as.vector(df_seropos$antiS_Jan2024[c(idx_male_50)]) / sum(weight_50) -> male_antiS_Jan2024[4, 4*(a-1)+b]
    weight_50%*%as.vector(df_seropos$antiN_Jan2024[c(idx_male_50)]) / sum(weight_50) -> male_antiN_Jan2024[4, 4*(a-1)+b]
    

    weight_adj_draws[a, b, c(idx_male_60)] %>%as.vector() -> weight_60
    # Mar2024
    weight_60%*%as.vector(df_seropos$antiS_Mar2024[c(idx_male_60)]) / sum(weight_60) -> male_antiS_Mar2024[5, 4*(a-1)+b]
    weight_60%*%as.vector(df_seropos$antiN_Mar2024[c(idx_male_60)]) / sum(weight_60) -> male_antiN_Mar2024[5, 4*(a-1)+b]
    # Jan2024
    weight_60%*%as.vector(df_seropos$antiS_Jan2024[c(idx_male_60)]) / sum(weight_60) -> male_antiS_Jan2024[5, 4*(a-1)+b]
    weight_60%*%as.vector(df_seropos$antiN_Jan2024[c(idx_male_60)]) / sum(weight_60) -> male_antiN_Jan2024[5, 4*(a-1)+b]
    
    weight_adj_draws[a, b, c(idx_male_70)] %>%as.vector() -> weight_70
    # Mar2024
    weight_70%*%as.vector(df_seropos$antiS_Mar2024[c(idx_male_70)]) / sum(weight_70) -> male_antiS_Mar2024[6, 4*(a-1)+b]
    weight_70%*%as.vector(df_seropos$antiN_Mar2024[c(idx_male_70)]) / sum(weight_70) -> male_antiN_Mar2024[6, 4*(a-1)+b]
    # Jan2024
    weight_70%*%as.vector(df_seropos$antiS_Jan2024[c(idx_male_70)]) / sum(weight_70) -> male_antiS_Jan2024[6, 4*(a-1)+b]
    weight_70%*%as.vector(df_seropos$antiN_Jan2024[c(idx_male_70)]) / sum(weight_70) -> male_antiN_Jan2024[6, 4*(a-1)+b]
    
  }
}  

male_antiS_Mar2024_summ = matrix(NA, nrow = 6, ncol = 4)%>%data.frame()
male_antiN_Mar2024_summ = matrix(NA, nrow = 6, ncol = 4)%>%data.frame()
male_antiS_Jan2024_summ = matrix(NA, nrow = 6, ncol = 4)%>%data.frame()
male_antiN_Jan2024_summ = matrix(NA, nrow = 6, ncol = 4)%>%data.frame()

male_antiS_Mar2024_summ[, 2:4] <- male_antiS_Mar2024%>%apply(MARGIN = 1, FUN = quantile, probs = c(.5, .025, .975))%>%t()
male_antiN_Mar2024_summ[, 2:4] <- male_antiN_Mar2024%>%apply(MARGIN = 1, FUN = quantile, probs = c(.5, .025, .975))%>%t()
male_antiS_Jan2024_summ[, 2:4] <- male_antiS_Jan2024%>%apply(MARGIN = 1, FUN = quantile, probs = c(.5, .025, .975))%>%t()
male_antiN_Jan2024_summ[, 2:4] <- male_antiN_Jan2024%>%apply(MARGIN = 1, FUN = quantile, probs = c(.5, .025, .975))%>%t()


male_antiN_Mar2024_summ[,1] <- c(paste0("ag" , 2:7*10))
male_antiS_Mar2024_summ[,1] <- c(paste0("ag" , 2:7*10))
male_antiN_Jan2024_summ[,1] <- c(paste0("ag" , 2:7*10))
male_antiS_Jan2024_summ[,1] <- c(paste0("ag" , 2:7*10))

colnames(male_antiN_Mar2024_summ) <- colnames(male_antiS_Mar2024_summ) <- c("age_group", "median", "lwr_95CrI", "upr_95CrI")
colnames(male_antiN_Jan2024_summ) <- colnames(male_antiS_Jan2024_summ) <- c("age_group", "median", "lwr_95CrI", "upr_95CrI")

male_antiN_Mar2024_summ
male_antiS_Mar2024_summ

male_antiN_Jan2024_summ
male_antiS_Jan2024_summ


###########################################################
df_analysis%>%
  filter(sex == 1 & age_group==2)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_female_20
df_analysis%>%
  filter(sex == 1 & age_group==3)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_female_30
df_analysis%>%
  filter(sex == 1 & age_group==4)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_female_40
df_analysis%>%
  filter(sex == 1 & age_group==5)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_female_50
df_analysis%>%
  filter(sex == 1 & age_group==6)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_female_60
df_analysis%>%
  filter(sex == 1 & age_group==7)%>%
  select(ID_new) %>%unlist()%>%as.vector() -> idx_female_70



######################

female_antiS_Mar2024 = matrix(NA, nrow = 6, ncol = niter*4)
female_antiN_Mar2024 = matrix(NA, nrow = 6, ncol = niter*4)

female_antiS_Jan2024 = matrix(NA, nrow = 6, ncol = niter*4)
female_antiN_Jan2024 = matrix(NA, nrow = 6, ncol = niter*4)

for(a in 1:niter){
  for(b in 1:4){
    # a=1
    # b=1
    weight_adj_draws[a, b, c(idx_female_20)] %>%as.vector() -> weight_20
    # Mar2024
    weight_20%*%as.vector(df_seropos$antiS_Mar2024[c(idx_female_20)]) / sum(weight_20) -> female_antiS_Mar2024[1, 4*(a-1)+b]
    weight_20%*%as.vector(df_seropos$antiN_Mar2024[c(idx_female_20)]) / sum(weight_20) -> female_antiN_Mar2024[1, 4*(a-1)+b]
    # Jan2024
    weight_20%*%as.vector(df_seropos$antiS_Jan2024[c(idx_female_20)]) / sum(weight_20) -> female_antiS_Jan2024[1, 4*(a-1)+b]
    weight_20%*%as.vector(df_seropos$antiN_Jan2024[c(idx_female_20)]) / sum(weight_20) -> female_antiN_Jan2024[1, 4*(a-1)+b]
    
    weight_adj_draws[a, b, c(idx_female_30)] %>%as.vector() -> weight_30
    # Mar2024
    weight_30%*%as.vector(df_seropos$antiS_Mar2024[c(idx_female_30)]) / sum(weight_30) -> female_antiS_Mar2024[2, 4*(a-1)+b]
    weight_30%*%as.vector(df_seropos$antiN_Mar2024[c(idx_female_30)]) / sum(weight_30) -> female_antiN_Mar2024[2, 4*(a-1)+b]
    # Jan2024
    weight_30%*%as.vector(df_seropos$antiS_Jan2024[c(idx_female_30)]) / sum(weight_30) -> female_antiS_Jan2024[2, 4*(a-1)+b]
    weight_30%*%as.vector(df_seropos$antiN_Jan2024[c(idx_female_30)]) / sum(weight_30) -> female_antiN_Jan2024[2, 4*(a-1)+b]
    
    weight_adj_draws[a, b, c(idx_female_40)] %>%as.vector() -> weight_40
    # Mar2024
    weight_40%*%as.vector(df_seropos$antiS_Mar2024[c(idx_female_40)]) / sum(weight_40) -> female_antiS_Mar2024[3, 4*(a-1)+b]
    weight_40%*%as.vector(df_seropos$antiN_Mar2024[c(idx_female_40)]) / sum(weight_40) -> female_antiN_Mar2024[3, 4*(a-1)+b]
    # Jan2024
    weight_40%*%as.vector(df_seropos$antiS_Jan2024[c(idx_female_40)]) / sum(weight_40) -> female_antiS_Jan2024[3, 4*(a-1)+b]
    weight_40%*%as.vector(df_seropos$antiN_Jan2024[c(idx_female_40)]) / sum(weight_40) -> female_antiN_Jan2024[3, 4*(a-1)+b]
    
    weight_adj_draws[a, b, c(idx_female_50)] %>%as.vector() -> weight_50
    # Mar2024
    weight_50%*%as.vector(df_seropos$antiS_Mar2024[c(idx_female_50)]) / sum(weight_50) -> female_antiS_Mar2024[4, 4*(a-1)+b]
    weight_50%*%as.vector(df_seropos$antiN_Mar2024[c(idx_female_50)]) / sum(weight_50) -> female_antiN_Mar2024[4, 4*(a-1)+b]
    # Jan2024
    weight_50%*%as.vector(df_seropos$antiS_Jan2024[c(idx_female_50)]) / sum(weight_50) -> female_antiS_Jan2024[4, 4*(a-1)+b]
    weight_50%*%as.vector(df_seropos$antiN_Jan2024[c(idx_female_50)]) / sum(weight_50) -> female_antiN_Jan2024[4, 4*(a-1)+b]
    
    
    weight_adj_draws[a, b, c(idx_female_60)] %>%as.vector() -> weight_60
    # Mar2024
    weight_60%*%as.vector(df_seropos$antiS_Mar2024[c(idx_female_60)]) / sum(weight_60) -> female_antiS_Mar2024[5, 4*(a-1)+b]
    weight_60%*%as.vector(df_seropos$antiN_Mar2024[c(idx_female_60)]) / sum(weight_60) -> female_antiN_Mar2024[5, 4*(a-1)+b]
    # Jan2024
    weight_60%*%as.vector(df_seropos$antiS_Jan2024[c(idx_female_60)]) / sum(weight_60) -> female_antiS_Jan2024[5, 4*(a-1)+b]
    weight_60%*%as.vector(df_seropos$antiN_Jan2024[c(idx_female_60)]) / sum(weight_60) -> female_antiN_Jan2024[5, 4*(a-1)+b]
    
    weight_adj_draws[a, b, c(idx_female_70)] %>%as.vector() -> weight_70
    # Mar2024
    weight_70%*%as.vector(df_seropos$antiS_Mar2024[c(idx_female_70)]) / sum(weight_70) -> female_antiS_Mar2024[6, 4*(a-1)+b]
    weight_70%*%as.vector(df_seropos$antiN_Mar2024[c(idx_female_70)]) / sum(weight_70) -> female_antiN_Mar2024[6, 4*(a-1)+b]
    # Jan2024
    weight_70%*%as.vector(df_seropos$antiS_Jan2024[c(idx_female_70)]) / sum(weight_70) -> female_antiS_Jan2024[6, 4*(a-1)+b]
    weight_70%*%as.vector(df_seropos$antiN_Jan2024[c(idx_female_70)]) / sum(weight_70) -> female_antiN_Jan2024[6, 4*(a-1)+b]
    
  }
}  

female_antiS_Mar2024_summ = matrix(NA, nrow = 6, ncol = 4)%>%data.frame()
female_antiN_Mar2024_summ = matrix(NA, nrow = 6, ncol = 4)%>%data.frame()
female_antiS_Jan2024_summ = matrix(NA, nrow = 6, ncol = 4)%>%data.frame()
female_antiN_Jan2024_summ = matrix(NA, nrow = 6, ncol = 4)%>%data.frame()

female_antiS_Mar2024_summ[, 2:4] <- female_antiS_Mar2024%>%apply(MARGIN = 1, FUN = quantile, probs = c(.5, .025, .975))%>%t()
female_antiN_Mar2024_summ[, 2:4] <- female_antiN_Mar2024%>%apply(MARGIN = 1, FUN = quantile, probs = c(.5, .025, .975))%>%t()
female_antiS_Jan2024_summ[, 2:4] <- female_antiS_Jan2024%>%apply(MARGIN = 1, FUN = quantile, probs = c(.5, .025, .975))%>%t()
female_antiN_Jan2024_summ[, 2:4] <- female_antiN_Jan2024%>%apply(MARGIN = 1, FUN = quantile, probs = c(.5, .025, .975))%>%t()


female_antiN_Mar2024_summ[,1] <- c(paste0("ag" , 2:7*10))
female_antiS_Mar2024_summ[,1] <- c(paste0("ag" , 2:7*10))
female_antiN_Jan2024_summ[,1] <- c(paste0("ag" , 2:7*10))
female_antiS_Jan2024_summ[,1] <- c(paste0("ag" , 2:7*10))

colnames(female_antiN_Mar2024_summ) <- colnames(female_antiS_Mar2024_summ) <- c("age_group", "median", "lwr_95CrI", "upr_95CrI")
colnames(female_antiN_Jan2024_summ) <- colnames(female_antiS_Jan2024_summ) <- c("age_group", "median", "lwr_95CrI", "upr_95CrI")

female_antiN_Mar2024_summ
female_antiS_Mar2024_summ

female_antiN_Jan2024_summ
female_antiS_Jan2024_summ


male_antiN_Mar2024_summ%<>%
  mutate(sex = "Male", 
         month = as.Date("2024/3/15"), 
         type = "anti-N")

male_antiN_Mar2024_summ%<>%
  mutate(obs_resid_serum_med = c(0.713, .752, .685, .586, .481, .382), 
         obs_resid_serum_lwr = c(0.641, .690, .621, .521, .411, .319), 
         obs_resid_serum_upr = c(0.777, .808, .745, .649, .551, .448)
         )%>%
  mutate(obs_donor_med = c(0.729, .711, .669, .598, .517, NA), 
         obs_donor_lwr = c(0.704, .687, .648, .575, .494, NA), 
         obs_donor_upr = c(0.752, .733, .690, .620, .540, NA)
  )

male_antiS_Mar2024_summ%<>%
  mutate(sex = "Male", 
         month = as.Date("2024/3/15"), 
         type = "anti-S")
male_antiS_Mar2024_summ%<>%
  mutate(obs_resid_serum_med = c(0.972, .968, .970, .971, .971, .996), 
         obs_resid_serum_lwr = c(0.937, .936, .939, .942, .938, .976), 
         obs_resid_serum_upr = c(0.991, .987, .988, .988, .989, 1.00)
  )%>%
  mutate(obs_donor_med = NA, 
         obs_donor_lwr = NA, 
         obs_donor_upr = NA
  )

  
male_antiN_Jan2024_summ%<>%
  mutate(sex = "Male", 
         month = as.Date("2024/1/15"), 
         type = "anti-N")
male_antiN_Jan2024_summ%<>%
  mutate(obs_resid_serum_med = c(0.718, .680, .652, .543, .435, .255), 
         obs_resid_serum_lwr = c(0.642, .611, .587, .478, .367, .198), 
         obs_resid_serum_upr = c(0.785, .743, .714, .606, .506, .318)
  )%>%
  mutate(obs_donor_med = c(0.674, .656, .603, .529, .459, NA), 
         obs_donor_lwr = c(0.648, .631, .581, .506, .436, NA), 
         obs_donor_upr = c(0.699, .680, .625, .552, .482, NA)
  )

male_antiS_Jan2024_summ%<>%
  mutate(sex = "Male", 
         month = as.Date("2024/1/15"), 
         type = "anti-S")
male_antiS_Jan2024_summ%<>%
  mutate(obs_resid_serum_med = c(0.957, .976, .974, .935, .967, .944), 
         obs_resid_serum_lwr = c(0.914, .944, .944, .897, .932, .905), 
         obs_resid_serum_upr = c(0.983, .992, .990, .963, .986, .971)
  )%>%
  mutate(obs_donor_med = NA, 
         obs_donor_lwr = NA, 
         obs_donor_upr = NA
  )

#####
female_antiN_Mar2024_summ%<>%
  mutate(sex = "Female", 
         month = as.Date("2024/3/15"), 
         type = "anti-N")
female_antiN_Mar2024_summ%<>%
  mutate(obs_resid_serum_med = c(0.742, .721, .734, .628, .559, .280), 
         obs_resid_serum_lwr = c(0.671, .655, .674, .566, .486, .228), 
         obs_resid_serum_upr = c(0.804, .781, .789, .688, .630, .336)
  )%>%
  mutate(obs_donor_med = c(0.723, .737, .667, .588, .515, NA), 
         obs_donor_lwr = c(0.697, .713, .646, .566, .492, NA), 
         obs_donor_upr = c(0.748, .759, .688, .610, .537, NA)
  )

female_antiS_Mar2024_summ%<>%
  mutate(sex = "Female", 
         month = as.Date("2024/3/15"), 
         type = "anti-S")

female_antiS_Mar2024_summ%<>%
  mutate(obs_resid_serum_med = c(0.994, .995, .983, .976, .959, .975), 
         obs_resid_serum_lwr = c(0.969, .974, .958, .949, .921, .949), 
         obs_resid_serum_upr = c(1.00, 1.00, .995, .991, .982, .990)
  ) %>%
  mutate(obs_donor_med = NA, 
         obs_donor_lwr = NA, 
         obs_donor_upr = NA
  )

female_antiN_Jan2024_summ%<>%
  mutate(sex = "Female", 
         month = as.Date("2024/1/15"), 
         type = "anti-N")
female_antiN_Jan2024_summ%<>%
  mutate(obs_resid_serum_med = c(0.707, .671, .589, .564, .472, .346), 
         obs_resid_serum_lwr = c(0.627, .603, .525, .501, .404, .289), 
         obs_resid_serum_upr = c(0.767, .735, .651, .626, .541, .406)
  )%>%
  mutate(obs_donor_med = c(0.684, .643, .617, .517, .508, NA), 
         obs_donor_lwr = c(0.657, .618, .595, .494, .486, NA), 
         obs_donor_upr = c(0.709, .667, .639, .539, .531, NA)
  )


female_antiS_Jan2024_summ%<>%
  mutate(sex = "Female", 
         month = as.Date("2024/1/15"), 
         type = "anti-S")
female_antiS_Jan2024_summ%<>%
  mutate(obs_resid_serum_med = c(0.983, .976, .980, .969, .977, .992), 
         obs_resid_serum_lwr = c(0.951, .945, .954, .940, .946, .973), 
         obs_resid_serum_upr = c(0.996, .992, .993, .986, .992, .999)
  )%>%
  mutate(obs_donor_med = NA, 
         obs_donor_lwr = NA, 
         obs_donor_upr = NA
  )


###
### PLOTS male
###
rbind(male_antiN_Jan2024_summ, male_antiN_Mar2024_summ, male_antiS_Jan2024_summ, male_antiS_Mar2024_summ) -> seropos_male_all

palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

plotN_male_March<- ggplot(seropos_male_all %>% filter(type == "anti-N", month == "2024-03-15"), 
       aes(x = age_group)) +
  geom_point(aes(y = median, color = "Present Study"), 
             size = 1, position = position_nudge(x = -0.2)) +
  geom_errorbar(aes(ymin = lwr_95CrI, ymax = upr_95CrI, color = "Present Study"), 
                width = 0.1, position = position_nudge(x = -0.2)) +
  
  geom_point(aes(y = obs_resid_serum_med, color = "Residual Serum Analysis"), 
             shape = 17, size = 1, position = position_nudge(x = 0)) +
  geom_errorbar(aes(ymin = obs_resid_serum_lwr, ymax = obs_resid_serum_upr, color = "Residual Serum Analysis"), 
                width = 0.1, position = position_nudge(x = 0)) +
  
  geom_point(aes(y = obs_donor_med, color = "Donated Blood Analysis"), 
             shape = 15, size = 1, position = position_nudge(x = 0.2)) +
  geom_errorbar(aes(ymin = obs_donor_lwr, ymax = obs_donor_upr, color = "Donated Blood Analysis"), 
                width = 0.1, position = position_nudge(x = 0.2)) +
  
  labs(
    title = "Past Infection: March 2024",
    x = "Age Group",
    y = "Value (with 95% CI)",
    color = "Source"
  ) +
  scale_color_manual(values = c("Present Study" = palette[2], 
                                "Residual Serum Analysis" = palette[5], 
                                "Donated Blood Analysis" = palette[4])) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.8, 0.85),
    legend.background = element_rect(fill = "white", color = "black")  
  )

plotN_male_Jan<- ggplot(seropos_male_all %>% filter(type == "anti-N", month == "2024-01-15"), 
                     aes(x = age_group)) +
  geom_point(aes(y = median, color = "Present Study"), 
             size = 1, position = position_nudge(x = -0.2)) +
  geom_errorbar(aes(ymin = lwr_95CrI, ymax = upr_95CrI, color = "Present Study"), 
                width = 0.1, position = position_nudge(x = -0.2)) +
  
  geom_point(aes(y = obs_resid_serum_med, color = "Residual Serum Analysis"), 
             shape = 17, size = 1, position = position_nudge(x = 0)) +
  geom_errorbar(aes(ymin = obs_resid_serum_lwr, ymax = obs_resid_serum_upr, color = "Residual Serum Analysis"), 
                width = 0.1, position = position_nudge(x = 0)) +
  
  geom_point(aes(y = obs_donor_med, color = "Donated Blood Analysis"), 
             shape = 15, size = 1, position = position_nudge(x = 0.2)) +
  geom_errorbar(aes(ymin = obs_donor_lwr, ymax = obs_donor_upr, color = "Donated Blood Analysis"), 
                width = 0.1, position = position_nudge(x = 0.2)) +
  
  labs(
    title = "Past Infection: Jan 2024",
    x = "Age Group",
    y = "Value (with 95% CI)",
    color = "Source"
  ) +
  scale_color_manual(values = c("Present Study" = palette[2], 
                                "Residual Serum Analysis" = palette[5], 
                                "Donated Blood Analysis" = palette[4])) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.8, 0.85), 
    legend.background = element_rect(fill = "white", color = "black") 
  )

###################################

plotS_male_March<- ggplot(seropos_male_all %>% filter(type == "anti-S", month == "2024-03-15"), 
                     aes(x = age_group)) +
  geom_point(aes(y = median, color = "Present Study"), 
             size = 1, position = position_nudge(x = -0.2)) +
  geom_errorbar(aes(ymin = lwr_95CrI, ymax = upr_95CrI, color = "Present Study"), 
                width = 0.1, position = position_nudge(x = -0.2)) +
  
  geom_point(aes(y = obs_resid_serum_med, color = "Residual Serum Analysis"), 
             shape = 17, size = 1, position = position_nudge(x = 0)) +
  geom_errorbar(aes(ymin = obs_resid_serum_lwr, ymax = obs_resid_serum_upr, color = "Residual Serum Analysis"), 
                width = 0.1, position = position_nudge(x = 0)) +
  labs(
    title = "Any Exposure: March 2024",
    x = "Age Group",
    y = "Value (with 95% CI)",
    color = "Source"
  ) +
  scale_color_manual(values = c("Present Study" = palette[2], 
                                "Residual Serum Analysis" = palette[5])
                     ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.8, 0.25),   
    legend.background = element_rect(fill = "white", color = "black") 
  )

plotS_male_Jan<- ggplot(seropos_male_all %>% filter(type == "anti-S", month == "2024-01-15"), 
                     aes(x = age_group)) +
  geom_point(aes(y = median, color = "Present Study"), 
             size = 1, position = position_nudge(x = -0.2)) +
  geom_errorbar(aes(ymin = lwr_95CrI, ymax = upr_95CrI, color = "Present Study"), 
                width = 0.1, position = position_nudge(x = -0.2)) +
  
  geom_point(aes(y = obs_resid_serum_med, color = "Residual Serum Analysis"), 
             shape = 17, size = 1, position = position_nudge(x = 0)) +
  geom_errorbar(aes(ymin = obs_resid_serum_lwr, ymax = obs_resid_serum_upr, color = "Residual Serum Analysis"), 
                width = 0.1, position = position_nudge(x = 0)) +
  labs(
    title = "Any Exposure: Jan 2024",
    x = "Age Group",
    y = "Value (with 95% CI)",
    color = "Source"
  ) +
  scale_color_manual(values = c("Present Study" = palette[2], 
                                "Residual Serum Analysis" = palette[5])
) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.8, 0.25), 
    legend.background = element_rect(fill = "white", color = "black") 
  )


gridExtra::grid.arrange(plotN_male_Jan, plotN_male_March, 
                        plotS_male_Jan, plotS_male_March, 
                        ncol=2)

tiff("./results_feb_2024/Fig_Sxxx_seropos_male.tiff", width = 18, height = 12, units = "in", res = 1000, compression = "lzw")
gridExtra::grid.arrange(plotN_male_Jan, plotN_male_March, 
                        plotS_male_Jan, plotS_male_March, 
                        ncol=2)
dev.off()




###
### PLOTS female
###
rbind(female_antiN_Jan2024_summ, female_antiN_Mar2024_summ, female_antiS_Jan2024_summ, female_antiS_Mar2024_summ)-> seropos_female_all

plotN_female_March<- ggplot(seropos_female_all %>% filter(type == "anti-N", month == "2024-03-15"), 
                          aes(x = age_group)) +
  geom_point(aes(y = median, color = "Present Study"), 
             size = 1, position = position_nudge(x = -0.2)) +
  geom_errorbar(aes(ymin = lwr_95CrI, ymax = upr_95CrI, color = "Present Study"), 
                width = 0.1, position = position_nudge(x = -0.2)) +
  
  geom_point(aes(y = obs_resid_serum_med, color = "Residual Serum Analysis"), 
             shape = 17, size = 1, position = position_nudge(x = 0)) +
  geom_errorbar(aes(ymin = obs_resid_serum_lwr, ymax = obs_resid_serum_upr, color = "Residual Serum Analysis"), 
                width = 0.1, position = position_nudge(x = 0)) +
  
  geom_point(aes(y = obs_donor_med, color = "Donated Blood Analysis"), 
             shape = 15, size = 1, position = position_nudge(x = 0.2)) +
  geom_errorbar(aes(ymin = obs_donor_lwr, ymax = obs_donor_upr, color = "Donated Blood Analysis"), 
                width = 0.1, position = position_nudge(x = 0.2)) +
  
  labs(
    title = "Past Infection: March 2024",
    x = "Age Group",
    y = "Value (with 95% CI)",
    color = "Source"
  ) +
  scale_color_manual(values = c("Present Study" = palette[2], 
                                "Residual Serum Analysis" = palette[5], 
                                "Donated Blood Analysis" = palette[4])) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.8, 0.85), 
    legend.background = element_rect(fill = "white", color = "black")  
  )

plotN_female_Jan<- ggplot(seropos_female_all %>% filter(type == "anti-N", month == "2024-01-15"), 
                        aes(x = age_group)) +
  geom_point(aes(y = median, color = "Present Study"), 
             size = 1, position = position_nudge(x = -0.2)) +
  geom_errorbar(aes(ymin = lwr_95CrI, ymax = upr_95CrI, color = "Present Study"), 
                width = 0.1, position = position_nudge(x = -0.2)) +
  
  geom_point(aes(y = obs_resid_serum_med, color = "Residual Serum Analysis"), 
             shape = 17, size = 1, position = position_nudge(x = 0)) +
  geom_errorbar(aes(ymin = obs_resid_serum_lwr, ymax = obs_resid_serum_upr, color = "Residual Serum Analysis"), 
                width = 0.1, position = position_nudge(x = 0)) +
  
  geom_point(aes(y = obs_donor_med, color = "Donated Blood Analysis"), 
             shape = 15, size = 1, position = position_nudge(x = 0.2)) +
  geom_errorbar(aes(ymin = obs_donor_lwr, ymax = obs_donor_upr, color = "Donated Blood Analysis"), 
                width = 0.1, position = position_nudge(x = 0.2)) +
  
  labs(
    title = "Past Infection: Jan 2024",
    x = "Age Group",
    y = "Value (with 95% CI)",
    color = "Source"
  ) +
  scale_color_manual(values = c("Present Study" = palette[2], 
                                "Residual Serum Analysis" = palette[5], 
                                "Donated Blood Analysis" = palette[4])) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.8, 0.85),  
    legend.background = element_rect(fill = "white", color = "black")  
  )

###################################

plotS_female_March<- ggplot(seropos_female_all %>% filter(type == "anti-S", month == "2024-03-15"), 
                          aes(x = age_group)) +
  geom_point(aes(y = median, color = "Present Study"), 
             size = 1, position = position_nudge(x = -0.2)) +
  geom_errorbar(aes(ymin = lwr_95CrI, ymax = upr_95CrI, color = "Present Study"), 
                width = 0.1, position = position_nudge(x = -0.2)) +
  
  geom_point(aes(y = obs_resid_serum_med, color = "Residual Serum Analysis"), 
             shape = 17, size = 1, position = position_nudge(x = 0)) +
  geom_errorbar(aes(ymin = obs_resid_serum_lwr, ymax = obs_resid_serum_upr, color = "Residual Serum Analysis"), 
                width = 0.1, position = position_nudge(x = 0)) +
  labs(
    title = "Any Exposure: March 2024",
    x = "Age Group",
    y = "Value (with 95% CI)",
    color = "Source"
  ) +
  scale_color_manual(values = c("Present Study" = palette[2], 
                                "Residual Serum Analysis" = palette[5])
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.8, 0.25),  
    legend.background = element_rect(fill = "white", color = "black")  
  )

plotS_female_Jan<- ggplot(seropos_female_all %>% filter(type == "anti-S", month == "2024-01-15"), 
                        aes(x = age_group)) +
  geom_point(aes(y = median, color = "Present Study"), 
             size = 1, position = position_nudge(x = -0.2)) +
  geom_errorbar(aes(ymin = lwr_95CrI, ymax = upr_95CrI, color = "Present Study"), 
                width = 0.1, position = position_nudge(x = -0.2)) +
  
  geom_point(aes(y = obs_resid_serum_med, color = "Residual Serum Analysis"), 
             shape = 17, size = 1, position = position_nudge(x = 0)) +
  geom_errorbar(aes(ymin = obs_resid_serum_lwr, ymax = obs_resid_serum_upr, color = "Residual Serum Analysis"), 
                width = 0.1, position = position_nudge(x = 0)) +
  labs(
    title = "Any Exposure: Jan 2024",
    x = "Age Group",
    y = "Value (with 95% CI)",
    color = "Source"
  ) +
  scale_color_manual(values = c("Present Study" = palette[2], 
                                "Residual Serum Analysis" = palette[5])
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.8, 0.25),
    legend.background = element_rect(fill = "white", color = "black")  
  )


gridExtra::grid.arrange(plotN_female_Jan, plotN_female_March, 
                        plotS_female_Jan, plotS_female_March, 
                        ncol=2)

tiff("./results_feb_2024/Fig_Sxxx_seropos_female.tiff", width = 18, height = 12, units = "in", res = 1000, compression = "lzw")
gridExtra::grid.arrange(plotN_female_Jan, plotN_female_March, 
                        plotS_female_Jan, plotS_female_March, 
                        ncol=2)
dev.off()

rbind(seropos_male_all, seropos_female_all)%>%readr::write_csv(file = "./results_feb_2024/seropos_data.csv")

