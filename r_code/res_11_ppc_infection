rm(list=ls())
library(posterior);library(tidyr);library(ggplot2)
df_analysis <- readr::read_csv(file = "./CODES/data_questionnaire.csv")%>%
  filter(infect_period2!=1)
df_analysis%<>%
  mutate(ID = 1:nrow(.))

load(file = "path_to_RData") # load cmdstan fit object

prob_infect_draws <- subset_draws(sampling_result_2, variable = "prob_infection")
prob_infect_mat <- as_draws_matrix(prob_infect_draws)

runif_mat <- matrix(runif(length(prob_infect_mat)), 
                    nrow = nrow(prob_infect_mat), 
                    ncol = ncol(prob_infect_mat))

set.seed(202402)
bin_outcome_mat <- (runif_mat < prob_infect_mat)*1L

#### by age and sex
ppc_res <- matrix(nrow = 12, ncol = 6) %>%
  data.frame()

colnames(ppc_res) <- c("age_group", "sex", "med", "lwr", "upr", "obs")


ppc_res_hist <- matrix(nrow = 12, ncol = 2 + 1 + 2 + 4000) %>%
  data.frame()

colnames(ppc_res_hist) <- c("age_group", "sex", "obs", "lwr", "upr", paste0("draw", 1:4000))



df_analysis$ID2 <- 1:nrow(df_analysis)

par(mfrow=c(3,4))
for(ag in 1:6){
  for(sex_idx in 0:1){
    row_idx = (ag-1)*2 + sex_idx + 1
    ppc_res[row_idx, 1] <- ag+1
    ppc_res[row_idx, 2] <- sex_idx
    
    ppc_res_hist[row_idx, 1] <- ag+1
    ppc_res_hist[row_idx, 2] <- sex_idx
    
    
    df_analysis%>%
      filter(age_group == ag+1 & sex == sex_idx) -> df_filtered
    
    bin_outcome_mat[, c(df_filtered$ID2)] -> ppc_mat_tmp
    
    ppc_mat_tmp %>% 
      apply(sum, MARGIN = 1) %>% 
      quantile(probs = c(.5, .025, .975)) -> ppc_res[row_idx, 3:5]
    
    ppc_mat_tmp %>% 
      apply(sum, MARGIN = 1) -> ppc_res_hist[row_idx, 5+1:4000]
    
    #### observed
    df_filtered$pos_feb2024 %>% sum() -> ppc_res[row_idx, 6]
    df_filtered$pos_feb2024 %>% sum() -> ppc_res_hist[row_idx, 3]
    
    ppc_res_hist[row_idx, 4:5]   <- ppc_mat_tmp %>% apply(1, sum) %>% quantile(c(.025, .975)) 
    
  }
}


ppc_res_hist_long <- ppc_res_hist %>%
  pivot_longer(
    cols = starts_with("draw"), 
    names_to = "draw_id", 
    values_to = "sim_count"
  )%>%
  mutate(
    age_label = case_when(
      age_group == 2 ~ "20-29", 
      age_group == 3 ~ "30-39", 
      age_group == 4 ~ "40-49", 
      age_group == 5 ~ "50-59", 
      age_group == 6 ~ "60-69", 
      age_group == 7 ~ "70+"
    ), 
    sex_label = case_when(
      sex ==0 ~ "Male", 
      sex ==1 ~ "Female"
    )
  )

ppc_res_hist_long%>%
  ggplot(aes(x = sim_count))+
  geom_histogram(binwidth = 1, fill = "skyblue", alpha = 0.8)+
  geom_vline(data = ppc_res_hist_long%>%
               distinct(age_label, sex_label, obs), 
             aes(xintercept = obs, color = "Observed"), 
             size = 1)+
  geom_vline(aes(xintercept = lwr, color = "95% PI"), linetype = "dashed") +
  geom_vline(aes(xintercept = upr, color = "95% PI"), linetype = "dashed") +
  geom_hline(aes(yintercept = 0, color = ""), 
             size = 0.25)+
  facet_wrap(~age_label + sex_label, scales = "free_y")+
  labs(
    title = "Posterior Predictive Check, COVID-19 infection in Feb 2024", 
    x = "Number of Infections", 
    y = "Frequency", 
    color = ""
  )+
  scale_color_manual(values = c("Observed" = "red"))+
  theme_bw() -> ppc_hist

tiff("path_tiff", width = 12, height = 6, units = "in", res = 1000, compression = "lzw")
ppc_hist
dev.off()


#### by prefecture ####
group_df_pref <- df_analysis %>%
  group_by(pref) %>%  
  summarise(
    id_list = list(ID2),
    n_g     = n(),
    obs_cnt = sum(pos_feb2024),
    .groups = "drop"
  )

summ_pref <- lapply(seq_len(nrow(group_df_pref)), function(i){
  ids   <- group_df_pref$id_list[[i]]
  sims  <- rowSums(bin_outcome_mat[, ids, drop = FALSE])   
  qs    <- quantile(sims, c(.025, .5, .975))
  sd_s  <- sd(sims)
  data.frame(
    prefecture = group_df_pref$pref[i],
    obs_cnt    = group_df_pref$obs_cnt[i],
    med_cnt    = qs[2],
    lwr_cnt    = qs[1],
    upr_cnt    = qs[3],
    sd_cnt     = sd_s
  )
})

ppc_pref_counts <- bind_rows(summ_pref) %>% 
  dplyr::left_join(dplyr::select(group_df_pref, pref, n_g),
                   by = c("prefecture" = "pref")) %>%
  dplyr::transmute(
    Prefecture = prefecture,
    N          = n_g,           
    Observed   = obs_cnt,
    Pred_med   = round(med_cnt),
    `95%PrI_lower`  = round(lwr_cnt),
    `95%PrI_lower`  = round(upr_cnt),
    `Covered`  = ifelse(Observed >= `95%PrI_lower` & Observed <= `95%PrI_lower`, "Yes", "No") 
    )

readr::write_csv(ppc_pref_counts, "path_csv")

