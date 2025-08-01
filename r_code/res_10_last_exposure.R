rm(list=ls())
library(posterior);library(tidyr);library(ggplot2)
df_analysis <- readr::read_csv(file = "./CODES/data_questionnaire.csv")%>%
  filter(infect_period2!=1)
df_analysis%<>%
  mutate(ID = 1:nrow(.))

load(file = "path_to_RData") # load cmdstan fit object

df_plot <- df_analysis %>%
  mutate(imm_type = case_when(
    imm_last == 1 & infect_period2 == 2 ~ "Infection (XBB)",
    imm_last == 1 & infect_period2 == 3 ~ "Infection (Pre-XBB Omicron)",
    imm_last == 1 & infect_period2 == 4 ~ "Infection (Pre-Omicron)",
    
    imm_last == 2 & vax_period2 == 1 ~ "Vaccine (XBB)",
    imm_last == 2 & vax_period2 == 2 ~ "Vaccine (Bivalent)",
    imm_last == 2 & vax_period2 == 3 ~ "Vaccine (Wuhan)",
    
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(imm_type))

plot_order <- c(
  "Infection (XBB)", 
  "Infection (Pre-XBB Omicron)",
  "Infection (Pre-Omicron)",
  "Vaccine (XBB)",
  "Vaccine (Bivalent)",
  "Vaccine (Wuhan)"
)

df_plot_ordered <- df_plot %>%
  mutate(imm_type = factor(imm_type, levels = plot_order))

df_plot_final <- df_plot_ordered %>%
  mutate(
    pos_feb2024_status = factor(
      pos_feb2024,
      levels = c(1, 0),
      labels = c("Positive (Feb 2024)", "Negative (Feb 2024)")
    )
  )

df_plot_final%>%
  ggplot(aes(x = time_since_imm)) +
  geom_histogram(aes(fill = imm_type), bins = 30, show.legend = FALSE) +
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  facet_grid(pos_feb2024_status ~ imm_type, , scales = "free_y") +
  labs(
    title = "Time since last immunization by exposure type and outcome",
    x = "Time since last immunization (days)",
    y = "Freq"
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 8)) -> plot_time_since_imm



tiff("path_tiff", width = 12, height = 6, units = "in", res = 1000, compression = "lzw")
plot_time_since_imm
dev.off()
