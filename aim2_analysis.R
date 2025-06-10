# Aim 2 analysis 
## Susan Hoffman 

#### Set up ####

library(ggplot2) 
library(tidyverse)
library(readr)

setwd("/Users/ssb214/Library/CloudStorage/OneDrive-EmoryUniversity/MyFiles/Dissertation/Aim 2/Output")
# setwd("C:/Users/SSBUCKE/OneDrive - Emory University/MyFiles/Dissertation/Aim 2/Output")

#### Loop plots #### 

#### p = 200 ####

##### Loading in simulation results #####

truth_200 <- read_csv("truth_200.csv", col_types = cols(...1 = col_skip()))

# list files in the output folder 
files_200 <- list.files(pattern = "^(HDMA200|HIMA200|MITM200).*\\.csv$")

# loading in the files and adding in the names 
data_files <- lapply(files_200, read.csv)

names(data_files) <- files_200

# combining the results from the different runs 
HDMA200 <- do.call(rbind, data_files[1:5]) 
HDMA200_sesp <- do.call(rbind, data_files[6:10]) 

HIMA200 <- do.call(rbind, data_files[11:15])
HIMA200_sesp <- do.call(rbind, data_files[16:20]) 

MITM200 <- do.call(rbind, data_files[21:25]) 
MITM200_sesp <- do.call(rbind, data_files[26:30]) 

# Visualization 

##### TIE #####

## Setting up the TIE results 

### HDMA

HDMA200_vis <- HDMA200 %>% 
  select(independent, med_beta, truth, n, indirect, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(TIE_HDMA = mean(indirect, na.rm = T),
         TIE_HDMA_min = min(indirect, na.rm = T),
         TIE_HDMA_max = max(indirect, na.rm = T)) %>% 
  distinct(group_id, TIE_HDMA, .keep_all = T)

HDMA200_vis_ind <- HDMA200_vis %>% 
  filter(independent == "yes") 

HDMA200_vis_cor <- HDMA200_vis %>% 
  filter(independent == "no")

### HIMA

HIMA200_vis <- HIMA200 %>% 
  select(independent, med_beta, truth, n, indirect, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(TIE_HIMA = mean(indirect, na.rm = T),
         TIE_HIMA_min = min(indirect, na.rm = T),
         TIE_HIMA_max = max(indirect, na.rm = T)) %>% 
  distinct(group_id, TIE_HIMA, .keep_all = T)

HIMA200_vis_ind <- HIMA200_vis %>% 
  filter(independent == "yes")

HIMA200_vis_cor <- HIMA200_vis %>% 
  filter(independent == "no")

### MITM

MITM200_vis <- MITM200 %>% 
  select(independent, med_beta, truth, n, tie, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(TIE_MITM = mean(tie, na.rm = T),
         TIE_MITM_min = min(tie, na.rm = T),
         TIE_MITM_max = max(tie, na.rm = T)) %>% 
  distinct(group_id, TIE_MITM, .keep_all = T)

MITM200_vis_ind <- MITM200_vis %>% 
  filter(independent == "yes")

MITM200_vis_cor <- MITM200_vis %>% 
  filter(independent == "no")

### Combined data 

vis_ind <- HDMA200_vis_ind %>% 
  select(independent, med_beta, truth, n, TIE_HDMA) %>% 
  left_join(HIMA200_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  left_join(MITM200_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  select(med_beta, truth, n, TIE_HDMA, TIE_HIMA, TIE_MITM, independent.x)

vis_cor <- HDMA200_vis_cor %>% 
  select(independent, med_beta, truth, n, TIE_HDMA) %>% 
  left_join(HIMA200_vis_cor, by = c("med_beta", "truth", "n")) %>% 
  left_join(MITM200_vis_cor, by = c("med_beta", "truth", "n")) %>% 
  select(med_beta, truth, n, TIE_HDMA, TIE_HIMA, TIE_MITM, independent.x)

vis_data <- vis_cor %>% 
  rbind(vis_ind) %>% 
  select(med_beta, truth, n, independent.x, TIE_HDMA, TIE_HIMA, TIE_MITM) %>% 
  rename("independent" = independent.x) %>% 
  left_join(truth_200, by = c("med_beta", "truth", "independent")) %>% 
  rename(
    "Mediator beta" = med_beta,
    "Number of mediators" = truth,
    "Independent" = independent,
    "HDMA" = TIE_HDMA,
    "HIMA" = TIE_HIMA,
    "MITM" = TIE_MITM,
    "Truth" = TIE 
  ) %>% 
  select(-c(p, cov_n, CIE))

# remotes::install_github("matherealize/looplot")
library(looplot)

nested_loop_plot(resdf = vis_data, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "TIE",
                 na_rm = F,
                 ylim = c(0, 0.3),
                 line_linetypes = c(1, 1, 1, 5),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 #colors = c("#1B9E77FF", "#D95F02FF", "#7570B3FF", "#E7298AFF" ),
                 colors = c("#1B9E77FF", "grey", "grey", "#E7298AFF" ),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.45, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.position = "top"
                   )
                 )
             )


##### CIE #####

### HDMA

HDMA200_vis <- HDMA200 %>% 
  select(independent, med_beta, truth, n, alpha_beta, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(CIE_HDMA = mean(alpha_beta, na.rm = T)) %>% 
  distinct(group_id, CIE_HDMA, .keep_all = T)

HDMA200_vis_ind <- HDMA200_vis %>% 
  filter(independent == "yes") 

HDMA200_vis_cor <- HDMA200_vis %>% 
  filter(independent == "no")

### HIMA

HIMA200_vis <- HIMA200 %>% 
  select(independent, med_beta, truth, n, alpha.beta, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(CIE_HIMA = mean(alpha.beta, na.rm = T)) %>% 
  distinct(group_id, CIE_HIMA, .keep_all = T)

HIMA200_vis_ind <- HIMA200_vis %>% 
  filter(independent == "yes")

HIMA200_vis_cor <- HIMA200_vis %>% 
  filter(independent == "no")

### MITM

MITM200_vis <- MITM200 %>% 
  select(independent, med_beta, truth, n, alphabeta, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(CIE_MITM = mean(alphabeta, na.rm = T)) %>% 
  distinct(group_id, CIE_MITM, .keep_all = T)

MITM200_vis_ind <- MITM200_vis %>% 
  filter(independent == "yes")

MITM200_vis_cor <- MITM200_vis %>% 
  filter(independent == "no")

### Combined data 

vis_ind <- HDMA200_vis_ind %>% 
  select(independent, med_beta, truth, n, CIE_HDMA) %>% 
  left_join(HIMA200_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  left_join(MITM200_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  select(med_beta, truth, n, CIE_HDMA, CIE_HIMA, CIE_MITM, independent.x)

## Correlated not meaningful 
# vis_cor <- HDMA200_vis_cor %>% 
#   select(independent, med_beta, truth, n, TIE_HDMA) %>% 
#   left_join(HIMA200_vis_cor, by = c("med_beta", "truth", "n")) %>% 
#   left_join(MITM200_vis_cor, by = c("med_beta", "truth", "n")) %>% 
#   select(med_beta, truth, n, TIE_HDMA, TIE_HIMA, TIE_MITM, independent.x)

vis_data <- vis_ind %>% 
  select(med_beta, truth, n, independent.x, CIE_HDMA, CIE_HIMA, CIE_MITM) %>% 
  rename("independent" = independent.x) %>% 
  left_join(truth_200, by = c("med_beta", "truth", "independent")) %>% 
  rename(
    "Mediator beta" = med_beta,
    "Number of mediators" = truth,
    "Independent" = independent,
    "HDMA" = CIE_HDMA,
    "HIMA" = CIE_HIMA,
    "MITM" = CIE_MITM,
    "Truth" = CIE 
  ) %>% 
  select(-c(p, cov_n, TIE))

# remotes::install_github("matherealize/looplot")
library(looplot)

nested_loop_plot(resdf = vis_data, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "CIE",
                 na_rm = F,
                 ylim = c(0, 0.15),
                 line_linetypes = c(1, 1, 1, 5),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 colors = scales::brewer_pal(palette = "Dark2"),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.35, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.position = "top"
                   )
                 )
             )


##### SE/SP #####

se_sp <- function(data){
  se <- numeric(nrow(data))
  sp <- numeric(nrow(data))
  
  for(i in 1:nrow(data)) {
    true_values <- paste0("M.", 1:data[i, "truth"])
    false_values <- paste0("M.", (data[i, "truth"] + 1):ncol(data)) 
    false_values <- false_values[false_values %in% names(data)]
    
    true_positives <- sum(data[i, true_values] == 1)
    false_negatives <- length(true_values) - true_positives
    
    true_negatives <- sum(data[i, false_values] == 0)
    false_positives <- sum(data[i, false_values] == 1)
    
    # Calculate se and sp for each row 
    se[i] <- true_positives/(true_positives+false_negatives)
    sp[i] <- true_negatives/(true_negatives+false_positives)
  }
  
  return(data.frame(
    Sensitivity = se, 
    Specificity = sp,
    sim_id = data$sim_id, 
    n = data$n, 
    p = data$p, 
    truth = data$truth, 
    independent = data$independent,
    cov_n = data$cov_n, 
    med_beta = data$med_beta
  ))
  
}

# Pulling Se and Sp 

HDMA_sesp <- se_sp(HDMA200_sesp)

HIMA_sesp <- se_sp(HIMA200_sesp)

MITM_sesp <- se_sp(MITM200_sesp)

# Removing those from NA runs 

rows_to_remove <- which(HDMA_sesp["Sensitivity"] == 0.00 &
                                HDMA_sesp["Specificity"] == 1.0000000) # 38
HDMA_sesp <- HDMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(HIMA_sesp["Sensitivity"] == 0.00 &
                          HIMA_sesp["Specificity"] == 1.0000000) # 546
HIMA_sesp <- HIMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(MITM_sesp["Sensitivity"] == 0.00 &
                          MITM_sesp["Specificity"] == 1.0000000) # 972
MITM_sesp <- MITM_sesp[-rows_to_remove, ]

HDMA_se <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Sensitivity = mean(Sensitivity, na.rm = T)) %>% 
  distinct(group_id, Sensitivity, .keep_all = T) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HDMA SE" = Sensitivity)

HDMA_sp <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Specificity = mean(Specificity, na.rm = T)) %>% 
  distinct(group_id, Specificity, .keep_all = T) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HDMA SP" = Specificity)

HIMA_se <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Sensitivity = mean(Sensitivity, na.rm = T)) %>% 
  distinct(group_id, Sensitivity, .keep_all = T) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HIMA SE" = Sensitivity)

HIMA_sp <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Specificity = mean(Specificity, na.rm = T)) %>% 
  distinct(group_id, Specificity, .keep_all = T) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HIMA SP" = Specificity)

MITM_se <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Sensitivity = mean(Sensitivity, na.rm = T)) %>% 
  distinct(group_id, Sensitivity, .keep_all = T) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("MITM SE" = Sensitivity)

MITM_sp <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Specificity = mean(Specificity, na.rm = T)) %>% 
  distinct(group_id, Specificity, .keep_all = T) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("MITM SP" = Specificity)

sesp_vis <- HDMA_se %>% 
  left_join(HDMA_sp, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(HIMA_se, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(HIMA_sp, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(MITM_se, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(MITM_sp, by = c("med_beta", "truth", "n", "independent")) %>% 
  select(med_beta, truth, n, independent, 
         `HDMA SE`, `HDMA SP`, 
         `HIMA SE`, `HIMA SP`,
         `MITM SE`, `MITM SP`) %>% 
  rename(
    "Mediator beta" = med_beta,
    "Number of mediators" = truth,
    "Independent" = independent
  )

se_vis <- sesp_vis %>% 
  select(`Mediator beta`, `Number of mediators`, n, Independent, 
         `HDMA SE`, 
         `HIMA SE`,
         `MITM SE`) %>% 
  rename("HIMA" = `HIMA SE`,
         "HDMA" = `HDMA SE`,
         "MITM" = `MITM SE`)

sp_vis <- sesp_vis %>% 
  select(`Mediator beta`, `Number of mediators`, n, Independent, 
         `HDMA SP`, 
         `HIMA SP`,
         `MITM SP`) %>% 
  rename("HIMA" = `HIMA SP`,
         "HDMA" = `HDMA SP`,
         "MITM" = `MITM SP`)

# Visualize 

library(looplot)

# Sensitivity 

nested_loop_plot(resdf = se_vis, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "Sensitivity",
                 na_rm = F,
                 ylim = c(0, 1.2),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 colors = scales::brewer_pal(palette = "Dark2"),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.45, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.position = "top"
                   )
                 )
)

# Specificity 

nested_loop_plot(resdf = sp_vis, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "Specificity",
                 na_rm = F,
                 ylim = c(0, 1.2),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 colors = scales::brewer_pal(palette = "Dark2"),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.45, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.position = "top"
                   )
                 )
)



#### p = 400 ####

##### Loading in simulation results #####

truth_400 <- read_csv("truth_400.csv", col_types = cols(...1 = col_skip()))

# list files in the output folder 
files_400 <- list.files(pattern = "^(HDMA400|HIMA400|MITM400).*\\.csv$")

# loading in the files and adding in the names 
data_files <- lapply(files_400, read.csv)

names(data_files) <- files_400

# combining the results from the different runs 
HDMA400 <- do.call(rbind, data_files[1:5]) 
HDMA400_sesp <- do.call(rbind, data_files[6:10]) 

HIMA400 <- do.call(rbind, data_files[11:15])
HIMA400_sesp <- do.call(rbind, data_files[16:20]) 

MITM400 <- do.call(rbind, data_files[21:25]) 
MITM400_sesp <- do.call(rbind, data_files[26:30]) 

# Visualization 

##### TIE #####

## Setting up the TIE results 

### HDMA

HDMA400_vis <- HDMA400 %>% 
  select(independent, med_beta, truth, n, indirect, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(TIE_HDMA = mean(indirect, na.rm = T)) %>% 
  distinct(group_id, TIE_HDMA, .keep_all = T)

HDMA400_vis_ind <- HDMA400_vis %>% 
  filter(independent == "yes") 

HDMA400_vis_cor <- HDMA400_vis %>% 
  filter(independent == "no")

### HIMA

HIMA400_vis <- HIMA400 %>% 
  select(independent, med_beta, truth, n, indirect, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(TIE_HIMA = mean(indirect, na.rm = T)) %>% 
  distinct(group_id, TIE_HIMA, .keep_all = T)

HIMA400_vis_ind <- HIMA400_vis %>% 
  filter(independent == "yes")

HIMA400_vis_cor <- HIMA400_vis %>% 
  filter(independent == "no")

### MITM

MITM400_vis <- MITM400 %>% 
  select(independent, med_beta, truth, n, tie, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(TIE_MITM = mean(tie, na.rm = T)) %>% 
  distinct(group_id, TIE_MITM, .keep_all = T)

MITM400_vis_ind <- MITM400_vis %>% 
  filter(independent == "yes")

MITM400_vis_cor <- MITM400_vis %>% 
  filter(independent == "no")

### Combined data 

vis_ind <- HDMA400_vis_ind %>% 
  select(independent, med_beta, truth, n, TIE_HDMA) %>% 
  left_join(HIMA400_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  left_join(MITM400_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  select(med_beta, truth, n, TIE_HDMA, TIE_HIMA, TIE_MITM, independent.x)

vis_cor <- HDMA400_vis_cor %>% 
  select(independent, med_beta, truth, n, TIE_HDMA) %>% 
  left_join(HIMA400_vis_cor, by = c("med_beta", "truth", "n")) %>% 
  left_join(MITM400_vis_cor, by = c("med_beta", "truth", "n")) %>% 
  select(med_beta, truth, n, TIE_HDMA, TIE_HIMA, TIE_MITM, independent.x)

vis_data <- vis_cor %>% 
  rbind(vis_ind) %>% 
  select(med_beta, truth, n, independent.x, TIE_HDMA, TIE_HIMA, TIE_MITM) %>% 
  rename("independent" = independent.x) %>% 
  left_join(truth_400, by = c("med_beta", "truth", "independent")) %>% 
  rename(
    "Mediator beta" = med_beta,
    "Number of mediators" = truth,
    "Independent" = independent,
    "HDMA" = TIE_HDMA,
    "HIMA" = TIE_HIMA,
    "MITM" = TIE_MITM,
    "Truth" = TIE 
  ) %>% 
  select(-c(p, cov_n, CIE))

# remotes::install_github("matherealize/looplot")
library(looplot)

nested_loop_plot(resdf = vis_data, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "TIE",
                 na_rm = F,
                 ylim = c(0, 0.5),
                 line_linetypes = c(1, 1, 1, 5),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 #colors = c("#1B9E77FF", "#D95F02FF", "#7570B3FF", "#E7298AFF" ),
                 colors = c("#1B9E77FF", "grey", "grey", "#E7298AFF" ),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.5, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.position = "top"
                   )
                 )
)

##### CIE #####

### HDMA

HDMA400_vis <- HDMA400 %>% 
  select(independent, med_beta, truth, n, alpha_beta, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(CIE_HDMA = mean(alpha_beta, na.rm = T)) %>% 
  distinct(group_id, CIE_HDMA, .keep_all = T)

HDMA400_vis_ind <- HDMA400_vis %>% 
  filter(independent == "yes") 

HDMA400_vis_cor <- HDMA400_vis %>% 
  filter(independent == "no")

### HIMA

HIMA400_vis <- HIMA400 %>% 
  select(independent, med_beta, truth, n, alpha.beta, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(CIE_HIMA = mean(alpha.beta, na.rm = T)) %>% 
  distinct(group_id, CIE_HIMA, .keep_all = T)

HIMA400_vis_ind <- HIMA400_vis %>% 
  filter(independent == "yes")

HIMA400_vis_cor <- HIMA400_vis %>% 
  filter(independent == "no")

### MITM

MITM400_vis <- MITM400 %>% 
  select(independent, med_beta, truth, n, alphabeta, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(CIE_MITM = mean(alphabeta, na.rm = T)) %>% 
  distinct(group_id, CIE_MITM, .keep_all = T)

MITM400_vis_ind <- MITM400_vis %>% 
  filter(independent == "yes")

MITM400_vis_cor <- MITM400_vis %>% 
  filter(independent == "no")

### Combined data 

vis_ind <- HDMA400_vis_ind %>% 
  select(independent, med_beta, truth, n, CIE_HDMA) %>% 
  left_join(HIMA400_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  left_join(MITM400_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  select(med_beta, truth, n, CIE_HDMA, CIE_HIMA, CIE_MITM, independent.x)

## Correlated not meaningful 
# vis_cor <- HDMA200_vis_cor %>% 
#   select(independent, med_beta, truth, n, TIE_HDMA) %>% 
#   left_join(HIMA200_vis_cor, by = c("med_beta", "truth", "n")) %>% 
#   left_join(MITM200_vis_cor, by = c("med_beta", "truth", "n")) %>% 
#   select(med_beta, truth, n, TIE_HDMA, TIE_HIMA, TIE_MITM, independent.x)

vis_data <- vis_ind %>% 
  select(med_beta, truth, n, independent.x, CIE_HDMA, CIE_HIMA, CIE_MITM) %>% 
  rename("independent" = independent.x) %>% 
  left_join(truth_400, by = c("med_beta", "truth", "independent")) %>% 
  rename(
    "Mediator beta" = med_beta,
    "Number of mediators" = truth,
    "Independent" = independent,
    "HDMA" = CIE_HDMA,
    "HIMA" = CIE_HIMA,
    "MITM" = CIE_MITM,
    "Truth" = CIE 
  ) %>% 
  select(-c(p, cov_n, TIE))

# remotes::install_github("matherealize/looplot")
library(looplot)

nested_loop_plot(resdf = vis_data, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "CIE",
                 na_rm = F,
                 ylim = c(0, 0.15),
                 line_linetypes = c(1, 1, 1, 5),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 colors = scales::brewer_pal(palette = "Dark2"),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.35, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.position = "top"
                   )
                 )
)


##### SE/SP #####

se_sp <- function(data){
  se <- numeric(nrow(data))
  sp <- numeric(nrow(data))
  
  for(i in 1:nrow(data)) {
    true_values <- paste0("M.", 1:data[i, "truth"])
    false_values <- paste0("M.", (data[i, "truth"] + 1):ncol(data)) 
    false_values <- false_values[false_values %in% names(data)]
    
    true_positives <- sum(data[i, true_values] == 1)
    false_negatives <- length(true_values) - true_positives
    
    true_negatives <- sum(data[i, false_values] == 0)
    false_positives <- sum(data[i, false_values] == 1)
    
    # Calculate se and sp for each row 
    se[i] <- true_positives/(true_positives+false_negatives)
    sp[i] <- true_negatives/(true_negatives+false_positives)
  }
  
  return(data.frame(
    Sensitivity = se, 
    Specificity = sp,
    sim_id = data$sim_id, 
    n = data$n, 
    p = data$p, 
    truth = data$truth, 
    independent = data$independent,
    cov_n = data$cov_n, 
    med_beta = data$med_beta
  ))
  
}

# Pulling Se and Sp 

HDMA_sesp <- se_sp(HDMA400_sesp)

HIMA_sesp <- se_sp(HIMA400_sesp)

MITM_sesp <- se_sp(MITM400_sesp)

# Removing those from NA runs 

rows_to_remove <- which(HDMA_sesp["Sensitivity"] == 0.00 &
                          HDMA_sesp["Specificity"] == 1.0000000) # 38
HDMA_sesp <- HDMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(HIMA_sesp["Sensitivity"] == 0.00 &
                          HIMA_sesp["Specificity"] == 1.0000000) # 546
HIMA_sesp <- HIMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(MITM_sesp["Sensitivity"] == 0.00 &
                          MITM_sesp["Specificity"] == 1.0000000) # 972
MITM_sesp <- MITM_sesp[-rows_to_remove, ]

HDMA_se <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Sensitivity = mean(Sensitivity, na.rm = T)) %>% 
  distinct(group_id, Sensitivity, .keep_all = T) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HDMA SE" = Sensitivity)

HDMA_sp <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Specificity = mean(Specificity, na.rm = T)) %>% 
  distinct(group_id, Specificity, .keep_all = T) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HDMA SP" = Specificity)

HIMA_se <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Sensitivity = mean(Sensitivity, na.rm = T)) %>% 
  distinct(group_id, Sensitivity, .keep_all = T) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HIMA SE" = Sensitivity)

HIMA_sp <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Specificity = mean(Specificity, na.rm = T)) %>% 
  distinct(group_id, Specificity, .keep_all = T) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HIMA SP" = Specificity)

MITM_se <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Sensitivity = mean(Sensitivity, na.rm = T)) %>% 
  distinct(group_id, Sensitivity, .keep_all = T) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("MITM SE" = Sensitivity)

MITM_sp <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Specificity = mean(Specificity, na.rm = T)) %>% 
  distinct(group_id, Specificity, .keep_all = T) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("MITM SP" = Specificity)

sesp_vis <- HDMA_se %>% 
  left_join(HDMA_sp, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(HIMA_se, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(HIMA_sp, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(MITM_se, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(MITM_sp, by = c("med_beta", "truth", "n", "independent")) %>% 
  select(med_beta, truth, n, independent, 
         `HDMA SE`, `HDMA SP`, 
         `HIMA SE`, `HIMA SP`,
         `MITM SE`, `MITM SP`) %>% 
  rename(
    "Mediator beta" = med_beta,
    "Number of mediators" = truth,
    "Independent" = independent
  )

se_vis <- sesp_vis %>% 
  select(`Mediator beta`, `Number of mediators`, n, Independent, 
         `HDMA SE`, 
         `HIMA SE`,
         `MITM SE`) %>% 
  rename("HIMA" = `HIMA SE`,
         "HDMA" = `HDMA SE`,
         "MITM" = `MITM SE`)

sp_vis <- sesp_vis %>% 
  select(`Mediator beta`, `Number of mediators`, n, Independent, 
         `HDMA SP`, 
         `HIMA SP`,
         `MITM SP`) %>% 
  rename("HIMA" = `HIMA SP`,
         "HDMA" = `HDMA SP`,
         "MITM" = `MITM SP`)

# Visualize 

library(looplot)

# Sensitivity 

nested_loop_plot(resdf = se_vis, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "Sensitivity",
                 na_rm = F,
                 ylim = c(0, 1.2),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 colors = scales::brewer_pal(palette = "Dark2"),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.45, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.position = "top"
                   )
                 )
)

# Specificity 

nested_loop_plot(resdf = sp_vis, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "Specificity",
                 na_rm = F,
                 ylim = c(0, 1.2),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 colors = scales::brewer_pal(palette = "Dark2"),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.45, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.position = "top"
                   )
                 )
)


#### p = 600 ####

##### Loading in simulation results #####

truth_600 <- read_csv("truth_600.csv", col_types = cols(...1 = col_skip()))

# list files in the output folder 
files_600 <- list.files(pattern = "^(HDMA600|HIMA600|MITM600).*\\.csv$")

# loading in the files and adding in the names 
data_files <- lapply(files_600, read.csv)

names(data_files) <- files_600

# combining the results from the different runs 
HDMA600 <- do.call(rbind, data_files[1:5]) 
HDMA600_sesp <- do.call(rbind, data_files[6:10]) 

HIMA600 <- do.call(rbind, data_files[11:15])
HIMA600_sesp <- do.call(rbind, data_files[16:20]) 

MITM600 <- do.call(rbind, data_files[21:25]) 
MITM600_sesp <- do.call(rbind, data_files[26:30]) 

# Visualization 

##### TIE #####

## Setting up the TIE results 

### HDMA

HDMA600_vis <- HDMA600 %>% 
  select(independent, med_beta, truth, n, indirect, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(TIE_HDMA = mean(indirect, na.rm = T)) %>% 
  distinct(group_id, TIE_HDMA, .keep_all = T)

HDMA600_vis_ind <- HDMA600_vis %>% 
  filter(independent == "yes") 

HDMA600_vis_cor <- HDMA600_vis %>% 
  filter(independent == "no")

### HIMA

HIMA600_vis <- HIMA600 %>% 
  select(independent, med_beta, truth, n, indirect, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(TIE_HIMA = mean(indirect, na.rm = T)) %>% 
  distinct(group_id, TIE_HIMA, .keep_all = T)

HIMA600_vis_ind <- HIMA600_vis %>% 
  filter(independent == "yes")

HIMA600_vis_cor <- HIMA600_vis %>% 
  filter(independent == "no")

### MITM

MITM600_vis <- MITM600 %>% 
  select(independent, med_beta, truth, n, tie, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(TIE_MITM = mean(tie, na.rm = T)) %>% 
  distinct(group_id, TIE_MITM, .keep_all = T)

MITM600_vis_ind <- MITM600_vis %>% 
  filter(independent == "yes")

MITM600_vis_cor <- MITM600_vis %>% 
  filter(independent == "no")

### Combined data 

vis_ind <- HDMA600_vis_ind %>% 
  select(independent, med_beta, truth, n, TIE_HDMA) %>% 
  left_join(HIMA600_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  left_join(MITM600_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  select(med_beta, truth, n, TIE_HDMA, TIE_HIMA, TIE_MITM, independent.x)

vis_cor <- HDMA600_vis_cor %>% 
  select(independent, med_beta, truth, n, TIE_HDMA) %>% 
  left_join(HIMA600_vis_cor, by = c("med_beta", "truth", "n")) %>% 
  left_join(MITM600_vis_cor, by = c("med_beta", "truth", "n")) %>% 
  select(med_beta, truth, n, TIE_HDMA, TIE_HIMA, TIE_MITM, independent.x)

vis_data <- vis_cor %>% 
  rbind(vis_ind) %>% 
  select(med_beta, truth, n, independent.x, TIE_HDMA, TIE_HIMA, TIE_MITM) %>% 
  rename("independent" = independent.x) %>% 
  left_join(truth_600, by = c("med_beta", "truth", "independent")) %>% 
  rename(
    "Mediator beta" = med_beta,
    "Number of mediators" = truth,
    "Independent" = independent,
    "HDMA" = TIE_HDMA,
    "HIMA" = TIE_HIMA,
    "MITM" = TIE_MITM,
    "Truth" = TIE 
  ) %>% 
  select(-c(p, cov_n, CIE))

# remotes::install_github("matherealize/looplot")
library(looplot)

nested_loop_plot(resdf = vis_data, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "TIE",
                 na_rm = F,
                 ylim = c(0, 6),
                 line_linetypes = c(1, 1, 1, 5),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 #colors = c("#1B9E77FF", "#D95F02FF", "#7570B3FF", "#E7298AFF" ),
                 colors = c("#1B9E77FF", "grey", "grey", "#E7298AFF" ),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.45, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.position = "top"
                   )
                 )
)


##### CIE #####

### HDMA

HDMA600_vis <- HDMA600 %>% 
  select(independent, med_beta, truth, n, alpha_beta, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(CIE_HDMA = mean(alpha_beta, na.rm = T)) %>% 
  distinct(group_id, CIE_HDMA, .keep_all = T)

HDMA600_vis_ind <- HDMA600_vis %>% 
  filter(independent == "yes") 

HDMA600_vis_cor <- HDMA600_vis %>% 
  filter(independent == "no")

### HIMA

HIMA600_vis <- HIMA600 %>% 
  select(independent, med_beta, truth, n, alpha.beta, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(CIE_HIMA = mean(alpha.beta, na.rm = T)) %>% 
  distinct(group_id, CIE_HIMA, .keep_all = T)

HIMA600_vis_ind <- HIMA600_vis %>% 
  filter(independent == "yes")

HIMA600_vis_cor <- HIMA600_vis %>% 
  filter(independent == "no")

### MITM

MITM600_vis <- MITM600 %>% 
  select(independent, med_beta, truth, n, alphabeta, sim_id) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(CIE_MITM = mean(alphabeta, na.rm = T)) %>% 
  distinct(group_id, CIE_MITM, .keep_all = T)

MITM600_vis_ind <- MITM600_vis %>% 
  filter(independent == "yes")

MITM600_vis_cor <- MITM600_vis %>% 
  filter(independent == "no")

### Combined data 

vis_ind <- HDMA600_vis_ind %>% 
  select(independent, med_beta, truth, n, CIE_HDMA) %>% 
  left_join(HIMA600_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  left_join(MITM600_vis_ind, by = c("med_beta", "truth", "n")) %>% 
  select(med_beta, truth, n, CIE_HDMA, CIE_HIMA, CIE_MITM, independent.x)

## Correlated not meaningful 
# vis_cor <- HDMA200_vis_cor %>% 
#   select(independent, med_beta, truth, n, TIE_HDMA) %>% 
#   left_join(HIMA200_vis_cor, by = c("med_beta", "truth", "n")) %>% 
#   left_join(MITM200_vis_cor, by = c("med_beta", "truth", "n")) %>% 
#   select(med_beta, truth, n, TIE_HDMA, TIE_HIMA, TIE_MITM, independent.x)

vis_data <- vis_ind %>% 
  select(med_beta, truth, n, independent.x, CIE_HDMA, CIE_HIMA, CIE_MITM) %>% 
  rename("independent" = independent.x) %>% 
  left_join(truth_600, by = c("med_beta", "truth", "independent")) %>% 
  rename(
    "Mediator beta" = med_beta,
    "Number of mediators" = truth,
    "Independent" = independent,
    "HDMA" = CIE_HDMA,
    "HIMA" = CIE_HIMA,
    "MITM" = CIE_MITM,
    "Truth" = CIE 
  ) %>% 
  select(-c(p, cov_n, TIE))

# remotes::install_github("matherealize/looplot")
library(looplot)

nested_loop_plot(resdf = vis_data, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "CIE",
                 na_rm = F,
                 ylim = c(0, 0.15),
                 line_linetypes = c(1, 1, 1, 5),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 colors = scales::brewer_pal(palette = "Dark2"),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.35, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.position = "top"
                   )
                 )
)


##### SE/SP #####

se_sp <- function(data){
  se <- numeric(nrow(data))
  sp <- numeric(nrow(data))
  
  for(i in 1:nrow(data)) {
    true_values <- paste0("M.", 1:data[i, "truth"])
    false_values <- paste0("M.", (data[i, "truth"] + 1):ncol(data)) 
    false_values <- false_values[false_values %in% names(data)]
    
    true_positives <- sum(data[i, true_values] == 1)
    false_negatives <- length(true_values) - true_positives
    
    true_negatives <- sum(data[i, false_values] == 0)
    false_positives <- sum(data[i, false_values] == 1)
    
    # Calculate se and sp for each row 
    se[i] <- true_positives/(true_positives+false_negatives)
    sp[i] <- true_negatives/(true_negatives+false_positives)
  }
  
  return(data.frame(
    Sensitivity = se, 
    Specificity = sp,
    sim_id = data$sim_id, 
    n = data$n, 
    p = data$p, 
    truth = data$truth, 
    independent = data$independent,
    cov_n = data$cov_n, 
    med_beta = data$med_beta
  ))
  
}

# Pulling Se and Sp 

HDMA_sesp <- se_sp(HDMA600_sesp)

HIMA_sesp <- se_sp(HIMA600_sesp)

MITM_sesp <- se_sp(MITM600_sesp)

# Removing those from NA runs 

rows_to_remove <- which(HDMA_sesp["Sensitivity"] == 0.00 &
                          HDMA_sesp["Specificity"] == 1.0000000) # 38
aHDMA_sesp <- HDMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(HIMA_sesp["Sensitivity"] == 0.00 &
                          HIMA_sesp["Specificity"] == 1.0000000) # 546
HIMA_sesp <- HIMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(MITM_sesp["Sensitivity"] == 0.00 &
                          MITM_sesp["Specificity"] == 1.0000000) # 972
MITM_sesp <- MITM_sesp[-rows_to_remove, ]

HDMA_se <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Sensitivity = mean(Sensitivity, na.rm = T)) %>% 
  distinct(group_id, Sensitivity, .keep_all = T) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HDMA SE" = Sensitivity)

HDMA_sp <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Specificity = mean(Specificity, na.rm = T)) %>% 
  distinct(group_id, Specificity, .keep_all = T) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HDMA SP" = Specificity)

HIMA_se <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Sensitivity = mean(Sensitivity, na.rm = T)) %>% 
  distinct(group_id, Sensitivity, .keep_all = T) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HIMA SE" = Sensitivity)

HIMA_sp <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Specificity = mean(Specificity, na.rm = T)) %>% 
  distinct(group_id, Specificity, .keep_all = T) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("HIMA SP" = Specificity)

MITM_se <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Sensitivity = mean(Sensitivity, na.rm = T)) %>% 
  distinct(group_id, Sensitivity, .keep_all = T) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("MITM SE" = Sensitivity)

MITM_sp <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  mutate(Specificity = mean(Specificity, na.rm = T)) %>% 
  distinct(group_id, Specificity, .keep_all = T) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) %>% 
  rename("MITM SP" = Specificity)

sesp_vis <- HDMA_se %>% 
  left_join(HDMA_sp, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(HIMA_se, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(HIMA_sp, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(MITM_se, by = c("med_beta", "truth", "n", "independent")) %>% 
  left_join(MITM_sp, by = c("med_beta", "truth", "n", "independent")) %>% 
  select(med_beta, truth, n, independent, 
         `HDMA SE`, `HDMA SP`, 
         `HIMA SE`, `HIMA SP`,
         `MITM SE`, `MITM SP`) %>% 
  rename(
    "Mediator beta" = med_beta,
    "Number of mediators" = truth,
    "Independent" = independent
  )

se_vis <- sesp_vis %>% 
  select(`Mediator beta`, `Number of mediators`, n, Independent, 
         `HDMA SE`, 
         `HIMA SE`,
         `MITM SE`) %>% 
  rename("HIMA" = `HIMA SE`,
         "HDMA" = `HDMA SE`,
         "MITM" = `MITM SE`)

sp_vis <- sesp_vis %>% 
  select(`Mediator beta`, `Number of mediators`, n, Independent, 
         `HDMA SP`, 
         `HIMA SP`,
         `MITM SP`) %>% 
  rename("HIMA" = `HIMA SP`,
         "HDMA" = `HDMA SP`,
         "MITM" = `MITM SP`)

# Visualize 

library(looplot)

# Sensitivity 

nested_loop_plot(resdf = se_vis, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "Sensitivity",
                 na_rm = F,
                 ylim = c(0, 1.2),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 colors = scales::brewer_pal(palette = "Dark2"),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.45, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.title = element_text(size = 14),
                     legend.position = "top"
                   )
                 )
)

# Specificity 

nested_loop_plot(resdf = sp_vis, 
                 x = "n", steps = "Number of mediators",
                 grid_rows = "Mediator beta", grid_cols = "Independent",
                 steps_y_base = -0.07, steps_y_height = 0.1,
                 x_name = "Sample Size", y_name = "Specificity",
                 na_rm = F,
                 ylim = c(0, 1.2),
                 spu_x_shift = 300,
                 hline_intercept = 0,
                 colors = scales::brewer_pal(palette = "Dark2"),
                 steps_values_annotate = TRUE, steps_annotation_size = 2,
                 y_expand_add = c(0.45, 0),
                 post_processing = list(
                   add_custom_theme = list(
                     strip.text.y = element_text(size = 14),
                     strip.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     axis.text.x = element_text(angle = -90),
                     legend.title = element_text(size = 14),
                     legend.position = "top"
                   )
                 )
)


#### TIE violin plots ####

group_names <- c(
  "1" = "n=200; B=0.1; Med=2%",
  "2" = "n=500; B=0.1; Med=2%",
  "3" = "n=1000; B=0.1; Med=2%",
  "4" = "n=200; B=0.1; Med=5%",
  "5" = "n=500; B=0.1; Med=5%",
  "6" = "n=1000; B=0.1; Med=5%",
  "7" = "n=200; B=0.1; Med=10%",
  "8" = "n=500; B=0.1; Med=10%",
  "9" = "n=1000; B=0.1; Med=10%",
  "10" = "n=200; B=0.3; Med=2%",
  "11" = "n=500; B=0.3; Med=2%",
  "12" = "n=1000; B=0.3; Med=2%",
  "13" = "n=200; B=0.3; Med=5%",
  "14" = "n=500; B=0.3; Med=5%",
  "15" = "n=1000; B=0.3; Med=5%",
  "16" = "n=200; B=0.3; Med=10%",
  "17" = "n=500; B=0.3; Med=10%",
  "18" = "n=1000; B=0.3; Med=10%",
  "19" = "n=200; B=0.1; Med=2%",
  "20" = "n=500; B=0.1; Med=2%",
  "21" = "n=1000; B=0.1; Med=2%",
  "22" = "n=200; B=0.1; Med=5%",
  "23" = "n=500; B=0.1; Med=5%",
  "24" = "n=1000; B=0.1; Med=5%",
  "25" = "n=200; B=0.1; Med=10%",
  "26" = "n=500; B=0.1; Med=10%",
  "27" = "n=1000; B=0.1; Med=10%",
  "28" = "n=200; B=0.3; Med=2%",
  "29" = "n=500; B=0.3; Med=2%",
  "30" = "n=1000; B=0.3; Med=2%",
  "31" = "n=200; B=0.3; Med=5%",
  "32" = "n=500; B=0.3; Med=5%",
  "33" = "n=1000; B=0.3; Med=5%",
  "34" = "n=200; B=0.3; Med=10%",
  "35" = "n=500; B=0.3; Med=10%",
  "36" = "n=1000; B=0.3; Med=10%"
)

#### p = 200 ####

truth_200 <- read_csv("truth_200.csv", col_types = cols(...1 = col_skip()))

# list files in the output folder 
files_200 <- list.files(pattern = "^(HDMA200|HIMA200|MITM200).*\\.csv$")

# loading in the files and adding in the names 
data_files <- lapply(files_200, read.csv)

names(data_files) <- files_200

# combining the results from the different runs 
HDMA200 <- do.call(rbind, data_files[1:5]) 
HDMA200_sesp <- do.call(rbind, data_files[6:10]) 

HIMA200 <- do.call(rbind, data_files[11:15])
HIMA200_sesp <- do.call(rbind, data_files[16:20]) 

MITM200 <- do.call(rbind, data_files[21:25]) 
MITM200_sesp <- do.call(rbind, data_files[26:30]) 

# Editing datasets for visualizations 

HDMA_density <- HDMA200 %>% 
  select(sim_id, indirect, n, p, independent, cov_n, truth, med_beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_200 %>% 
              select(p, truth, independent, cov_n, med_beta, TIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta"))

HIMA_density <- HIMA200 %>% 
  select(sim_id, indirect, n, p, independent, cov_n, truth, med_beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_200 %>% 
              select(p, truth, independent, cov_n, med_beta, TIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta"))

MITM_density <- MITM200 %>% 
  select(sim_id, tie, n, p, independent, cov_n, truth, med_beta) %>% 
  rename(indirect = tie) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_200 %>% 
              select(p, truth, independent, cov_n, med_beta, TIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta"))

## Combining the datasets 
HDMA_density <- HDMA_density %>% 
  mutate(method = "HDMA")
MITM_density <- MITM_density %>% 
  mutate(method = "MITM")
HIMA_density <- HIMA_density %>% 
  mutate(method = "HIMA")

combined_density <- bind_rows(HDMA_density, MITM_density, HIMA_density)

# Density plot 
# ggplot(combined_density, aes(x = indirect, fill = method, color = method)) +
#   geom_density(alpha = 0.4, position = "identity") +  # Use alpha for transparency
#   geom_vline(aes(xintercept = mean(TIE)), color = "red", linetype = "dashed") +
#   labs(title = "Density of Indirect Estimates by Method",
#        x = "Indirect Estimate",
#        y = "Density") +
#   facet_wrap(~group_id, scales = "free") +  # Facet by group_id
#   scale_fill_manual(values = c("HDMA" = "blue", "MITM" = "green", "HIMA" = "purple"), name = "Method") +
#   scale_color_manual(values = c("HDMA" = "blue", "MITM" = "green", "HIMA" = "purple"), name = "Method") +
#   theme_minimal()

# violin plot 

ggplot(combined_density, aes(x = factor(group_id), y = indirect, fill = method)) +
  geom_violin(alpha = 0.7) +
  geom_hline(aes(yintercept = TIE), color = "red", linetype = "dashed") +  # Use the truth (TIE) for each group_id
  labs(title = "TIE, p=200",
       x = "Simulation scenario",
       y = "TIE") +
  facet_wrap(~group_id, scales = "free") +  
  scale_x_discrete(labels = group_names) +  # Add custom labels
  scale_fill_manual(values = c("HDMA" = "blue", "MITM" = "green", "HIMA" = "purple"), name = "Method") +
  theme_minimal()
#### 1-18 are independent == no

# Table data

summary_table_200 <- combined_density %>%
  group_by(group_id, method) %>%
  summarise(
    mean_indirect = mean(indirect, na.rm = TRUE),
    sd_indirect = sd(indirect, na.rm = TRUE),
    tie_value = unique(TIE),  # Assuming `tie` is consistent within group_id
    bias = mean_indirect - tie_value
  ) %>%
  # Create the bias column with tie value in parentheses
  mutate(
    bias = paste0(round(bias, 3), " (", tie_value, ")"),
    mean_sd_indirect = paste0(round(mean_indirect, 3), " (", round(sd_indirect, 3), ")"),
    group_id = factor(group_id, levels = names(group_names), labels = group_names)  # Replace IDs with names
  ) %>%
  select(group_id, method, mean_sd_indirect, bias)


#### p = 400 ####

truth_400 <- read_csv("truth_400.csv", col_types = cols(...1 = col_skip()))

# list files in the output folder 
files_400 <- list.files(pattern = "^(HDMA400|HIMA400|MITM400).*\\.csv$")

# loading in the files and adding in the names 
data_files <- lapply(files_400, read.csv)

names(data_files) <- files_400

# combining the results from the different runs 
HDMA400 <- do.call(rbind, data_files[1:5]) 
HDMA400_sesp <- do.call(rbind, data_files[6:10]) 

HIMA400 <- do.call(rbind, data_files[11:15])
HIMA400_sesp <- do.call(rbind, data_files[16:20]) 

MITM400 <- do.call(rbind, data_files[21:25]) 
MITM400_sesp <- do.call(rbind, data_files[26:30]) 


# Editing datasets for visualizations 

HDMA_density <- HDMA400 %>% 
  select(sim_id, indirect, n, p, independent, cov_n, truth, med_beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_400 %>% 
              select(p, truth, independent, cov_n, med_beta, TIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta"))

HIMA_density <- HIMA400 %>% 
  select(sim_id, indirect, n, p, independent, cov_n, truth, med_beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_400 %>% 
              select(p, truth, independent, cov_n, med_beta, TIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta"))

MITM_density <- MITM400 %>% 
  select(sim_id, tie, n, p, independent, cov_n, truth, med_beta) %>% 
  rename(indirect = tie) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_400 %>% 
              select(p, truth, independent, cov_n, med_beta, TIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta"))

## Combining the datasets 
HDMA_density <- HDMA_density %>% 
  mutate(method = "HDMA")
MITM_density <- MITM_density %>% 
  mutate(method = "MITM")
HIMA_density <- HIMA_density %>% 
  mutate(method = "HIMA")

combined_density <- bind_rows(HDMA_density, MITM_density, HIMA_density)

# violin plot 

ggplot(combined_density, aes(x = factor(group_id), y = indirect, fill = method)) +
  geom_violin(alpha = 0.7) +
  geom_hline(aes(yintercept = TIE), color = "red", linetype = "dashed") +  # Use the truth (TIE) for each group_id
  labs(title = "TIE, p=400",
       x = "Simulation scenario",
       y = "TIE") +
  facet_wrap(~group_id, scales = "free") +  
  scale_x_discrete(labels = group_names) +  # Add custom labels
  scale_fill_manual(values = c("HDMA" = "blue", "MITM" = "green", "HIMA" = "purple"), name = "Method") +
  theme_minimal()

# Table data

summary_table_400 <- combined_density %>%
  group_by(group_id, method) %>%
  summarise(
    mean_indirect = mean(indirect, na.rm = TRUE),
    sd_indirect = sd(indirect, na.rm = TRUE),
    tie_value = unique(TIE),  # Assuming `tie` is consistent within group_id
    bias = mean_indirect - tie_value
  ) %>%
  # Create the bias column with tie value in parentheses
  mutate(
    bias = paste0(round(bias, 3), " (", round(tie_value, 3), ")"),
    mean_sd_indirect = paste0(round(mean_indirect, 3), " (", round(sd_indirect, 3), ")"),
    group_id = factor(group_id, levels = names(group_names), labels = group_names)  # Replace IDs with names
  ) %>%
  select(group_id, method, mean_sd_indirect, bias)

#### p = 600 ####

truth_600 <- read_csv("truth_600.csv", col_types = cols(...1 = col_skip()))

# list files in the output folder 
files_600 <- list.files(pattern = "^(HDMA600|HIMA600|MITM600).*\\.csv$")

# loading in the files and adding in the names 
data_files <- lapply(files_600, read.csv)

names(data_files) <- files_600

# combining the results from the different runs 
HDMA600 <- do.call(rbind, data_files[1:5]) 
HDMA600_sesp <- do.call(rbind, data_files[6:10]) 

HIMA600 <- do.call(rbind, data_files[11:15])
HIMA600_sesp <- do.call(rbind, data_files[16:20]) 

MITM600 <- do.call(rbind, data_files[21:25]) 
MITM600_sesp <- do.call(rbind, data_files[26:30]) 

# Editing datasets for visualizations 

HDMA_density <- HDMA600 %>% 
  select(sim_id, indirect, n, p, independent, cov_n, truth, med_beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_600 %>% 
              select(p, truth, independent, cov_n, med_beta, TIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta"))

HIMA_density <- HIMA600 %>% 
  select(sim_id, indirect, n, p, independent, cov_n, truth, med_beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_600 %>% 
              select(p, truth, independent, cov_n, med_beta, TIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta"))

MITM_density <- MITM600 %>% 
  select(sim_id, tie, n, p, independent, cov_n, truth, med_beta) %>% 
  rename(indirect = tie) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_600 %>% 
              select(p, truth, independent, cov_n, med_beta, TIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta"))

## Combining the datasets 
HDMA_density <- HDMA_density %>% 
  mutate(method = "HDMA")
MITM_density <- MITM_density %>% 
  mutate(method = "MITM")
HIMA_density <- HIMA_density %>% 
  mutate(method = "HIMA")

combined_density <- bind_rows(HDMA_density, MITM_density, HIMA_density)

# violin plot 

ggplot(combined_density, aes(x = factor(group_id), y = indirect, fill = method)) +
  geom_violin(alpha = 0.7) +
  geom_hline(aes(yintercept = TIE), color = "red", linetype = "dashed") +  # Use the truth (TIE) for each group_id
  labs(title = "TIE, p=600",
       x = "Simulation scenario",
       y = "TIE") +
  facet_wrap(~group_id, scales = "free") +  
  scale_x_discrete(labels = group_names) +  # Add custom labels
  scale_fill_manual(values = c("HDMA" = "blue", "MITM" = "green", "HIMA" = "purple"), name = "Method") +
  theme_minimal()

# Table data

summary_table_600 <- combined_density %>%
  group_by(group_id, method) %>%
  summarise(
    mean_indirect = mean(indirect, na.rm = TRUE),
    sd_indirect = sd(indirect, na.rm = TRUE),
    tie_value = unique(TIE),  # Assuming `tie` is consistent within group_id
    bias = mean_indirect - tie_value
  ) %>%
  # Create the bias column with tie value in parentheses
  mutate(
    bias = paste0(round(bias, 3), " (", tie_value, ")"),
    mean_sd_indirect = paste0(round(mean_indirect, 3), " (", round(sd_indirect, 3), ")"),
    group_id = factor(group_id, levels = names(group_names), labels = group_names)  # Replace IDs with names
  ) %>%
  select(group_id, method, mean_sd_indirect, bias)


##### Combined table #####

group_names <- c(
  "1" = "n=200; B=0.1; Med=2%; Ind=no",
  "2" = "n=500; B=0.1; Med=2%; Ind=no",
  "3" = "n=1000; B=0.1; Med=2%; Ind=no",
  "4" = "n=200; B=0.1; Med=5%; Ind=no",
  "5" = "n=500; B=0.1; Med=5%; Ind=no",
  "6" = "n=1000; B=0.1; Med=5%; Ind=no",
  "7" = "n=200; B=0.1; Med=10%; Ind=no",
  "8" = "n=500; B=0.1; Med=10%; Ind=no",
  "9" = "n=1000; B=0.1; Med=10%; Ind=no",
  "10" = "n=200; B=0.3; Med=2%; Ind=no",
  "11" = "n=500; B=0.3; Med=2%; Ind=no",
  "12" = "n=1000; B=0.3; Med=2%; Ind=no",
  "13" = "n=200; B=0.3; Med=5%; Ind=no",
  "14" = "n=500; B=0.3; Med=5%; Ind=no",
  "15" = "n=1000; B=0.3; Med=5%; Ind=no",
  "16" = "n=200; B=0.3; Med=10%; Ind=no",
  "17" = "n=500; B=0.3; Med=10%; Ind=no",
  "18" = "n=1000; B=0.3; Med=10%; Ind=no",
  "19" = "n=200; B=0.1; Med=2%; Ind=yes",
  "20" = "n=500; B=0.1; Med=2%; Ind=yes",
  "21" = "n=1000; B=0.1; Med=2%; Ind=yes",
  "22" = "n=200; B=0.1; Med=5%; Ind=yes",
  "23" = "n=500; B=0.1; Med=5%; Ind=yes",
  "24" = "n=1000; B=0.1; Med=5%; Ind=yes",
  "25" = "n=200; B=0.1; Med=10%; Ind=yes",
  "26" = "n=500; B=0.1; Med=10%; Ind=yes",
  "27" = "n=1000; B=0.1; Med=10%; Ind=yes",
  "28" = "n=200; B=0.3; Med=2%; Ind=yes",
  "29" = "n=500; B=0.3; Med=2%; Ind=yes",
  "30" = "n=1000; B=0.3; Med=2%; Ind=yes",
  "31" = "n=200; B=0.3; Med=5%; Ind=yes",
  "32" = "n=500; B=0.3; Med=5%; Ind=yes",
  "33" = "n=1000; B=0.3; Med=5%; Ind=yes",
  "34" = "n=200; B=0.3; Med=10%; Ind=yes",
  "35" = "n=500; B=0.3; Med=10%; Ind=yes",
  "36" = "n=1000; B=0.3; Med=10%; Ind=yes"
)

## pulling from the summary tables created for each P 

summary_table_combined <- summary_table_200 %>% 
  rename(c(meansd_200 = mean_sd_indirect,
           bias_200 = bias)) %>% 
  cbind(summary_table_400[,3:4]) %>% 
  rename(c(meansd_400 = mean_sd_indirect,
           bias_400 = bias)) %>% 
  cbind(summary_table_600[,3:4]) %>% 
  rename(c(meansd_600 = mean_sd_indirect,
           bias_600 = bias))

tie <- flextable(summary_table_combined) %>%
  set_header_labels(
    group_id = "Simulation Scenario",
    method = "Method",
    meansd_200 = "Mean (SD)",
    bias_200 = "Bias (Truth)",
    meansd_400 = "Mean (SD)",
    bias_400 = "Bias (Truth)",
    meansd_600 = "Mean (SD)",
    bias_600 = "Bias (Truth)"
  ) %>%
  add_header_row(
    values = c("", "", "p=200", "p=400", "p=600"),
    colwidths = c(1, 1, 2, 2, 2)  # Indicating how many columns each header spans
  ) %>%
  merge_v(j = "group_id") %>%
  theme_vanilla() %>%
  align(align = "center", part = "all") %>%
  set_table_properties(width = 0.5, layout = "autofit") %>%
  autofit()

tie

save_as_docx(tie, path = "~/Downloads/tie_table.docx")


#### CIE violin plots ####

group_names <- c(
  "1" = "n=200; B=0.1; Med=2%",
  "2" = "n=500; B=0.1; Med=2%",
  "3" = "n=1000; B=0.1; Med=2%",
  "4" = "n=200; B=0.1; Med=5%",
  "5" = "n=500; B=0.1; Med=5%",
  "6" = "n=1000; B=0.1; Med=5%",
  "7" = "n=200; B=0.1; Med=10%",
  "8" = "n=500; B=0.1; Med=10%",
  "9" = "n=1000; B=0.1; Med=10%",
  "10" = "n=200; B=0.3; Med=2%",
  "11" = "n=500; B=0.3; Med=2%",
  "12" = "n=1000; B=0.3; Med=2%",
  "13" = "n=200; B=0.3; Med=5%",
  "14" = "n=500; B=0.3; Med=5%",
  "15" = "n=1000; B=0.3; Med=5%",
  "16" = "n=200; B=0.3; Med=10%",
  "17" = "n=500; B=0.3; Med=10%",
  "18" = "n=1000; B=0.3; Med=10%",
  "19" = "n=200; B=0.1; Med=2%",
  "20" = "n=500; B=0.1; Med=2%",
  "21" = "n=1000; B=0.1; Med=2%",
  "22" = "n=200; B=0.1; Med=5%",
  "23" = "n=500; B=0.1; Med=5%",
  "24" = "n=1000; B=0.1; Med=5%",
  "25" = "n=200; B=0.1; Med=10%",
  "26" = "n=500; B=0.1; Med=10%",
  "27" = "n=1000; B=0.1; Med=10%",
  "28" = "n=200; B=0.3; Med=2%",
  "29" = "n=500; B=0.3; Med=2%",
  "30" = "n=1000; B=0.3; Med=2%",
  "31" = "n=200; B=0.3; Med=5%",
  "32" = "n=500; B=0.3; Med=5%",
  "33" = "n=1000; B=0.3; Med=5%",
  "34" = "n=200; B=0.3; Med=10%",
  "35" = "n=500; B=0.3; Med=10%",
  "36" = "n=1000; B=0.3; Med=10%"
)

#### p = 200 ####

truth_200 <- read_csv("truth_200.csv", col_types = cols(...1 = col_skip()))

# list files in the output folder 
files_200 <- list.files(pattern = "^(HDMA200|HIMA200|MITM200).*\\.csv$")

# loading in the files and adding in the names 
data_files <- lapply(files_200, read.csv)

names(data_files) <- files_200

# combining the results from the different runs 
HDMA200 <- do.call(rbind, data_files[1:5]) 
HDMA200_sesp <- do.call(rbind, data_files[6:10]) 

HIMA200 <- do.call(rbind, data_files[11:15])
HIMA200_sesp <- do.call(rbind, data_files[16:20]) 

MITM200 <- do.call(rbind, data_files[21:25]) 
MITM200_sesp <- do.call(rbind, data_files[26:30]) 

# Editing datasets for visualizations 

HDMA_density <- HDMA200 %>% 
  select(sim_id, alpha_beta, n, p, independent, cov_n, truth, med_beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_200 %>% 
              select(p, truth, independent, cov_n, med_beta, CIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta")) %>% 
  filter(!is.na(CIE))

HIMA_density <- HIMA200 %>% 
  select(sim_id, alpha.beta, n, p, independent, cov_n, truth, med_beta) %>% 
  rename(alpha_beta = alpha.beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_200 %>% 
              select(p, truth, independent, cov_n, med_beta, CIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta")) %>% 
  filter(!is.na(CIE))

MITM_density <- MITM200 %>% 
  select(sim_id, alphabeta, n, p, independent, cov_n, truth, med_beta) %>% 
  rename(alpha_beta = alphabeta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_200 %>% 
              select(p, truth, independent, cov_n, med_beta, CIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta")) %>% 
  filter(!is.na(CIE))

## Combining the datasets 
HDMA_density <- HDMA_density %>% 
  mutate(method = "HDMA")
MITM_density <- MITM_density %>% 
  mutate(method = "MITM")
HIMA_density <- HIMA_density %>% 
  mutate(method = "HIMA")

combined_density <- bind_rows(HDMA_density, MITM_density, HIMA_density)


# violin plot 

ggplot(combined_density, aes(x = factor(group_id), y = alpha_beta, fill = method)) +
  geom_violin(alpha = 0.7) +
  geom_hline(aes(yintercept = CIE), color = "red", linetype = "dashed") +
  labs(title = "CIE, p=200",
       x = "Simulation scenario",
       y = "CIE") +
  facet_wrap(~group_id, scales = "free") +  
  scale_x_discrete(labels = group_names) +# Add custom labels
  scale_fill_manual(values = c("HDMA" = "blue", "MITM" = "green", "HIMA" = "purple"), name = "Method") +
  theme_minimal()

# Table data

summary_table_200 <- combined_density %>%
  group_by(group_id, method) %>%
  summarise(
    mean_cie = mean(alpha_beta, na.rm = TRUE),
    sd_cie = sd(alpha_beta, na.rm = TRUE),
    cie_value = unique(CIE),  # Assuming `tie` is consistent within group_id
    bias = mean_cie - cie_value
  ) %>%
  # Create the bias column with tie value in parentheses
  mutate(
    bias = paste0(round(bias, 3), " (", cie_value, ")"),
    mean_sd_cie = paste0(round(mean_cie, 3), " (", round(sd_cie, 3), ")"),
    group_id = factor(group_id, levels = names(group_names), labels = group_names)  # Replace IDs with names
  ) %>%
  select(group_id, method, mean_sd_cie, bias)

#### p = 400 ####

truth_400 <- read_csv("truth_400.csv", col_types = cols(...1 = col_skip()))

# list files in the output folder 
files_400 <- list.files(pattern = "^(HDMA400|HIMA400|MITM400).*\\.csv$")

# loading in the files and adding in the names 
data_files <- lapply(files_400, read.csv)

names(data_files) <- files_400

# combining the results from the different runs 
HDMA400 <- do.call(rbind, data_files[1:5]) 
HDMA400_sesp <- do.call(rbind, data_files[6:10]) 

HIMA400 <- do.call(rbind, data_files[11:15])
HIMA400_sesp <- do.call(rbind, data_files[16:20]) 

MITM400 <- do.call(rbind, data_files[21:25]) 
MITM400_sesp <- do.call(rbind, data_files[26:30]) 

# Editing datasets for visualizations 

HDMA_density <- HDMA400 %>% 
  select(sim_id, alpha_beta, n, p, independent, cov_n, truth, med_beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_400 %>% 
              select(p, truth, independent, cov_n, med_beta, CIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta")) %>% 
  filter(!is.na(CIE))

HIMA_density <- HIMA400 %>% 
  select(sim_id, alpha.beta, n, p, independent, cov_n, truth, med_beta) %>% 
  rename(alpha_beta = alpha.beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_400 %>% 
              select(p, truth, independent, cov_n, med_beta, CIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta")) %>% 
  filter(!is.na(CIE))

MITM_density <- MITM400 %>% 
  select(sim_id, alphabeta, n, p, independent, cov_n, truth, med_beta) %>% 
  rename(alpha_beta = alphabeta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_400 %>% 
              select(p, truth, independent, cov_n, med_beta, CIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta")) %>% 
  filter(!is.na(CIE))

## Combining the datasets 
HDMA_density <- HDMA_density %>% 
  mutate(method = "HDMA")
MITM_density <- MITM_density %>% 
  mutate(method = "MITM")
HIMA_density <- HIMA_density %>% 
  mutate(method = "HIMA")

combined_density <- bind_rows(HDMA_density, MITM_density, HIMA_density)


# violin plot 

ggplot(combined_density, aes(x = factor(group_id), y = alpha_beta, fill = method)) +
  geom_violin(alpha = 0.7) +
  geom_hline(aes(yintercept = CIE), color = "red", linetype = "dashed") +
  labs(title = "CIE, p=400",
       x = "Simulation scenario",
       y = "CIE") +
  facet_wrap(~group_id, scales = "free") +  
  scale_x_discrete(labels = group_names) +# Add custom labels
  scale_fill_manual(values = c("HDMA" = "blue", "MITM" = "green", "HIMA" = "purple"), name = "Method") +
  theme_minimal()

# Table data 

summary_table_400 <- combined_density %>%
  group_by(group_id, method) %>%
  summarise(
    mean_cie = mean(alpha_beta, na.rm = TRUE),
    sd_cie = sd(alpha_beta, na.rm = TRUE),
    cie_value = unique(CIE),  # Assuming `tie` is consistent within group_id
    bias = mean_cie - cie_value
  ) %>%
  # Create the bias column with tie value in parentheses
  mutate(
    bias = paste0(round(bias, 3), " (", cie_value, ")"),
    mean_sd_cie = paste0(round(mean_cie, 3), " (", round(sd_cie, 3), ")"),
    group_id = factor(group_id, levels = names(group_names), labels = group_names)  # Replace IDs with names
  ) %>%
  select(group_id, method, mean_sd_cie, bias)


#### p = 600 ####

truth_600 <- read_csv("truth_600.csv", col_types = cols(...1 = col_skip()))

# list files in the output folder 
files_600 <- list.files(pattern = "^(HDMA600|HIMA600|MITM600).*\\.csv$")

# loading in the files and adding in the names 
data_files <- lapply(files_600, read.csv)

names(data_files) <- files_600

# combining the results from the different runs 
HDMA600 <- do.call(rbind, data_files[1:5]) 
HDMA600_sesp <- do.call(rbind, data_files[6:10]) 

HIMA600 <- do.call(rbind, data_files[11:15])
HIMA600_sesp <- do.call(rbind, data_files[16:20]) 

MITM600 <- do.call(rbind, data_files[21:25]) 
MITM600_sesp <- do.call(rbind, data_files[26:30]) 

# Editing datasets for visualizations 

HDMA_density <- HDMA600 %>% 
  select(sim_id, alpha_beta, n, p, independent, cov_n, truth, med_beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_600 %>% 
              select(p, truth, independent, cov_n, med_beta, CIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta")) %>% 
  filter(!is.na(CIE))

HIMA_density <- HIMA600 %>% 
  select(sim_id, alpha.beta, n, p, independent, cov_n, truth, med_beta) %>% 
  rename(alpha_beta = alpha.beta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_600 %>% 
              select(p, truth, independent, cov_n, med_beta, CIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta")) %>% 
  filter(!is.na(CIE))

MITM_density <- MITM600 %>% 
  select(sim_id, alphabeta, n, p, independent, cov_n, truth, med_beta) %>% 
  rename(alpha_beta = alphabeta) %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>%
  left_join(truth_600 %>% 
              select(p, truth, independent, cov_n, med_beta, CIE), 
            by = c("p", "truth", "independent", "cov_n", "med_beta")) %>% 
  filter(!is.na(CIE))

## Combining the datasets 
HDMA_density <- HDMA_density %>% 
  mutate(method = "HDMA")
MITM_density <- MITM_density %>% 
  mutate(method = "MITM")
HIMA_density <- HIMA_density %>% 
  mutate(method = "HIMA")

combined_density <- bind_rows(HDMA_density, MITM_density, HIMA_density)


# violin plot 

ggplot(combined_density, aes(x = factor(group_id), y = alpha_beta, fill = method)) +
  geom_violin(alpha = 0.7) +
  geom_hline(aes(yintercept = CIE), color = "red", linetype = "dashed") +
  labs(title = "CIE, p=600",
       x = "Simulation scenario",
       y = "CIE") +
  facet_wrap(~group_id, scales = "free") +  
  scale_x_discrete(labels = group_names) +# Add custom labels
  scale_fill_manual(values = c("HDMA" = "blue", "MITM" = "green", "HIMA" = "purple"), name = "Method") +
  theme_minimal()

# Table data 

summary_table_600 <- combined_density %>%
  group_by(group_id, method) %>%
  summarise(
    mean_cie = mean(alpha_beta, na.rm = TRUE),
    sd_cie = sd(alpha_beta, na.rm = TRUE),
    cie_value = unique(CIE),  # Assuming `tie` is consistent within group_id
    bias = mean_cie - cie_value
  ) %>%
  # Create the bias column with tie value in parentheses
  mutate(
    bias = paste0(round(bias, 3), " (", cie_value, ")"),
    mean_sd_cie = paste0(round(mean_cie, 3), " (", round(sd_cie, 3), ")"),
    group_id = factor(group_id, levels = names(group_names), labels = group_names)  # Replace IDs with names
  ) %>%
  select(group_id, method, mean_sd_cie, bias)


##### Combined table #####

group_names <- c(
  "1" = "n=200; B=0.1; Med=2%; Ind=no",
  "2" = "n=500; B=0.1; Med=2%; Ind=no",
  "3" = "n=1000; B=0.1; Med=2%; Ind=no",
  "4" = "n=200; B=0.1; Med=5%; Ind=no",
  "5" = "n=500; B=0.1; Med=5%; Ind=no",
  "6" = "n=1000; B=0.1; Med=5%; Ind=no",
  "7" = "n=200; B=0.1; Med=10%; Ind=no",
  "8" = "n=500; B=0.1; Med=10%; Ind=no",
  "9" = "n=1000; B=0.1; Med=10%; Ind=no",
  "10" = "n=200; B=0.3; Med=2%; Ind=no",
  "11" = "n=500; B=0.3; Med=2%; Ind=no",
  "12" = "n=1000; B=0.3; Med=2%; Ind=no",
  "13" = "n=200; B=0.3; Med=5%; Ind=no",
  "14" = "n=500; B=0.3; Med=5%; Ind=no",
  "15" = "n=1000; B=0.3; Med=5%; Ind=no",
  "16" = "n=200; B=0.3; Med=10%; Ind=no",
  "17" = "n=500; B=0.3; Med=10%; Ind=no",
  "18" = "n=1000; B=0.3; Med=10%; Ind=no",
  "19" = "n=200; B=0.1; Med=2%; Ind=yes",
  "20" = "n=500; B=0.1; Med=2%; Ind=yes",
  "21" = "n=1000; B=0.1; Med=2%; Ind=yes",
  "22" = "n=200; B=0.1; Med=5%; Ind=yes",
  "23" = "n=500; B=0.1; Med=5%; Ind=yes",
  "24" = "n=1000; B=0.1; Med=5%; Ind=yes",
  "25" = "n=200; B=0.1; Med=10%; Ind=yes",
  "26" = "n=500; B=0.1; Med=10%; Ind=yes",
  "27" = "n=1000; B=0.1; Med=10%; Ind=yes",
  "28" = "n=200; B=0.3; Med=2%; Ind=yes",
  "29" = "n=500; B=0.3; Med=2%; Ind=yes",
  "30" = "n=1000; B=0.3; Med=2%; Ind=yes",
  "31" = "n=200; B=0.3; Med=5%; Ind=yes",
  "32" = "n=500; B=0.3; Med=5%; Ind=yes",
  "33" = "n=1000; B=0.3; Med=5%; Ind=yes",
  "34" = "n=200; B=0.3; Med=10%; Ind=yes",
  "35" = "n=500; B=0.3; Med=10%; Ind=yes",
  "36" = "n=1000; B=0.3; Med=10%; Ind=yes"
)

## pulling from the summary tables created for each P 

summary_table_combined <- summary_table_200 %>% 
  rename(c(meansd_200 = mean_sd_cie,
           bias_200 = bias)) %>% 
  cbind(summary_table_400[,3:4]) %>% 
  rename(c(meansd_400 = mean_sd_cie,
           bias_400 = bias)) %>% 
  cbind(summary_table_600[,3:4]) %>% 
  rename(c(meansd_600 = mean_sd_cie,
           bias_600 = bias))

cie <- flextable(summary_table_combined) %>%
  set_header_labels(
    group_id = "Simulation Scenario",
    method = "Method",
    meansd_200 = "Mean (SD)",
    bias_200 = "Bias (Truth)",
    meansd_400 = "Mean (SD)",
    bias_400 = "Bias (Truth)",
    meansd_600 = "Mean (SD)",
    bias_600 = "Bias (Truth)"
  ) %>%
  add_header_row(
    values = c("", "", "p=200", "p=400", "p=600"),
    colwidths = c(1, 1, 2, 2, 2)  # Indicating how many columns each header spans
  ) %>%
  merge_v(j = "group_id") %>%
  theme_vanilla() %>%
  align(align = "center", part = "all") %>%
  set_table_properties(width = 0.5, layout = "autofit") %>%
  autofit()

cie

save_as_docx(cie, path = "~/Downloads/cie_table.docx")




#### Se/Sp #### 


group_names <- c(
  "1" = "1. n=200; B=0.1; Med=2%",
  "2" = "2. n=500; B=0.1; Med=2%",
  "3" = "3. n=1000; B=0.1; Med=2%",
  "4" = "4. n=200; B=0.1; Med=5%",
  "5" = "5. n=500; B=0.1; Med=5%",
  "6" = "6. n=1000; B=0.1; Med=5%",
  "7" = "7. n=200; B=0.1; Med=10%",
  "8" = "8. n=500; B=0.1; Med=10%",
  "9" = "9. n=1000; B=0.1; Med=10%",
  "10" = "10. n=200; B=0.3; Med=2%",
  "11" = "11. n=500; B=0.3; Med=2%",
  "12" = "12. n=1000; B=0.3; Med=2%",
  "13" = "13. n=200; B=0.3; Med=5%",
  "14" = "14. n=500; B=0.3; Med=5%",
  "15" = "15. n=1000; B=0.3; Med=5%",
  "16" = "16. n=200; B=0.3; Med=10%",
  "17" = "17. n=500; B=0.3; Med=10%",
  "18" = "18. n=1000; B=0.3; Med=10%",
  "19" = "19. n=200; B=0.1; Med=2%",
  "20" = "20. n=500; B=0.1; Med=2%",
  "21" = "21. n=1000; B=0.1; Med=2%",
  "22" = "22. n=200; B=0.1; Med=5%",
  "23" = "23. n=500; B=0.1; Med=5%",
  "24" = "24. n=1000; B=0.1; Med=5%",
  "25" = "25. n=200; B=0.1; Med=10%",
  "26" = "26. n=500; B=0.1; Med=10%",
  "27" = "27. n=1000; B=0.1; Med=10%",
  "28" = "28. n=200; B=0.3; Med=2%",
  "29" = "29. n=500; B=0.3; Med=2%",
  "30" = "30. n=1000; B=0.3; Med=2%",
  "31" = "31. n=200; B=0.3; Med=5%",
  "32" = "32. n=500; B=0.3; Med=5%",
  "33" = "33. n=1000; B=0.3; Med=5%",
  "34" = "34. n=200; B=0.3; Med=10%",
  "35" = "35. n=500; B=0.3; Med=10%",
  "36" = "36. n=1000; B=0.3; Med=10%"
)

#### p = 200 ####


se_sp <- function(data){
  se <- numeric(nrow(data))
  sp <- numeric(nrow(data))
  
  for(i in 1:nrow(data)) {
    true_values <- paste0("M.", 1:data[i, "truth"])
    false_values <- paste0("M.", (data[i, "truth"] + 1):ncol(data)) 
    false_values <- false_values[false_values %in% names(data)]
    
    true_positives <- sum(data[i, true_values] == 1)
    false_negatives <- length(true_values) - true_positives
    
    true_negatives <- sum(data[i, false_values] == 0)
    false_positives <- sum(data[i, false_values] == 1)
    
    # Calculate se and sp for each row 
    se[i] <- true_positives/(true_positives+false_negatives)
    sp[i] <- true_negatives/(true_negatives+false_positives)
  }
  
  return(data.frame(
    Sensitivity = se, 
    Specificity = sp,
    sim_id = data$sim_id, 
    n = data$n, 
    p = data$p, 
    truth = data$truth, 
    independent = data$independent,
    cov_n = data$cov_n, 
    med_beta = data$med_beta
  ))
  
}

# Pulling Se and Sp 

HDMA_sesp <- se_sp(HDMA200_sesp)

HIMA_sesp <- se_sp(HIMA200_sesp)

MITM_sesp <- se_sp(MITM200_sesp)

# Removing those from NA runs 

rows_to_remove <- which(HDMA_sesp["Sensitivity"] == 0.00 &
                          HDMA_sesp["Specificity"] == 1.0000000) # 38
HDMA_sesp <- HDMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(HIMA_sesp["Sensitivity"] == 0.00 &
                          HIMA_sesp["Specificity"] == 1.0000000) # 546
HIMA_sesp <- HIMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(MITM_sesp["Sensitivity"] == 0.00 &
                          MITM_sesp["Specificity"] == 1.0000000) # 972
MITM_sesp <- MITM_sesp[-rows_to_remove, ]

HDMA_se <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

HDMA_sp <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

HIMA_se <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>%
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

HIMA_sp <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

MITM_se <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>%
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

MITM_sp <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

## Se visualization 

combined_data_se <- bind_rows(
  HDMA_se %>% mutate(method = "HDMA"),
  HIMA_se %>% mutate(method = "HIMA"),
  MITM_se %>% mutate(method = "MITM")
)

####### Boxplot option 
# ggplot(combined_data_se, aes(x = factor(group_id), y = Sensitivity, fill = factor(group_id))) +
#   geom_boxplot() +
#   labs(title = "Sensitivity, p=200", 
#        x = "Simulation scenario", y = "Sensitivity") +
#   scale_fill_discrete(labels = group_names, name = "Group Details") +
#   facet_wrap(~ method, scales = "free_y") +  # Separate facet for each method
#   theme_minimal() +
#   theme(legend.position = "bottom")

summary_data <- combined_data_se %>%
  group_by(group_id, method) %>%
  summarise(mean_sensitivity = mean(Sensitivity, na.rm = TRUE),
            sd_sensitivity = sd(Sensitivity, na.rm = TRUE))

ggplot(summary_data, aes(x = factor(group_id), y = mean_sensitivity, group = method, color = method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_sensitivity - sd_sensitivity, 
                    ymax = mean_sensitivity + sd_sensitivity), width = 0.2) +
  labs(title = "Sensitivity, p=200", 
       x = "Simulation scenario", y = "Sensitivity") +
  scale_color_discrete(name = "Method") +
  scale_x_discrete(labels = group_names) +  # Apply custom labels here
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, hjust = 1))


## Sp visualization 

combined_data_sp <- bind_rows(
  HDMA_sp %>% mutate(method = "HDMA"),
  HIMA_sp %>% mutate(method = "HIMA"),
  MITM_sp %>% mutate(method = "MITM")
)

summary_data <- combined_data_sp %>%
  group_by(group_id, method) %>%
  summarise(mean_specificity = mean(Specificity, na.rm = TRUE),
            sd_specificity = sd(Specificity, na.rm = TRUE))

ggplot(summary_data, aes(x = factor(group_id), y = mean_specificity, group = method, color = method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_specificity - sd_specificity, 
                    ymax = mean_specificity + sd_specificity), width = 0.2) +
  labs(title = "Specificity, p=200", 
       x = "Simulation scenario", y = "Specificity") +
  scale_color_discrete(name = "Method") +
  scale_x_discrete(labels = group_names) +  # Apply custom labels here
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, hjust = 1))

#### Se/Sp Table 

table_data <- combined_data_se %>% 
  cbind(combined_data_sp$Specificity) %>% 
  rename(Specificity = `...11`)


# table data 
summary_table_200 <- table_data %>%
  group_by(group_id, method) %>%
  summarise(
    mean_sensitivity = mean(Sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(Sensitivity, na.rm = TRUE),
    mean_specificty = mean(Specificity, na.rm = TRUE),
    sd_specificity = sd(Specificity, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    mean_sensitivity = round(mean_sensitivity, 2),
    sd_sensitivity = round(sd_sensitivity, 2),
    mean_specificty = round(mean_specificty, 2),
    sd_specificity = round(sd_specificity, 2),
    group_id = factor(group_id, levels = names(group_names), labels = group_names)  # Replace IDs with names
  )

# p=200 se table 
sesp_table <- flextable(summary_table) %>%
  set_header_labels(
    group_id = "Simulation Scenario",
    method = "Method",
    mean_sensitivity = "Mean Sensitivity",
    sd_sensitivity = "SD Sensitivity",
    mean_specificty = "Mean Specificity",
    sd_specificity = "SD Specificity"
  ) %>%
  set_table_properties(width = 0.5, layout = "autofit") %>%
  theme_vanilla() %>%
  align(align = "center", part = "all") %>%  # Align all text to center
  merge_v(j = "group_id")

print(sesp_table, preview = "html")

save_as_docx(sesp_table, path = "~/Downloads/p200_sesp_table.docx")

#### p = 400 ####

se_sp <- function(data){
  se <- numeric(nrow(data))
  sp <- numeric(nrow(data))
  
  for(i in 1:nrow(data)) {
    true_values <- paste0("M.", 1:data[i, "truth"])
    false_values <- paste0("M.", (data[i, "truth"] + 1):ncol(data)) 
    false_values <- false_values[false_values %in% names(data)]
    
    true_positives <- sum(data[i, true_values] == 1)
    false_negatives <- length(true_values) - true_positives
    
    true_negatives <- sum(data[i, false_values] == 0)
    false_positives <- sum(data[i, false_values] == 1)
    
    # Calculate se and sp for each row 
    se[i] <- true_positives/(true_positives+false_negatives)
    sp[i] <- true_negatives/(true_negatives+false_positives)
  }
  
  return(data.frame(
    Sensitivity = se, 
    Specificity = sp,
    sim_id = data$sim_id, 
    n = data$n, 
    p = data$p, 
    truth = data$truth, 
    independent = data$independent,
    cov_n = data$cov_n, 
    med_beta = data$med_beta
  ))
  
}

# Pulling Se and Sp 

HDMA_sesp <- se_sp(HDMA400_sesp)

HIMA_sesp <- se_sp(HIMA400_sesp)

MITM_sesp <- se_sp(MITM400_sesp)

# Removing those from NA runs 

rows_to_remove <- which(HDMA_sesp["Sensitivity"] == 0.00 &
                          HDMA_sesp["Specificity"] == 1.0000000) # 38
HDMA_sesp <- HDMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(HIMA_sesp["Sensitivity"] == 0.00 &
                          HIMA_sesp["Specificity"] == 1.0000000) # 546
HIMA_sesp <- HIMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(MITM_sesp["Sensitivity"] == 0.00 &
                          MITM_sesp["Specificity"] == 1.0000000) # 972
MITM_sesp <- MITM_sesp[-rows_to_remove, ]

HDMA_se <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

HDMA_sp <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

HIMA_se <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

HIMA_sp <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

MITM_se <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

MITM_sp <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

## Se visualization 

combined_data_se <- bind_rows(
  HDMA_se %>% mutate(method = "HDMA"),
  HIMA_se %>% mutate(method = "HIMA"),
  MITM_se %>% mutate(method = "MITM")
)

summary_data <- combined_data_se %>%
  group_by(group_id, method) %>%
  summarise(mean_sensitivity = mean(Sensitivity, na.rm = TRUE),
            sd_sensitivity = sd(Sensitivity, na.rm = TRUE))

ggplot(summary_data, aes(x = factor(group_id), y = mean_sensitivity, group = method, color = method)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_sensitivity - sd_sensitivity, 
                    ymax = mean_sensitivity + sd_sensitivity), width = 0.2) +
  labs(title = "Sensitivity, p=400", 
       x = "Simulation scenario", y = "Sensitivity") +
  scale_color_discrete(name = "Method") +
  scale_x_discrete(labels = group_names) +  # Apply custom labels here
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))


## Sp visualization 

combined_data_sp <- bind_rows(
  HDMA_sp %>% mutate(method = "HDMA"),
  HIMA_sp %>% mutate(method = "HIMA"),
  MITM_sp %>% mutate(method = "MITM")
)

summary_data <- combined_data_sp %>%
  group_by(group_id, method) %>%
  summarise(mean_specificity = mean(Specificity, na.rm = TRUE),
            sd_specificity = sd(Specificity, na.rm = TRUE))

ggplot(summary_data, aes(x = factor(group_id), y = mean_specificity, group = method, color = method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_specificity - sd_specificity, 
                    ymax = mean_specificity + sd_specificity), width = 0.2) +
  labs(title = "Specificity, p=400", 
       x = "Simulation scenario", y = "Specificity") +
  scale_color_discrete(name = "Method") +
  scale_x_discrete(labels = group_names) +  # Apply custom labels here
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))

#### Se/Sp Table 

table_data <- combined_data_se %>% 
  cbind(combined_data_sp$Specificity) %>% 
  rename(Specificity = `...11`)


# table data 
summary_table_400 <- table_data %>%
  group_by(group_id, method) %>%
  summarise(
    mean_sensitivity = mean(Sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(Sensitivity, na.rm = TRUE),
    mean_specificty = mean(Specificity, na.rm = TRUE),
    sd_specificity = sd(Specificity, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    mean_sensitivity = round(mean_sensitivity, 2),
    sd_sensitivity = round(sd_sensitivity, 2),
    mean_specificty = round(mean_specificty, 2),
    sd_specificity = round(sd_specificity, 2),
    group_id = factor(group_id, levels = names(group_names), labels = group_names)  # Replace IDs with names
  )

# p=400 se table 
sesp_table <- flextable(summary_table) %>%
  set_header_labels(
    group_id = "Simulation Scenario",
    method = "Method",
    mean_sensitivity = "Mean Sensitivity",
    sd_sensitivity = "SD Sensitivity",
    mean_specificty = "Mean Specificity",
    sd_specificity = "SD Specificity"
  ) %>%
  set_table_properties(width = 0.5, layout = "autofit") %>%
  theme_vanilla() %>%
  align(align = "center", part = "all") %>%  # Align all text to center
  merge_v(j = "group_id")

print(sesp_table, preview = "html")

save_as_docx(sesp_table, path = "~/Downloads/p400_sesp_table.docx")

#### p = 600 ####

se_sp <- function(data){
  se <- numeric(nrow(data))
  sp <- numeric(nrow(data))
  
  for(i in 1:nrow(data)) {
    true_values <- paste0("M.", 1:data[i, "truth"])
    false_values <- paste0("M.", (data[i, "truth"] + 1):ncol(data)) 
    false_values <- false_values[false_values %in% names(data)]
    
    true_positives <- sum(data[i, true_values] == 1)
    false_negatives <- length(true_values) - true_positives
    
    true_negatives <- sum(data[i, false_values] == 0)
    false_positives <- sum(data[i, false_values] == 1)
    
    # Calculate se and sp for each row 
    se[i] <- true_positives/(true_positives+false_negatives)
    sp[i] <- true_negatives/(true_negatives+false_positives)
  }
  
  return(data.frame(
    Sensitivity = se, 
    Specificity = sp,
    sim_id = data$sim_id, 
    n = data$n, 
    p = data$p, 
    truth = data$truth, 
    independent = data$independent,
    cov_n = data$cov_n, 
    med_beta = data$med_beta
  ))
  
}

# Pulling Se and Sp 

HDMA_sesp <- se_sp(HDMA600_sesp)

HIMA_sesp <- se_sp(HIMA600_sesp)

MITM_sesp <- se_sp(MITM600_sesp)

# Removing those from NA runs 

rows_to_remove <- which(HDMA_sesp["Sensitivity"] == 0.00 &
                          HDMA_sesp["Specificity"] == 1.0000000) # 38
aHDMA_sesp <- HDMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(HIMA_sesp["Sensitivity"] == 0.00 &
                          HIMA_sesp["Specificity"] == 1.0000000) # 546
HIMA_sesp <- HIMA_sesp[-rows_to_remove, ]

rows_to_remove <- which(MITM_sesp["Sensitivity"] == 0.00 &
                          MITM_sesp["Specificity"] == 1.0000000) # 972
MITM_sesp <- MITM_sesp[-rows_to_remove, ]

HDMA_se <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

HDMA_sp <- HDMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

HIMA_se <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

HIMA_sp <- HIMA_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

MITM_se <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Specificity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

MITM_sp <- MITM_sesp %>% 
  group_by(independent, med_beta, truth, n) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id, sim_id, .keep_all = T) %>% 
  group_by(group_id) %>% 
  select(-Sensitivity) %>% 
  mutate(independent = ifelse(independent == 1, "yes", "no")) 

## Se visualization 

combined_data_se <- bind_rows(
  HDMA_se %>% mutate(method = "HDMA"),
  HIMA_se %>% mutate(method = "HIMA"),
  MITM_se %>% mutate(method = "MITM")
)

summary_data <- combined_data_se %>%
  group_by(group_id, method) %>%
  summarise(mean_sensitivity = mean(Sensitivity, na.rm = TRUE),
            sd_sensitivity = sd(Sensitivity, na.rm = TRUE))

ggplot(summary_data, aes(x = factor(group_id), y = mean_sensitivity, group = method, color = method)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_sensitivity - sd_sensitivity, 
                    ymax = mean_sensitivity + sd_sensitivity), width = 0.2) +
  labs(title = "Sensitivity, p=600", 
       x = "Simulation scenario", y = "Sensitivity") +
  scale_color_discrete(name = "Method") +
  scale_x_discrete(labels = group_names) +  # Apply custom labels here
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))


## Sp visualization 

combined_data_sp <- bind_rows(
  HDMA_sp %>% mutate(method = "HDMA"),
  HIMA_sp %>% mutate(method = "HIMA"),
  MITM_sp %>% mutate(method = "MITM")
)

summary_data <- combined_data_sp %>%
  group_by(group_id, method) %>%
  summarise(mean_specificity = mean(Specificity, na.rm = TRUE),
            sd_specificity = sd(Specificity, na.rm = TRUE))

ggplot(summary_data, aes(x = factor(group_id), y = mean_specificity, group = method, color = method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_specificity - sd_specificity, 
                    ymax = mean_specificity + sd_specificity), width = 0.2) +
  labs(title = "Specificity, p=600", 
       x = "Simulation scenario", y = "Specificity") +
  scale_color_discrete(name = "Method") +
  scale_x_discrete(labels = group_names) +  # Apply custom labels here
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))

#### Se/Sp Table 

table_data <- combined_data_se %>% 
  cbind(combined_data_sp$Specificity) %>% 
  rename(Specificity = `...11`)


# table data 
summary_table_600 <- table_data %>%
  group_by(group_id, method) %>%
  summarise(
    mean_sensitivity = mean(Sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(Sensitivity, na.rm = TRUE),
    mean_specificty = mean(Specificity, na.rm = TRUE),
    sd_specificity = sd(Specificity, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    mean_sensitivity = round(mean_sensitivity, 2),
    sd_sensitivity = round(sd_sensitivity, 2),
    mean_specificty = round(mean_specificty, 2),
    sd_specificity = round(sd_specificity, 2),
    group_id = factor(group_id, levels = names(group_names), labels = group_names)  # Replace IDs with names
  )

# p=600 se table 
sesp_table <- flextable(summary_table) %>%
  set_header_labels(
    group_id = "Simulation Scenario",
    method = "Method",
    mean_sensitivity = "Mean Sensitivity",
    sd_sensitivity = "SD Sensitivity",
    mean_specificty = "Mean Specificity",
    sd_specificity = "SD Specificity"
  ) %>%
  set_table_properties(width = 0.5, layout = "autofit") %>%
  theme_vanilla() %>%
  align(align = "center", part = "all") %>%  # Align all text to center
  merge_v(j = "group_id")

print(sesp_table, preview = "html")

save_as_docx(sesp_table, path = "~/Downloads/p600_sesp_table.docx")


##### Combined Se/Sp table #####

# pulled the data from the summary_table_X00 files created above 

group_names <- c(
  "1" = "n=200; B=0.1; Med=2%; Ind=no",
  "2" = "n=500; B=0.1; Med=2%; Ind=no",
  "3" = "n=1000; B=0.1; Med=2%; Ind=no",
  "4" = "n=200; B=0.1; Med=5%; Ind=no",
  "5" = "n=500; B=0.1; Med=5%; Ind=no",
  "6" = "n=1000; B=0.1; Med=5%; Ind=no",
  "7" = "n=200; B=0.1; Med=10%; Ind=no",
  "8" = "n=500; B=0.1; Med=10%; Ind=no",
  "9" = "n=1000; B=0.1; Med=10%; Ind=no",
  "10" = "n=200; B=0.3; Med=2%; Ind=no",
  "11" = "n=500; B=0.3; Med=2%; Ind=no",
  "12" = "n=1000; B=0.3; Med=2%; Ind=no",
  "13" = "n=200; B=0.3; Med=5%; Ind=no",
  "14" = "n=500; B=0.3; Med=5%; Ind=no",
  "15" = "n=1000; B=0.3; Med=5%; Ind=no",
  "16" = "n=200; B=0.3; Med=10%; Ind=no",
  "17" = "n=500; B=0.3; Med=10%; Ind=no",
  "18" = "n=1000; B=0.3; Med=10%; Ind=no",
  "19" = "n=200; B=0.1; Med=2%; Ind=yes",
  "20" = "n=500; B=0.1; Med=2%; Ind=yes",
  "21" = "n=1000; B=0.1; Med=2%; Ind=yes",
  "22" = "n=200; B=0.1; Med=5%; Ind=yes",
  "23" = "n=500; B=0.1; Med=5%; Ind=yes",
  "24" = "n=1000; B=0.1; Med=5%; Ind=yes",
  "25" = "n=200; B=0.1; Med=10%; Ind=yes",
  "26" = "n=500; B=0.1; Med=10%; Ind=yes",
  "27" = "n=1000; B=0.1; Med=10%; Ind=yes",
  "28" = "n=200; B=0.3; Med=2%; Ind=yes",
  "29" = "n=500; B=0.3; Med=2%; Ind=yes",
  "30" = "n=1000; B=0.3; Med=2%; Ind=yes",
  "31" = "n=200; B=0.3; Med=5%; Ind=yes",
  "32" = "n=500; B=0.3; Med=5%; Ind=yes",
  "33" = "n=1000; B=0.3; Med=5%; Ind=yes",
  "34" = "n=200; B=0.3; Med=10%; Ind=yes",
  "35" = "n=500; B=0.3; Med=10%; Ind=yes",
  "36" = "n=1000; B=0.3; Med=10%; Ind=yes"
)

## pulling from the summary tables created for each P 

summary_table_200f <- summary_table_200 %>% 
  mutate(
    sesd_200 = paste0(round(mean_sensitivity, 3), " (", round(sd_sensitivity, 3), ")"),
    spsd_200 = paste0(round(mean_specificty, 3), " (", round(sd_specificity, 3), ")"),
  ) %>% 
  select(group_id, method, sesd_200, spsd_200)

summary_table_400f <- summary_table_400 %>% 
  mutate(
    sesd_400 = paste0(round(mean_sensitivity, 3), " (", round(sd_sensitivity, 3), ")"),
    spsd_400 = paste0(round(mean_specificty, 3), " (", round(sd_specificity, 3), ")"),
  ) %>% 
  select(group_id, method, sesd_400, spsd_400)

summary_table_600f <- summary_table_600 %>% 
  mutate(
    sesd_600 = paste0(round(mean_sensitivity, 3), " (", round(sd_sensitivity, 3), ")"),
    spsd_600 = paste0(round(mean_specificty, 3), " (", round(sd_specificity, 3), ")"),
  ) %>% 
  select(group_id, method, sesd_600, spsd_600)

summary_table_combined <- summary_table_200f %>% 
  cbind(summary_table_400f[,3:4]) %>% 
  cbind(summary_table_600f[,3:4]) 

sesp <- flextable(summary_table_combined) %>%
  set_header_labels(
    group_id = "Simulation Scenario",
    method = "Method",
    sesd_200 = "Sensitivity (SD)",
    spsd_200 = "Specificity (SD)",
    sesd_400 = "Sensitivity (SD)",
    spsd_400 = "Specificity (SD)",
    sesd_600 = "Sensitivity (SD)",
    spsd_600 = "Specificity (SD)"
  ) %>%
  add_header_row(
    values = c("", "", "p=200", "p=400", "p=600"),
    colwidths = c(1, 1, 2, 2, 2)  # Indicating how many columns each header spans
  ) %>%
  merge_v(j = "group_id") %>%
  theme_vanilla() %>%
  align(align = "center", part = "all") %>%
  set_table_properties(width = 0.5, layout = "autofit") %>%
  autofit()

sesp

save_as_docx(tie, path = "~/Downloads/sesp_table.docx")


