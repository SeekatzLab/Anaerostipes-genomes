### Figure 8 and Supplementary Figure 7: 

# Authors:
# - all figures: DB

## Input files:
## CM3.tsv
## PBM.tsv
## BMAA.tsv
## YCFA.tsv
## sugar_invitro.xlsb.xlsx


## Output files:
## Figure 8A pdf
## Figure 8B pdf
## Figure 8C pdf
## Supplementary Figure 7 pdf

## libraries (common to most)
library(stringr)
library(stringi)
library(readxl)
library(tidyverse)
library(ggplot2)
library(readr)
library(gplots)
library(tidyr)
library(plyr)
library(RColorBrewer)
library(dplyr)
library(janitor)
library(growthcurver)

############################### 

setwd("~/Data/")

ids <- c("02_2", "02_5", "02_91A", "06_28", "03_47", "06_29")

## The following was used to create the original CBM3. 
##  creating a df to incorporate all PBM1 plates
# Initialize empty lists and data frame
# meta_list_PBM1 <- list()
# data_list_PBM1 <- list()
# CM_all_PBM1 <- data.frame()
# 
# for (id in ids) {
#   
#   # Format parts
#   #id_label <- gsub("_", "", id)       # For CM object naming like 0291
#   #file_suffix <- gsub("_", "", id)    # For file suffix like 0291A
#   
#   # File paths
#   meta_file_PBM1 <- paste0("./PBM1_plates/ANA_PBM1_meta_", id, ".xlsx")
#   data_file_PBM1 <- paste0("./PBM1_plates/ANA_PBM1_data_", id, ".xlsx")
#   
#   # Load metadata
#   meta_PBM1 <- read_excel(meta_file_PBM1, col_names = TRUE)
#   
#   # Load and process main data
#   data_PBM1 <- read_excel(data_file_PBM1, skip = 2, col_names = TRUE, n_max = 97) %>%
#     dplyr::rename("time" = 1, "temperature" = 2) %>%
#     select(-temperature) %>%
#     mutate(time = as.numeric(gsub("s", "", time)) / 60 / 60) %>%
#     mutate_if(is.numeric, round, 2) %>%
#     column_to_rownames(var = "time") %>%
#     sweep(., 2, colMeans(.[1:3,])) %>%
#     rownames_to_column(var = "time") %>%
#     mutate(time = as.numeric(time)) %>%
#     pivot_longer(cols = !time, names_to = "wellID", values_to = "OD600") %>%
#     merge(meta_PBM1, ., by = "wellID")
#   
#   # Create summarized CM table
#   CM_PBM1 <- data_PBM1 %>%
#     group_by(group, media, sugar, species, strain, condition, time) %>%
#     dplyr::summarise(
#       mean = mean(OD600, na.rm = TRUE),
#       n = length(OD600),
#       sd = sd(OD600, na.rm = TRUE),
#       se = sd / sqrt(n),
#     ) %>%
#     mutate(sample_id = id)
#   
#   # Store meta and data in lists
#   meta_list_PBM1[[id]] <- meta_PBM1
#   data_list_PBM1[[id]] <- data_PBM1
#   
#   # Append to combined CM_all
#   CM_all_PBM1 <- bind_rows(CM_all_PBM1, CM_PBM1)
# }

## Creating a df for all PBM2 plates

# meta_list_PBM2 <- list()
# data_list_PBM2 <- list()
# CM_all_PBM2 <- data.frame()
# 
# for (id in ids) {
#   
#   # Format parts
#   #id_label <- gsub("_", "", id)       # For CM object naming like 0291
#   #file_suffix <- gsub("_", "", id)    # For file suffix like 0291A
#   
#   # File paths
#   meta_file_PBM2 <- paste0("./PBM2_plates/ANA_PBM2_meta_", id, ".xlsx")
#   data_file_PBM2 <- paste0("./PBM2_plates/ANA_PBM2_data_", id, ".xlsx")
#   
#   # Load metadata
#   meta_PBM2 <- read_excel(meta_file_PBM2, col_names = TRUE)
#   
#   # Load and process main data
#   data_PBM2 <- read_excel(data_file_PBM2, skip = 2, col_names = TRUE, n_max = 97) %>%
#     dplyr::rename("time" = 1, "temperature" = 2) %>%
#     select(-temperature) %>%
#     mutate(time = as.numeric(gsub("s", "", time)) / 60 / 60) %>%
#     mutate_if(is.numeric, round, 2) %>%
#     column_to_rownames(var = "time") %>%
#     sweep(., 2, colMeans(.[1:3,])) %>%
#     rownames_to_column(var = "time") %>%
#     mutate(time = as.numeric(time)) %>%
#     pivot_longer(cols = !time, names_to = "wellID", values_to = "OD600") %>%
#     merge(meta_PBM2, ., by = "wellID")
#   
#   # Create summarized CM table
#   CM_PBM2 <- data_PBM2 %>%
#     group_by(group, media, sugar, species, strain, condition, time) %>%
#     dplyr::summarise(
#       mean = mean(OD600, na.rm = TRUE),
#       n = length(OD600),
#       sd = sd(OD600, na.rm = TRUE),
#       se = sd / sqrt(n),
#     ) %>%
#     mutate(sample_id = id)
#   
#   # Store meta and data in lists
#   meta_list_PBM2[[id]] <- meta_PBM2
#   data_list_PBM2[[id]] <- data_PBM2
#   
#   # Append to combined CM_all
#   CM_all_PBM2 <- bind_rows(CM_all_PBM2, CM_PBM2)
# }
# 
# CM2 <- rbind(CM_all_PBM1, CM_all_PBM2)

# CM3 <- subset(CM2, CM2$sugar %in% c("chondroitin_sodium_sulfate", "sucrose", "d_trehalose", "arabinogalactan", "GOS_betaglucan", "lentinan", "chitosan", "xyloglucan"))
# 
# write_tsv(CM3, "~/Paper4_Anaerostipes/Resubmission/Resubmission/CM3.tsv")
CM3 <- read_tsv("/Data/figure8_supplementary_figure7/CM3.tsv")
CM3_2 <- CM3
CM3_2$mean[CM3_2$mean < 0] <- 0

## Lactate plates were made separately so have different inputs

# For the current L-lactate results, you'll need the BMAA plate 2 data and the YCFA plate 1 data sheets. I will go and check to see if I can recover the rest of the timepoints for the BMAA_plate2_data (I guess we have an issue with converting the Tecan data to an excel sheet?). If not, I can also re-run these conditions on the D-lactate plate. 

setwd("~/Data/")

# YCFA_1_meta <- read_excel("./lactate_YCFA_plate1_meta.xlsx")
# BMAA_1_meta <- read_excel("./lactate_BMAA_plate1_meta.xlsx")
# BMAA_2_meta <- read_excel("./lactate_BMAA_plate2_meta.xlsx")
# YCFA_2_meta <- read_excel("./lactate_combo_plate1_meta.xlsx")
# BMAA_3_meta <- read_excel("./lactate_BMAA_upd_conc_meta.xlsx") 
# 
# YCFA_1_data <- read_excel("./lactate_YCFA_plate1_data.xlsx", skip = 2, col_names = TRUE, n_max = 97) %>%
#   dplyr::rename("time" = 1, "temperature" = 2) %>%
#   select(-temperature) %>%
#   mutate(time = as.numeric(gsub("s", "", time)) / 60 / 60) %>%
#   mutate_if(is.numeric, round, 2) %>%
#   column_to_rownames(var = "time") %>%
#   sweep(., 2, colMeans(.[1:3,])) %>%
#   rownames_to_column(var = "time") %>%
#   mutate(time = as.numeric(time)) %>%
#   pivot_longer(cols = !time, names_to = "wellID", values_to = "OD600") %>%
#   merge(YCFA_1_meta, ., by = "wellID")
# 
# BMAA_1_data <- read_excel("./lactate_BMAA_plate1_data.xlsx", skip = 2, col_names = TRUE, n_max = 97) %>%
#   dplyr::rename("time" = 1, "temperature" = 2) %>%
#   select(-temperature) %>%
#   mutate(time = as.numeric(gsub("s", "", time)) / 60 / 60) %>%
#   mutate_if(is.numeric, round, 2) %>%
#   column_to_rownames(var = "time") %>%
#   sweep(., 2, colMeans(.[1:3,])) %>%
#   rownames_to_column(var = "time") %>%
#   mutate(time = as.numeric(time)) %>%
#   pivot_longer(cols = !time, names_to = "wellID", values_to = "OD600") %>%
#   merge(BMAA_1_meta, ., by = "wellID")
# 
# BMAA_2_data <- read_excel("./lactate_BMAA_plate2_data.xlsx", skip = 2, col_names = TRUE, n_max = 97) %>%
#   dplyr::rename("time" = 1, "temperature" = 2) %>%
#   select(-temperature) %>%
#   mutate(time = as.numeric(gsub("s", "", time)) / 60 / 60) %>%
#   mutate_if(is.numeric, round, 2) %>%
#   column_to_rownames(var = "time") %>%
#   sweep(., 2, colMeans(.[1:3,])) %>%
#   rownames_to_column(var = "time") %>%
#   mutate(time = as.numeric(time)) %>%
#   pivot_longer(cols = !time, names_to = "wellID", values_to = "OD600") %>%
#   merge(BMAA_2_meta, ., by = "wellID")
# 
# YCFA_2_data <- read_excel("./lactate_combo_plate1_data.xlsx", skip = 2, col_names = TRUE, n_max = 97) %>%
#   dplyr::rename("time" = 1, "temperature" = 2) %>%
#   select(-temperature) %>%
#   mutate(time = as.numeric(gsub("s", "", time)) / 60 / 60) %>%
#   mutate_if(is.numeric, round, 2) %>%
#   column_to_rownames(var = "time") %>%
#   sweep(., 2, colMeans(.[1:3,])) %>%
#   rownames_to_column(var = "time") %>%
#   mutate(time = as.numeric(time)) %>%
#   pivot_longer(cols = !time, names_to = "wellID", values_to = "OD600") %>%
#   merge(YCFA_2_meta, ., by = "wellID")
# 
# BMAA_3_data <- read_excel("./lactate_BMAA_upd_conc_data.xlsx", skip = 2, col_names = TRUE, n_max = 97) %>%
#   dplyr::rename("time" = 1, "temperature" = 2) %>%
#   select(-temperature) %>%
#   mutate(time = as.numeric(gsub("s", "", time)) / 60 / 60) %>%
#   mutate_if(is.numeric, round, 2) %>%
#   column_to_rownames(var = "time") %>%
#   sweep(., 2, colMeans(.[1:3,])) %>%
#   rownames_to_column(var = "time") %>%
#   mutate(time = as.numeric(time)) %>%
#   pivot_longer(cols = !time, names_to = "wellID", values_to = "OD600") %>%
#   merge(BMAA_3_meta, ., by = "wellID")
# 
# 
# # BMAA_nobac, BMAA_bac, BMAA_glu, BMAA_d_lac, BMAA_l_lac, BMAA_glu_l_lac, BMAA_glu_d_lac
# BMAA1 <- subset(BMAA_3_data, BMAA_3_data$group %in% c("BMAA", "BMAA_l_lac", "BMAA_d_lac"))
# BM1 <- subset(YCFA_2_data, YCFA_2_data$group %in% c("BMAA_glu_l_lac", "BMAA_glu_d_lac"))
# BM2 <- subset(BMAA_2_data, BMAA_2_data$group %in% c("BMAA_glu"))
# BMAA <- rbind(BMAA1, BM1, BM2)
# write_tsv(BMAA, "./BMAA.tsv")

BMAA <- read_tsv("~/Data/figure_8_Supplementary_figure_7/BMAA.tsv")
BMAA_2 <- BMAA %>% group_by(group, media, sugar, species, strain, condition, time) %>%
  dplyr::summarise(
    mean = mean(OD600, na.rm = TRUE),
    n = length(OD600),
    sd = sd(OD600, na.rm = TRUE),
    se = sd / sqrt(n),
  )

BMAA_2 %>% ggplot(., aes(x=time, y=mean, colour=strain)) +
  geom_line() + 
  facet_wrap(~sugar) + # this line is great if you want to separate the graphs by strain / sugar / etc 
  geom_point(size=0.5) +
  ggtitle("Sugar utilization of each strain") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width =0.01, alpha=0.5) +
  expand_limits(y=0) +     # Expand y range
  theme_bw() +
  scale_color_manual(values = c("CM02_02" = "darkorchid4", "CM02_05" = "mediumorchid3", "CM02_91A" = "goldenrod3", "CM03_47" = "mediumpurple1", "CM06_28" = "lightgoldenrod", "CM06_29" = "plum1")) +
  xlab("time(hrs)") +
  ylab("growth (OD600)") +
  scale_x_continuous(breaks = seq(0, max(BMAA_2$time), by = 4))

# YCFA_nobac, YCFA_bac, YCFA_L-lac, YCFA_D-lac, YCFA_glu, YCFA_glu_l-lac, YCFA_glu_d-lac
# YCFA <- subset(YCFA_2_data, YCFA_2_data$group %in% c("YCFA_nobac", "YCFA_d_lac", "YCFA_l_lac"))
# YC2 <- subset(BMAA_2_data, BMAA_2_data$group=="YCFA")
# YCFA <- rbind(YCFA, YC2)
# write_tsv(YCFA, "./YCFA.tsv")
YCFA <- read_tsv("~/Data/figure8_supplementary_figure7/YCFA.tsv")
YCFA_2 <- YCFA %>% group_by(group, media, sugar, species, strain, condition, time) %>%
  dplyr::summarise(
    mean = mean(OD600, na.rm = TRUE),
    n = length(OD600),
    sd = sd(OD600, na.rm = TRUE),
    se = sd / sqrt(n),
  )

YCFA_2 %>% ggplot(., aes(x=time, y=mean, colour=strain)) +
  geom_line() + 
  facet_wrap(~sugar) + # this line is great if you want to separate the graphs by strain / sugar / etc 
  geom_point(size=0.5) +
  ggtitle("Sugar utilization of each strain") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width =0.01, alpha=0.5) +
  expand_limits(y=0) +     # Expand y range
  theme_bw() +
  scale_color_manual(values = c("CM02_02" = "darkorchid4", "CM02_05" = "mediumorchid3", "CM02_91A" = "goldenrod3", "CM03_47" = "mediumpurple1", "CM06_28" = "lightgoldenrod", "CM06_29" = "plum1")) +
  xlab("time(hrs)") +
  ylab("growth (OD600)") +
  scale_x_continuous(breaks = seq(0, max(YCFA_2$time), by = 4))

## Supplementary figure 7

CM3_3 <- CM3_2[,-12]
YCFA_6 <- subset(YCFA_2, YCFA_2$group %in% c("YCFA", "YCFA_d_lac", "YCFA_l_lac"))
YCFA_6$sugar[YCFA_6$sugar=="d_lactate"] <- "YCFA_d_lactate"
YCFA_6$sugar[YCFA_6$sugar=="l_lactate"] <- "YCFA_l_lactate"
YCFA_6$sugar[YCFA_6$sugar=="NA"] <- "YCFA"
BMAA_2$sugar[BMAA_2$sugar=="NA"] <- "BMAA"
BMAA_2$sugar[BMAA_2$sugar=="d_lactate"] <- "BMAA_d_lactate"
BMAA_2$sugar[BMAA_2$sugar=="l_lactate"] <- "BMAA_l_lactate"
BMAA_2$sugar[BMAA_2$sugar=="glucose"] <- "BMAA_glucose"
BMAA_2$sugar[BMAA_2$sugar=="glu_d_lactate"] <- "BMAA_glucose_d_lactate"
BMAA_2$sugar[BMAA_2$sugar=="glu_l_lactate"] <- "BMAA_glucose_l_lactate"
INVIT <- rbind(CM3_3, YCFA_6, BMAA_2)
INVIT$sugar <- factor(INVIT$sugar, levels = c("chondroitin_sodium_sulfate", "d_trehalose", "sucrose", "GOS_betaglucan", "arabinogalactan", "chitosan", "lentinan", "xyloglucan", "YCFA", "YCFA_d_lactate", "YCFA_l_lactate", "BMAA", "BMAA_l_lactate", "BMAA_d_lactate", "BMAA_glucose", "BMAA_glucose_d_lactate", "BMAA_glucose_l_lactate"))

## supplementary figure 7

INVIT %>% ggplot(., aes(x=time, y=mean, colour=strain)) +
  geom_line() + 
  facet_wrap(~sugar) + # this line is great if you want to separate the graphs by strain / sugar / etc 
  geom_point(size=0.5) +
  ggtitle("Sugar utilization of each strain") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width =0.01, alpha=0.5) +
  expand_limits(y=0) +     # Expand y range
  theme_bw() +
  scale_color_manual(values = c("CM02_02" = "darkorchid4", "CM02_05" = "mediumorchid3", "CM02_91A" = "goldenrod3", "CM03_47" = "mediumpurple1", "CM06_28" = "lightgoldenrod", "CM06_29" = "plum1")) +
  xlab("time(hrs)") +
  ylab("growth (OD600)") +
  scale_x_continuous(breaks = seq(0, max(INVIT$time), by = 4))


## Figure 8 with growthcurver

YCFA_3 <- subset(YCFA, YCFA$sugar %in% c("d_lactate", "l_lactate", "NA"))
YCFA_4 <- select(YCFA_3, c(5,8,1,15,19,20)) %>% 
  pivot_wider(id_cols= c("strain", "time"), names_from = "sugar", values_from = "OD600")
colnames(YCFA_4) <- c("strain", "time", "d_lactate", "l_lactate", "no_sugar")
YCFA_5 <- YCFA_4 %>%
  mutate(
    d_lactate = str_remove_all(d_lactate, "c\\(|\\)"),  # remove c(...)
    d_lactate = str_remove_all(d_lactate, '"'),         # remove quotes
    l_lactate = str_remove_all(l_lactate, "c\\(|\\)"),
    l_lactate = str_remove_all(l_lactate, '"'),
    no_sugar = str_remove_all(no_sugar, "c\\(|\\)"),
    no_sugar = str_remove_all(no_sugar, '"')
  ) %>%
  separate(d_lactate, into = c("d_lac1","d_lac2","d_lac3"), sep = ",", convert = TRUE) %>%
  separate(l_lactate, into = c("l_lac1","l_lac2","l_lac3"), sep = ",", convert = TRUE) %>% 
  separate(no_sugar, into = c("ns1","ns2","ns3"), sep = ",", convert = TRUE)
YCFA_5[,c(3,6,9)] <- as.numeric(YCFA_5[,c(3,6,9)])
YCFA_5$d_lac1 <- as.numeric(YCFA_5$d_lac1)
YCFA_5$l_lac1 <- as.numeric(YCFA_5$l_lac1)
YCFA_5$ns1 <- as.numeric(YCFA_5$ns1)

results_list <- list()

for (strain_id in strain_list) {
  
  temp_df <- YCFA_5 %>%
    filter(strain == strain_id) %>%
    select(-1)  # assuming first column is time or index you want removed
  
  if (nrow(temp_df) > 0) {  # only run if subset has data
    gcvr <- SummarizeGrowthByPlate(temp_df)
    results_list[[strain_id]] <- gcvr
  } else {
    warning(paste("No data for strain:", strain_id))
    results_list[[strain_id]] <- NULL
  }
}
final_results_lac_ycfa <- dplyr::bind_rows(results_list, .id = "strain")

final_results_lac_ycfa$sugar[final_results_lac_ycfa$sample %in% c("d_lac1", "d_lac2", "d_lac3")] <- "d_lactate"
final_results_lac_ycfa$sugar[final_results_lac_ycfa$sample %in% c("l_lac1", "l_lac2", "l_lac3")] <- "l_lactate"

## stats for YCFA
d_lactate <- subset(final_results_lac_ycfa, final_results_lac_ycfa$sugar=="d_lactate")
aov_dlac <- aov(auc_l~strain, d_lactate)
summary(aov_dlac)
# Df Sum Sq Mean Sq F value  Pr(>F)    
# strain       5  435.6   87.12     108 1.5e-09 ***
#   Residuals   12    9.7    0.81                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_dlac, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = d_lactate)
# 
# $strain
#                     diff         lwr        upr     p adj
# CM02_05-CM02_02   14.5636435  12.1004839 17.0268030 0.0000000***
# CM02_91A-CM02_02   2.4903390   0.0271795  4.9534986 0.0469884*
# CM03_47-CM02_02    6.8264576   4.3632980  9.2896171 0.0000090***
# CM06_28-CM02_02   11.2975168   8.8343573 13.7606764 0.0000000***
# CM06_29-CM02_02    6.5841592   4.1209997  9.0473188 0.0000132***
# CM02_91A-CM02_05 -12.0733044 -14.5364640 -9.6101449 0.0000000***
# CM03_47-CM02_05   -7.7371859 -10.2003454 -5.2740264 0.0000024***
# CM06_28-CM02_05   -3.2661266  -5.7292862 -0.8029671 0.0079208**
# CM06_29-CM02_05   -7.9794842 -10.4426438 -5.5163247 0.0000017***
# CM03_47-CM02_91A   4.3361185   1.8729590  6.7992781 0.0007768***
# CM06_28-CM02_91A   8.8071778   6.3440183 11.2703373 0.0000006***
# CM06_29-CM02_91A   4.0938202   1.6306607  6.5569797 0.0012873**
# CM06_28-CM03_47    4.4710593   2.0078998  6.9342188 0.0005897***
# CM06_29-CM03_47   -0.2422983  -2.7054578  2.2208612 0.9993282
# CM06_29-CM06_28   -4.7133576  -7.1765171 -2.2501981 0.0003632***
# CM02_02, CM02_05, CM02_91A, CM06_28 are ***

l_lactate <- subset(final_results_lac_ycfa, final_results_lac_ycfa$sugar=="l_lactate")
aov_llac <- aov(auc_l~strain, l_lactate)
summary(aov_llac)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# strain       5  322.4   64.47    88.8 4.71e-09 ***
#   Residuals   12    8.7    0.73                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_llac, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = l_lactate)
# 
# $strain
# diff         lwr        upr     p adj
# CM02_05-CM02_02  10.5481241   8.2112958 12.8849524 0.0000000***
# CM02_91A-CM02_02  2.1278368  -0.2089915  4.4646651 0.0822955
# CM03_47-CM02_02   9.7393896   7.4025613 12.0762179 0.0000001***
# CM06_28-CM02_02   8.3680041   6.0311758 10.7048324 0.0000006***
# CM06_29-CM02_02   1.7332399  -0.6035884  4.0700682 0.2010577
# CM02_91A-CM02_05 -8.4202873 -10.7571156 -6.0834590 0.0000005***
# CM03_47-CM02_05  -0.8087345  -3.1455628  1.5280938 0.8458620
# CM06_28-CM02_05  -2.1801200  -4.5169483  0.1567083 0.0727275
# CM06_29-CM02_05  -8.8148842 -11.1517125 -6.4780559 0.0000003***
# CM03_47-CM02_91A  7.6115528   5.2747245  9.9483811 0.0000016***
# CM06_28-CM02_91A  6.2401673   3.9033390  8.5769956 0.0000133***
# CM06_29-CM02_91A -0.3945969  -2.7314252  1.9422314 0.9914700
# CM06_28-CM03_47  -1.3713855  -3.7082138  0.9654428 0.4097750
# CM06_29-CM03_47  -8.0061497 -10.3429780 -5.6693214 0.0000009***
# CM06_29-CM06_28  -6.6347642  -8.9715925 -4.2979359 0.0000070***
# CM02_02, CM02_05, CM02_91A, CM06_29 is ***


BMAA2 <- select(BMAA, c(5,8,1,15,19,20)) %>% 
  pivot_wider(id_cols= c("strain", "time"), names_from = "sugar", values_from = "OD600")
colnames(BMAA2) <- c("strain", "time", "l_lactate", "no_sugar", "d_lactate", "glu_d_lactate", "glu_l_lactate", "glucose")
BMAA3 <- BMAA2 %>%
  mutate(
    d_lactate = str_remove_all(d_lactate, "c\\(|\\)"),  # remove c(...)
    d_lactate = str_remove_all(d_lactate, '"'),         # remove quotes
    l_lactate = str_remove_all(l_lactate, "c\\(|\\)"),
    l_lactate = str_remove_all(l_lactate, '"'),
    no_sugar = str_remove_all(no_sugar, "c\\(|\\)"),
    no_sugar = str_remove_all(no_sugar, '"'),
    glu_l_lactate = str_remove_all(glu_l_lactate, "c\\(|\\)"),
    glu_l_lactate = str_remove_all(glu_l_lactate, '"'),
    glu_d_lactate = str_remove_all(glu_d_lactate, "c\\(|\\)"),
    glu_d_lactate = str_remove_all(glu_d_lactate, '"'),
    glucose = str_remove_all(glucose, "c\\(|\\)"),
    glucose = str_remove_all(glucose, '"')
  ) %>%
  separate(d_lactate, into = c("d_lac1","d_lac2","d_lac3"), sep = ",", convert = TRUE) %>%
  separate(l_lactate, into = c("l_lac1","l_lac2","l_lac3"), sep = ",", convert = TRUE) %>% 
  separate(no_sugar, into = c("no_sugar1","no_sugar2","no_sugar3"), sep = ",", convert = TRUE) %>%
  separate(glu_l_lactate, into = c("glu_l_lac1","glu_l_lac2","glu_l_lac3"), sep = ",", convert = TRUE) %>%
  separate(glu_d_lactate, into = c("glu_d_lac1","glu_d_lac2","glu_d_lac3"), sep = ",", convert = TRUE) %>%
  separate(glucose, into = c("glu_1","glu_2","glu_3"), sep = ",", convert = TRUE)

BMAA3$glu_1 <- as.numeric(BMAA3$glu_1)

results_list_bmaa <- list()

for (strain_id in strain_list) {
  
  temp_df <- BMAA3 %>%
    filter(strain == strain_id) %>%
    select(-1)  # assuming first column is time or index you want removed
  
  if (nrow(temp_df) > 0) {  # only run if subset has data
    gcvr <- SummarizeGrowthByPlate(temp_df)
    results_list_bmaa[[strain_id]] <- gcvr
  } else {
    warning(paste("No data for strain:", strain_id))
    results_list_bmaa[[strain_id]] <- NULL
  }
}
final_results_lac_bmaa <- dplyr::bind_rows(results_list_bmaa, .id = "strain")

final_results_lac_bmaa$sugar[final_results_lac_bmaa$sample %in% c("d_lac1", "d_lac2", "d_lac3")] <- "d_lactate"
final_results_lac_bmaa$sugar[final_results_lac_bmaa$sample %in% c("l_lac1", "l_lac2", "l_lac3")] <- "l_lactate"
final_results_lac_bmaa$sugar[final_results_lac_bmaa$sample %in% c("no_sugar1", "no_sugar2", "no_sugar3")] <- "no_sugar"
final_results_lac_bmaa$sugar[final_results_lac_bmaa$sample %in% c("glu_l_lac1", "glu_l_lac2", "glu_l_lac3")] <- "glucose_l_lactate"
final_results_lac_bmaa$sugar[final_results_lac_bmaa$sample %in% c("glu_d_lac1", "glu_d_lac2", "glu_d_lac3")] <- "glucose_d_lactate"
final_results_lac_bmaa$sugar[final_results_lac_bmaa$sample %in% c("glu_1", "glu_2", "glu_3")] <- "glucose"

## stats for BMAA
d_lactate <- subset(final_results_lac_bmaa, final_results_lac_bmaa$sugar=="d_lactate")
aov_dlac <- aov(auc_l~strain, d_lactate)
summary(aov_dlac)
# Df Sum Sq Mean Sq F value Pr(>F)   
# strain       5 0.4903 0.09806   7.679 0.0019 **
#   Residuals   12 0.1532 0.01277                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_dlac, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = d_lactate)
# 
# $strain
# diff         lwr         upr     p adj
# CM02_05-CM02_02  -0.193784922 -0.50370672  0.11613688 0.3481122
# CM02_91A-CM02_02 -0.181315190 -0.49123699  0.12860661 0.4128490
# CM03_47-CM02_02   0.263040770 -0.04688103  0.57296257 0.1152105
# CM06_28-CM02_02  -0.123598582 -0.43352038  0.18632322 0.7592975
# CM06_29-CM02_02  -0.200691798 -0.51061360  0.10923000 0.3152499
# CM02_91A-CM02_05  0.012469732 -0.29745207  0.32239153 0.9999917
# CM03_47-CM02_05   0.456825692  0.14690389  0.76674749 0.0035014**
# CM06_28-CM02_05   0.070186340 -0.23973546  0.38010814 0.9692496
# CM06_29-CM02_05  -0.006906876 -0.31682867  0.30301492 0.9999996
# CM03_47-CM02_91A  0.444355960  0.13443416  0.75427776 0.0043616**
# CM06_28-CM02_91A  0.057716609 -0.25220519  0.36763841 0.9867755
# CM06_29-CM02_91A -0.019376608 -0.32929841  0.29054519 0.9999267
# CM06_28-CM03_47  -0.386639352 -0.69656115 -0.07671755 0.0123054*
# CM06_29-CM03_47  -0.463732568 -0.77365437 -0.15381077 0.0031027**
# CM06_29-CM06_28  -0.077093217 -0.38701502  0.23282858 0.9547593
# C03_47is **

l_lactate <- subset(final_results_lac_bmaa, final_results_lac_bmaa$sugar=="l_lactate")
aov_llac <- aov(auc_l~strain, l_lactate)
summary(aov_llac)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# strain       5  7.549  1.5097   6.962 0.00287 **
#   Residuals   12  2.602  0.2169                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_llac, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = l_lactate)
# 
# $strain
# diff        lwr        upr     p adj
# CM02_05-CM02_02  -0.128892338 -1.4060308  1.1482461 0.9992398
# CM02_91A-CM02_02 -0.003552413 -1.2806909  1.2735861 1.0000000
# CM03_47-CM02_02   1.591545939  0.3144075  2.8686844 0.0124000*
# CM06_28-CM02_02  -0.209608764 -1.4867473  1.0675297 0.9925029
# CM06_29-CM02_02  -0.280977532 -1.5581160  0.9961610 0.9727733
# CM02_91A-CM02_05  0.125339925 -1.1517986  1.4024784 0.9993357
# CM03_47-CM02_05   1.720438277  0.4432998  2.9975768 0.0070410**
# CM06_28-CM02_05  -0.080716426 -1.3578549  1.1964221 0.9999227
# CM06_29-CM02_05  -0.152085194 -1.4292237  1.1250533 0.9983202
# CM03_47-CM02_91A  1.595098351  0.3179599  2.8722368 0.0122070*
# CM06_28-CM02_91A -0.206056351 -1.4831948  1.0710821 0.9930647
# CM06_29-CM02_91A -0.277425120 -1.5545636  0.9997134 0.9742017
# CM06_28-CM03_47  -1.801154703 -3.0782932 -0.5240162 0.0049617**
# CM06_29-CM03_47  -1.872523471 -3.1496620 -0.5953850 0.0036535**
# CM06_29-CM06_28  -0.071368768 -1.3485073  1.2057697 0.9999579
#CM 03_47 is **
no_sugar <- subset(final_results_lac_bmaa, final_results_lac_bmaa$sugar=="no_sugar")
aov_nos <- aov(auc_l~strain, no_sugar)
summary(aov_nos)
# Df Sum Sq Mean Sq F value Pr(>F)  
# strain       5  1.664  0.3328    2.67 0.0759 .
# Residuals   12  1.496  0.1246                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
glu_l_lac <- subset(final_results_lac_bmaa, final_results_lac_bmaa$sugar=="glucose_l_lactate")
aov_glullac <- aov(auc_l~strain, glu_l_lac)
summary(aov_glullac)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# strain       5  6.227  1.2453   9.836 0.000634 ***
#   Residuals   12  1.519  0.1266                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_glullac, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = glu_l_lac)
# 
# $strain
# diff        lwr        upr     p adj
# CM02_05-CM02_02   0.232189619 -0.7436832  1.2080625 0.9622576
# CM02_91A-CM02_02 -0.089700677 -1.0655735  0.8861722 0.9995163
# CM03_47-CM02_02   1.625503646  0.6496308  2.6013765 0.0012629**
# CM06_28-CM02_02   0.240917056 -0.7349558  1.2167899 0.9561243
# CM06_29-CM02_02   0.021625643 -0.9542472  0.9974985 0.9999996
# CM02_91A-CM02_05 -0.321890296 -1.2977632  0.6539826 0.8690724
# CM03_47-CM02_05   1.393314027  0.4174412  2.3691869 0.0045076**
# CM06_28-CM02_05   0.008727437 -0.9671454  0.9846003 1.0000000
# CM06_29-CM02_05  -0.210563976 -1.1864368  0.7653089 0.9749274
# CM03_47-CM02_91A  1.715204323  0.7393315  2.6910772 0.0007878***
# CM06_28-CM02_91A  0.330617733 -0.6452551  1.3064906 0.8565134
# CM06_29-CM02_91A  0.111326320 -0.8645465  1.0871992 0.9986310
# CM06_28-CM03_47  -1.384586590 -2.3604594 -0.4087137 0.0047346**
# CM06_29-CM03_47  -1.603878003 -2.5797509 -0.6280051 0.0014175**
# CM06_29-CM06_28  -0.219291413 -1.1951643  0.7565814 0.9702341
# CM03_47 is ***
glu_d_lac <- subset(final_results_lac_bmaa, final_results_lac_bmaa$sugar=="glucose_d_lactate")
aov_gludlac <- aov(auc_l~strain, glu_d_lac)
summary(aov_gludlac)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# strain       5  9.190  1.8381   77.09 1.07e-08 ***
#   Residuals   12  0.286  0.0238                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_gludlac, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = glu_d_lac)
# 
# $strain
# diff        lwr        upr     p adj
# CM02_05-CM02_02  -0.29776421 -0.7212564  0.1257280 0.2430871
# CM02_91A-CM02_02 -0.19614769 -0.6196399  0.2273445 0.6388805
# CM03_47-CM02_02   1.68345535  1.2599631  2.1069476 0.0000002***
# CM06_28-CM02_02  -0.26533963 -0.6888318  0.1581526 0.3461637
# CM06_29-CM02_02  -0.30834833 -0.7318405  0.1151439 0.2151336
# CM02_91A-CM02_05  0.10161652 -0.3218757  0.5251087 0.9609261
# CM03_47-CM02_05   1.98121956  1.5577273  2.4047118 0.0000000***
# CM06_28-CM02_05   0.03242458 -0.3910676  0.4559168 0.9998016
# CM06_29-CM02_05  -0.01058412 -0.4340763  0.4129081 0.9999992
# CM03_47-CM02_91A  1.87960304  1.4561108  2.3030953 0.0000000***
# CM06_28-CM02_91A -0.06919194 -0.4926842  0.3543003 0.9926555
# CM06_29-CM02_91A -0.11220064 -0.5356929  0.3112916 0.9418407
# CM06_28-CM03_47  -1.94879498 -2.3722872 -1.5253028 0.0000000***
# CM06_29-CM03_47  -1.99180368 -2.4152959 -1.5683115 0.0000000***
# CM06_29-CM06_28  -0.04300870 -0.4665009  0.3804835 0.9992164
# CM03_47 is ***
glu <- subset(final_results_lac_bmaa, final_results_lac_bmaa$sugar=="glucose")
aov_glu <- aov(auc_l~strain, glu)
summary(aov_glu)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# strain       5  350.0   69.99   163.7 1.31e-10 ***
#   Residuals   12    5.1    0.43                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_glu, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = glu)
# 
# $strain
# diff         lwr         upr     p adj
# CM02_05-CM02_02   -6.891893  -8.6849512  -5.0988353 0.0000003***
# CM02_91A-CM02_02 -12.009318 -13.8023762 -10.2162603 0.0000000***
# CM03_47-CM02_02  -11.031253 -12.8243111  -9.2381952 0.0000000***
# CM06_28-CM02_02   -2.722964  -4.5160219  -0.9299061 0.0027505**
# CM06_29-CM02_02   -9.835764 -11.6288221  -8.0427062 0.0000000***
# CM02_91A-CM02_05  -5.117425  -6.9104829  -3.3243670 0.0000066***
# CM03_47-CM02_05   -4.139360  -5.9324178  -2.3463020 0.0000592***
# CM06_28-CM02_05    4.168929   2.3758713   5.9619872 0.0000551***
# CM06_29-CM02_05   -2.943871  -4.7369288  -1.1508129 0.0014303**
# CM03_47-CM02_91A   0.978065  -0.8149929   2.7711230 0.4824940
# CM06_28-CM02_91A   9.286354   7.4932963  11.0794122 0.0000000***
# CM06_29-CM02_91A   2.173554   0.3804961   3.9666120 0.0150249*
# CM06_28-CM03_47    8.308289   6.5152312  10.1013471 0.0000000***
# CM06_29-CM03_47    1.195489  -0.5975689   2.9885470 0.2887621
# CM06_29-CM06_28   -7.112800  -8.9058581  -5.3197422 0.0000002***
# CM02_2, 02_5, 06_28, 06_29 is ***


## Negatives for BMAA and YCFA; this section determines the upper and lower limits for no, strong and weak growth
BMYC_neg <- subset(BMAA, BMAA$sugar=="NA")
YC_neg <- subset(YCFA, YCFA$sugar=="NA")
BMYC_neg <- rbind(BMYC_neg, YC_neg)
BMYC_neg <- subset(BMYC_neg, BMYC_neg$strain %in% c("CM02_02","CM02_05","CM06_29","CM02_91A","CM06_28","CM03_47"))
BMYC_neg2 <- select(BMYC_neg, c(5,4, 1, 15,19,20)) %>% 
  mutate(wellID_strain = paste0(media, "_", sugar)) %>% 
  #select(-1,-2) %>% 
  pivot_wider(id_cols= c("strain", "time"), names_from = "wellID_strain", values_from = "OD600")

BMYC_neg3 <- BMYC_neg2 %>%
  mutate(
    BMAA_NA = str_remove_all(BMAA_NA, "c\\(|\\)"),  # remove c(...)
    BMAA_NA = str_remove_all(BMAA_NA, '"'),         # remove quotes
    YCFA_NA = str_remove_all(YCFA_NA, "c\\(|\\)"),
    YCFA_NA = str_remove_all(YCFA_NA, '"'),
  ) %>%
  separate(BMAA_NA, into = c("bmaa1","bmaa2","bmaa3"), sep = ",", convert = TRUE) %>%
  separate(YCFA_NA, into = c("ycfa1","ycfa2","ycfa3"), sep = ",", convert = TRUE)


results_list_neg_8b <- list()
for (strain_id in strain_list) {
  
  temp_df_neg <- BMYC_neg3 %>%
    filter(strain == strain_id) %>%
    select(-1)  # assuming first column is time or index you want removed
  
  if (nrow(temp_df_neg) > 0) {  # only run if subset has data
    gcvr_neg_8b <- SummarizeGrowthByPlate(temp_df_neg)
    results_list_neg_8b[[strain_id]] <- gcvr_neg_8b
  } else {
    warning(paste("No data for strain:", strain_id))
    results_list_neg_8b[[strain_id]] <- NULL
  }
}
final_results_neg_8b <- dplyr::bind_rows(results_list_neg_8b, .id = "strain")
neg_mean_bmaa_8b <- mean(final_results_neg_8b$auc_l[final_results_neg_8b$sample %in% c("bmaa1","bmaa2","bmaa3")], na.rm = TRUE)
neg_sd_bmaa_8b <- sd(final_results_neg_8b$auc_l[final_results_neg_8b$sample %in% c("bmaa1","bmaa2","bmaa3")], na.rm = TRUE)
cutoff_bmaa_8b <- neg_mean_bmaa_8b + 3 * neg_sd_bmaa_8b
lower_cutoff_bmaa_8b <- neg_mean_bmaa_8b+neg_sd_bmaa_8b

neg_mean_ycfa_8b <- mean(final_results_neg_8b$auc_l[final_results_neg_8b$sample %in% c("ycfa1","ycfa2","ycfa3")], na.rm = TRUE)
neg_sd_ycfa_8b <- sd(final_results_neg_8b$auc_l[final_results_neg_8b$sample %in% c("ycfa1","ycfa2","ycfa3")], na.rm = TRUE)
upper_cutoff_ycfa_8b <- neg_mean_ycfa_8b + 3 * neg_sd_ycfa_8b
lower_cutoff_ycfa_8b <- neg_mean_ycfa_8b+neg_sd_ycfa_8b

## figure 8 A

# data_list_PBM1_fixed <- lapply(data_list_PBM1, function(df) {
#   df$date <- as.character(df$date)
#   df$strain_origin_date <- as.character(df$strain_origin_date)
#   return(df)
# })
# PBM1 <- bind_rows(data_list_PBM1_fixed)
# 
# data_list_PBM2_fixed <- lapply(data_list_PBM2, function(df) {
#   df$date <- as.character(df$date)
#   df$strain_origin_date <- as.character(df$strain_origin_date)
#   return(df)
# })
# PBM2 <- bind_rows(data_list_PBM2_fixed)
# 
# PBM <- rbind(PBM1, PBM2)
PBM <- read_tsv("~/Data/figure_8_Supplementary_figure_7/PBM.tsv")
PBM_sugar <- subset(PBM, PBM$sugar %in% c("chondroitin_sodium_sulfate", "sucrose", "d_trehalose", "arabinogalactan", "GOS_betaglucan", "lentinan", "chitosan", "xyloglucan"))

PBM_sugar2 <- select(PBM_sugar, c(5,1,15,19,20)) %>% 
  mutate(wellID_strain = paste0(wellID, "_", sugar)) #%>% 
#select(-1) %>% 
#pivot_wider(id_cols= c("strain", "time"), names_from = "wellID_strain", values_from = "OD600")

PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="E03_d_trehalose"] <- 1
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="E07_d_trehalose"] <- 2
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="E11_d_trehalose"] <- 3
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="G01_chondroitin_sodium_sulfate"] <- 1
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="G05_chondroitin_sodium_sulfate"] <- 2
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="G09_chondroitin_sodium_sulfate"] <- 3
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="H01_sucrose"] <- 1
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="H05_sucrose"] <- 2
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="H09_sucrose"] <- 3
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="D01_xyloglucan"] <- 1
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="D05_xyloglucan"] <- 2
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="D09_xyloglucan"] <- 3
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="D02_GOS_betaglucan"] <- 1
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="D06_GOS_betaglucan"] <- 2
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="D10_GOS_betaglucan"] <- 3
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="G03_chitosan"] <- 1
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="G07_chitosan"] <- 2
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="G11_chitosan"] <- 3
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="H02_arabinogalactan"] <- 1
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="H06_arabinogalactan"] <- 2
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="H10_arabinogalactan"] <- 3
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="H03_lentinan"] <- 1
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="H07_lentinan"] <- 2
PBM_sugar2$replicate[PBM_sugar2$wellID_strain=="H11_lentinan"] <- 3


## summarize by plate; this works

PBM_sugar3 <- select(PBM_sugar, c(5,1,15,19,20)) %>% 
  mutate(wellID_strain = paste0(wellID, "_", sugar)) %>% 
  select(-1) %>% 
  pivot_wider(id_cols= c("strain", "time"), names_from = "wellID_strain", values_from = "OD600")

## check
PBM_02_2 <- subset(PBM_sugar3, PBM_sugar3$strain=="CM02_02") %>% 
  select(-1)

gcvr_02_2 <- SummarizeGrowthByPlate(PBM_02_2)


strain_list <- unique(PBM_sugar3$strain)

# Create a list to store results
results_list <- list()

for (strain_id in strain_list) {
  
  temp_df <- PBM_sugar3 %>%
    filter(strain == strain_id) %>%
    select(-1)  # assuming first column is time or index you want removed
  
  if (nrow(temp_df) > 0) {  # only run if subset has data
    gcvr <- SummarizeGrowthByPlate(temp_df)
    results_list[[strain_id]] <- gcvr
  } else {
    warning(paste("No data for strain:", strain_id))
    results_list[[strain_id]] <- NULL
  }
}
final_results <- dplyr::bind_rows(results_list, .id = "strain")
PBM_sugar_small <- PBM_sugar2[,c(6,1)] %>% distinct()
gcvr_results <- inner_join(final_results, PBM_sugar_small, by = c("sample"="wellID_strain"))

## stats for PBM1 and PBM2 sugars
arabinogalactan <- subset(gcvr_results, gcvr_results$sugar=="arabinogalactan")
aov_ara <- aov(auc_l~strain, arabinogalactan)
summary(aov_ara)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# strain       5  6.003   1.201      15 8.36e-05 ***
#   Residuals   12  0.960   0.080                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_ara, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = arabinogalactan)
# 
# $strain
#                       diff         lwr          upr     p adj
# CM02_05-CM02_02  -0.58323581 -1.35912033  0.192648708 0.1912017
# CM02_91A-CM02_02 -1.35670363 -2.13258815 -0.580819118 0.0008247***
# CM03_47-CM02_02   0.16378217 -0.61210234  0.939666690 0.9771669
# CM06_28-CM02_02  -0.75777560 -1.53366011  0.018108919 0.0569994
# CM06_29-CM02_02  -1.28387820 -2.05976272 -0.507993687 0.0013371**
# CM02_91A-CM02_05 -0.77346783 -1.54935234  0.002416692 0.0508835
# CM03_47-CM02_05   0.74701798 -0.02886653  1.522902499 0.0615956
# CM06_28-CM02_05  -0.17453979 -0.95042431  0.601344729 0.9700989
# CM06_29-CM02_05  -0.70064239 -1.47652691  0.075242123 0.0857797
# CM03_47-CM02_91A  1.52048581  0.74460129  2.296370325 0.0002907***
# CM06_28-CM02_91A  0.59892804 -0.17695648  1.374812554 0.1725041
# CM06_29-CM02_91A  0.07282543 -0.70305909  0.848709948 0.9994647
# CM06_28-CM03_47  -0.92155777 -1.69744229 -0.145673254 0.0172571*
# CM06_29-CM03_47  -1.44766038 -2.22354489 -0.671775860 0.0004587***
# CM06_29-CM06_28  -0.52610261 -1.30198712  0.249781911 0.2738887

chitosan <- subset(gcvr_results, gcvr_results$sugar=="chitosan")
aov_chi <- aov(auc_l~strain, chitosan)
summary(aov_chi)
# Df Sum Sq Mean Sq F value Pr(>F)  
# strain       5  15.08   3.017   2.741 0.0707 .
# Residuals   12  13.21   1.100                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1    
chondroitin_sodium_sulfate <- subset(gcvr_results, gcvr_results$sugar=="chondroitin_sodium_sulfate")
aov_cho <- aov(auc_l~strain, chondroitin_sodium_sulfate)
summary(aov_cho)
# Df Sum Sq Mean Sq F value Pr(>F)
# strain       5  4.225  0.8450   1.591  0.236
# Residuals   12  6.374  0.5312      
d_trehalose <- subset(gcvr_results, gcvr_results$sugar=="d_trehalose")
aov_tre <- aov(auc_l~strain, d_trehalose)
summary(aov_tre)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# strain       5  296.7   59.34   79.14 9.18e-09 ***
#   Residuals   12    9.0    0.75                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_tre, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = d_trehalose)
# 
# $strain
#                     diff        lwr        upr     p adj
# CM02_05-CM02_02   -9.9159303 -12.290737 -7.5411237 0.0000001***
# CM02_91A-CM02_02  -2.6736280  -5.048435 -0.2988214 0.0245240*
# CM03_47-CM02_02  -10.1907865 -12.565593 -7.8159799 0.0000001***
# CM06_28-CM02_02   -3.6809607  -6.055767 -1.3061540 0.0023241**
# CM06_29-CM02_02   -9.6822558 -12.057062 -7.3074492 0.0000001***
# CM02_91A-CM02_05   7.2423023   4.867496  9.6171089 0.0000033***
# CM03_47-CM02_05   -0.2748562  -2.649663  2.0999504 0.9985333
# CM06_28-CM02_05    6.2349697   3.860163  8.6097763 0.0000159***
# CM06_29-CM02_05    0.2336745  -2.141132  2.6084811 0.9993273
# CM03_47-CM02_91A  -7.5171585  -9.891965 -5.1423519 0.0000022***
# CM06_28-CM02_91A  -1.0073326  -3.382139  1.3674740 0.7130687
# CM06_29-CM02_91A  -7.0086278  -9.383434 -4.6338212 0.0000046***
# CM06_28-CM03_47    6.5098259   4.135019  8.8846325 0.0000101***
# CM06_29-CM03_47    0.5085307  -1.866276  2.8833373 0.9757261
# CM06_29-CM06_28   -6.0012952  -8.376102 -3.6264885 0.0000236***
GOS <- subset(gcvr_results, gcvr_results$sugar=="GOS_betaglucan")
aov_GOS <- aov(auc_l~strain, GOS)
summary(aov_GOS)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# strain       5  31.91   6.382   25.66 5.12e-06 ***
#   Residuals   12   2.98   0.249                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_GOS, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = GOS)
# 
# $strain
# diff       lwr        upr     p adj
# CM02_05-CM02_02  -0.6884816 -2.056184  0.6792205 0.5613699
# CM02_91A-CM02_02 -0.8255395 -2.193242  0.5421626 0.3821532
# CM03_47-CM02_02   2.4720881  1.104386  3.8397902 0.0006128***
# CM06_28-CM02_02  -1.2020361 -2.569738  0.1656659 0.0978868
# CM06_29-CM02_02  -1.5731315 -2.940834 -0.2054294 0.0213546*
# CM02_91A-CM02_05 -0.1370579 -1.504760  1.2306442 0.9992653
# CM03_47-CM02_05   3.1605697  1.792868  4.5282718 0.0000586***
# CM06_28-CM02_05  -0.5135545 -1.881257  0.8541475 0.7994556
# CM06_29-CM02_05  -0.8846499 -2.252352  0.4830522 0.3163075
# CM03_47-CM02_91A  3.2976276  1.929926  4.6653297 0.0000382***
# CM06_28-CM02_91A -0.3764967 -1.744199  0.9912054 0.9325026
# CM06_29-CM02_91A -0.7475920 -2.115294  0.6201100 0.4804324
# CM06_28-CM03_47  -3.6741242 -5.041826 -2.3064222 0.0000125***
# CM06_29-CM03_47  -4.0452196 -5.412922 -2.6775175 0.0000045***
# CM06_29-CM06_28  -0.3710954 -1.738797  0.9966067 0.9361721
lentinan <- subset(gcvr_results, gcvr_results$sugar=="lentinan")
aov_len <- aov(auc_l~strain, lentinan)
summary(aov_len)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# strain       5  6.991  1.3982   7.728 0.00185 **
# Residuals   12  2.171  0.1809                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_len, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = lentinan)
# 
# $strain
# diff        lwr         upr     p adj
# CM02_05-CM02_02  -0.8131234 -1.9796381  0.35339128 0.2502833
# CM02_91A-CM02_02 -1.4719197 -2.6384343 -0.30540498 0.0113542*
# CM03_47-CM02_02   0.1036258 -1.0628888  1.27014050 0.9995902
# CM06_28-CM02_02  -1.2245150 -2.3910297 -0.05800038 0.0377652*
# CM06_29-CM02_02  -1.2935527 -2.4600674 -0.12703807 0.0269980*
# CM02_91A-CM02_05 -0.6587963 -1.8253109  0.50771840 0.4478867
# CM03_47-CM02_05   0.9167492 -0.2497654  2.08326389 0.1605325
# CM06_28-CM02_05  -0.4113917 -1.5779063  0.75512301 0.8359345
# CM06_29-CM02_05  -0.4804294 -1.6469440  0.68608531 0.7358041
# CM03_47-CM02_91A  1.5755455  0.4090308  2.74206015 0.0069036**
# CM06_28-CM02_91A  0.2474046 -0.9191101  1.41391927 0.9767015
# CM06_29-CM02_91A  0.1783669 -0.9881478  1.34488158 0.9945785
# CM06_28-CM03_47  -1.3281409 -2.4946555 -0.16162621 0.0228143*
# CM06_29-CM03_47  -1.3971786 -2.5636932 -0.23066391 0.0163076*
# CM06_29-CM06_28  -0.0690377 -1.2355524  1.09747697 0.9999441
sucrose <- subset(gcvr_results, gcvr_results$sugar=="sucrose")
aov_suc <- aov(auc_l~strain, sucrose)
summary(aov_suc)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# strain       5  832.6  166.52   57.91 5.53e-08 ***
#   Residuals   12   34.5    2.88                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_suc, "strain")
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = sucrose)
# 
# $strain
# diff         lwr        upr     p adj
# CM02_05-CM02_02    6.5381915   1.8876343  11.188749 0.0050837**
# CM02_91A-CM02_02 -10.7598587 -15.4104159  -6.109301 0.0000579***
# CM03_47-CM02_02   -8.6905770 -13.3411343  -4.040020 0.0004522***
# CM06_28-CM02_02    4.4997816  -0.1507756   9.150339 0.0599709
# CM06_29-CM02_02   -8.3804562 -13.0310135  -3.729899 0.0006297***
# CM02_91A-CM02_05 -17.2980502 -21.9486075 -12.647493 0.0000004***
# CM03_47-CM02_05  -15.2287685 -19.8793258 -10.578211 0.0000015***
# CM06_28-CM02_05   -2.0384099  -6.6889672   2.612147 0.6864748
# CM06_29-CM02_05  -14.9186477 -19.5692050 -10.268090 0.0000019***
# CM03_47-CM02_91A   2.0692817  -2.5812756   6.719839 0.6738455
# CM06_28-CM02_91A  15.2596403  10.6090830  19.910198 0.0000015***
# CM06_29-CM02_91A   2.3794025  -2.2711548   7.029960 0.5456288
# CM06_28-CM03_47   13.1903586   8.5398014  17.840916 0.0000071***
# CM06_29-CM03_47    0.3101208  -4.3404364   4.960678 0.9998993
# CM06_29-CM06_28  -12.8802378 -17.5307951  -8.229681 0.0000091***
xyloglucan <- subset(gcvr_results, gcvr_results$sugar=="xyloglucan")
aov_xyl <- aov(auc_l~strain, xyloglucan)
summary(aov_xyl)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# strain       5  11.98  2.3964   23.57 8.08e-06 ***
#   Residuals   12   1.22  0.1017                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov_xyl)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = auc_l ~ strain, data = xyloglucan)
# 
# $strain
# diff        lwr        upr     p adj
# CM02_05-CM02_02  -0.38270251 -1.2572950  0.4918899 0.6878720
# CM02_91A-CM02_02 -0.05992418 -0.9345166  0.8146683 0.9998850
# CM03_47-CM02_02   1.74061000  0.8660176  2.6152025 0.0002513***
# CM06_28-CM02_02  -0.46524515 -1.3398376  0.4093473 0.5073962
# CM06_29-CM02_02  -0.77357852 -1.6481710  0.1010139 0.0949307
# CM02_91A-CM02_05  0.32277833 -0.5518141  1.1973708 0.8100354
# CM03_47-CM02_05   2.12331252  1.2487201  2.9979050 0.0000356***
# CM06_28-CM02_05  -0.08254264 -0.9571351  0.7920498 0.9994502
# CM06_29-CM02_05  -0.39087601 -1.2654685  0.4837164 0.6700832
# CM03_47-CM02_91A  1.80053419  0.9259417  2.6751266 0.0001820***
# CM06_28-CM02_91A -0.40532096 -1.2799134  0.4692715 0.6383547
# CM06_29-CM02_91A -0.71365433 -1.5882468  0.1609381 0.1371568
# CM06_28-CM03_47  -2.20585515 -3.0804476 -1.3312627 0.0000241***
# CM06_29-CM03_47  -2.51418852 -3.3887810 -1.6395961 0.0000061***
# CM06_29-CM06_28  -0.30833337 -1.1829258  0.5662591 0.8361212


PBM_avg <- gcvr_results %>% 
  group_by(strain, sugar) %>% 
  dplyr::summarise(mean_auc_l = mean(auc_l), mean_auc_e = mean(auc_e), mean_t = mean(t_gen), 
                   n_auc_l = length(auc_l),
                   sd_auc_l = sd(auc_l, na.rm = TRUE),
                   se_auc_l = sd_auc_l / sqrt(n_auc_l))

## Negatives for PBM plates; this section determines the upper and lower limits for no, strong and weak growth
PBM_neg <- subset(PBM, PBM$group=="negative")

PBM_neg2 <- select(PBM_neg, c(5,8, 1, 15,19,20)) %>% 
  mutate(wellID_strain = paste0(expID, "_", wellID)) %>% 
  select(-1,-2) %>% 
  pivot_wider(id_cols= c("strain", "time"), names_from = "wellID_strain", values_from = "OD600")

results_list_neg <- list()
for (strain_id in strain_list) {
  
  temp_df_neg <- PBM_neg2 %>%
    filter(strain == strain_id) %>%
    select(-1)  # assuming first column is time or index you want removed
  
  if (nrow(temp_df_neg) > 0) {  # only run if subset has data
    gcvr_neg <- SummarizeGrowthByPlate(temp_df_neg)
    results_list_neg[[strain_id]] <- gcvr_neg
  } else {
    warning(paste("No data for strain:", strain_id))
    results_list_neg[[strain_id]] <- NULL
  }
}
final_results_neg <- dplyr::bind_rows(results_list_neg, .id = "strain")
neg_mean <- mean(final_results_neg$auc_l, na.rm = TRUE)
neg_sd <- sd(final_results_neg$auc_l, na.rm = TRUE)
cutoff <- neg_mean + 3 * neg_sd
cutoff
neg_mean
neg_sd
upper_cutoff <- cutoff
lower_cutoff <- neg_mean+neg_sd
upper_cutoff
lower_cutoff


## numbers determined from the negatives
PBM_avg$gstat[PBM_avg$mean_auc_l < 0.6672] <- "no_growth"
PBM_avg$gstat[PBM_avg$mean_auc_l > 0.6672 & PBM_avg$mean_auc_l < 1.65822] <- "weak_growth"
PBM_avg$gstat[PBM_avg$mean_auc_l > 1.65822] <- "strong_growth"

## figure 8A
ggplot(PBM_avg, aes(x = strain, y= sugar, fill = gstat)) + geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  coord_fixed() +
  scale_fill_manual(values = c("grey70", "navy", "orchid")) +
  scale_x_discrete(limits = c("CM02_91A", "CM06_28", "CM02_02", "CM02_05", "CM03_47", "CM06_29"))


## fig 8B

final_results_lac_bmaa$media <- "bmaa"
final_results_lac_ycfa$media <- "ycfa"
fig8b_gcvr <- rbind(final_results_lac_bmaa, final_results_lac_ycfa)
fig8b_gcvr$sugar[fig8b_gcvr$sample %in% c("ns1", "ns2", "ns3")] <- "no_sugar"
fig8b_avg <- fig8b_gcvr %>% 
  group_by(strain, sugar, media) %>% 
  dplyr::summarise(mean_auc_l = mean(auc_l), mean_auc_e = mean(auc_e), mean_t = mean(t_gen), 
                   n_auc_l = length(auc_l),
                   sd_auc_l = sd(auc_l, na.rm = TRUE),
                   se_auc_l = sd_auc_l / sqrt(n_auc_l))

fig8b_avg_bmaa <- subset(fig8b_avg, fig8b_avg$media == "bmaa")
fig8b_avg_bmaa$gstat[fig8b_avg_bmaa$mean_auc_l < 0.9243027] <- "no_growth"
fig8b_avg_bmaa$gstat[fig8b_avg_bmaa$mean_auc_l > 0.9243027 & fig8b_avg_bmaa$mean_auc_l < 1.767432] <- "weak_growth"
fig8b_avg_bmaa$gstat[fig8b_avg_bmaa$mean_auc_l > 1.767432] <- "strong_growth"

fig8b_avg_ycfa <- subset(fig8b_avg, fig8b_avg$media == "ycfa")
fig8b_avg_ycfa$gstat[fig8b_avg_ycfa$mean_auc_l < 15.3859399] <- "no_growth"
fig8b_avg_ycfa$gstat[fig8b_avg_ycfa$mean_auc_l > 15.3859399 & fig8b_avg_ycfa$mean_auc_l < 21.017107] <- "weak_growth"
fig8b_avg_ycfa$gstat[fig8b_avg_ycfa$mean_auc_l > 21.017107] <- "strong_growth"

fig8b_avg <- rbind(fig8b_avg_bmaa,fig8b_avg_ycfa)
fig8b_avg$sugar_media <- paste0(fig8b_avg$media, "_", fig8b_avg$sugar)
fig8b_avg2 <- subset(fig8b_avg, !(fig8b_avg$sugar %in% c("no_sugar")))
## figure 8B
ggplot(fig8b_avg2, aes(x = strain, y= sugar_media, fill = gstat)) + geom_tile() +
  #facet_wrap(~media) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  coord_fixed() +
  scale_fill_manual(values = c("grey70", "navy", "orchid")) +
  scale_x_discrete(limits = c("CM02_91A", "CM06_28", "CM02_02", "CM02_05", "CM03_47", "CM06_29"))

## figure 8C

sugar_invitro <- read_excel("~/Data/figure_8_Supplementary_figure_7/sugar_invitro.xlsb.xlsx")
sugar_invitro$Actual[sugar_invitro$Actual==1] <- 2
sugar_invitro$total <- sugar_invitro$Predicted + sugar_invitro$Actual

ggplot(sugar_invitro, aes(y=Sugar, x=Strain, fill=as.factor(total))) + geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  coord_fixed() +
  scale_fill_manual(values = c("0" ="grey70", "1"="slateblue", "2"="green3", "3"="darkgreen")) +
  scale_x_discrete(limits = c("CM02_91", "CM06_28", "CM02_2", "CM02_5", "CM03_47", "CM06_29"))

