# Figure 4:

# Authors:
  # - all figures: DB

## Input files:
  ## pcoa_anaero126.tsv
  ## anaero_meta_final.txt
  ## cog-24.def.tab.txt
  ## functl_enr_newcluster.txt
  ## cazy_anaero126.tsv

## Output files:
  ## Figure 4A pdf
  ## Figure 4B pdf
  ## Figure 4C pdf
  ## Figure 4D pdf
  ## Figure 4E pdf

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

################################
# Figure 4A

## specific libraries
library(vegan)
library(ape)
library(factoextra)
library(cluster)
library(ggforce)
#install.packages('devtools')
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)


# curated_list_anaero <- read_tsv("/anaero_final_list.txt") ## this is just a list of the genome names
# file_list <- as.vector(curated_list_anaero$list)
# ## create an empty list
# list.data.pcoa<-list()
# ### file path to .tsv from prokka
# filepath_pcoa = "/prokka_tsv/"
# 
# for(i in (1:length(file_list))){
#   #print(i) ## checks if list will be created correctly
#   list.data.pcoa[[i]] <- read_tsv(paste0(filepath_pcoa, file_list[i], '.tsv'))
# }
# 
# ### combining all the .tsv from a list to one large data frame and simultaneously adding their species name
# 
# all_pcoa <- list.data.pcoa[[1]]
# all_pcoa$genome_name <- "CM01_24_S247" ## whatever the first name is on YOUR curated list
# 
# for (i in 2:126) {
#   all_pcoa <- rbind.fill(all_pcoa, list.data.pcoa[[i]])
#   all_pcoa <- mutate(all_pcoa, genome_name = replace(genome_name, is.na(genome_name), file_list[[i]]))
# }
# 
# write_tsv(all_pcoa, file = "/Data/figure3/pcoa_anaero126.tsv")

all_pcoa <- read_tsv("/Data/figure4/pcoa_anaero126.tsv")
prokka_cog <- subset(all_pcoa, !is.na(all_pcoa$product))
pcoa_cog <- select(prokka_cog, product, genome_name) %>% distinct()
pcoa_cog$prab <- 1
## assigned relative abundance to the ones present as 1, now assign the ones absent as 0
pcoa_cog2 <- pcoa_cog %>% pivot_wider(names_from = product, values_from = prab, values_fill = 0)
pcoa_cog3 <- pcoa_cog2 
pcoa_cog3 <- pcoa_cog3[,-1]
rownames(pcoa_cog3) <- pcoa_cog2$genome_name

pcoa_cog4 <- as.matrix(pcoa_cog3)
## vegdist or avedist to generate distance 
set.seed(12345)
bray_dist_pcoa <- vegdist(pcoa_cog4, method = "jaccard")
## pcoa calculation
pcoa_vals <- pcoa(bray_dist_pcoa)
pcoa_axes <- pcoa_vals$vectors[,1:2]
pcoa_axes2 <- as.data.frame(pcoa_axes) %>% rownames_to_column("genome_name")
anaero_meta <- read_tsv("/Data/figure4/anaero_meta_final.txt")

pcoa_axes_meta <- inner_join(pcoa_axes2, anaero_meta, by = c('genome_name' = 'list'))
## calculate proportion variation 
pcoa_eigen <- as.data.frame(pcoa_vals$values$Eigenvalues)
prop_var <- (pcoa_eigen$`pcoa_vals$values$Eigenvalues`/sum(pcoa_eigen$`pcoa_vals$values$Eigenvalues`)) * 100
## assign dark boundary to only the ones associated with this study
pcoa_axes_meta$isolate_id <- "0"
pcoa_axes_meta$isolate_id[grep("^CM", pcoa_axes_meta$genome_name)] <- "this_study"
pcoa_axes_meta$isolate_id[pcoa_axes_meta$isolate_id == "0"] <- "other"

## clustering with PAM 
# library(factoextra)
# library(cluster)
# library(ggforce)

kmed <- pam(pcoa_axes, k=3)
sil <- silhouette(kmed$clustering, dist(pcoa_axes))
mean(sil[, 'sil_width'])
kmed_cluster <- as.data.frame(kmed$clustering) %>% rownames_to_column("genome_name")
colnames(kmed_cluster) <- c("genome_name", "clust")
pcoa_axes_meta2 <- inner_join(pcoa_axes_meta, kmed_cluster, by = "genome_name")
## final figure
ggplot(pcoa_axes_meta2, aes(x=Axis.1, y=Axis.2,fill=species, color=isolate_id)) + geom_point(shape = 21, size=5) + 
  theme_bw() +
  scale_fill_manual(breaks = pcoa_axes_meta2$species, values =  pcoa_axes_meta2$color_species) +
  scale_color_manual(breaks = c('this_study', 'other'), values = c("black", rgb(0, 0, 0, alpha=0))) +
  #geom_mark_ellipse(data= pcoa_axes_meta2, aes(fill = as.factor(clust))) + ## used it to demarcate the 3 clusters
  theme(panel.grid = element_blank()) + xlab("Axis.1 (63.2%)") + ylab("Axis.2 (10.47%)") +
  theme(text = element_text(size = 10), panel.background = element_rect(colour = "black", size=0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed()
## save as pdf
 
## Statistics: permanova

adonis2(bray_dist_pcoa ~ species, data=pcoa_axes_meta, method = "jaccard")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray_dist_pcoa ~ species, data = pcoa_axes_meta, method = "jaccard")
# Df SumOfSqs      R2      F Pr(>F)    
# species   11  2.16001 0.87028 69.531  0.001 ***
#   Residual 114  0.32195 0.12972                  
# Total    125  2.48196 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

bd <- betadisper(d = bray_dist_pcoa, group = pcoa_axes_meta$species, type="median")
permutest(bd, pairwise = TRUE)
anova(bd)
# Analysis of Variance Table
# 
# Response: Distances
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Groups     11 0.086286 0.0078442   11.16 8.514e-14 ***
#   Residuals 114 0.080127 0.0007029                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pairwise.adonis2(bray_dist_pcoa ~ species, data=pcoa_axes_meta, method = "jaccard")

# $parent_call
# [1] "bray_dist_pcoa ~ species , strata = Null , permutations 999"
# 
# $hadrus_vs_caccae
# Df SumOfSqs      R2      F Pr(>F)    
# species   1  0.95073 0.76259 311.57  0.001 ***
#   Residual 97  0.29598 0.23741                  
# Total    98  1.24671 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $hadrus_vs_sp000508985
# Df SumOfSqs     R2      F Pr(>F)    
# species   1  0.75454 0.7369 260.48  0.001 ***
#   Residual 93  0.26940 0.2631                  
# Total    94  1.02394 1.0000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $hadrus_vs_butyraticus
# Df SumOfSqs      R2      F Pr(>F)    
# species   1  0.17591 0.39773 56.132  0.001 ***
#   Residual 85  0.26637 0.60227                  
# Total    86  0.44228 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $hadrus_vs_sp018381315
# Df SumOfSqs      R2      F Pr(>F)  
# species   1  0.07527 0.22463 24.045  0.015 *
#   Residual 83  0.25981 0.77537                
# Total    84  0.33508 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $hadrus_vs_excrementavium
# Df SumOfSqs      R2      F Pr(>F)   
# species   1  0.15197 0.36578 48.446  0.002 **
#   Residual 84  0.26350 0.63422                 
# Total    85  0.41547 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $hadrus_vs_avistercoris
# Df SumOfSqs      R2      F Pr(>F)   
# species   1  0.05740 0.18095 18.336  0.004 **
#   Residual 83  0.25981 0.81905                 
# Total    84  0.31721 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $hadrus_vs_sp900066705
# Df SumOfSqs      R2      F Pr(>F)   
# species   1 0.036752 0.12393 11.883  0.002 **
#   Residual 84 0.259811 0.87607                 
# Total    85 0.296563 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $hadrus_vs_hadrus_A
# Df SumOfSqs      R2      F Pr(>F)    
# species   1  0.05402 0.16882 17.468  0.001 ***
#   Residual 86  0.26594 0.83118                  
# Total    87  0.31996 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $hadrus_vs_rhamnosivorans
# Df SumOfSqs      R2      F Pr(>F)  
# species   1  0.07747 0.22969 24.749  0.012 *
#   Residual 83  0.25981 0.77031                
# Total    84  0.33728 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $hadrus_vs_sp018918155
# Df SumOfSqs      R2      F Pr(>F)  
# species   1  0.06655 0.20392 21.261  0.012 *
#   Residual 83  0.25981 0.79608                
# Total    84  0.32636 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $hadrus_vs_sp001940315
# Df SumOfSqs      R2      F Pr(>F)  
# species   1  0.05807 0.18268 18.552  0.013 *
#   Residual 83  0.25981 0.81732                
# Total    84  0.31788 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $caccae_vs_sp000508985
# Df SumOfSqs     R2      F Pr(>F)    
# species   1 0.083030 0.6447 43.549  0.001 ***
#   Residual 24 0.045758 0.3553                  
# Total    25 0.128788 1.0000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $caccae_vs_butyraticus
# Df SumOfSqs      R2      F Pr(>F)   
# species   1 0.186733 0.81377 69.914  0.002 **
#   Residual 16 0.042734 0.18623                 
# Total    17 0.229467 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $caccae_vs_sp018381315
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.012130 0.25112 4.6947  0.123
# Residual 14 0.036172 0.74888              
# Total    15 0.048302 1.00000              
# 
# $caccae_vs_excrementavium
# Df SumOfSqs      R2      F Pr(>F)   
# species   1 0.148851 0.78878 56.016  0.005 **
#   Residual 15 0.039859 0.21122                 
# Total    16 0.188710 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $caccae_vs_avistercoris
# Df SumOfSqs      R2     F Pr(>F)  
# species   1 0.070924 0.66225 27.45  0.071 .
# Residual 14 0.036172 0.33775               
# Total    15 0.107096 1.00000               
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $caccae_vs_sp900066705
# Df SumOfSqs      R2      F Pr(>F)   
# species   1 0.160938 0.81649 66.739  0.009 **
#   Residual 15 0.036172 0.18351                 
# Total    16 0.197111 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $caccae_vs_hadrus_A
# Df SumOfSqs      R2     F Pr(>F)    
# species   1  0.28244 0.86973 113.5  0.001 ***
#   Residual 17  0.04231 0.13027                 
# Total    18  0.32475 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $caccae_vs_rhamnosivorans
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.012195 0.25214 4.7201  0.128
# Residual 14 0.036172 0.74786              
# Total    15 0.048368 1.00000              
# 
# $caccae_vs_sp018918155
# Df SumOfSqs      R2      F Pr(>F)  
# species   1 0.084353 0.69988 32.648  0.062 .
# Residual 14 0.036172 0.30012                
# Total    15 0.120525 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $caccae_vs_sp001940315
# Df SumOfSqs      R2      F Pr(>F)  
# species   1 0.082748 0.69583 32.026  0.059 .
# Residual 14 0.036172 0.30417                
# Total    15 0.118920 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $sp000508985_vs_butyraticus
# Df SumOfSqs      R2     F Pr(>F)   
# species   1 0.174126 0.91513 129.4  0.003 **
#   Residual 12 0.016148 0.08487                
# Total    13 0.190274 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $sp000508985_vs_sp018381315
# Df  SumOfSqs      R2      F Pr(>F)  
# species   1 0.0203843 0.68015 21.265  0.089 .
# Residual 10 0.0095858 0.31985                
# Total    11 0.0299701 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $sp000508985_vs_excrementavium
# Df SumOfSqs      R2     F Pr(>F)  
# species   1 0.144074 0.91565 119.4  0.017 *
#   Residual 11 0.013273 0.08435               
# Total    12 0.157347 1.00000               
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $sp000508985_vs_avistercoris
# Df SumOfSqs      R2      F Pr(>F)  
# species   1 0.067787 0.87611 70.716  0.084 .
# Residual 10 0.009586 0.12389                
# Total    11 0.077373 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $sp000508985_vs_sp900066705
# Df SumOfSqs      R2      F Pr(>F)   
# species   1 0.161957 0.94412 185.85  0.008 **
#   Residual 11 0.009586 0.05588                 
# Total    12 0.171543 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $sp000508985_vs_hadrus_A
# Df SumOfSqs     R2      F Pr(>F)    
# species   1 0.272696 0.9455 225.52  0.001 ***
#   Residual 13 0.015719 0.0545                  
# Total    14 0.288415 1.0000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $sp000508985_vs_rhamnosivorans
# Df SumOfSqs      R2      F Pr(>F)  
# species   1 0.023834 0.71317 24.864  0.076 .
# Residual 10 0.009586 0.28683                
# Total    11 0.033420 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $sp000508985_vs_sp018918155
# Df SumOfSqs      R2      F Pr(>F)  
# species   1 0.083887 0.89745 87.511  0.095 .
# Residual 10 0.009586 0.10255                
# Total    11 0.093472 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $sp000508985_vs_sp001940315
# Df SumOfSqs      R2      F Pr(>F)  
# species   1 0.080178 0.89321 83.642  0.092 .
# Residual 10 0.009586 0.10679                
# Total    11 0.089764 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $butyraticus_vs_sp018381315
# Df SumOfSqs      R2     F Pr(>F)
# species   1 0.057843 0.89811 17.63   0.25
# Residual  2 0.006562 0.10189             
# Total     3 0.064405 1.00000             
# 
# $butyraticus_vs_excrementavium
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.043958 0.81093 12.867    0.1
# Residual  3 0.010249 0.18907              
# Total     4 0.054207 1.00000              
# 
# $butyraticus_vs_avistercoris
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.010495 0.61528 3.1986   0.25
# Residual  2 0.006562 0.38472              
# Total     3 0.017056 1.00000              
# 
# $butyraticus_vs_sp900066705
# Df SumOfSqs     R2      F Pr(>F)
# species   1 0.075670 0.9202 34.595    0.1
# Residual  3 0.006562 0.0798              
# Total     4 0.082232 1.0000              
# 
# $butyraticus_vs_hadrus_A
# Df SumOfSqs      R2      F Pr(>F)  
# species   1 0.124376 0.90738 48.984  0.025 *
#   Residual  5 0.012695 0.09262                
# Total     6 0.137071 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $butyraticus_vs_rhamnosivorans
# Df SumOfSqs      R2    F Pr(>F)
# species   1 0.060041 0.90148 18.3   0.25
# Residual  2 0.006562 0.09852            
# Total     3 0.066602 1.00000            
# 
# $butyraticus_vs_sp018918155
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.062954 0.90561 19.188   0.25
# Residual  2 0.006562 0.09439              
# Total     3 0.069516 1.00000              
# 
# $butyraticus_vs_sp001940315
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.040267 0.85987 12.273   0.25
# Residual  2 0.006562 0.14013              
# Total     3 0.046829 1.00000              
# 
# $sp018381315_vs_excrementavium
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.058709 0.94091 15.923 0.3333
# Residual  1 0.003687 0.05909              
# Total     2 0.062396 1.00000              
# 
# $sp018381315_vs_avistercoris
# Df SumOfSqs R2 F Pr(>F)
# Model     1 0.036854  1         
# Residual  0 0.000000  0         
# Total     1 0.036854  1         
# 
# $sp018381315_vs_sp900066705
# Df SumOfSqs R2          F Pr(>F)
# species   1 0.066844  1 2.3597e+16 0.3333
# Residual  1 0.000000  0                  
# Total     2 0.066844  1                  
# 
# $sp018381315_vs_hadrus_A
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.071992 0.92149 35.212    0.2
# Residual  3 0.006134 0.07851              
# Total     4 0.078125 1.00000              
# 
# $sp018381315_vs_rhamnosivorans
# Df  SumOfSqs R2 F Pr(>F)
# Model     1 0.0064961  1         
# Residual  0 0.0000000  0         
# Total     1 0.0064961  1         
# 
# $sp018381315_vs_sp018918155
# Df SumOfSqs R2 F Pr(>F)
# Model     1 0.044113  1         
# Residual  0 0.000000  0         
# Total     1 0.044113  1         
# 
# $sp018381315_vs_sp001940315
# Df SumOfSqs R2 F Pr(>F)
# Model     1 0.042535  1         
# Residual  0 0.000000  0         
# Total     1 0.042535  1         
# 
# $excrementavium_vs_avistercoris
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.027619 0.88223 7.4909 0.3333
# Residual  1 0.003687 0.11777              
# Total     2 0.031306 1.00000              
# 
# $excrementavium_vs_sp900066705
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.091380 0.96122 49.568 0.3333
# Residual  2 0.003687 0.03878              
# Total     3 0.095067 1.00000              
# 
# $excrementavium_vs_hadrus_A
# Df SumOfSqs      R2      F  Pr(>F)  
# species   1 0.126334 0.92787 51.457 0.06667 .
# Residual  4 0.009821 0.07213                 
# Total     5 0.136155 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $excrementavium_vs_rhamnosivorans
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.060747 0.94278 16.476 0.3333
# Residual  1 0.003687 0.05722              
# Total     2 0.064434 1.00000              
# 
# $excrementavium_vs_sp018918155
# Df SumOfSqs     R2      F Pr(>F)
# species   1 0.063722 0.9453 17.283 0.3333
# Residual  1 0.003687 0.0547              
# Total     2 0.067409 1.0000              
# 
# $excrementavium_vs_sp001940315
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.052529 0.93441 14.247 0.3333
# Residual  1 0.003687 0.06559              
# Total     2 0.056216 1.00000              
# 
# $avistercoris_vs_sp900066705
# Df SumOfSqs R2          F Pr(>F)
# species   1  0.04307  1 3.0409e+16 0.3333
# Residual  1  0.00000  0                  
# Total     2  0.04307  1                  
# 
# $avistercoris_vs_hadrus_A
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.054268 0.89846 26.544    0.2
# Residual  3 0.006134 0.10154              
# Total     4 0.060402 1.00000              
# 
# $avistercoris_vs_rhamnosivorans
# Df SumOfSqs R2 F Pr(>F)
# Model     1 0.040678  1         
# Residual  0 0.000000  0         
# Total     1 0.040678  1         
# 
# $avistercoris_vs_sp018918155
# Df SumOfSqs R2 F Pr(>F)
# Model     1 0.039675  1         
# Residual  0 0.000000  0         
# Total     1 0.039675  1         
# 
# $avistercoris_vs_sp001940315
# Df SumOfSqs R2 F Pr(>F)
# Model     1 0.026147  1         
# Residual  0 0.000000  0         
# Total     1 0.026147  1         
# 
# $sp900066705_vs_hadrus_A
# Df  SumOfSqs      R2    F  Pr(>F)  
# species   1 0.0220803 0.78261 14.4 0.06667 .
# Residual  4 0.0061335 0.21739               
# Total     5 0.0282138 1.00000               
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $sp900066705_vs_rhamnosivorans
# Df SumOfSqs R2          F Pr(>F)
# species   1 0.061596  1 1.5991e+31 0.3333
# Residual  1 0.000000  0                  
# Total     2 0.061596  1                  
# 
# $sp900066705_vs_sp018918155
# Df SumOfSqs R2           F Pr(>F)
# species   1 0.040687  1 -2.8726e+16      1
# Residual  1 0.000000  0                   
# Total     2 0.040687  1                   
# 
# $sp900066705_vs_sp001940315
# Df SumOfSqs R2          F Pr(>F)
# species   1 0.043999  1 3.1064e+16 0.3333
# Residual  1 0.000000  0                  
# Total     2 0.043999  1                  
# 
# $hadrus_A_vs_rhamnosivorans
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.072337 0.92184 35.381    0.2
# Residual  3 0.006134 0.07816              
# Total     4 0.078470 1.00000              
# 
# $hadrus_A_vs_sp018918155
# Df SumOfSqs      R2      F Pr(>F)
# species   1 0.052083 0.89464 25.475    0.2
# Residual  3 0.006134 0.10536              
# Total     4 0.058217 1.00000              
# 
# $hadrus_A_vs_sp001940315
# Df SumOfSqs     R2    F Pr(>F)
# species   1 0.049272 0.8893 24.1    0.2
# Residual  3 0.006134 0.1107            
# Total     4 0.055405 1.0000            
# 
# $rhamnosivorans_vs_sp018918155
# Df SumOfSqs R2 F Pr(>F)
# Model     1    0.047  1         
# Residual  0    0.000  0         
# Total     1    0.047  1         
# 
# $rhamnosivorans_vs_sp001940315
# Df SumOfSqs R2 F Pr(>F)
# Model     1 0.047762  1         
# Residual  0 0.000000  0         
# Total     1 0.047762  1         
# 
# $sp018918155_vs_sp001940315
# Df SumOfSqs R2 F Pr(>F)
# Model     1 0.036514  1         
# Residual  0 0.000000  0         
# Total     1 0.036514  1         
# 
# attr(,"class")
# [1] "pwadstrata" "list"      

##########################################

# Figure 4B:

## specific libraries
library(ggvenn)

## input files
all_pcoa <- read_tsv("/Data/figure4/pcoa_anaero126.tsv")
gene_pres <- all_pcoa %>% select(genome_name, product, COG) %>% distinct()
#na <- subset(all_pcoa, is.na(all_pcoa$product))
gene_pres2 <- subset(gene_pres, !is.na(gene_pres$product))
gene_pres3 <- gene_pres 
gene_pres3$prab <- 1
gene_pres4 <- inner_join(gene_pres3, anaero_meta, by = c('genome_name' = 'list'))

tally_gene <- gene_pres4 %>% group_by(product) %>% tally()

gene_hadrus <- as.vector(unique(gene_pres4$product[gene_pres4$cluster_sp=="hadrus"]))
gene_caccae <- as.vector(unique(gene_pres4$product[gene_pres4$cluster_sp=="caccae"]))
gene_chicken <- (as.vector(unique(gene_pres4$product[gene_pres4$cluster_sp=="chicken"])))

venn <- list(hadrus = gene_hadrus, caccae = gene_caccae, chicken = gene_chicken)
## final figure
ggvenn(venn)
## save as pdf

###############################

# Figure 4C:

cog <- read_tsv("/Data/figure4/cog-24.def.tab.txt", col_names = c("COG", "COG_functional_category", "COG_name", "Gene_name", "Functional_pathway", "PubMed_ID", "PDB_ID"))
# all_pcoa <- read_tsv("/Data/figure3/pcoa_anaero126.tsv")
# gene_pres <- all_pcoa %>% select(genome_name, product, COG) %>% distinct()
# #na <- subset(all_pcoa, is.na(all_pcoa$product))
# gene_pres2 <- subset(gene_pres, !is.na(gene_pres$product))
# gene_pres3 <- gene_pres 
# gene_pres3$prab <- 1
# gene_pres4 <- inner_join(gene_pres3, anaero_meta, by = c('genome_name' = 'list'))
## gene_pres4 is from Figure 3B
gene_pres5 <- subset(gene_pres4, !is.na(gene_pres4$COG)) 
gene_hadrus2 <- subset(gene_pres5, gene_pres5$cluster_paper == "A")
gene_caccae2 <- subset(gene_pres5, gene_pres5$cluster_paper == "B")
gene_chicken2 <- subset(gene_pres5, gene_pres5$cluster_paper == "C")
core <- as.data.frame(Reduce(intersect, list(gene_hadrus2$COG, gene_caccae2$COG, gene_chicken2$COG)))
colnames(core) <- "COG"
core2 <- inner_join(core, cog, by = "COG")
core3 <- core2 %>% separate_rows(COG_functional_category, sep = "") %>% distinct()
core4 <- subset(core3, !(core3$COG_functional_category==""))
core4$group <- "core"
summ_core <- core4 %>% group_by(COG_functional_category,group) %>% tally()
summ_core$perc <- (summ_core$n/sum(summ_core$n))*100

a <- intersect(gene_hadrus2$COG, gene_caccae2$COG)
b <- intersect(gene_hadrus2$COG, gene_chicken2$COG)
gene_hadrus_unique <- subset(gene_hadrus2, !(gene_hadrus2$COG %in% a) & !(gene_hadrus2$COG %in% b)) 
gene_hadrus_unique <- as.data.frame(unique(gene_hadrus_unique$COG))
colnames(gene_hadrus_unique) <- "COG"
c <- intersect(gene_caccae2$COG, gene_chicken2$COG)
gene_caccae_unique <- subset(gene_caccae2, !(gene_caccae2$COG %in% a) & !(gene_caccae2$COG %in% c))
gene_caccae_unique <- as.data.frame(unique(gene_caccae_unique$COG))
colnames(gene_caccae_unique) <- "COG"
gene_chicken_unique <- subset(gene_chicken2, !(gene_chicken2$COG %in% b) & !(gene_chicken2$COG %in% c))
gene_chicken_unique <- as.data.frame(unique(gene_chicken_unique$COG))
colnames(gene_chicken_unique) <- "COG"

hadrus_unique <- inner_join(gene_hadrus_unique, cog, by = "COG")
hadrus_unique2 <- hadrus_unique %>% separate_rows(COG_functional_category, sep = "") %>% distinct()
hadrus_unique3 <- subset(hadrus_unique2, !(hadrus_unique2$COG_functional_category==""))
hadrus_unique3$group <- "hadrus"
summ_hadrus <- hadrus_unique3 %>% group_by(COG_functional_category, group) %>% tally()
summ_hadrus$perc <- (summ_hadrus$n/sum(summ_hadrus$n))*100

caccae_unique <- inner_join(gene_caccae_unique, cog, by = "COG")
caccae_unique2 <- caccae_unique %>% separate_rows(COG_functional_category, sep = "") %>% distinct()
caccae_unique3 <- subset(caccae_unique2, !(caccae_unique2$COG_functional_category==""))
caccae_unique3$group <- "caccae"
summ_caccae <- caccae_unique3 %>% group_by(COG_functional_category, group) %>% tally()
summ_caccae$perc <- (summ_caccae$n/sum(summ_caccae$n))*100

chicken_unique <- inner_join(gene_chicken_unique, cog, by = "COG")
chicken_unique2 <- chicken_unique %>% separate_rows(COG_functional_category, sep = "") %>% distinct()
chicken_unique3 <- subset(chicken_unique2, !(chicken_unique2$COG_functional_category==""))
chicken_unique3$group <- "chicken"
summ_chicken <- chicken_unique3 %>% group_by(COG_functional_category,group) %>% tally()
summ_chicken$perc <- (summ_chicken$n/sum(summ_chicken$n))*100

final_summ <- rbind(summ_core, summ_caccae, summ_hadrus, summ_chicken)

color_pal <- as.data.frame(unique(final_summ$COG_functional_category))
colnames(color_pal) <- "newcog"
color_pal$color[color_pal$newcog == "R" | color_pal$newcog == "S" | color_pal$newcog == "W" | color_pal$newcog == "Not_Assigned"] <- "grey67"
color_pal$color[color_pal$newcog == "J" | color_pal$newcog == "K" | color_pal$newcog == "L"] <- "#A71B4B"
color_pal$color[color_pal$newcog == "D" | color_pal$newcog == "V" | color_pal$newcog == "T" | color_pal$newcog == "M" | color_pal$newcog == "N" | color_pal$newcog == "O" | color_pal$newcog == "U"] <- "#D04939"
color_pal$color[color_pal$newcog == "C"] <- "#EB7803"
color_pal$color[color_pal$newcog == "G"] <- "#F9BC53"
color_pal$color[color_pal$newcog == "E"] <- "#FEF1A6"
color_pal$color[color_pal$newcog == "F"] <- "#E2F8B5"
color_pal$color[color_pal$newcog == "H"] <- "#9CE5AD"
color_pal$color[color_pal$newcog == "I"] <- "#43CBB1"
color_pal$color[color_pal$newcog == "P"] <- "#00AAB6"
color_pal$color[color_pal$newcog == "Q"] <- "#0080B2"
color_pal$color[color_pal$newcog == "X"] <- "#584B9F"

ord <- c("Not_Assigned", "R", "S", "W", "J", "K", "L", "D", "V", "T", "M", "N", "O", "U", "C", "G", "E", "F", "H", "I", "P", "Q", "X", "")
final_summ2 <- final_summ %>% arrange(factor(COG_functional_category, levels = ord))
final_summ3 <- inner_join(final_summ2, color_pal, by = c("COG_functional_category" = "newcog"))
color_pal <- color_pal %>% arrange(factor(newcog, levels = ord))
## final figure
ggplot(final_summ3, aes(x=factor(group, levels = c("core", "hadrus", "caccae", "chicken")), y=perc, fill=factor(COG_functional_category, levels = ord))) + geom_col() +
  theme_bw() +
  scale_fill_manual(breaks = final_summ3$COG_functional_category, values = final_summ3$color)
## saved as pdf


###################################

#Figure 4D:

funenr <- read_tsv("/Data/figure4/functl_enr_newcluster.txt")
funenr2 <- subset(funenr, funenr$adjusted_q_value < 0.005)
funenr2$func <- funenr2$`function`
funenr3 <- funenr2 %>% select(func, p_caccae, p_hadrus, p_chicken)
funenr4 <- funenr3 %>% pivot_longer(cols = !func, names_to = "group", values_to = "fraction")
## final figure
ggplot(funenr4, aes(x=factor(group, levels = c("p_hadrus", "p_caccae", "p_chicken")), y=func, fill=(fraction))) +geom_raster() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(high = "#0c7b74", low = "white") +
  coord_fixed()
## save as pdf


###################################

# Figure 4E:

## specific libraries
library(pheatmap)

# curated_list_anaero <- read_tsv("/anaero_final_list.txt")
# file_list <- as.vector(curated_list_anaero$list)
# ### create an empty list which will be populated with all the overview.txt
# list.data.cazy<-list()
# ### file path to overview.txt from run_dbcan
# filepath_cazy = "~/Library/CloudStorage/Box-Box/Disha/Projects/paper4_anaero/anaero126/dbcan4/dbcan_"
# 
# for(i in (1:length(file_list))){
#   #print(i)
#   list.data.cazy[[i]] <- read_tsv(paste0(filepath_cazy, file_list[i], '/overview.txt'))
# }
# 
# ### combining all the overview.txt from a list to one large data frame and simultaneously adding their species name
# 
# all_cazy <- list.data.cazy[[1]]
# all_cazy$genome_name <- "CM01_24_S247"
# 
# for (i in 2:126) {
#   all_cazy <- rbind.fill(all_cazy, list.data.cazy[[i]])
#   all_cazy <- mutate(all_cazy, genome_name = replace(genome_name, is.na(genome_name), file_list[[i]]))
# }
# 
# write_tsv(all_cazy, file = "~/Library/CloudStorage/Box-Box/Disha/Projects/paper4_anaero/figures/cazy_anaero126.tsv")

all_cazy <- read_tsv("/Data/figure4/cazy_anaero126.tsv")
all_cazy$newgenome_name <- all_cazy$genome_name
all_cazy$newgenome_name <- sub("^([^_]*_[^_]*)_.*$", "\\1", all_cazy$genome_name)
all_cazy$newgenome_name <- gsub("\\..*","",all_cazy$newgenome_name)

### removing brackets
super_dbcan <- all_cazy
super_dbcan$HMMER <- gsub("\\s*\\([^\\)]+\\)","",as.character(super_dbcan$HMMER))

### removing the blank spaces and the '-' from the 3 tools hmmer, hotpep and diamond so it doesnot interfere in the coalesce stage

super_dbcan$HMMER <- replace(super_dbcan$HMMER, grepl("^\\s*$", super_dbcan$HMMER) == TRUE, NA)
super_dbcan$dbCAN_sub <- replace(super_dbcan$dbCAN_sub, grepl("^\\s*$", super_dbcan$dbCAN_sub) == TRUE, NA)
super_dbcan$DIAMOND <- replace(super_dbcan$DIAMOND, grepl("^\\s*$", super_dbcan$DIAMOND) == TRUE, NA)
super_dbcan$HMMER <- replace(super_dbcan$HMMER, grepl("^\\-*$", super_dbcan$HMMER) == TRUE, NA)
super_dbcan$dbCAN_sub <- replace(super_dbcan$dbCAN_sub, grepl("^\\-*$", super_dbcan$dbCAN_sub) == TRUE, NA)
super_dbcan$DIAMOND <- replace(super_dbcan$DIAMOND, grepl("^\\-*$", super_dbcan$DIAMOND) == TRUE, NA)
super_dbcan <- super_dbcan %>% mutate(dbCAN_sub = gsub("_e[0-9]+", "", dbCAN_sub))
## coalesce
super_dbcan <- mutate(super_dbcan, cazyfam= coalesce(HMMER, dbCAN_sub, DIAMOND))
## splitting rows at cazyfam
super_dbcan2 <- separate_rows(super_dbcan, cazyfam, sep = "\\+")
### extract the first 2 letters of the cazyfam to identify the type of cazy
super_dbcan2$cazyfaminit <- substr(super_dbcan2$cazyfam, 1, 2)
### assigning category to cazyfam
super_dbcan2$cazyfamcategory[super_dbcan2$cazyfaminit=="GT"] <- "GlycosylTransferases"
super_dbcan2$cazyfamcategory[super_dbcan2$cazyfaminit=="GH"] <- "Glycoside Hydrolases"
super_dbcan2$cazyfamcategory[super_dbcan2$cazyfaminit=="CB"] <- "Carbohydrate-Binding Modules"
super_dbcan2$cazyfamcategory[super_dbcan2$cazyfaminit=="CE"] <- "Carbohydrate Esterases"
super_dbcan2$cazyfamcategory[super_dbcan2$cazyfaminit=="PL"] <- "Polysaccharide Lyases"
super_dbcan2$cazyfamcategory[super_dbcan2$cazyfaminit=="AA"] <- "Auxiliary Activities"
super_dbcan3 <- super_dbcan2 %>% group_by(genome_name, newgenome_name, cazyfaminit) %>% tally()
anaero_meta <- read_tsv("/Data/figure2/anaero_meta_final.txt")
### merge the all_meta and super_dbcan5
merge_meta <- inner_join(super_dbcan2, anaero_meta, by = c("newgenome_name" = "genome_name"))
merge_meta$cazyfam <- sub("_.*", "", merge_meta$cazyfam)
### now the entire file is ready to be counted and plotted; ordering all tables for correct coloring
summary_dbcan_species <- merge_meta %>% group_by(newgenome_name, cazyfam, cazyfamcategory, species, color_species, cluster_sp) %>% tally() %>% arrange(species)
### modifying the summary_dbcan into numric matrix
sum_mat <- summary_dbcan_species[,c(1,2,4,7)]
sum_mat1 <- sum_mat %>% pivot_wider(names_from = cazyfam, values_from = n, values_fill = 0) 
sum_mat2 <- sum_mat1 %>% column_to_rownames("newgenome_name")
sum_mat2_2 <- sum_mat2[,-1]
sum_mat3 <- t(sum_mat2_2)
ordered.list <- as.character(sort(unique(summary_dbcan_species$cazyfam)))
sum_mat3_ordered <- sum_mat3[ order(match(rownames(sum_mat3), rev(ordered.list))), ]
### for annotation_col
cluster_group2 <- merge_meta[,c(8,13,22)] %>% distinct() %>% arrange(species, newgenome_name)
cluster_group3 <- cluster_group2 %>% column_to_rownames("newgenome_name")
ann_colors <- list(species = c(amylophilus = "#0C5650", avistercoris = "#6E792B", butyraticus = "#D09D07", caccae = "#EDB700", excrementavium = "#7A3036", faecalis = "#dc252b", 
                               hadrus = "#E56635", hominis = "#f59d95", rhamnosivorans = "#CC5F79", sp018381315 = "#b3b3b3", sp018918155 = "#b3b3b3", sp900066705 = "#b3b3b3"),
                   cluster_sp = c(hadrus = "#0000ff", chicken = "#00ff00", caccae = "#ffff00"))

## for heatmap colors/breaks:
cols <- brewer.pal(7, "Blues")
# to visualize colors:
heat_colors <- c("#9E9E9E", cols[2:7])
heat_breaks <- c(0, 0.99, 1.99, 4.99, 9.99, 19.99, 30, 37)
## final figure
pheatmap(sum_mat3_ordered, cluster_rows = FALSE, annotation_col= cluster_group3, annotation_colors = ann_colors,
         cutree_cols = 3, border_color = "black", breaks = heat_breaks, color = heat_colors, 
         legend_breaks = heat_breaks, clustering_method = "complete", fontsize = 3)
## saved as pdf







