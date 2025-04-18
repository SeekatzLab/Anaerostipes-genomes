## Supplementary figure 5
# Authors:
# - all figures: DB
## Libraries
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

# Figure S5
## PathoFact predictions
# curated_list_anaero <- read_tsv("/anaero_final_list") ## list of genomes
# file_list <- as.vector(curated_list_anaero$list)
# list.data.path<-list()
# filepath_path = "//Pathofact/path_final/"
# for(i in (1:length(file_list))){
#   #print(i)
#   list.data.path[[i]] <- read_tsv(paste0(filepath_path, 'Pathofact_',file_list[i], '_predictions.tsv'))
# }
# all_path <- list.data.path[[1]]
# all_path$genome_name <- "CM01_24_S247" ## whatever your first genome name is
# for (i in 2:126) {
#   all_path <- rbind.fill(all_path, list.data.path[[i]])
#   all_path <- mutate(all_path, genome_name = replace(genome_name, is.na(genome_name), file_list[[i]]))
# }
# write_tsv(all_path, file = "/Data/sfigure5/pathofact_anaero126.tsv") ## not used in the following ana,ysis; email corresponding author for this file
# 
# ## Toxin prediciton report from pathofact
# list.data.tox<-list()
# filepath_path = "/path_final/"
# for(i in (1:length(file_list))){
#   #print(i)
#   list.data.tox[[i]] <- read_tsv(paste0(filepath_path, 'Toxin_prediction_',file_list[i], '_report.tsv'))
# }
# all_tox <- list.data.tox[[1]]
# all_tox$genome_name <- "CM01_24_S247" ## whatever your first genomes name is
# for (i in 2:126) {
#   all_tox <- rbind.fill(all_tox, list.data.tox[[i]])
#   all_tox <- mutate(all_tox, genome_name = replace(genome_name, is.na(genome_name), file_list[[i]]))
# }
# write_tsv(all_tox, file = "/Data/sfigure5/tox_anaero126.tsv") 
# 
# ## Toxin library 
# list.data.toxlib<-list()
# filepath_path = "/path_final/"
# for(i in (1:length(file_list))){
#   #print(i)
#   list.data.toxlib[[i]] <- read_tsv(paste0(filepath_path, 'Toxin_gene_library_',file_list[i], '_report.tsv'))
# }
# all_toxlib <- list.data.toxlib[[1]]
# all_toxlib$genome_name <- "CM01_24_S247" ## whatever your first genome name is
# for (i in 2:126) {
#   all_toxlib <- rbind.fill(all_toxlib, list.data.toxlib[[i]])
#   all_toxlib <- mutate(all_toxlib, genome_name = replace(genome_name, is.na(genome_name), file_list[[i]]))
# }
# write_tsv(all_toxlib, file = "/Data/sfigure5/toxlib_anaero126.tsv")

all_tox <- read_tsv("/Data/sfigure5/tox_anaero126.tsv")
all_toxlib <- read_tsv("/Data/sfigure5/toxlib_anaero126.tsv")
all_tox$genorf <- paste0(all_tox$genome_name,'_', all_tox$ORF)
all_toxlib$genorf <- paste0(all_toxlib$genome_name,'_', all_toxlib$ORF)
full_tox <- inner_join(all_toxlib, all_tox, by = "genorf")
full_tox2 <- full_tox %>% select(genome_name.x, Significance_evalue, Toxin_confidence_level, Toxin_prediction, Description) %>% distinct()
full_tox3 <- inner_join(full_tox2, anaero_meta, by = c("genome_name.x" = "list"))
full_tox4 <- subset(full_tox3, full_tox3$species == "hadrus" | full_tox3$species == "hominis" | full_tox3$species == "caccae")
full_tox4$prab <- 1
full_tox7 <- full_tox4 %>% arrange(factor(species, levels = c("caccae", "hominis", "hadrus")))
## final figure
ggplot(full_tox7, aes(x=genome_name, y = Description, color = factor(prab))) + geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle =90), panel.grid = element_blank()) +
  scale_color_manual(values = c("#96296E")) +
  scale_x_discrete(limits = unique(full_tox7$genome_name))
## saved as pdf
