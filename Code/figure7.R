# Figure 7:

# Authors:
# - all figures: DB

## Input files:
  ## amr_anaero126.tsv
  ## anaero_meta_final.txt
  ## AN8ANPAN_gene_clusters_summary.txt
  ## spo.txt
  ## RAxML_bipartitionsBranchLabels.germination
  ## gerAB_anaero126.tsv

## Output files:
  ## Figure 7A pdf
  ## Figure 7C pdf

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

#########################################

# Figure 7A:

# list.data.amr<-list()
# filepath_path = "/path_final/"
# 
# for(i in (1:length(file_list))){
#   #print(i)
#   list.data.amr[[i]] <- read_tsv(paste0(filepath_path, 'AMR_MGE_prediction_',file_list[i], '_report.tsv'))
# }
# all_amr <- list.data.amr[[1]]
# all_amr$genome_name <- "CM01_24_S247" ## whatever is first on the list
# 
# for (i in 2:126) {
#   all_amr <- rbind.fill(all_amr, list.data.amr[[i]])
#   all_amr <- mutate(all_amr, genome_name = replace(genome_name, is.na(genome_name), file_list[[i]]))
# }
# 
# write_tsv(all_amr, file = "/Data/figure7/amr_anaero126.tsv")

all_amr <- read_tsv("/Data/figure7/amr_anaero126.tsv")
full_amr2 <- all_amr %>% select(c(5:10, 12)) %>% distinct()
anaero_meta <- read_tsv("/Data/figure7/anaero_meta_final.txt")
full_amr3 <- inner_join(full_amr2, anaero_meta, by = c("genome_name" = "list"))
full_amr4 <- subset(full_amr3, full_amr3$species == "hadrus" | full_amr3$species == "hominis" | full_amr3$species == "caccae")
full_amr4$prab <- 1
full_amr5 <- full_amr4 %>% separate_rows(c(AMR_sub_class), sep = "\\;.*") %>% distinct()
full_amr6 <- subset(full_amr5, !(full_amr5$AMR_sub_class == "-" | full_amr5$AMR_sub_class==""))
full_amr7 <- full_amr6 %>% arrange(factor(species, levels = c("caccae", "hominis", "hadrus")))
## part of final figure
ggplot(full_amr7, aes(x=genome_name.y, y = AMR_category, color = factor(prab))) + geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle =90), panel.grid = element_blank()) +
  scale_color_manual(values = c("#96296E")) +
  scale_x_discrete(limits = unique(full_amr7$genome_name.y))
## saved as pdf 7"x8"

## sporulation and germination
anaero126_pan <- read_tsv(file = "/Data/figure7/AN8ANPAN_gene_clusters_summary.txt")
anaero126_pan2 <- anaero126_pan[,-27]
anaero126_pan2$pangenome_type[anaero126_pan2$num_genomes_gene_cluster_has_hits == 1] <- "singletons"
anaero126_pan2$pangenome_type[anaero126_pan2$num_genomes_gene_cluster_has_hits ==126 ] <- "core"
anaero126_pan2$pangenome_type[is.na(anaero126_pan2$pangenome_type) ] <- "soft-shell"
anaero126_pan3 <- anaero126_pan2 %>% separate_rows(c(COG20_CATEGORY_ACC, COG20_CATEGORY), sep = "\\|.*") %>% distinct()

anaero_meta <- read_tsv(file = "/Data/figure7/anaero_meta_final.txt") 
rows_with_spore <- grep("spore", anaero126_pan3$COG20_FUNCTION, ignore.case = TRUE)
spore_A <- anaero126_pan3[rows_with_spore, ]
rows_with_sporulation <- grep("sporulat", anaero126_pan3$COG20_FUNCTION, ignore.case = TRUE)
spore_B <- anaero126_pan3[rows_with_sporulation, ]
spore <- rbind(spore_A, spore_B)
spore2 <- spore %>% select(genome_name, COG20_FUNCTION) %>% distinct()
spo <- spore2 %>% group_by(COG20_FUNCTION) %>% tally()
write_tsv(spo, "/Data/figure7/spo.txt")
## added a new column to define "cat" sporulation and germination and genes I dont want to focus on
spo <- read_tsv("/Data/figure7/spo.txt")
spore2$prab <- 1
spore3 <- spore2 %>% pivot_wider(names_from = "genome_name", values_from = "prab", values_fill = 0)
spore4 <- spore3 %>% pivot_longer(cols = !COG20_FUNCTION, names_to = "genome_name", values_to = "presence_absence")
spore5 <- inner_join(spore4, anaero_meta, by= c("genome_name"="genome_name2")) %>% 
  inner_join(., spo)
spore6 <- subset(spore5, spore5$species %in% c("caccae", "hominis", "hadrus"))
spore7 <- subset(spore6, !(spore6$cat == "rem"))
spore8 <- spore7 %>% arrange(factor(species, levels = c("caccae", "hominis", "hadrus"))) %>% arrange(cat)
## 2nd part of final figure
ggplot(spore8, aes(x= genome_name.y, y= COG20_FUNCTION, color= as.factor(presence_absence))) + geom_point() +
  theme_bw() +
  theme(text = element_text(size = 6), axis.text.x = element_text(angle = 90), panel.grid = element_blank()) +
  scale_x_discrete(limits = unique(spore8$genome_name.y)) +
  scale_y_discrete(limits = unique(spore8$COG20_FUNCTION)) +
  scale_color_manual(values = c("white", "#96296E"))
## saved as pdf
## combined the two final figure parts in Adobe illustrator

#################################################################

# Figure 7B: figures from phase contrast in Adobe Illustrator


################################################################

# Figure 7C:

## tree specific libraries
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)

## input
tree_gerAB <- read.raxml("/Data/figure7/RAxML_bipartitionsBranchLabels.germination")
anaero_meta <- read_tsv(file = "/Data/figure7/anaero_meta_final.txt") 
gerAB <- read_tsv("/Data/figure7/gerAB_anaero126.tsv") ## email corresponding author top find out hpow this file was created
gerAB_id <- select(gerAB, qseqid, genome_name) %>% distinct()
gerAB_meta <- inner_join(gerAB_id, anaero_meta, by = c("genome_name" = "list"))

gerAB_meta$cluster_sp_color[gerAB_meta$cluster_sp=="hadrus"] <- "#0000ff"
gerAB_meta$cluster_sp_color[gerAB_meta$cluster_sp=="caccae"] <- "#ffff00"
gerAB_meta$cluster_sp_color[gerAB_meta$cluster_sp=="chicken"] <- "#00ff00"

gerAB_meta_cluster <- gerAB_meta[,c(1,13,14)]
colnames(gerAB_meta_cluster) <- c("ID", "cluster", "cluster_color")
gerAB_meta_cluster$value <- 1
gerAB_meta_sp <- gerAB_meta[,c(1,3,4)]
colnames(gerAB_meta_sp) <- c("ID", "species", "species_color")
gerAB_meta_sp$value <- 1
## final figure
ggtree(tree_gerAB, layout = "fan", open.angle = 10) + 
  geom_tree() +
  new_scale_fill() +
  geom_fruit(data=gerAB_meta_sp, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=species),
             pwidth = 0.25, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = gerAB_meta_sp$species, values = gerAB_meta_sp$species_color) +
  new_scale_fill() +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=5.5),
        legend.spacing.y = unit(0.02, "cm")) +
  geom_fruit(data=gerAB_meta_cluster, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=cluster),
             pwidth = 0.25, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = gerAB_meta_cluster$cluster, values = gerAB_meta_cluster$cluster_color) +
  new_scale_fill()
## saved as pdf

