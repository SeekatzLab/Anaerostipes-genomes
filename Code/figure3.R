# Figure 3:

# Authors:
# - all figures: DB

## Input files:
  ## RAxML_bipartitionsBranchLabels.core_snp_final
  ## anaero_meta_final.txt
  ## ANIb_percentage_identity.txt
  ## AN8ANPAN_gene_clusters_summary.txt

## Output files:
  ## Figure 3A pdf
  ## Figure 3B pdf
  ## Figure 3D pdf

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


###################################

# Figure 3A
## tree specific libraries
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
## input files
tree_anaero <- read.raxml("/Data/figure3/RAxML_bipartitionsBranchLabels.core_snp_final")
anaero_meta <- read_tsv(file = "/Data/figure3/anaero_meta_final.txt") 
## color dataframes
data_species <- anaero_meta[,c(11,2,3)]
colnames(data_species) <- c("ID", "species", "color_species")
data_species$value <- 1

data_host <-  anaero_meta[,c(11,4,5)]       
colnames(data_host) <- c("ID", "origin", "color_origin")
data_host$value <- 1
data_host$color_origin[data_host$origin=="human"] <- "#b5c99a"
data_host$color_origin[data_host$origin=="chicken"] <- "#915330"
data_host$color_origin[data_host$origin=="pig"] <- "#c2926b"
data_host$color_origin[data_host$origin=="mice"] <- "#feecad"
data_host$color_origin[data_host$origin=="macaque"] <- "#ffca61"

data_sample <- anaero_meta[,c(11,7,6)]        
colnames(data_sample) <- c("ID", "sample", "color_sample")
data_sample$value <- 1
data_sample$color_sample[data_sample$sample == "GC"] <- "gray90"
## final figure
ggtree(tree_anaero, layout = "fan", yscale_mapping = data_species$speceies) + 
  geom_tree() +
  new_scale_fill() +
  geom_fruit(data=data_species, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=species),
             pwidth = 0.1, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = data_species$species, values = data_species$color_species) +
  new_scale_fill() +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=5.5),
        legend.spacing.y = unit(0.02, "cm")) +
  geom_fruit(data=data_host, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=origin),
             pwidth = 0.1, offset = 0, grid.params=list(size=0.0)) +
  scale_fill_manual(name = data_host$origin, breaks = data_host$origin, values = data_host$color_origin) +
  new_scale_fill() +
  geom_fruit(data=data_sample, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=sample),
             pwidth = 0.1, offset = 0, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = data_sample$sample, values = data_sample$color_sample)
## save as pdf

###########################

# Figure 3B:

ani_anaero126 <- read_tsv("/Data/figure3/ANIb_percentage_identity.txt")
ani_anaero126_2 <- ani_anaero126 %>% column_to_rownames("key") %>% as.matrix()
ani_anaero126_3 <- ani_anaero126 %>% pivot_longer(cols = !key, names_to = "key2", values_to = "ani")
anaero_meta <- read_tsv(file = "/Data/figure3/anaero_meta_final.txt") 

meta_ani <- inner_join(ani_anaero126_3, anaero_meta, by = c("key" = "genome_name2"))
table(anaero_meta$species)

sp <- c("caccae", "hominis", "rhamnosivorans", "sp018381315", "avistercoris", "faecalis", "butyraticus", "excrementavium", "amylophilus", "hadrus","sp018918155","sp900066705")
meta_ani2 <- meta_ani %>% arrange(factor(species, levels = sp))
## final figure
ggplot(meta_ani2, aes(x=key, y=key2, fill=ani)) + geom_raster() +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())+
  scale_x_discrete(limits = meta_ani2$key) +
  scale_y_discrete(limits = meta_ani2$key) +
  scale_fill_gradient2(low="white", high = "red", mid = "white", midpoint = 0.8)
## save as pdf


############################################

# Figure 3C: figure from Anvi'o


####################################

# Figure 3D:

anaero126_pan <- read_tsv(file = "/Data/figure3/AN8ANPAN_gene_clusters_summary.txt")
anaero126_pan2 <- anaero126_pan[,-27]

anaero126_pan2$pangenome_type[anaero126_pan2$num_genomes_gene_cluster_has_hits == 1] <- "singletons"
anaero126_pan2$pangenome_type[anaero126_pan2$num_genomes_gene_cluster_has_hits ==126 ] <- "core"
anaero126_pan2$pangenome_type[is.na(anaero126_pan2$pangenome_type) ] <- "soft-shell"

num_genes <- anaero126_pan2 %>% group_by(gene_cluster_id, pangenome_type) %>% tally()
anaero126_pan3 <- anaero126_pan2 %>% separate_rows(c(COG20_CATEGORY_ACC, COG20_CATEGORY), sep = "\\|.*") %>% distinct()

summ_anvio <- anaero126_pan3 %>% group_by(pangenome_type, COG20_CATEGORY_ACC, COG20_CATEGORY) %>% tally()

summ_anvio1 <- summ_anvio
summ_anvio1$COG20_CATEGORY_ACC[is.na(summ_anvio1$COG20_CATEGORY_ACC)] <- "Not_Assigned"
summ_anvio1$COG20_CATEGORY[is.na(summ_anvio1$COG20_CATEGORY)] <- "Not_Assigned"

summ_anvio1 <- subset(summ_anvio1, !(summ_anvio1$COG20_CATEGORY==""))
summ_anvio1$newcog <- summ_anvio1$COG20_CATEGORY_ACC

core_summary <- subset(summ_anvio1, summ_anvio1$pangenome_type == "core")
#core_summary <- core_summary %>% group_by(bin_name, COG20_CATEGORY, newcog) %>% tally()
core_summary$perc <- round((core_summary$n/sum(core_summary$n)) *100, digits = 6)

accsin_summary <- subset(summ_anvio1, summ_anvio1$pangenome_type == "singletons")
#accsin_summary <- accsin_summary %>% group_by(bin_name, newcog) %>% tally()
accsin_summary$perc <- round((accsin_summary$n/sum(accsin_summary$n)) *100, digits = 6)

soft_summary <- subset(summ_anvio1, summ_anvio1$pangenome_type == "soft-shell")
#soft_summary <- soft_summary %>% group_by(bin_name, newcog) %>% tally()
soft_summary$perc <- round((soft_summary$n/sum(soft_summary$n)) *100, digits = 6)

final_summary <- rbind(core_summary, accsin_summary, soft_summary)


color_pal <- as.data.frame(unique(final_summary$newcog))
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
#color_pal$color[color_pal$newcog == ""] <- "black"

fin_sum <- inner_join(final_summary, color_pal, by = "newcog")
ord <- c("Not_Assigned", "R", "S", "W", "J", "K", "L", "D", "V", "T", "M", "N", "O", "U", "C", "G", "E", "F", "H", "I", "P", "Q", "X", "")
fin_sum1 <- fin_sum %>% arrange(factor(newcog, levels = ord))
## final figure
ggplot(fin_sum1, aes(x=factor(pangenome_type, levels = c('core', 'soft-shell', 'singletons')), y=perc, fill =factor(newcog, levels = ord))) + geom_col(width = 0.4) +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", size=0.5), text = element_text(size = 5), legend.key.size = unit(5, "mm"), axis.line = element_blank()) +
  scale_fill_manual(breaks = fin_sum1$newcog, values = fin_sum1$color)
## save as pdf





