## Supplementary figure 2
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

## Figure S2A

## tree specific libraries
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
#library(ggplot2)
library(ggnewscale)

## input:
tree_anaero <- read.raxml("/Data/sfigure2/RAxML_bipartitionsBranchLabels.core_snp_final")
anaero_meta <- read_tsv(file = "/Data/sfigure2/anaero_meta_final.txt") 
## creating color dfs
data_species <- anaero_meta[,c(11,2,3)]
colnames(data_species) <- c("ID", "species", "color_species")
data_species$value <- 1
## there is nothing worng with the origin color, I just want it to match the original figure
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
## figure
ggtree(tree_anaero, layout = "rectangular", yscale_mapping = data_species$speceies) + 
  geom_tree() +
  geom_text(aes(label = bootstrap), hjust = 1, vjust = -0.4, size = 3, color="red") +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.5) +
  #geom_tiplab(colour = 'red', size = 2) +
  #geom_hilight(node=200, fill="springgreen", alpha=.3) +
  new_scale_fill() +
  #coord_fixed() +
  geom_fruit(data=data_species, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=species),
             pwidth = 0.1, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = data_species$species, values = data_species$color_species) +
  new_scale_fill() +
  #scale_color_manual(values = data_species$color, breaks = data_species$species) +
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

## Figure S2B

### input file
aai_nwk <- read.newick("/Data/sfigure2/aai_anaero.nwk")
anero_meta

data_species <- anaero_meta[,c(1,2,3)]
colnames(data_species) <- c("ID", "species", "color_species")
data_species$value <- 1

data_host <-  anaero_meta[,c(1,4,5)]          #read_tsv("~/Library/CloudStorage/Box-Box/Disha/Projects/paper4_anaero/data_host.txt")
colnames(data_host) <- c("ID", "origin", "color_origin")
data_host$value <- 1
data_host$color_origin[data_host$origin=="human"] <- "#b5c99a"
data_host$color_origin[data_host$origin=="chicken"] <- "#915330"
data_host$color_origin[data_host$origin=="pig"] <- "#c2926b"
data_host$color_origin[data_host$origin=="mice"] <- "#feecad"
data_host$color_origin[data_host$origin=="macaque"] <- "#ffca61"

data_sample <- anaero_meta[,c(1,7,6)]               #read_tsv("~/Library/CloudStorage/Box-Box/Disha/Projects/paper4_anaero/anaero_meta_final.txt")
colnames(data_sample) <- c("ID", "sample", "color_sample")
data_sample$value <- 1
data_sample$color_sample[data_sample$sample == "GC"] <- "gray90"

ggtree(aai_nwk, layout = "fan", yscale_mapping = data_species$species) + 
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

## Figure S2C
##aai box

aai_anaero126 <- read_tsv("/Data/sfigure2/aai_anaero.tsv")
aai_anaero126_2 <- aai_anaero126 %>% select(`Label 1`, `Label 2`, AAI)
anaero_meta <- read_tsv(file = "/Data/sfigure2/anaero_meta_final.txt") 
meta_aai <- inner_join(aai_anaero126_2, anaero_meta, by = c("Label 1" = "list"))
table(anaero_meta$species)
sp <- c("caccae", "hominis", "rhamnosivorans", "sp018381315", "faecalis", "butyraticus", "excrementavium", "avistercoris","sp018918155", "amylophilus", "hadrus","sp900066705")
meta_aai2 <- meta_aai %>% arrange(factor(species, levels = sp))

ggplot(meta_aai2, aes(x=`Label 1`, y=`Label 2`, fill=AAI)) + geom_raster() +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())+
  scale_x_discrete(limits = meta_aai2$`Label 1`) +
  scale_y_discrete(limits = meta_aai2$`Label 1`) +
  scale_fill_gradient2(low="white", high = "red", midpoint = 80)

