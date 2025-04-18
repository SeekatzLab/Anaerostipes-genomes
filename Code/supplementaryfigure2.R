# Supplementary Figure 1:

# Authors
  ## Sfigure 1: DB

# Input files
  ## RAxML_bipartitionsBranchLabels.core_snp_final
  ## anaero_meta_final.txt
  ## RAxML_bestTree.final126_refined.tre


# Output files
  ## supplementary figure 1A pdf
  ## supplementary figure 1B pdf

## Library

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

#####################

## Supplementary Figure 1A:

## more libraries for trees:
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)

tree_anaero <- read.raxml("/Data/supplementaryfigure1/RAxML_bipartitionsBranchLabels.core_snp_final")
anaero_meta <- read_tsv(file = "/Data/supplementaryfigure1/anaero_meta_final.txt") 
## data coloring
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
ggtree(tree_anaero, layout = "rectangular", yscale_mapping = data_species$speceies) + 
  geom_tree() +
  geom_text(aes(label = bootstrap), hjust = 1, vjust = -0.4, size = 3) +
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

#################

# Supplementary figuer 1B: phylophlan input

tree_anaero <- read.tree("/Data/supplementaeryfigure1/RAxML_bestTree.final126_refined.tre")
anaero_meta <- read_tsv(file = "/Data/supplementaryfigure1/anaero_meta_final.txt") 
## data coloring
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
ggtree(tree_anaero, layout = "fan", yscale_mapping = data_species$species) + 
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






