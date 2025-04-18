## Supplementary figure 4
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

# Figure S4A
## the pangenome is from Anvio
had_pan <- read_tsv(file = "/Data/sfigure4/AN8HADPAN_gene_clusters_summary.txt")
had_pan2 <- had_pan[,-27]
had_pan2$pangenome_type[had_pan2$num_genomes_gene_cluster_has_hits == 1] <- "singletons"
had_pan2$pangenome_type[had_pan2$num_genomes_gene_cluster_has_hits ==84 ] <- "core"
had_pan2$pangenome_type[is.na(had_pan2$pangenome_type) ] <- "soft-shell"
had_pan3 <- had_pan2 %>% separate_rows(c(COG20_CATEGORY_ACC, COG20_CATEGORY), sep = "\\|.*") %>% distinct()
summ_anvio <- had_pan3 %>% group_by(pangenome_type, COG20_CATEGORY_ACC, COG20_CATEGORY) %>% tally()
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
## saved as pdf

############################################

# Figure S4B
## pangenome from anvio
cac_pan <- read_tsv(file = "/Data/sfigure4/AN8CACPAN_gene_clusters_summary.txt")
cac_pan2 <- cac_pan[,-27]
cac_pan2$pangenome_type[cac_pan2$num_genomes_gene_cluster_has_hits == 1] <- "singletons"
cac_pan2$pangenome_type[cac_pan2$num_genomes_gene_cluster_has_hits == 15] <- "core"
cac_pan2$pangenome_type[is.na(cac_pan2$pangenome_type) ] <- "soft-shell"
cac_pan3 <- cac_pan2 %>% separate_rows(c(COG20_CATEGORY_ACC, COG20_CATEGORY), sep = "\\|.*") %>% distinct()
summ_anvio <- cac_pan3 %>% group_by(pangenome_type, COG20_CATEGORY_ACC, COG20_CATEGORY) %>% tally()
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
fin_sum1_cac <- fin_sum %>% arrange(factor(newcog, levels = ord))
## final figure
ggplot(fin_sum1_cac, aes(x=factor(pangenome_type, levels = c('core', 'soft-shell', 'singletons')), y=perc, fill =factor(newcog, levels = ord))) + geom_col(width = 0.4) +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", size=0.5), text = element_text(size = 5), legend.key.size = unit(5, "mm"), axis.line = element_blank()) +
  scale_fill_manual(breaks = fin_sum1$newcog, values = fin_sum1$color)
## saved as pdf

##########################################################

# Figure S4C
## pangenome from anvio

hom_pan <- read_tsv(file = "/Data/sfigure4/AN8HOMPAN_gene_clusters_summary.txt")
hom_pan2 <- hom_pan[,-27]

hom_pan2$pangenome_type[hom_pan2$num_genomes_gene_cluster_has_hits == 1] <- "singletons"
hom_pan2$pangenome_type[hom_pan2$num_genomes_gene_cluster_has_hits == 11] <- "core"
hom_pan2$pangenome_type[is.na(hom_pan2$pangenome_type) ] <- "soft-shell"

hom_pan3 <- hom_pan2 %>% separate_rows(c(COG20_CATEGORY_ACC, COG20_CATEGORY), sep = "\\|.*") %>% distinct()
summ_anvio <- hom_pan3 %>% group_by(pangenome_type, COG20_CATEGORY_ACC, COG20_CATEGORY) %>% tally()
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
fin_sum1_hom <- fin_sum %>% arrange(factor(newcog, levels = ord))
## final figure
ggplot(fin_sum1_hom, aes(x=factor(pangenome_type, levels = c('core', 'soft-shell', 'singletons')), y=perc, fill =factor(newcog, levels = ord))) + geom_col(width = 0.4) +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", size=0.5), text = element_text(size = 5), legend.key.size = unit(5, "mm"), axis.line = element_blank()) +
  scale_fill_manual(breaks = fin_sum1$newcog, values = fin_sum1$color)
## saved as pdf


