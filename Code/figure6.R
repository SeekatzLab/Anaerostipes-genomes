# Figure 6:

# Authors:
# - all figures: DB

## Input files:
  ## gene_presence_absence_caccae.Rtab; number_of_genes_in_pan_genome_caccae.Rtab; number_of_new_genes_caccae.Rtab
  ## gene_presence_absence_hominis.Rtab; number_of_genes_in_pan_genome_hominis.Rtab; number_of_new_genes_hominis.Rtab
  ## gene_presence_absence_hadrus.Rtab; number_of_genes_in_pan_genome_hadrus.Rtab; number_of_new_genes_hadrus.Rtab
  ## anaero_meta_final.txt
  ## funenr_kmod_3sp.txt
  ## sugar_anaero126.tsv
  ## functl_enr_newcluster.txt
  ## cazy_anaero126.tsv

## Output files:
  ## Figure 6A pdf
  ## Figure 6B pdf
  ## Figure 6C pdf

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

# Figure 6A:
## specific libraries
library(micropan)

## there are 3 panels in figure 6A

## panel 1: caccae

gene1 <- read.table("/Data/figure5/gene_presence_absence_caccae.Rtab")
gene2 <- gene1
colnames(gene2) <- gene1[1, ]
rownames(gene2) <- gene1$V1
gene3 <- gene2[-1, ]
gene4 <- gene3[ , -1]
gene5 <- t(gene4)
set.seed(1234)
h.est <- heaps(gene5, n.perm = 500) 
h.est

# Intercept        alpha 
# 660.3467709   0.3891369 

mydata = read_tsv("/Data/figure6/number_of_genes_in_pan_genome_caccae.Rtab", col_names = FALSE)
data2 <- pivot_longer(mydata, cols = everything(), names_to = "no_of_genomes", values_to = "no_of_genes")
data2$no_of_genomes <- gsub("X", "", data2$no_of_genomes)
data2$no_of_genomes <- as.numeric(data2$no_of_genomes)
data3 <- data2 %>% group_by(no_of_genomes) %>% summarise_at(vars(no_of_genes), list(mean = mean, sd = sd)) %>% as.data.frame()
new_genes = read.table(file = "/Data/figure6/number_of_new_genes_caccae.Rtab")
colnames(new_genes) <- 1:15
data_new <- as.data.frame(colMeans(new_genes))
data_new$no_of_genomes <- 1:15
colnames(data_new) <- c('avg_no_of_new_genes', 'no_of_genomes')

full_new <- inner_join(data3[,c(1:2)], data_new, by = "no_of_genomes")
full_new2 <- full_new %>% pivot_longer(cols = -no_of_genomes, names_to = "type", values_to = "no_genes")
## final figure 5A panel 1
ggplot(full_new2, aes(x=no_of_genomes, y= no_genes, color=type)) + geom_point() +
  theme_bw() +
  geom_text(x = 5, y= 3000, label= "N=3005.2619n^0.33545", size = 2) +
  geom_text(x= 6, y = 4000, label = "alpha = 0.38913 (Heap's law), p-value < 2e-16", size = 2) +
  theme(axis.text.x = element_text(angle = 90), panel.grid = element_blank()) +
  scale_color_manual(values = c("goldenrod2", "navy")) +
  scale_x_continuous(breaks = seq(1,126,1)) +
  scale_y_continuous(breaks = seq(10, 32000, 1000))
## save as pdf

## panel 2: hominis

gene1 <- read.table("/Data/figure6/gene_presence_absence_hominis.Rtab")
gene2 <- gene1
colnames(gene2) <- gene1[1, ]
rownames(gene2) <- gene1$V1
gene3 <- gene2[-1, ]
gene4 <- gene3[ , -1]
gene5 <- t(gene4)
set.seed(1234)
h.est <- heaps(gene5, n.perm = 500) 
h.est

# Intercept      alpha 
# 836.307693   1.081467 
mydata = read_tsv("/Data/figure6/number_of_genes_in_pan_genome_hominis.Rtab", col_names = FALSE)
data2 <- pivot_longer(mydata, cols = everything(), names_to = "no_of_genomes", values_to = "no_of_genes")
data2$no_of_genomes <- gsub("X", "", data2$no_of_genomes)
data2$no_of_genomes <- as.numeric(data2$no_of_genomes)
data3 <- data2 %>% group_by(no_of_genomes) %>% summarise_at(vars(no_of_genes), list(mean = mean, sd = sd)) %>% as.data.frame()

new_genes = read.table(file = "/Data/figure6/number_of_new_genes_hominis.Rtab")
colnames(new_genes) <- 1:11
data_new <- as.data.frame(colMeans(new_genes))
data_new$no_of_genomes <- 1:11
colnames(data_new) <- c('avg_no_of_new_genes', 'no_of_genomes')

full_new <- inner_join(data3[,c(1:2)], data_new, by = "no_of_genomes")
full_new2 <- full_new %>% pivot_longer(cols = -no_of_genomes, names_to = "type", values_to = "no_genes")
## final figure 5A panel 2
ggplot(full_new2, aes(x=no_of_genomes, y= no_genes, color=type)) + geom_point() +
  theme_bw() +
  geom_text(x = 5, y= 3000, label= "N=3005.2619n^0.33545", size = 2) +
  geom_text(x= 6, y = 1000, label = "alpha = 1.081467 (Heap's law), p-value < 2e-16", size = 2) +
  theme(axis.text.x = element_text(angle = 90), panel.grid = element_blank()) +
  scale_color_manual(values = c("goldenrod2", "navy")) +
  scale_x_continuous(breaks = seq(1,126,1)) +
  scale_y_continuous(breaks = seq(10, 32000, 1000))
## saved as pdf

## panel 3: hadrus

gene1 <- read.table("/Data/figure6/gene_presence_absence_hadrus.Rtab")
gene2 <- gene1
colnames(gene2) <- gene1[1, ]
rownames(gene2) <- gene1$V1
gene3 <- gene2[-1, ]
gene4 <- gene3[ , -1]
gene5 <- t(gene4)
set.seed(1234)
h.est.hadrus <- heaps(gene5, n.perm = 500) 
h.est.hadrus

# Intercept        alpha 
# 1396.4785112    0.7288171

mydata = read_tsv("/Data/figure6/number_of_genes_in_pan_genome_hadrus.Rtab", col_names = FALSE)
data2 <- pivot_longer(mydata, cols = everything(), names_to = "no_of_genomes", values_to = "no_of_genes")
data2$no_of_genomes <- gsub("X", "", data2$no_of_genomes)
data2$no_of_genomes <- as.numeric(data2$no_of_genomes)
data3 <- data2 %>% group_by(no_of_genomes) %>% summarise_at(vars(no_of_genes), list(mean = mean, sd = sd)) %>% as.data.frame()
new_genes = read.table(file = "/Data/figure6/number_of_new_genes_hadrus.Rtab")
colnames(new_genes) <- 1:84
data_new <- as.data.frame(colMeans(new_genes))
data_new$no_of_genomes <- 1:84
colnames(data_new) <- c('avg_no_of_new_genes', 'no_of_genomes')

full_new <- inner_join(data3[,c(1:2)], data_new, by = "no_of_genomes")
full_new2 <- full_new %>% pivot_longer(cols = -no_of_genomes, names_to = "type", values_to = "no_genes")
## figure 5A panel 3
ggplot(full_new2, aes(x=no_of_genomes, y= no_genes, color=type)) + geom_point() +
  theme_bw() + theme(legend.position = "none") + 
  geom_text(x = 50, y= 6000, label= "N=3166.6451n^0.344588", size = 2) +
  geom_text(x= 60, y = 5000, label = "alpha = 0.7288171 (Heap's law), p-value < 2e-16", size = 2) +
  theme(axis.text.x = element_text(angle = 90), panel.grid = element_blank()) +
  scale_color_manual(values = c("goldenrod2", "navy")) +
  scale_x_continuous(breaks = seq(1,84,6)) +
  scale_y_continuous(breaks = seq(10, 32000, 1000))
## save as pdf

################################################################

# Figure 5B:

kmod_3sp <- read_tsv("/Data/figure6/funenr_kmod_3sp.txt")

kmod_3sp_2 <- subset(kmod_3sp, kmod_3sp$adjusted_q_value < 0.001)

kmod_3sp_3 <- kmod_3sp_2[, c(7,8,10,12)] %>% pivot_longer(cols = !`function`, names_to = "species", values_to = "fraction")
## final figure
ggplot(kmod_3sp_3, aes(x=factor(species, levels = c("p_caccae", "p_hominis", "p_hadrus")), y=`function`, fill=fraction)) + geom_tile(color="black") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
  scale_fill_gradient2(high = "#0c5650", low = "white") +
  coord_fixed()
## saved as pdf

###################################################################

# Figure 6C:

# curated_list_anaero <- read_tsv("/anaero_final_list.txt")
# file_list <- as.vector(curated_list_anaero$list)
# ### create an empty list which will be populated
# list.data.sugar<-list()
# ### file path to dbsub.out from run_dbcan
# filepath_cazy = "/dbcan4/dbcan_"
# 
# for(i in (1:length(file_list))){
#   #print(i)
#   list.data.sugar[[i]] <- read_tsv(paste0(filepath_cazy, file_list[i], '/dbsub.out'))
# }
# 
# all_sugar <- list.data.sugar[[1]]
# all_sugar$genome_name <- "CM01_24_S247" ## whichever one is first on your list 
# 
# for (i in 2:126) {
#   all_sugar <- rbind.fill(all_sugar, list.data.sugar[[i]])
#   all_sugar <- mutate(all_sugar, genome_name = replace(genome_name, is.na(genome_name), file_list[[i]]))
# }

# write_tsv(all_sugar, file = "/Data/figure5/sugar_anaero126.tsv")

all_sugar <- read_tsv(file="/Data/figure6/sugar_anaero126.tsv")
### modifying all_dbcan_cinn; shortening the genome name for easier use later on
all_sugar$newgenome_name <- sub("^([^_]*_[^_]*)_.*$", "\\1", all_sugar$genome_name)
all_sugar$newgenome_name <- gsub("\\..*","",all_sugar$newgenome_name)
### removing brackets
all_sugar2 <- all_sugar[,c(1,2,4,6,14,15)]
## subset the ones with sugars
all_sugar3 <- all_sugar2 %>% separate_longer_delim(Substrate, delim = ", ")
all_sugar4 <- subset(all_sugar3, !(all_sugar3$Substrate == "-"))
all_sugar4$prab <- 1
all_sugar5 <- all_sugar4[,c(3,6,7)] %>% distinct() %>% 
  pivot_wider(names_from = Substrate, values_from = prab, values_fill = 0)
all_sugar6 <- all_sugar5 %>% pivot_longer(cols = !(newgenome_name), names_to = "Substrate", values_to = "prab")
anaero_meta <- read_tsv("/Data/figure6/anaero_final_meta.txt")
merge_sugar <- inner_join(all_sugar6, anaero_meta, by = c('newgenome_name' = 'genome_name')) %>% arrange(species)
merge_sugar2 <- subset(merge_sugar, merge_sugar$species %in% c("hadrus", "caccae", "hominis")) 
merge_sugar2 <- merge_sugar2 %>% arrange(factor(species, levels = c("caccae", "hominis", "hadrus")))
## final figure
ggplot(merge_sugar2, aes(y=newgenome_name, x=Substrate, fill=factor(prab))) + geom_tile(color="black") +
  coord_fixed()+
  theme_bw() +
  theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(limits = unique(merge_sugar2$newgenome_name))+
  scale_fill_manual(values = c("gray90","#96296E"))
## saved as pdf






