# Figure 5

## Input files
  ###anaero_meta_final.txt -- contains genome metadata
  ###AN8ANPAN_gene_clusters_summary.txt -- from Anvio'8 saved from anvi-display-pan
  ###cob.txt -- added a new column in the cob df to assign pathways to cobalamin
  ###acetate_anaero126.tsv
  ###butyrate_anaero126.tsv
  ###propionate_anaero126.tsv
  ###bsh_anaero126_2.tsv
  ###RAxML_bipartitionsBranchLabels.butbutrax
  ###RAxML_bipartitionsBranchLabels.bukbutrax
  ###RAxML_bipartitionsBranchLabels.bshrax

## Output files
  ### figure 5A pdf
  ### figure 5B pdf
  ### figure 5C pdf
  ### figure 5D pdf
  ### figure 5E pdf

### Libraries
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

## Figure 3A: Cobalamin synthesis across clusters
### Input
anaero_meta <- read_tsv("/Data/figure5/anaero_meta_final.txt")
anaero126_pan <- read_tsv("/Data/figure5/AN8ANPAN_gene_clusters_summary.txt")

na_static <- anaero126_pan[is.na(anaero126_pan$bin_name),]
## should be 0 obs
anaero126_pan2 <- anaero126_pan[,-27]
anaero126_pan2$pangenome_type[anaero126_pan2$num_genomes_gene_cluster_has_hits == 1] <- "singletons"
anaero126_pan2$pangenome_type[anaero126_pan2$num_genomes_gene_cluster_has_hits ==126 ] <- "core"
anaero126_pan2$pangenome_type[is.na(anaero126_pan2$pangenome_type) ] <- "soft-shell"
anaero126_pan3 <- anaero126_pan2 %>% separate_rows(c(COG20_FUNCTION_ACC, COG20_FUNCTION), sep = "\\|.*") %>% distinct()
### collect cobalamin related genes
rows_with_cob <- grep("cob", anaero126_pan3$COG20_FUNCTION, ignore.case = TRUE)
cob_A <- anaero126_pan3[rows_with_cob, ]
rows_with_b12 <- grep("B12", anaero126_pan3$COG20_PATHWAY, ignore.case = TRUE)
cob_B <- anaero126_pan3[rows_with_b12, ]
cob <- rbind(cob_A, cob_B) ### final cobalamin collection
### create presence absence column
cob2 <- cob %>% select(genome_name, COG20_FUNCTION) %>% distinct()
cob2$prab <- 1
cob3 <- cob2 %>% pivot_wider(names_from = "genome_name", values_from = "prab", values_fill = 0)
cob4 <- cob3 %>% pivot_longer(cols = !COG20_FUNCTION, names_to = "genome_name", values_to = "presence_absence")
cob5 <- inner_join(cob4, anaero_meta, by= c("genome_name"="genome_name2"))
cob6 <- cob5 %>% arrange(factor(cluster_paper, levels = c("A", "B", "C")))
cob7 <- cob6 %>% group_by(cluster_paper, species, COG20_FUNCTION) %>% dplyr::summarise(fraction=mean(presence_absence))
### order for species in accordance to cluster_paper
sp_ord <- c("amylophilus", "hadrus","sp018918155","sp900066705", "caccae", "hominis", "rhamnosivorans", "sp018381315", "avistercoris", "butyraticus", "excrementavium",  "faecalis")

cobal <- cob2 %>% group_by(COG20_FUNCTION) %>% tally()
cobal$cob_path <- "a"
cobal$cob_path[grep("Cob", cob$COG20_FUNCTION, ignore.case = FALSE)] <- "aerobic"
### write out cob.txt so you can add cob_path classifcaiotn of pathways (aerobic, anaerobic etc) in Excel
write_tsv(cob, "/Data/figure4/cob.txt")
cobalamin2 <- read_tsv("/Data/figure5/cob.txt")
cobalamin3 <- subset(cobalamin2, !(cobalamin2$cob_path=="rem")) ## I wanted to discuss only a few genes
cobalamin4 <- inner_join(cobalamin3, cob7, by = "COG20_FUNCTION")

cobalamin5 <- cobalamin4 %>% arrange(cob_path)
### final figure
ggplot(cobalamin5, aes(x= factor(species, levels = sp_ord), y= COG20_FUNCTION, fill= fraction)) + geom_tile(linetype = 1, color="black") +
  theme_bw() +
  theme(text = element_text(size = 9), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = unique(cobalamin5$COG20_FUNCTION)) +
  scale_fill_gradient2(low = "white", high = "#0c5650", mid = "green", midpoint = 0.5) + coord_fixed()
## save as pdf

########

## Figure 5B: SCFA genes across different species by clusters

### how to get all the diamond.tsv together to create single dataframes
# curated_list_anaero <- read_tsv("/anaero_final_list.txt") ---> list of genome names matching the file names in *_diamond.tsv
# file_list <- as.vector(curated_list_anaero$list)
# ### create an empty list which will be populated with all the overview.txt
# list.data.ace<-list()
# filepath_scfa = "/scfa/"
# ace <- as.vector(c("eutD_ace", "tdcD_ace"))
# list.data.ace <- vector("list", length(ace))
# 
# for(x in seq_along(ace)){
#   list.data.ace[[x]] <- vector("list", length(file_list))
#   for(i in seq_along(file_list)){
#     #print(i)
#     file_path <- file.path(filepath_scfa, ace[x], paste0(file_list[i], '_', ace[x], '_diamond.tsv'))
#     list.data.ace[[x]][[i]] <- read_tsv(file_path, col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
#   }
# }
# 
# 
# acetate <- data.frame()
# 
# 
# for (x in seq_along(ace)) {
#   for (i in seq_along(file_list)) {
#     # Add genome_name and ace_name columns to the current data frame
#     list.data.ace[[x]][[i]] <- mutate(list.data.ace[[x]][[i]], 
#                                       genome_name = file_list[i], 
#                                       ace_name = ace[x])
#     
#     # Combine the current data frame with all_sugar
#     acetate <- rbind.fill(acetate, list.data.ace[[x]][[i]])
#   }
# }
# 
# names(acetate)[names(acetate)=="ace_name"] <- "gene_name"
# acetate$comp <- "acetate"
# write_tsv(acetate, file = "/Data/figure5/acetate_anaero126.tsv")
# 
# but <- as.vector(c("abfD_but", "ato_but", "buk_but", "but_but", "cro_but", "gcd_but", "kal_but"))
# list.data.but<-list()
# list.data.but <- vector("list", length(but))
# 
# for(x in seq_along(but)){
#   list.data.but[[x]] <- vector("list", length(file_list))
#   for(i in seq_along(file_list)){
#     #print(i)
#     file_path.but <- file.path(filepath_scfa, but[x], paste0(file_list[i], '_', but[x], '_diamond.tsv'))
#     list.data.but[[x]][[i]] <- read_tsv(file_path.but, col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
#   }
# }
# 
# 
# 
# butyrate <- data.frame()
# 
# 
# for (x in seq_along(but)) {
#   for (i in seq_along(file_list)) {
#     # Add genome_name and ace_name columns to the current data frame
#     list.data.but[[x]][[i]] <- mutate(list.data.but[[x]][[i]], 
#                                       genome_name = file_list[i], 
#                                       gene_name = but[x])
#     
#     # Combine the current data frame with all_sugar
#     butyrate <- rbind.fill(butyrate, list.data.but[[x]][[i]])
#   }
# }
# butyrate$comp <- "butyrate"
# write_tsv(butyrate, file = "/Data/figure5/butyrate_anaero126.tsv")
# 
# 
# pro <- as.vector(c("epimce_prop","lcda_prop","mmdA_prop","mut_prop","pduC_prop","pduP_prop","ygfH_prop"))
# list.data.pro<-list()
# list.data.pro <- vector("list", length(pro))
# 
# for(x in seq_along(pro)){
#   list.data.pro[[x]] <- vector("list", length(file_list))
#   for(i in seq_along(file_list)){
#     #print(i)
#     file_path.pro <- file.path(filepath_scfa, pro[x], paste0(file_list[i], '_', pro[x], '_diamond.tsv'))
#     list.data.pro[[x]][[i]] <- read_tsv(file_path.pro, col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
#   }
# }
# 
# 
# propionate <- data.frame()
# 
# 
# for (x in seq_along(pro)) {
#   for (i in seq_along(file_list)) {
#     # Add genome_name and ace_name columns to the current data frame
#     list.data.pro[[x]][[i]] <- mutate(list.data.pro[[x]][[i]], 
#                                       genome_name = file_list[i], 
#                                       gene_name = pro[x])
#     
#     # Combine the current data frame with all_sugar
#     propionate <- rbind.fill(propionate, list.data.pro[[x]][[i]])
#   }
# }
# propionate$comp <- "propionate"
# write_tsv(propionate, file = "/Data/figure5/propionate_anaero126.tsv")
# 
# file_list <- anaero_meta$genome_name2
# #file_list <- as.vector(curated_list_anaero$list)
# ### create an empty list which will be populated with all the overview.txt
# list.data.bsh<-list()
# ### file path to overview.txt from run_dbcan
# filepath_bsh = "~/Library/CloudStorage/Box-Box/Disha (dbhatta@clemson.edu)/Projects/paper4_anaero/diamond_new_bsh/"
# 
# for(i in (1:length(file_list))){
#   #print(i)
#   list.data.bsh[[i]] <- read_tsv(paste0(filepath_bsh, file_list[i], '_diamond.tsv'), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
# }
# 
# ### combining all the overview.txt from a list to one large data frame and simultaneously adding their species name
# 
# bsh <- list.data.bsh[[1]]
# bsh$genome_name <- "CM01_24_S247"
# 
# for (i in 2:126) {
#   bsh <- rbind.fill(bsh, list.data.bsh[[i]])
#   bsh <- mutate(bsh, genome_name = replace(genome_name, is.na(genome_name), file_list[[i]]))
# }
# 
# bsh$gene_name <- "bsh"
# bsh$comp <- "bsh"
# write_tsv(bsh, file = "/Data/figure5/bsh_anaero126_2.tsv")

### Analysis starts here
### Input files from Data 
acetate <- read_tsv("/Data/figure5/acetate_anaero126.tsv")
butyrate <- read_tsv("/Data/figure5/butyrate_anaero126.tsv")
propionate <- read_tsv("/Data/figure5/propionate_anaero126.tsv")
bsh <- read_tsv("/Data/figure5/bsh_anaero126_2.tsv")

## combine
full_scfa <- rbind(acetate, butyrate, propionate, bsh)

full_scfa2_1 <- inner_join(full_scfa[!(full_scfa$gene_name=="bsh"),], anaero_meta, by = c("genome_name" = "list")) %>% select(-genome_name2)
full_scfa2_2 <- inner_join(full_scfa[(full_scfa$gene_name=="bsh"),], anaero_meta, by = c("genome_name" = "genome_name2")) %>% select(-list)
full_scfa2 <- rbind(full_scfa2_1, full_scfa2_2)
full_scfa3 <- full_scfa2 %>% group_by(genome_name.y, gene_name, comp, species) %>% tally()
## presence or absence is denoted by 1 or 0 respectively
full_scfa3$prab <- 1
full_scfa4 <- full_scfa3[,c(1,2,6)] %>% pivot_wider(names_from = "gene_name", values_from = "prab", values_fill = 0)
## adding back genes that didnot have any diamond output
full_scfa4$ygfH_prop <- 0
full_scfa4$pduC_prop <- 0
full_scfa4$pduP_prop <- 0

full_scfa5 <- full_scfa4 %>% pivot_longer(cols = !genome_name.y, values_to = "prab", names_to = "gene_name")
full_scfa6 <- inner_join(full_scfa5, anaero_meta, by = c("genome_name.y" = "genome_name"))

full_scfa7 <- full_scfa6 %>% group_by(species, gene_name) %>% dplyr::summarise(num = sum(prab))
## just checking
numbers_sp <- anaero_meta %>% group_by(species) %>% tally()

full_scfa8 <- inner_join(full_scfa7, numbers_sp, by = "species")
full_scfa8$fraction <- full_scfa8$num/full_scfa8$n

full_scfa8$comp[full_scfa8$gene_name == "bsh"] <- "bsh"
full_scfa8$comp[full_scfa8$gene_name %in% c("eutD_ace", "tdcD_ace")] <- "acetate"
full_scfa8$comp[full_scfa8$gene_name %in% c("abfD_but", "ato_but", "buk_but", "but_but", "cro_but", "gcd_but", "kal_but")] <- "butyrate"
full_scfa8$comp[full_scfa8$gene_name %in% c("epimce_prop","lcda_prop","mmdA_prop","mut_prop","pduC_prop","pduP_prop","ygfH_prop")] <- "propionate"

ord_comp <- c("acetate", "butyrate", "propionate", "bsh")
full_scfa8 <- full_scfa8 %>% arrange(factor(comp, levels = ord_comp))
sp_ord <- c("amylophilus","hadrus","sp018918155","sp900066705", "caccae", "hominis", "rhamnosivorans", "sp018381315", "avistercoris", "butyraticus", "excrementavium", "faecalis")
## final figure
ggplot(full_scfa8, aes(x= factor(species, levels = sp_ord), y=gene_name, fill=fraction)) +geom_tile(color="black", linetype=1) + 
  coord_fixed() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = unique(full_scfa8$gene_name)) +
  #facet_wrap(~comp, ncol = 1, scales = "free_y") +
  scale_fill_gradient2(high = "#0c5650", low = "white", mid = "green", midpoint = 0.5)
## save as pdf

########

## Figure 3C: Butyrate genes but across all species tree

### libraries 
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)

### tree with characteristics input files
tree_but <- read.raxml("/Data/figure5/RAxML_bipartitionsBranchLabels.butbutrax")
anaero_meta <- read_tsv(file = "/Data/figure5/anaero_meta_final.txt") 
butyrate <- read_tsv("/Data/figure5/butyrate_anaero126.tsv")
### making the df for the boxes around the tree
but_id <- butyrate[butyrate$gene_name=="but_but",] %>% select(qseqid, genome_name) %>% distinct()
but_meta <- inner_join(but_id, anaero_meta, by = c("genome_name" = "list"))

but_meta$cluster_sp_color[but_meta$cluster_sp=="hadrus"] <- "#0000ff"
but_meta$cluster_sp_color[but_meta$cluster_sp=="caccae"] <- "#ffff00"
but_meta$cluster_sp_color[but_meta$cluster_sp=="chicken"] <- "#00ff00"

but_meta_cluster <- but_meta[,c(1,13,14)]
colnames(but_meta_cluster) <- c("ID", "cluster", "cluster_color")
but_meta_cluster$value <- 1

but_meta_sp <- but_meta[,c(1,3,4)]
colnames(but_meta_sp) <- c("ID", "species", "species_color")
but_meta_sp$value <- 1
### tree
ggtree(tree_but, layout = "fan", open.angle = 10) + 
  geom_tree() +
  new_scale_fill() +
  geom_fruit(data=but_meta_sp, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=species),
             pwidth = 0.25, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = but_meta_sp$species, values = but_meta_sp$species_color) +
  new_scale_fill() +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=5.5),
        legend.spacing.y = unit(0.02, "cm")) +
  geom_fruit(data=but_meta_cluster, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=cluster),
             pwidth = 0.25, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = but_meta_cluster$cluster, values = but_meta_cluster$cluster_color) +
  new_scale_fill()
## save as pdf

##############

## Figure 3D: Butyrate genes buk across all species tree

### libraries 
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)

### tree with characteristics input files
tree_buk <- read.raxml("/Data/figure5/RAxML_bipartitionsBranchLabels.bukbutrax")
anaero_meta <- read_tsv(file = "/Data/figure5/anaero_meta_final.txt") 
butyrate <- read_tsv("/Data/figure5/butyrate_anaero126.tsv")
### making the df for the boxes around the tree
buk_id <- butyrate[butyrate$gene_name=="buk_but",] %>% select(qseqid, genome_name) %>% distinct()
buk_meta <- inner_join(buk_id, anaero_meta, by = c("genome_name" = "list"))

buk_meta$cluster_sp_color[buk_meta$cluster_sp=="hadrus"] <- "#0000ff"
buk_meta$cluster_sp_color[buk_meta$cluster_sp=="caccae"] <- "#ffff00"
buk_meta$cluster_sp_color[buk_meta$cluster_sp=="chicken"] <- "#00ff00"

buk_meta_cluster <- buk_meta[,c(1,13,14)]
colnames(buk_meta_cluster) <- c("ID", "cluster", "cluster_color")
buk_meta_cluster$value <- 1
buk_meta_sp <- buk_meta[,c(1,3,4)]
colnames(buk_meta_sp) <- c("ID", "species", "species_color")
buk_meta_sp$value <- 1
## tree
ggtree(tree_buk, layout = "fan", open.angle = 10) + 
  geom_tree() +
  new_scale_fill() +
  geom_fruit(data=buk_meta_sp, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=species),
             pwidth = 0.25, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = buk_meta_sp$species, values = buk_meta_sp$species_color) +
  new_scale_fill() +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=5.5),
        legend.spacing.y = unit(0.02, "cm")) +
  geom_fruit(data=buk_meta_cluster, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=cluster),
             pwidth = 0.25, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = buk_meta_cluster$cluster, values = buk_meta_cluster$cluster_color) +
  new_scale_fill()
## save as pdf

##############################

## Figure 3E: BSH genes across all genomes

tree_bsh <- read.raxml("/Data/figure5/RAxML_bipartitionsBranchLabels.bshrax")
anaero_meta <- read_tsv(file = "/Data/figure5/anaero_meta_final.txt") 
bsh <- read_tsv("/Data/figure5/bsh_anaero126_2.tsv")
bsh_id <- bsh %>% select(qseqid, genome_name) %>% distinct()
bsh_meta <- inner_join(bsh_id, anaero_meta, by = c("genome_name" = "list"))
bsh_meta$cluster_sp_color[bsh_meta$cluster_paper=="A"] <- "#3a53a4"
bsh_meta$cluster_sp_color[bsh_meta$cluster_paper=="B"] <- "#f6eb16"
bsh_meta$cluster_sp_color[bsh_meta$cluster_paper=="C"] <- "#6abd45"
bsh_meta_cluster <- bsh_meta[,c(1,14,15)]
colnames(bsh_meta_cluster) <- c("ID", "cluster", "cluster_color")
bsh_meta_cluster$value <- 1
bsh_meta_sp <- bsh_meta[,c(1,3,4)]
colnames(bsh_meta_sp) <- c("ID", "species", "species_color")
bsh_meta_sp$value <- 1
## final figure
ggtree(tree_bsh, layout = "fan", open.angle = 10) + 
  geom_tree() +
  new_scale_fill() +
  geom_fruit(data=bsh_meta_sp, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=species),
             pwidth = 0.25, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = bsh_meta_sp$species, values = bsh_meta_sp$species_color) +
  new_scale_fill() +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=5.5),
        legend.spacing.y = unit(0.02, "cm")) +
  geom_fruit(data=bsh_meta_cluster, geom=geom_col,
             mapping=aes(y=ID, x=value, fill=cluster),
             pwidth = 0.25, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = bsh_meta_cluster$cluster, values = bsh_meta_cluster$cluster_color) +
  new_scale_fill() #+
geom_fruit(data=data_sample, geom=geom_col,
           mapping=aes(y=ID, x=value, fill=sample),
           pwidth = 0.25, offset = 0.1, grid.params=list(size=0.0)) +
  scale_fill_manual(breaks = data_sample$sample, values = data_sample$color_sample)
## save as pdf


