# Figure 2

# Authors: 
  # -Figure 2A, 2B, 2C, 2D: MS
  # -Figure 2E: DB

## Input files:
  ## custom_taxonomy.tax
  ## custom_fasta.fasta
  ## anprevclost.filtered.unique.precluster.denovo.vsearch.pick.custom_taxonomy.wang.tax.summary
  ## 16S_metadata.xlsx
  ## 
  ## 

## Output files:
  ## prevelance_final.txt
  ## Figure 2D pdf


###################################

# Figure 1: Made in adobe illustrator


##################################

# Figure 2A:

library(stringr)
library(stringi)
library(readxl)
library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(plyr)
library(janitor)

# files needed:  ----------------------------------------------------------

CGM_Master <- read_excel("CGM_Master.csv", col_names = T) 


# Fig. 2A - families per sample - organize the data -----------------------------------------
detach(package:plyr) #group_by() or count() won't work

# select only the columns you need
totals <- CGM_Master %>%
  dplyr::select(SampleID, FecalSampleID, IsolateNumber, Isolated_On, Isolation_colony_morphology, 
                Phylum, Class, Order, Family, Genus, Species) 

# subset the data by fecal sample
CM01_isolates <- subset(totals, FecalSampleID == "CM001")
CM02_isolates <- subset(totals, FecalSampleID == "CM002")
CM03_isolates <- subset(totals, FecalSampleID == "CM003")
CM06_isolates <- subset(totals, FecalSampleID == "CM006")
CM07_isolates <- subset(totals, FecalSampleID == "CM007")
CM08_isolates <- subset(totals, FecalSampleID == "CM008")
CM09_isolates <- subset(totals, FecalSampleID == "CM009")
CM10_isolates <- subset(totals, FecalSampleID == "CM010")
CM11_isolates <- subset(totals, FecalSampleID == "CM011")

## count the number of families isolated per sample
cm01_fam <- CM01_isolates %>%
  group_by(FecalSampleID, Phylum, Family) %>% 
  count(Family)
cm02_fam <- CM02_isolates %>%
  group_by(FecalSampleID, Phylum, Family) %>% 
  count(Family)
cm03_fam <- CM03_isolates %>%
  group_by(FecalSampleID, Phylum, Family) %>% 
  count(Family)
cm06_fam <- CM06_isolates %>%
  group_by(FecalSampleID, Phylum, Family) %>% 
  count(Family)
cm07_fam <- CM07_isolates %>%
  group_by(FecalSampleID, Phylum, Family) %>% 
  count(Family)
cm08_fam <- CM08_isolates %>%
  group_by(FecalSampleID, Phylum, Family) %>% 
  count(Family)
cm09_fam <- CM09_isolates %>%
  group_by(FecalSampleID, Phylum, Family) %>% 
  count(Family)
cm10_fam <- CM10_isolates %>%
  group_by(FecalSampleID, Phylum, Family) %>% 
  count(Family)
cm11_fam <- CM11_isolates %>%
  group_by(FecalSampleID, Phylum, Family) %>% 
  count(Family)

# join the sample dataframes together and clean up NAs
fam_counts <- rbind(cm01_fam, cm02_fam, cm03_fam, cm06_fam, cm07_fam,
                    cm08_fam, cm09_fam, cm10_fam, cm11_fam)

fam_counts_fixed <- fam_counts %>% 
  mutate(Family = if_else(str_detect(Family, "N/A"),"NA", Family)) %>%
  mutate(Phylum = if_else(str_detect(Phylum, "N/A"),"NA", Phylum)) 


# Fig. 2A - set up the coloring scheme ----------------------------------------------

# create a new dataframe with just the Family and Phylum
df <- fam_counts_fixed %>%
  select(Family, Phylum)

# what phyla do you have? 
summary(as.factor(fam_counts_fixed$Phylum))

# ordering the phyla by least abundant to most abundant
fam_counts_fixed$Phylum <- factor(fam_counts_fixed$Phylum, levels = c( "NA", "Gemmatimonadetes","Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes"))
df <- fam_counts_fixed[order(fam_counts_fixed$Phylum),]

# coloring at the family level
# the first line removes any duplicated lines, so that each family only appears once
fam_color <- df[!duplicated(df$Family), ]
table(fam_color$Phylum)

## setting the colors for each phylum
# the numbers = how many unique families w/in the phyla 

actino <- colorRampPalette(c("tan", "brown", "orangered4"))(n=3)

pro <- colorRampPalette(c("gold", "orange2"))(n=2)

firm <- colorRampPalette(c("midnightblue","dodgerblue4","deepskyblue4",
                           "skyblue3","skyblue", "lightskyblue3","deepskyblue1","dodgerblue1","steelblue1","royalblue1",
                           "#4911EA","#4911AA","#7911AA","#8911AA", "#7939A0"))(n=15)

bac <- colorRampPalette(c("#6AB100", "#3AB100","#9EFD80"))(n=3)

gemma <- c("aquamarine2")

N_A <- c("grey") 

Color <- c(N_A, gemma, pro, bac, actino, firm)

# assigns the set colors to each row, based on the phylum
fam_color$fam_color <- Color

# setting the Phylum column as a factor, with a specific set of variables 
fam_color$Phylum <- factor(fam_color$Phylum, levels = c("NA", "Gemmatimonadetes","Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes"))

# orders the families in the specified phyla order
fam_color <- fam_color[order(fam_color$Phylum),]

# now at the phylum level, removes duplicated phlya
phyl_color <- df[!duplicated(df$Phylum), ] 

# what phyla do you have?
phyl_color <- phyl_color %>% select(Phylum) 
table(phyl_color$Phylum)

# one color assigned for each unique phylum
act<-c("brown")

pr<-c("orange2")

fir<-c("lightskyblue3")

ba<-c("forestgreen")

gemma <- c("aquamarine2")

na <- c("grey")

Color2 <- c(na, gemma, pr, ba, act, fir)

phyl_color$phyl_color <- Color2

# joins the assigned phyla and family colors into 1 dataframe to use for graphing 
tax_colors <- inner_join(phyl_color, fam_color, by=c("Phylum"))

# order the families for graphing
fam_ordered <- fam_counts_fixed %>%
  mutate(Family = factor(Family, levels = fam_color$Family)) 

fam_ordered$Family <- factor(fam_ordered$Family) 

# to add the total number of isolates per sample ID
total_n <- fam_ordered %>%
  group_by(FecalSampleID) %>%
  summarize(total_n = sum(n))

# join the total numbers to the ordered dataframe
fam_ordered <- fam_ordered %>%
  left_join(total_n, by = "FecalSampleID")


# Fig. 2A - graphing -------------------------------------------------------

# graphing the number of families per fecal sample
fam_ordered %>% 
  ggplot(aes(x = FecalSampleID, y = n, fill = factor(Family), group = factor(Phylum, levels = c("NA", "Gemmatimonadetes","Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes")))) +
  geom_col() +
  labs(title = "Frequency of Families Isolated per Sample", 
       x = "Sample ID", y = "Frequency", 
       fill = "Families") +
  theme_classic() +
  scale_x_discrete(labels = function(x) {
    paste(x, "\nn =", total_n$total_n[match(x, total_n$FecalSampleID)]) 
  }) +   # Add total n values to the x-axis labels
  scale_fill_manual(breaks = tax_colors$Family.y, values = tax_colors$fam_color) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")


#################################

# Figure 2B:

## filter the master sheet by successful sanger IDs
# gives total number of isolates ID'd by sanger sequencing 
sanger_IDs <- isolates %>%
  filter(SangerSeq_Results %in% c("success", "failed multispecies", "multispecies", "failed, success on repeat")) %>%
  dplyr::select(SampleID, FecalSampleID, IsolateNumber, Isolated_On, Isolation_colony_morphology, 
                Phylum, Class, Order, Family, Genus, Species) #cleaning up the df a bit


# make a dataframe with the info (column names) you need
cgm <- select(sanger_IDs, SampleID, FecalSampleID, IsolateNumber, Isolated_On, Phylum, Class, Order, Family, Genus, Species)

# unite the genus and species columns into one 
cgm2<- unite(data = cgm, col = "gen_sp", sep = "_", Genus, Species)

# fix the NAs
cgm3<- subset(cgm2, !(cgm2$gen_sp == "NA_NA"))

#count the number of times a species is isolated
species001 <- cgm3 %>% filter(FecalSampleID %in% c("CM001")) %>%
  group_by(FecalSampleID, Phylum, Family, gen_sp) %>% 
  count(gen_sp)
species002 <- cgm3 %>% filter(FecalSampleID %in% c("CM002")) %>%
  group_by(FecalSampleID, Phylum, Family, gen_sp) %>% 
  count(gen_sp)
species003 <- cgm3 %>% filter(FecalSampleID %in% c("CM003")) %>%
  group_by(FecalSampleID, Phylum, Family, gen_sp) %>% 
  count(gen_sp)
species006 <- cgm3 %>% filter(FecalSampleID %in% c("CM006")) %>%
  group_by(FecalSampleID, Phylum, Family, gen_sp) %>% 
  count(gen_sp)
species007 <- cgm3 %>% filter(FecalSampleID %in% c("CM007")) %>%
  group_by(FecalSampleID, Phylum, Family, gen_sp) %>% 
  count(gen_sp)
species008 <- cgm3 %>% filter(FecalSampleID %in% c("CM008")) %>%
  group_by(FecalSampleID, Phylum, Family, gen_sp) %>% 
  count(gen_sp)
species009 <- cgm3 %>% filter(FecalSampleID %in% c("CM009")) %>%
  group_by(FecalSampleID, Phylum, Family, gen_sp) %>% 
  count(gen_sp)
species010 <- cgm3 %>% filter(FecalSampleID %in% c("CM010")) %>%
  group_by(FecalSampleID, Phylum, Family, gen_sp) %>% 
  count(gen_sp)
species011 <- cgm3 %>% filter(FecalSampleID %in% c("CM011")) %>%
  group_by(FecalSampleID, Phylum, Family, gen_sp) %>% 
  count(gen_sp)

# group all of the fecal samples together
unique_counts <- rbind(species001, species002, species003, species006, species007,
                       species008, species009, species010, species011)


# calculate the total number of isolates per sample ID
total_n <- unique_counts %>%
  group_by(FecalSampleID) %>%
  summarise(total_n = sum(n, na.rm = TRUE))

# join the unique species counts with the total count dataframe and calculate percent
counts <- unique_counts %>%
  left_join(total_n, by = "FecalSampleID")%>%
  mutate(percent_unique_species = (n / total_n) * 100)


# Fig. 2B - set up the coloring scheme ----------------------------

# note - want to color by family (ex. all firmicutes are blue)

# make a new dataframe 
df <- counts %>%
  select(Family, Phylum)

# ordering the phyla
counts$Phylum <- factor(counts$Phylum, levels = c( "NA", "Gemmatimonadetes","Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes"))
df <- counts[order(counts$Phylum),]


# coloring at the family level
# the first line removes any duplicated lines, so that each family only appears once
fam_color_unique <- df[!duplicated(df$gen_sp), ]
table(fam_color_unique$Phylum)

## setting the colors for the phyla
# the numbers = how many unique families w/in the phyla 

actino_u <- colorRampPalette(c("tan", "salmon4", "brown", "orangered4", "darkred"))(n=6)

pro_u <- colorRampPalette(c("yellow","yellow2","gold", "darkgoldenrod3","goldenrod2","orange2"))(n=7)

firm_u <- colorRampPalette(c("midnightblue","royalblue4","dodgerblue4","deepskyblue4","steelblue4","#1211AC","#0000CD",
                             "skyblue3","#1291AA","skyblue", "lightskyblue3","deepskyblue1","dodgerblue1","steelblue1","royalblue1",
                             "#4911EA","#4911AA","#7911AA","#8911AA", "#7939A0","plum4"))(n=127)

bac_u <- colorRampPalette(c("#1A6700","#1A4700","#6AB100", "#9AB100","#3AB100","#9EFD80","#7EFD00"))(n=14)

gemma_u <- c("aquamarine2")

N_A_u <- c("grey") 

Color <- c(N_A_u, gemma_u, pro_u, bac_u, actino_u, firm_u)

#assigns the set colors to each row, based on the phylum
fam_color_unique$fam_color_unique <- Color

#setting the Phylum column as a factor, with a specific set of variables (the list)
fam_color_unique$Phylum <- factor(fam_color_unique$Phylum, levels = c("NA", "Gemmatimonadetes","Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes"))

#orders the families in the specified phyla order
fam_color_unique <- fam_color_unique[order(fam_color_unique$Phylum),]

# now at the phylum level, removes duplicated phlya
phyl_color_unique <- df[!duplicated(df$Phylum), ] 

phyl_color_unique <- phyl_color_unique %>% select(Phylum) 
table(phyl_color_unique$Phylum)

# one color assigned for each unique phylum
act<-c("brown")

pr<-c("orange2")

fir<-c("lightskyblue3")

ba<-c("forestgreen")

gemma <- c("aquamarine2")

na <- c("grey")

Color2 <- c(na, gemma, pr, ba, act, fir)

phyl_color_unique$phyl_color_unique <- Color2

# joins the assigned phyla and family colors into 1 dataframe to use for graphing 
tax_colors_unique <- inner_join(phyl_color_unique, fam_color_unique, by=c("Phylum"))

##order the families 
fam_ordered_unique <- counts %>%
  mutate(gen_sp = factor(gen_sp, levels = fam_color_unique$gen_sp)) 

fam_ordered_unique$gen_sp <- factor(fam_ordered_unique$gen_sp) 


# Fig. 2B - graphing -------------------------------------------------------

# gives the distribution of unique species across the samples
# with coloring scheme - by unique gen_sp
fam_ordered_unique %>% 
  ggplot(aes(x = FecalSampleID, y = percent_unique_species, fill = factor(gen_sp), group = factor(Phylum, levels = c("NA", "Gemmatimonadetes","Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes")))) +
  geom_col() +
  labs(title = "Unique Species Distribution Across Samples", 
       x = "Sample ID", y = "Species Abundance (%)", 
       fill = "Genus_Species") +
  theme_classic() +
  scale_fill_manual(breaks = tax_colors_unique$gen_sp.y, values = tax_colors_unique$fam_color_unique) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


####################

# Figure 2C

## filter the master sheet by successful sanger IDs
# gives total # of isolates ID'd by sanger sequencing 
sanger_IDs <- isolates %>%
  filter(SangerSeq_Results %in% c("success", "failed multispecies", "multispecies", "failed, success on repeat")) %>%
  dplyr::select(SampleID, FecalSampleID, IsolateNumber, Isolated_On, Isolation_colony_morphology, 
                Phylum, Class, Order, Family, Genus, Species) #cleaning up the df a bit


# subset the samples by sample ID
CM01<- subset(sanger_IDs, FecalSampleID == "CM001")
CM02<- subset(sanger_IDs, FecalSampleID == "CM002")
CM03<- subset(sanger_IDs, FecalSampleID == "CM003")
CM06<- subset(sanger_IDs, FecalSampleID == "CM006")
CM07<- subset(sanger_IDs, FecalSampleID == "CM007")
CM08<- subset(sanger_IDs, FecalSampleID== "CM008")
CM09<- subset(sanger_IDs, FecalSampleID == "CM009")
CM10<- subset(sanger_IDs, FecalSampleID == "CM010")
CM11<- subset(sanger_IDs, FecalSampleID == "CM011")

# count the number of times each media type is used to initially isolate a colony
cm01_media <- CM01 %>%
  group_by(FecalSampleID, Isolated_On) %>%
  count(Isolated_On)
cm02_media <- CM02 %>%
  group_by(FecalSampleID, Isolated_On) %>%
  count(Isolated_On)
cm03_media <- CM03 %>%
  group_by(FecalSampleID, Isolated_On) %>%
  count(Isolated_On)
cm06_media <- CM06 %>%
  group_by(FecalSampleID, Isolated_On) %>%
  count(Isolated_On)
cm07_media <- CM07 %>%
  group_by(FecalSampleID, Isolated_On) %>%
  count(Isolated_On)
cm08_media <- CM08 %>%
  group_by(FecalSampleID, Isolated_On) %>%
  count(Isolated_On)
cm09_media <- CM09 %>%
  group_by(FecalSampleID, Isolated_On) %>%
  count(Isolated_On)
cm10_media <- CM10 %>%
  group_by(FecalSampleID, Isolated_On) %>%
  count(Isolated_On)
cm11_media <- CM11 %>%
  group_by(FecalSampleID, Isolated_On) %>%
  count(Isolated_On)

# group all of the samples together
media_counts <- rbind(cm01_media, cm02_media, cm03_media, cm06_media, cm07_media,
                      cm08_media, cm09_media, cm10_media, cm11_media) 


# Fig. 2C - graphing ----------------------------------------------------------------

# to add the total number of isolates per sample ID
total_n_media <- media_counts %>%
  group_by(FecalSampleID) %>%
  summarize(total_n = sum(n))

# join the isolate totals with the media totals
media_counts <- media_counts %>%
  left_join(total_n_media, by = "FecalSampleID")

# color palette used for graphing
media <- brewer.pal(6, "Dark2")


media_counts %>% 
  ggplot(aes(x = FecalSampleID, y = n, fill = Isolated_On)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Frequency of Isolates per Media Type", 
       x = "Sample ID", y = "Number of Identified Isolates") +
  theme_classic() +
  scale_x_discrete(labels = function(x) {
    paste(x, "\nn =", total_n_media$total_n[match(x, total_n_media$FecalSampleID)]) 
  }) +   # Add total n values to the x-axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_fill_manual(values = media)




################################
# Figure 2D:

## HPCC: mothur to create a custom classifier and applying it on the datasets

## Datasets

## download/cp all datasets
### select refers to those samples selected for healthy adult humans/controls
afolayan_select
bian_select      
daquinan
dwiyanto
kameoka_select
macpherson_select
sprockett_select
tang_select

### Creating a custom classifier
### Input files are custom_taxonomy and custom_fasta; make sure they are in the correct format
## Used multiple sequences for same species to give better confidence
## first create custom_taxonomy.txt; CM07_96 was dropped because 16S was only 878 nt long
## This is custom_taxonomy:
```
CM07_153_Anaerobutyricum_hallii	Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Anaerobutyricum;Anaerobutyricum_hallii;
CM07_176_Anaerobutyricum_hallii	Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Anaerobutyricum;Anaerobutyricum_hallii;
CM07_23_Anaerobutyricum_hallii	Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Anaerobutyricum;Anaerobutyricum_hallii;
etc...
```

## create the custom classifier using mothur:

align.seqs(fasta=custom_fasta.fasta, reference=silva.seed_v138_2.align, processors=16)
summary.seqs(fasta=custom_fasta.align)
# let's also align to v4 region
pcr.seqs(fasta=silva.seed_v138_2.align, start=13862, end=23444, keepdots=F, processors=16)
system(mv silva.seed_v138_2.pcr.align silva.v4.fasta)
summary.seqs(fasta=silva.v4.fasta, processors=16)
# then, align the silva v4 to our alignment
align.seqs(fasta=custom_fasta.fasta, reference=silva.v4.fasta, processors=16)
summary.seqs(fasta=custom_fasta.align)
# checked classification of selected sequences, these two not required 
unique.seqs(fasta=custom_fasta.align, output=names)
classify.seqs(fasta=custom_fasta.unique.align, count=custom_fasta.count_table, reference=trainset18_062020.rdp.fasta, taxonomy=trainset18_062020.rdp.tax, cutoff=80)
## I had 17 unique seqs and all the things are called at the correct tax in rdp

## Prep the datasets; Did it in batch to save time

##slurm script:
```
#!/bin/sh
#SBATCH --job-name=16S_prev
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=750G
#SBATCH --time=24:00:00
#SBATCH --constraint=interconnect_hdr

cd /pathtodatasets/16S_prev/datasets/ 
  
module load mothur/1.48.0

mothur mbatch.txt
```
## mbatch.txt:
```
make.contigs(file=anprevclost.files, processors=16)
summary.seqs(fasta=anprevclost.trim.contigs.fasta, processors=16)
screen.seqs(fasta=anprevclost.trim.contigs.fasta, count=anprevclost.contigs.count_table, maxambig=0, maxlength=275, processors=16)
unique.seqs(fasta=anprevclost.trim.contigs.good.fasta, count=anprevclost.contigs.good.count_table, processors=16)
summary.seqs(fasta=anprevclost.trim.contigs.good.unique.fasta, count=anprevclost.trim.contigs.good.count_table, processors=16)
align.seqs(fasta=anprevclost.trim.contigs.good.unique.fasta, reference=silva.seed_v138_2.align, processors=16)
summary.seqs(fasta=anprevclost.trim.contigs.good.unique.align, count=anprevclost.trim.contigs.good.count_table, processors=16)
screen.seqs(fasta=anprevclost.trim.contigs.good.unique.align, count=anprevclost.trim.contigs.good.count_table, summary=anprevclost.trim.contigs.good.unique.summary, start=13862, end=23444, maxhomop=8, processors=16)
summary.seqs(fasta=anprevclost.trim.contigs.good.unique.good.align, count=anprevclost.trim.contigs.good.good.count_table, processors=16)
filter.seqs(fasta=anprevclost.trim.contigs.good.unique.good.align, vertical=T, trump=., processors=16)
system(mv anprevclost.trim.contigs.good.unique.good.filter.fasta anprevclost.filtered.fasta)
system(mv anprevclost.trim.contigs.good.good.count_table anprevclost.filtered.count_table)
unique.seqs(fasta=anprevclost.filtered.fasta, count=anprevclost.filtered.count_table, processors=16)
pre.cluster(fasta=anprevclost.filtered.unique.fasta, count=anprevclost.filtered.unique.count_table, diffs=1, processors=4)
chimera.vsearch(fasta=anprevclost.filtered.unique.precluster.fasta, count=anprevclost.filtered.unique.precluster.count_table, dereplicate=t, processors=4)
summary.seqs(fasta=current, count=current, processors=16)
classify.seqs(fasta=anprevclost.filtered.unique.precluster.denovo.vsearch.fasta, count=anprevclost.filtered.unique.precluster.denovo.vsearch.count_table, reference=trainset18_062020.rdp.fasta, taxonomy=trainset18_062020.rdp.tax, cutoff=80)
remove.lineage(fasta=anprevclost.filtered.unique.precluster.denovo.vsearch.fasta, count=anprevclost.filtered.unique.precluster.denovo.vsearch.count_table, taxonomy=anprevclost.filtered.unique.precluster.denovo.vsearch.rdp.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
```

# then, align .fasta file of 16S rRNA sample files to custom classifier files:
# using interactive I did the followiung:
classify.seqs(fasta=anprevclost.filtered.unique.precluster.denovo.vsearch.pick.fasta, count=anprevclost.filtered.unique.precluster.denovo.vsearch.pick.count_table, reference=custom_fasta.align , taxonomy=custom_taxonomy.tax, cutoff=80)
## this where I get the final file for RStudio

# Rstudio

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

## input files

tax_sum <- read_tsv(file="/Data/figure2/anprevclost.filtered.unique.precluster.denovo.vsearch.pick.custom_taxonomy.wang.tax.summary")

meta_16S <- read_excel(path="/Data/figure2/16S_metadata.xlsx")
a <- colnames(tax_sum[,c(6:550)])
meta2 <- subset(meta_16S, meta_16S$SRA %in% a)
tax2 <- subset(tax_sum, tax_sum$taxlevel==7)
tax3 <- tax2 %>% select(3,5)
tax5 <- tax2[,c(3,6:550)]
tax5_2 <- tax5 %>% pivot_longer(cols = !taxon, names_to = "SRA", values_to = "relab") %>% filter(relab>0)
## filter SRA greater than 1500
tax6 <- tax5 %>% column_to_rownames("taxon") %>% as.matrix() 
tax6_2 <- as.data.frame(colSums(tax6)) %>% rownames_to_column("SRA")
tax6_2_2 <- subset(tax6_2, tax6_2$colsums >1500)

tax5_2 <- subset(tax5_2, tax5_2$SRA %in% c(tax6_2_2$SRA))
tax5_2$prev <- 1
tax5_3 <- tax5_2 %>% group_by(taxon) %>% dplyr::summarise(sum_prev=sum(prev))
tax5_3$prev_percent <- (tax5_3$sum_prev/418)*100
prevalent <- subset(tax5_3, !(tax5_3$taxon %in% c("unknown_unclassified", "Dorea_longicatena", "Lachnospiraceae_unclassified")))
write_tsv(prevalent, "/Data/figure2/prevelance_final.txt")

tax6 <- tax5 %>% column_to_rownames("taxon") %>% as.matrix() 
tax6_2 <- as.data.frame(colSums(tax6)) %>% rownames_to_column("SRA") 
colnames(tax6_2) <- c("SRA", "colsums")
tax6_2_2 <- subset(tax6_2, tax6_2$colsums >1500)
tax6_3 <- inner_join(tax5_2, tax6_2_2, by ="SRA")
tax6_4 <- tax6_3
tax6_4$relab_per <- (tax6_4$relab/tax6_4$colsums)*100
tax6_4$logra <- log10(tax6_4$relab_per)
tax6_5 <- inner_join(tax6_4, meta2, by="SRA")
tax6_6 <- subset(tax6_5, !(tax6_5$taxon %in% c("unknown_unclassified", "Dorea_longicatena", "Lachnospiraceae_unclassified")))
## final figure
ggplot(tax6_6, aes(y=taxon, x=logra, alpha = 0.5)) + geom_point(position = position_jitter(width = 0.7), color = "grey67") + geom_boxplot(outlier.shape = NA) +
  theme_bw()
## save as pdf











