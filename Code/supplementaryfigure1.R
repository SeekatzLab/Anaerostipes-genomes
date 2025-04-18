### Supplementary Figure 1
## Author: Meagan Seesengood
# 3.7.2025


# libraries needed: -------------------------------------------------------

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


# files needed: -----------------------------------------------------------

CGM_Master <- read_excel("CGM_Master.csv", col_names = T) 


# organize the data -------------------------------------------------------

detach(package:plyr)

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

# count the number of times a species is isolated
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


# group all of the samples together
gen_counts <- rbind(species001, species002, species003, species006, species007,
                    species008, species009, species010, species011)


# want to add the total number of colonies picked per sample to the x-axis 
colonies_picked <- isolates %>%
  select(FecalSampleID, Genus, Species) %>%
  group_by(FecalSampleID) %>%
  tally() %>%
  slice(-10) # remove the NAs

# tally the number of unique species per fecal sample
total_unique <- gen_counts %>%
  group_by(FecalSampleID) %>%
  summarise(total_n = sum(n, na.rm = TRUE))

# join unique species count with total colonies picked
total_unique <- total_unique %>%
  left_join(colonies_picked, by = "FecalSampleID")

# rename columns
total_unique <- total_unique %>%
  rename(number_colonies_picked = n, unique_sp_n = total_n) 


# graphing ----------------------------------------------------------------

total_unique %>%
  ggplot(aes(x = FecalSampleID, y = unique_sp_n)) +
  geom_col(fill = "dodgerblue4") +
  labs(title = "Unique Species Picked", 
       x = "Colonies Picked", y = "Total Number of Unique Species") +
  scale_x_discrete(labels = function(x) {
    paste(x, "\nn =", total_unique$number_colonies_picked[match(x, total_unique$FecalSampleID)]) 
  }) +   # Add total n values to the x-axis labels
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")
