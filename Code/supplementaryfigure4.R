## Supplementary figure 3
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

# Figure S3A
## specific library
library(micropan)

mydata = read_tsv("/Data/sfigure3/number_of_genes_in_pan_genome.txt")
data2 <- pivot_longer(mydata, cols = everything(), names_to = "no_of_genomes", values_to = "no_of_genes")
data2$no_of_genomes <- as.numeric(data2$no_of_genomes)
### calculating stats
model <- lm(log(data2$no_of_genes) ~ log(data2$no_of_genomes))
summary(model)

# Call:
#   lm(formula = log(data2$no_of_genes) ~ log(data2$no_of_genomes))
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.57576 -0.03266  0.00488  0.05374  0.24759 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              8.281835   0.011757   704.4   <2e-16 ***
#   log(data2$no_of_genomes) 0.429857   0.002958   145.3   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09811 on 1258 degrees of freedom
# Multiple R-squared:  0.9438,	Adjusted R-squared:  0.9437 
# F-statistic: 2.111e+04 on 1 and 1258 DF,  p-value: < 2.2e-16

gene1 <- read.table("/Data/sfigure3/gene_presence_absence.Rtab")
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
# 2599.3157645    0.6954943 

data3 <- data2 %>% group_by(no_of_genomes) %>% summarise_at(vars(no_of_genes), list(mean = mean, sd = sd)) %>% as.data.frame()
data3$mino3 <- data3$mean - data3$sd
data3$maxo3 <- data3$mean + data3$sd
## final figure
ggplot(data3, aes(x=(no_of_genomes), y=mean)) + 
  geom_ribbon(aes(ymin = (mean - sd), ymax = (mean + sd)), alpha = .5, fill = "#E8A419", color = "transparent") +
  geom_line(size = 0.5, color = "#E8A419") +
  #geom_boxplot(outlier.shape = NA, fill= "#E8A419") +
  #geom_point() +
  # geom_text(x = 600, y = 1, label = power_eqn(DD), parse = TRUE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=0.3), text = element_text(size = 5), axis.line = element_blank()) +
  scale_x_continuous(breaks = seq(1, length(data3$no_of_genomes), by =10)) +
  geom_text(x = 20, y= 28000, label= "N=3951.4386n^0.4298", size = 2) +
  geom_text(x= 30, y = 30000, label = "alpha = 0.695 (Heap's law), p-value < 2e-16", size = 2) +
  ggtitle("Total number of genes vs number of genomes")
## save as pdf

## Figure S3B
new_genes = read.table(file = "/Data/sfigure3/number_of_new_genes.Rtab")
colnames(new_genes) <- 1:126
data_new <- as.data.frame(colMeans(new_genes))
data_new$no_of_genomes <- 1:126
colnames(data_new) <- c('avg_no_of_new_genes', 'no_of_genomes')

full_new <- inner_join(data3[,c(1:2)], data_new, by = "no_of_genomes")
full_new2 <- full_new %>% pivot_longer(cols = -no_of_genomes, names_to = "type", values_to = "no_genes")

# Fig S3C
gene_presence_absence <- read_csv("/Data/sfigure3/gene_presence_absence.csv")

colnames(gene_presence_absence)[colnames(gene_presence_absence) == "No. isolates"] <- "no_isolates" 
gene_presence_absence$`Avg sequences per isolate`[gene_presence_absence$`Avg sequences per isolate` == 1.00] <- 1
gene_presence_absence1 <- subset(gene_presence_absence, gene_presence_absence$`Avg sequences per isolate` < 1.02)

gene_count <- gene_presence_absence1 %>% group_by(no_isolates) %>% tally()
## final figure
ggplot(gene_count, aes(x=no_isolates, y=n)) + geom_col(fill = "goldenrod") + theme_classic() +
  scale_x_continuous(breaks = seq(1, 126, by=10)) +
  #scale_x_discrete(breaks = seq(1, length(gene_count$no_isolates), by =10)) +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=0.5), text = element_text(size = 8), axis.line = element_blank()) +
  ggtitle("Gene frequency vs number of genomes") +
  scale_y_continuous(breaks = seq(0,12000, 1000)) +
  xlab("No. of genomes") +
  ylab("No. of genes") +
  geom_text(x = 126, y= 4000, label = "core", size =3)
## saved as pdf

