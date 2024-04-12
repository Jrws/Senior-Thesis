### Plot normalized read counts for a gene of interest between tissues and conditions
##### Set working directory to the folder containing this file
setwd("C:/Users/jonny/OneDrive/Desktop/Code/Research/Matthews Transcriptomics")

#BiocManager::install("ggplot2")
#BiocManager::install("tidyverse")
library(ggplot2)
library(tidyverse)

# Read keys and normalized read counts for all tissues
a_sampleTable <- read.csv(file="Antennae/AntennaeKeys.csv", header=TRUE)
b_sampleTable <- read.csv(file="Brain/BrainKeys.csv", header=TRUE)
h_sampleTable <- read.csv(file="Hindlegs/HindlegsKeys.csv", header=TRUE)
r_sampleTable <- read.csv(file="Rostrums/RostrumsKeys.csv", header=TRUE)

antennae_gene_list <- read.csv("Antennae/Antennae_Deseq2_normalized.csv",header = TRUE, row.names=1,  fileEncoding = "UTF-8-BOM")
brain_gene_list <- read.csv("Brain/Brain_Deseq2_normalized.csv",header = TRUE, row.names=1,  fileEncoding = "UTF-8-BOM")
hindlegs_gene_list <- read.csv("Hindlegs/Hindlegs_Deseq2_normalized.csv",header = TRUE, row.names=1,  fileEncoding = "UTF-8-BOM")
rostrums_gene_list <- read.csv("Rostrums/Rostrums_Deseq2_normalized.csv",header = TRUE, row.names=1,  fileEncoding = "UTF-8-BOM")

######### Input gene of interest here! #########
GeneOfInterest <- "LOC5571255" # Potassium channel

# Plot only data for gene of interest
antennae_gene_plot<-antennae_gene_list[GeneOfInterest,]
brain_gene_plot<-brain_gene_list[GeneOfInterest,]
hindlegs_gene_plot<-hindlegs_gene_list[GeneOfInterest,]
rostrums_gene_plot<-rostrums_gene_list[GeneOfInterest,]

a_T_gene_plot<-data.frame(t(antennae_gene_plot))
b_T_gene_plot<-data.frame(t(brain_gene_plot))
h_T_gene_plot<-data.frame(t(hindlegs_gene_plot))
r_T_gene_plot<-data.frame(t(rostrums_gene_plot))

# Rename conditions
a_Condition<-factor(a_sampleTable$Condition)
b_Condition<-factor(b_sampleTable$Condition)
h_Condition<-factor(h_sampleTable$Condition)
r_Condition<-factor(r_sampleTable$Condition)
levels(a_Condition) <- c("Bloodfed","Non-bloodfed")
levels(b_Condition) <- c("Bloodfed","Non-bloodfed")
levels(h_Condition) <- c("Bloodfed","Non-bloodfed")
levels(r_Condition) <- c("Bloodfed","Non-bloodfed")

a_T_gene_plot <- cbind(a_T_gene_plot, a_Condition)
b_T_gene_plot <- cbind(b_T_gene_plot, b_Condition)
h_T_gene_plot <- cbind(h_T_gene_plot, h_Condition)
r_T_gene_plot <- cbind(r_T_gene_plot, r_Condition)

# Calculate average and standard deviation for each tissue + condititon
a_avg <- do.call(data.frame, aggregate(a_T_gene_plot[,GeneOfInterest], list(a_T_gene_plot$a_Condition), FUN=function(x) c(avg = mean(x), std = sd(x))))
b_avg <- do.call(data.frame, aggregate(b_T_gene_plot[,GeneOfInterest], list(b_T_gene_plot$b_Condition), FUN=function(x) c(avg = mean(x), std = sd(x))))
h_avg <- do.call(data.frame, aggregate(h_T_gene_plot[,GeneOfInterest], list(h_T_gene_plot$h_Condition), FUN=function(x) c(avg = mean(x), std = sd(x))))
r_avg <- do.call(data.frame, aggregate(r_T_gene_plot[,GeneOfInterest], list(r_T_gene_plot$r_Condition), FUN=function(x) c(avg = mean(x), std = sd(x))))

colnames(a_avg) = c('Condition', 'Average', 'StdDev')
colnames(b_avg) = c('Condition', 'Average', 'StdDev')
colnames(h_avg) = c('Condition', 'Average', 'StdDev')
colnames(r_avg) = c('Condition', 'Average', 'StdDev')

a_avg = mutate(a_avg, Tissue="Antennae")
b_avg = mutate(b_avg, Tissue="Brain")
h_avg = mutate(h_avg, Tissue="Hindlegs")
r_avg = mutate(r_avg, Tissue="Rostrums")

# Merge averages + standard deviations into summary table
avg_read_by_tissue = list(a_avg, b_avg, h_avg, r_avg) %>% reduce(full_join)
avg_read_by_tissue

# Plot average normalized read counts for gene of interest by tissue for each bloodfeed condition
p1<-ggplot(avg_read_by_tissue, aes(fill=Tissue, x=Condition, y=Average)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=Average-StdDev, ymax=Average+StdDev), position=position_dodge(0.9), alpha=0.8, width=0.4) +
  labs(title=paste("Average expression levels of", "Ir21a"),x="Condition", y = "Normalized Read Counts") + theme(plot.title = element_text(hjust = 0.5))
p1

# Plot average normalized read counts for gene of interest by bloodfeed condition for each tissue
p2<-ggplot(avg_read_by_tissue, aes(fill=Condition, x=Tissue, y=Average)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=Average-StdDev, ymax=Average+StdDev), position=position_dodge(0.9), alpha=0.8, width=0.4) +
  labs(title=paste("Average expression levels of", "Ir21a"),x="Tissue", y = "Normalized Read Counts") + theme(plot.title = element_text(hjust = 0.5))
p2
