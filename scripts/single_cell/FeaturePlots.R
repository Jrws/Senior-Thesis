suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(sctransform)
  library(ComplexHeatmap)
  library(future)
  library(scales)
  library(data.table)
  library(ggplot2)
})

theme_set(theme_cowplot())

projectFolder="C:/Users/jonny/OneDrive/Desktop/Code/Research/scRNASeq"
setwd(projectFolder)

#Load one of the following datasets
load(".RData")
#load("AntennaAllCells.RData")
#load("MaxPalpNeurons.RData") 
#load("MaxPalpAllCells.RData")

#Name dataset of interest "dataset"
dataset <- AntennaNeurons
#dataset <- AntennaAllCells
#dataset <- MaxPalpNeurons
#dataset <- MaxPalpAllCells

############################# Pseudo bulk (may be required for code further on)

counts <- dataset@assays$RNA@counts                  #Pull raw counts from Seurat object
pseudo_bulk <- rowSums(counts)                                             #Sum values across cells

pseudo_bulk <- as.data.frame(pseudo_bulk)                                  #Convert to dataframe
pseudo_bulk <- tibble::rownames_to_column(pseudo_bulk, "Gene")             #Pull out rownames as column

pseudo_bulk <- pseudo_bulk[!grepl("LOC*", pseudo_bulk$Gene),]              #get rid of LOC & MC genes (for just chemoreceptor analysis)                     #
pseudo_bulk <- pseudo_bulk[!grepl("MT*", pseudo_bulk$Gene),]                    
pseudo_bulk <- pseudo_bulk[order(pseudo_bulk[,2],decreasing=TRUE),]        #Sort genes in decreasing order 

pseudo_bulk.ORs <- pseudo_bulk[grepl("Or*", pseudo_bulk$Gene),]            #Separate ORs, IRs, GRs into three variables
pseudo_bulk.IRs <- pseudo_bulk[grepl("Ir*", pseudo_bulk$Gene),]
pseudo_bulk.GRs <- pseudo_bulk[grepl("Gr*", pseudo_bulk$Gene),]

rm(counts)

################################################################### Feature Plots (Figure S9H, Figure S13F)


geneName1 <- toString("LOC110678282") #repo, glial marker
name1 <- toString("repo (Glial)")

geneName2 <- toString("LOC5564305") #grh, epithelial marker
name2 <- toString("grh (Epithelial)")
#Neural markers
geneName3 <- toString("LOC5565901") #syt1, neuronal marker
name3 <- toString("syt1 (Neuronal)")

geneName4 <- toString("LOC5575210") #nompC, mechanosensory channel
name4 <- toString("nompC (Mechanosensory)")

Dsc1_gene <- toString("LOC5575962") # Dsc1
Dsc1_name <- toString("DSC1")

Rdl_gene <- toString("LOC5570466")  # Rdl
Rdl_name <- toString("Rdl")

Orco_gene <- toString("Orco")       # Orco
Orco_name <- toString("Orco")

Para_gene <- toString("LOC5567355") # Para
Para_name <- toString("Para")

Potassium_channel <- toString("LOC5568798")
K_chan_name <- toString("Shaw")

Calcium_channel <- toString("LOC5564339")  #Cacophony ortholog Ca2+ channel
Ca_chan_name <- toString("LOC5564339")

# geneName5 <- toString("LOC5564848") #CadN
# geneName6 <- toString("LOC5570381") #brp
# geneName7 <- toString("LOC5570204") #elav

features <- c(Dsc1_gene, geneName1, geneName2, geneName3, geneName4, Orco_gene, Rdl_gene, Para_gene, Potassium_channel, Calcium_channel)
titles <- c(Dsc1_name, name1, name2, name3, name4, Orco_name, Rdl_name, Para_name, K_chan_name, Ca_chan_name)

p <- FeaturePlot(dataset,  reduction = 'tsne', features = features, pt.size = .1, slot = "data", order=TRUE, combine = FALSE)
p <- mapply(function(x, t) x + ggtitle(t) + theme(plot.title = element_text(size = 10, face = "bold")), p, titles, SIMPLIFY = FALSE)
wrap_plots(p, ncol=5)

FeaturePlot(dataset, reduction = 'tsne', features = c(Dsc1_gene), pt.size = .1, slot = 'data', order=TRUE)

# Vln.filename <- toString("FeaturePlot.pdf")
# ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 3000, height = 2000)

#genes <- paste(toString(name1),toString(name2),toString(name3),toString(name4), Dsc1_name, Rdl_name, sep="_")
#Vln.filename <- paste("FeaturePlot_",toString(genes),".pdf", sep="")

###################### Feature Plots (Figure S9H, Figure S14H, )

nrow(pseudo_bulk.ORs) 
nrow(pseudo_bulk.IRs) 
nrow(pseudo_bulk.GRs) 

print(pseudo_bulk.ORs)

#Feature Plots for ORs
for (i in 1:20) {
  #i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.ORs$Gene[c])
  geneName2 <- toString(pseudo_bulk.ORs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.ORs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.ORs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.ORs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.ORs$Gene[c+5])
  FeaturePlot(dataset,  reduction = 'tsne', features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", order=TRUE, ncol = 3)
  # genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  # Vln.filename <- paste("FeaturePlot_",toString(genes),".pdf", sep="")
  # ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 3000, height = 2000) #width 2000, height 1000
}
dev.off()

#Feature Plots for IRs
for (i in 1:20) {
  #i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.IRs$Gene[c])
  geneName2 <- toString(pseudo_bulk.IRs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.IRs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.IRs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.IRs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.IRs$Gene[c+5])
  FeaturePlot(dataset,  reduction = 'tsne', features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", order=TRUE, ncol = 3)
  # genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  # Vln.filename <- paste("FeaturePlot_",toString(genes),".pdf", sep="")
  # ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 3000, height = 2000) #width 2000, height 1000
}
dev.off()

#Feature Plots for GRs
for (i in 1:20) {
  i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.GRs$Gene[c])
  geneName2 <- toString(pseudo_bulk.GRs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.GRs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.GRs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.GRs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.GRs$Gene[c+5])
  FeaturePlot(dataset,  reduction = 'tsne', features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", order=TRUE, ncol = 3)
  # genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  # Vln.filename <- paste("FeaturePlot_",toString(genes),".pdf", sep="")
  # ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 3000, height = 2000) #width 2000, height 1000
}
dev.off()

################################################################### Neural Marker Dot Plot (Figure S9I, Figure S13H)

DotPlot(dataset, features = features) + scale_y_discrete(limits = rev)

ggsave("neuron_markers.pdf", limitsize = FALSE, units = "px", width = 2200, height = 4000)


################################################################### Neuron Cluster Marker Dot Plot (Figure S14E)

receptorL <- c("Orco", "Ir25a", "Ir76b", "Gr2", "Gr3", "Or8", "Or49", "LOC5575210")
DotPlot(dataset, features = receptorL, dot.scale = 25) + scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
# Vln.filename <- paste0("DotPlot_MaxPalp.pdf")
# ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 4000, height = 2500)

#######################################

heatmap_features <- c(Dsc1_gene, pseudo_bulk.ORs$Gene)
# 
DoHeatmap(dataset, features = heatmap_features)

Or31_features <- c(Dsc1_gene, "Orco", "Or31")
DoHeatmap(dataset, features = Or31_features)

# DoHeatmap(dataset, features = c(Dsc1_gene, Orco_gene))

#######################################

# Get cells for specific ORs (not working?)

all_features <- c(Dsc1_gene, pseudo_bulk.ORs$Gene)

### create list w/ genes of interest
genelist <- c(Rdl_gene, Dsc1_gene)
Rdl_gene

### subset w/ pasted output
subset_SeuratObj <- subset(dataset, LOC5570466 > 0)
#FeaturePlot(subset_SeuratObj, features=c('Or31', Dsc1_gene))

#DotPlot(subset_SeuratObj, features = genelist) + scale_y_discrete(limits = rev)
DoHeatmap(subset_SeuratObj, features = genelist)

#######################################
KH6 <- "LOC5571255"
Shaw <- "LOC5568798"

Shab <- "LOC5572202"
Shaker <- "LOC5572028"
KQT2 <- "LOC5571969"

TrpA1 <- "LOC5571938"

FeaturePlot(pos_cells, features = c(Dsc1, Shaw), blend=TRUE)
#FeaturePlot(pos_cells, features = c(Dsc1, KH6), cols=c('gray','blue','red'), blend=TRUE)
