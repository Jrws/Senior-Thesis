# Navigate to BAM files, save data
setwd("C:/Users/jonny/OneDrive/Desktop/Code/Research/Knockout/BAM")
BAM_files <- list.files()
BAM_files <- BAM_files[order(c(4,5,6,1,2,3))]

# External Aedes aegypti genome annotation file
annot <- "C:/Users/jonny/OneDrive/Desktop/Code/Research/ref/a_aegypti/genomic.gtf"

# Load Subread package with featureCounts
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
  BiocManager::install("Rsubread") 
  install.packages('Rsubread')  
library("Rsubread")

# featureCounts to count reads mapping across genome
fc <- featureCounts(files=BAM_files,
                    annot.ext=annot,  # external annotation file
                    isGTFAnnotationFile=TRUE,  # GTF format for annotations
                    isPairedEnd=TRUE,
                    countReadPairs=TRUE,     # Counts pairs of reads rather than individual reads, on by default
                    GTF.featureType="exon",  # Coding regions only, exon is default already
                    GTF.attrType="gene_id")  # attribute type already gene_id by default

setwd("C:/Users/jonny/OneDrive/Desktop/Code/Research/Knockout/DiffExp")
write.csv(fc$counts, file="Subreads_counts.csv")

################################################################################

#### five major steps for DESeq2
#1. dds <- makeExampleDESeqDataSet()
#2. dds <- DESeq(dds)
#3. plotMA(dds)
#4. res <- results(dds)
#4. plotMA(res)

# Load DESeq2 package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")

# Load keys for samples
sampleTable <- read.csv(file="sample_keys.csv", header=TRUE)

# Prepare DESeq2 parameters
counts <- fc$counts
sample_cols <- sampleTable
sample_cols$condition <- factor(sample_cols$condition)

# Running differential expression testing using DESeq2
DESeqMat <- DESeqDataSetFromMatrix(countData=counts,
                                 colData=sample_cols,
                                 design=~condition)
dds <- DESeq(DESeqMat)

# Calling results without any arguments will extract the estimated log2 fold changes
# and p values for the last variable in the design formula. If there are more than 2
# levels for this variable, results will extract the results table for a comparison
# of the last level over the first level. This comparison is printed at the top
# of the output: dex trt vs untrt.

res <- results(dds, contrast=c('condition','DSC1_KO', 'WT'))

# Shows meanings of columns of results
mcols(res, use.names=TRUE)

# Summarizes results of DESeq2
summary(res)

# Statistically significantly different results with p-adj < 0.05
res.05 <- results(dds, contrast=c('condition','DSC1_KO', 'WT'), alpha=0.05)

table(res.05$padj < 0.05)  # number of significantly different expressions
table(is.na(res.05$padj))  # number lacking data
dim(res.05)  # number of genes

# Data wrangling/cleaning
res.05.NA <- res.05[is.na(res.05$padj),]       # rows with p-adj == NA

res.05.rmNA <- res.05[!(is.na(res.05$padj)),]  # remove NA in p-adj
res.05.rmNA <- res.05.rmNA[res.05.rmNA$padj < 0.05,]  # filter by sufficient p-adj

dim(res.05.rmNA)  # number of remaining genes not including NA p-adj and p-adj >= 0.05

dim(res.05.rmNA[res.05.rmNA$log2FoldChange < 0,]) # knockout downregulated
upreg <- res.05.rmNA[res.05.rmNA$log2FoldChange > 0,] # kncokout upregulated
upreg.orderpadj <- upreg[order(upreg$padj),]
head(upreg.orderpadj)

res.05.rmNA.orderpadj <- res.05.rmNA[order(res.05.rmNA$padj),]  # sort by p-adj
res.05.rmNA.orderLog2 <- res.05.rmNA[order(res.05.rmNA$log2FoldChange),]  # sort by log2FC

head(res.05.rmNA.orderpadj)
head(res.05.rmNA.orderLog2)
write.csv(res.05.rmNA.orderpadj, file="Knockout_DESeq2.csv")

# lots of transferrins? Iron-binding proteins

################################################################################

# Adding annotations to DE gene list
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install("AnnotationDbi")
  BiocManager::install("AnnotationHub")
library (AnnotationDbi) 
library(AnnotationHub)

# Bug with dbplyr compatibility, downgrade to fix
# devtools::install_version("dbplyr", version = "2.3.4")

# Get Aedes aegypti annotations
ah <- AnnotationHub()
ah <- subset(ah, species == "Aedes aegypti")
# AH10366  | hom.Aedes_aegypti.inp8.sqlite
# AH112940 | org.Aedes_aegypti.eg.sqlite 
org.Aaeg.eg.db <- ah[["AH112940"]]
columns(org.Aaeg.eg.db)

# Get gene IDs for RNA-Seq dataset
Annotated_res_geneIDs <- read.csv(file="Knockout_DESeq2.csv")
row.names(Annotated_res_geneIDs) <- substring(Annotated_res_geneIDs[,1], first = 4)
head(Annotated_res_geneIDs)

# We can use the mapIds function to add individual columns to our results table.
# We provide the row names of our results table as a key, and specify that
# keytype=ENSEMBL. The column argument tells the mapIds function which information
# we want, and the multiVals argument tells the function what to do if there are
# multiple possible values for a single input value. Here we ask to just give us
# back the first one that occurs in the database. To add the gene symbol and
# Entrez ID, we call mapIds twice.

Annotated_res_geneIDs$REFSEQ <- mapIds(org.Aaeg.eg.db,
                                       keys=row.names(Annotated_res_geneIDs),
                                       column="REFSEQ",
                                       keytype="GID",
                                       multiVals="first")
#y$genes$symbol <- res$symbol
Annotated_res_geneIDs$GENENAME <- mapIds(org.Aaeg.eg.db,
                                         keys=row.names(Annotated_res_geneIDs),
                                         column="GENENAME",
                                         keytype="GID",
                                         multiVals="first")
#y$genes$symbol <- res$symbol
Annotated_res_geneIDs$GO <- mapIds(org.Aaeg.eg.db,
                                   keys=row.names(Annotated_res_geneIDs),
                                   column="GO",
                                   keytype="GID",
                                   multiVals="first")

write.csv(Annotated_res_geneIDs, file="Knockout_DESeq2_annotated.csv")

################################################################################

# DESeq2 Normalized Reads, estimate library size factors

dds <- estimateSizeFactors(dds)
sizeFactors(dds)  # view size factors for each sample

normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="Knockout_DESeq2_normalized.csv")

################################################################################

# MA Plots for gene data
dim(res)  # view number of genes again

# No labeled genes
DESeq2::plotMA(res, alpha=0.05, ylim=c(-8,8), colSig = "red",cex=.7, xlab="mean expression", ylab="log2(FC)")

# Labeled genes: top two lowest p-adj genes
plotMA(res, ylim=c(-8,8),colSig = "red",cex=.7)
  with(res.05.rmNA.orderpadj[c(1,2), ], {
    points(baseMean, log2FoldChange, col="blue", cex=2, lwd=2)  # circle
    text(baseMean, log2FoldChange, row.names(res.05.rmNA.orderpadj[c(1,2),])  , pos=4, col="blue")  #lable: pos=1 bottom,   pos=2  left, 3 up and 4 right
  }
)

# Labeled genes: top two lowest Log2FC
plotMA(res, ylim=c(-8,8),colSig = "red",cex=.7)
  with(res.05.rmNA.orderLog2[c(1,2), ], {
    points(baseMean, log2FoldChange, col="blue", cex=2, lwd=2)  # circle
    text(baseMean, log2FoldChange, row.names(res.05.rmNA.orderLog2[c(1,2),])  , pos=c(4,4), col="red")  
  }
)

### Look for genes of interest
genes_of_interest <- c("LOC5575962","LOC5567058", "LOC5567355", "LOC5570466")
                      # Dsc1,        Orco,         Para,         Rdl
res[genes_of_interest,]

plotMA(res, ylim=c(-8,8),colSig = "red",cex=.7)
  with(res[genes_of_interest,], {
    points(baseMean, log2FoldChange, col="blue", cex=2, lwd=2)  # circle
    text(baseMean, log2FoldChange, genes_of_interest, pos=c(4,4), font=6, col="blue")
    #text(x=0.8e4, y=-2.7, "LOC5580178",col="blue", font=2, cex=1)  #manually locate the text by (x,y) coordinates
  }
)

################################################################################

# Gene expression pattern plots
BiocManager::install("ggplot2")
library(ggplot2)

# Load gene list
gene_list <- read.csv("Knockout_DESeq2_normalized.csv", header=TRUE, row.names=1, fileEncoding="UTF-8-BOM")
head(gene_list)
dim(gene_list)

GeneOfInterest="LOC5567355"  # transferrin LOC5579417 LOC5577029
#GeneOfInterest=""   # Para

gene_list[GeneOfInterest,]
gene_plot<-gene_list[GeneOfInterest,]
T_gene_plot<-data.frame(t(gene_plot))
condition<-factor(sampleTable$condition)
condition <- relevel(condition, "WT")
T_gene_plot <- cbind(T_gene_plot, condition)
T_gene_plot

# Box plot
p1<-ggplot(T_gene_plot, aes(x=condition, y=get(GeneOfInterest), fill = condition)) + geom_dotplot(binaxis='y', stackdir='center') + scale_x_discrete(limits = levels(condition)) +
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), 
                 geom="crossbar", width=0.5) + scale_fill_manual(values=c("red", "blue")) +
  labs(title=paste("Expression levels of", GeneOfInterest),x="Conditions", y = "Normalized Read Counts") + theme(plot.title = element_text(hjust = 0.5))
p1

# Error bars as standard deviation
avg_sd <- do.call(data.frame, aggregate(T_gene_plot[,GeneOfInterest], list(T_gene_plot$condition), FUN=function(x) c(avg = mean(x), std = sd(x))))
colnames(avg_sd) <- c("Genotype","Average","StdDev")
avg_sd

p2<-ggplot(avg_sd, aes(fill=Genotype, x=Genotype, y=Average)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=Average-StdDev, ymax=Average+StdDev), alpha=0.8, width=0.4) +
  labs(title=paste("Average expression levels of", GeneOfInterest),x="Condition", y = "Normalized Read Counts") + theme(plot.title = element_text(hjust = 0.5))
p2

################################################################################

##### PCA

### Option 1
vsd <- vst(dds)
pcaData <- plotPCA(vsd, intgroup="condition", returnData = TRUE)
pcaData
plotPCA(vsd, intgroup="condition")  # interesting group 

### Option 2
t.ReadTable<-t(counts)
pca.obj <- prcomp(t.ReadTable)
#library(ggplot2)  # restart the Rstudio  if there is a version error with DESeq2 
samples <- sampleTable$condition
samples
dtp <- data.frame('Samples'=samples, pca.obj$x[,1:3]) # the first two components are selected (NB: you can also select 3 for 3D plottings or 3+)
ggplot(data = dtp) + 
  geom_point(aes(x=PC1,y=PC2,col=Samples),size = 3) + 
  theme_minimal() 

### 3D
install.packages("plotly")
library(plotly)
plot_ly(dtp, x=pca.obj$x[,1], y=pca.obj$x[,2], z=pca.obj$x[,3], type="scatter3d", mode="markers", color= dtp$Samples, colors=c("#48bf2a",'#BF382A', '#0C4B8E' ))

dtp
################################################################################################################
