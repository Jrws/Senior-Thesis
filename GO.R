# Functional Enrichment Analysis using TopGO
# Note: depends on some objects in environment from DiffExp.R
setwd("C:/Users/jonny/OneDrive/Desktop/Code/Research/Knockout/FuncEnrich")

load("C:/Users/jonny/OneDrive/Desktop/Code/Research/Knockout/DiffExp/DE.Rdata")

save.image(file="C:/Users/jonny/OneDrive/Desktop/Code/Research/Knockout/DiffExp/DE.Rdata")

########################### (1) Data preparation ###############################

##### Installation of topGO, AnnotationForge #####
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")
BiocManager::install("AnnotationForge")

library("topGO")
library("AnnotationForge")

install.packages('rjson')
library('rjson')

##### Preparing gene lists #####
res.rmNA <- res[!is.na(res$padj),]    # remove genes with NA padj values

all_genes <- setNames(res.rmNA$padj, rownames(res.rmNA))  # named vector of genes, padj values
DE_genes <- function(padj) {  # filter function for differentially expressed genes with padj < 0.05
  return(padj < 0.05)
}

length(all_genes)   # num genes in gene universe
length(DE_genes)    # num differentially expressed genes

write.csv(res.rmNA,file='allRes.csv')


##### Create new gene annotations object from AnnotationHub #####
# reload Aedes Aegypti annotations database
org.Aaeg.eg.db <- ah[["AH112940"]]

# Data wrangling to format dataframe, may be unnecessary
allAnnotations <- read.csv(file='allRes.csv', header=TRUE, row.names=1, fileEncoding="UTF-8-BOM")
allAnnotations <- tibble::rownames_to_column(allAnnotations, var='X')
row.names(allAnnotations) <- substring(allAnnotations[,1], first = 4)

allAnnotations$REFSEQ <- mapIds(org.Aaeg.eg.db,
                               keys=row.names(allAnnotations),
                               column="REFSEQ",
                               keytype="GID",
                               multiVals="first")
#y$genes$symbol <- res$symbol
allAnnotations$GENENAME <- mapIds(org.Aaeg.eg.db,
                                 keys=row.names(allAnnotations),
                                 column="GENENAME",
                                 keytype="GID",
                                 multiVals="first")

# create new column for list of all GO annotations from org.Aaeg.eg.db object
allAnnotations$GO <- mapIds(org.Aaeg.eg.db,
                            keys=row.names(allAnnotations),
                            column="GO",
                            keytype="GID",
                            multiVals="list")

# Create named list of character vectors of all genes and all their GO annotations
allgene2GO <- setNames(allAnnotations$GO, allAnnotations$X)

# Number of GO annotations
length(allAnnotations$GO[!is.na(allAnnotations$GO)])
length(unique(allAnnotations$GO[!is.na(allAnnotations$GO)]))  # unique GO annotations

### Didn't work, no such column gene_id
### org.Aaeg.eg.sqlite package formation
# library(AnnotationForge)
# library(AnnotationHub)
# file.copy(AnnotationHub::cache(ah["AH112940"]), "./org.Aaeg.eg.sqlite")
# seed <- new("AnnDbPkgSeed", Package='org.Aaeg.eg.db', Version="0.0.1", Author="Jonathan Shi",Maintainer="Jonathan Shi <jonathan.shi@duke.edu>",PkgTemplate="NOSCHEMA.DB",AnnObjPrefix="org.Aaeg.eg",organism="Aedes aegypti",species="Aedes aegypti",biocViews="annotation",manufacturerUrl="none",manufacturer="none",chipName="none")
# makeAnnDbPkg(seed, "org.Aaeg.eg.sqlite")
# 
# # install built org.Aaeg.eg.sqlite package for use in topGO
# install.packages("org.Aaeg.eg.db", type="source", repos=NULL)
# library("org.Aaeg.eg.db")

# Use NCBI RefSeq Annotations
allgene2GO <- fromJSON(file='ncbi_annot.json')

##### Create topGO data objects #####
CC_GOdata <- new("topGOdata",
                 ontology="CC",          # CC (cellular component), MF (molecular function), or BP (biological process)
                 allGenes=all_genes,     # all genes and their padj-values 
                 geneSel=DE_genes,       # differentially expressed genes, padj < 0.05
                 nodeSize=10,            # Number of genes needed for a particular GO annotation for graph
                 annot=annFUN.gene2GO,   # Map each gene to a GO annotation
                 gene2GO=allgene2GO)     # Custom annotations from AnnotationHub GO annotations

BP_GOdata <- new("topGOdata",
              ontology="BP",          # CC (cellular component), MF (molecular function), or BP (biological process)
              allGenes=all_genes,     # all genes and their padj-values 
              geneSel=DE_genes,       # differentially expressed genes, padj < 0.05
              nodeSize=10,            # Number of genes needed for a particular GO annotation for graph
              annot=annFUN.gene2GO,   # Map each gene to a GO annotation
              gene2GO=allgene2GO)     # Custom annotations from AnnotationHub GO annotations

MF_GOdata <- new("topGOdata",
                 ontology="MF",          # CC (cellular component), MF (molecular function), or BP (biological process)
                 allGenes=all_genes,     # all genes and their padj-values 
                 geneSel=DE_genes,       # differentially expressed genes, padj < 0.05
                 nodeSize=10,            # Number of genes needed for a particular GO annotation for graph
                 annot=annFUN.gene2GO,   # Map each gene to a GO annotation
                 gene2GO=allgene2GO)     # Custom annotations from AnnotationHub GO annotations

CC_GOdata
BP_GOdata
MF_GOdata

########################### (2) Enrichment Tests ###############################

### Fisher Tests
CC_resultFisher <- runTest(CC_GOdata, algorithm="classic", statistic="fisher")
BP_resultFisher <- runTest(BP_GOdata, algorithm="classic", statistic="fisher")
MF_resultFisher <- runTest(MF_GOdata, algorithm="classic", statistic="fisher")

CC_resultFisher
BP_resultFisher
MF_resultFisher

### Kolmogorov-Smirnov Test: classic
CC_resultKS_classic <- runTest(CC_GOdata, algorithm="classic", statistic="ks")
BP_resultKS_classic <- runTest(BP_GOdata, algorithm="classic", statistic="ks")
MF_resultKS_classic <- runTest(MF_GOdata, algorithm="classic", statistic="ks")

CC_resultKS_classic
BP_resultKS_classic
MF_resultKS_classic

### Kolmogorov-Smirnov Test: elim
CC_resultKS_elim <- runTest(CC_GOdata, algorithm="elim", statistic="ks")
BP_resultKS_elim <- runTest(BP_GOdata, algorithm="elim", statistic="ks")
MF_resultKS_elim <- runTest(MF_GOdata, algorithm="elim", statistic="ks")

CC_resultKS_elim
BP_resultKS_elim
MF_resultKS_elim

### Summary tables of statistical tests for GO terms

CC_allRes <- GenTable(CC_GOdata, classicFisher=CC_resultFisher,
                      classicKS=CC_resultKS_classic, elimKS=CC_resultKS_elim,
                      orderBy="elimKS", ranksOf="classicFisher", topNodes=length(usedGO(CC_GOdata)))
CC_allRes

BP_allRes <- GenTable(BP_GOdata, classicFisher=BP_resultFisher,
                      classicKS=BP_resultKS_classic, elimKS=BP_resultKS_elim,
                      orderBy="elimKS", ranksOf="classicFisher", topNodes=length(usedGO(BP_GOdata)))
BP_allRes

MF_allRes <- GenTable(MF_GOdata, classicFisher=MF_resultFisher,
                      classicKS=MF_resultKS_classic, elimKS=MF_resultKS_elim,
                      orderBy="elimKS", ranksOf="classicFisher", topNodes=length(usedGO(MF_GOdata)))
MF_allRes


write.csv(CC_allRes, 'CC_allRes.csv')
write.csv(BP_allRes, 'BP_allRes.csv')
write.csv(MF_allRes, 'MF_allRes.csv')

### Check differences between classic and elim methods of K-S test

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

# Plot classic vs elim p-values in a scatterplot
CC_pValue.classic <- score(CC_resultKS_classic)
CC_pValue.elim <- score(CC_resultKS_elim)[names(CC_pValue.classic)]
CC_gstat <- termStat(CC_GOdata, names(CC_pValue.classic))
CC_gSize <- CC_gstat$Annotated / max(CC_gstat$Annotated) * 4
CC_gCol <- colMap(CC_gstat$Significant)
plot(CC_pValue.classic, CC_pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = CC_gSize, col = CC_gCol)

# Find terms where elim p-value less conservative than classic p-value
CC_sel.go <- names(CC_pValue.classic)[CC_pValue.elim < CC_pValue.classic]
cbind(termStat(CC_GOdata, CC_sel.go),
      elim = CC_pValue.elim[CC_sel.go],
      classic = CC_pValue.classic[CC_sel.go]) # None

# BP
BP_pValue.classic <- score(BP_resultKS_classic)
BP_pValue.elim <- score(BP_resultKS_elim)[names(BP_pValue.classic)]
BP_gstat <- termStat(BP_GOdata, names(BP_pValue.classic))
BP_gSize <- BP_gstat$Annotated / max(BP_gstat$Annotated) * 4
BP_gCol <- colMap(BP_gstat$Significant)
plot(BP_pValue.classic, BP_pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = BP_gSize, col = BP_gCol)

BP_sel.go <- names(BP_pValue.classic)[BP_pValue.elim < BP_pValue.classic]
cbind(termStat(BP_GOdata, BP_sel.go),
      elim = BP_pValue.elim[BP_sel.go],
      classic = BP_pValue.classic[BP_sel.go]) # None

# MF
MF_pValue.classic <- score(MF_resultKS_classic)
MF_pValue.elim <- score(MF_resultKS_elim)[names(MF_pValue.classic)]
MF_gstat <- termStat(MF_GOdata, names(MF_pValue.classic))
MF_gSize <- MF_gstat$Annotated / max(MF_gstat$Annotated) * 4
MF_gCol <- colMap(MF_gstat$Significant)
plot(MF_pValue.classic, MF_pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = MF_gSize, col = MF_gCol)

MF_sel.go <- names(MF_pValue.classic)[MF_pValue.elim < MF_pValue.classic]
cbind(termStat(MF_GOdata, MF_sel.go),
      elim = MF_pValue.elim[MF_sel.go],
      classic = MF_pValue.classic[MF_sel.go]) # None

# Can list GO terms and numbers of genes where elim p-value < classic p-value, which is more expected

BiocManager::install("Rgraphviz")

### GO graph
#   - Significant nodes = rectangles
par(bg='white', cex=0.7)  # doesn't work
showSigOfNodes(CC_GOdata,
               score(CC_resultKS_elim),
               firstSigNodes=5,   # subgraph induced by 5 most significant GO terms
               useInfo = 'all')
printGraph(CC_GOdata, CC_resultKS_elim, 5, useInfo='all', pdfSW=TRUE)

par(bg='white', cex=0.5)
showSigOfNodes(BP_GOdata,
               score(BP_resultKS_elim),
               firstSigNodes=5,
               useInfo = 'all')

par(bg='white', cex=0.8)  # doesn't work?
showSigOfNodes(MF_GOdata,
               score(MF_resultKS_elim),
               firstSigNodes=5,
               useInfo = 'all')
