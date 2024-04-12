# Generating a clustering heatmap of differentially expressed genes
setwd("C:/Users/jonny/OneDrive/Desktop/Code/Research/Knockout/FuncEnrich")

load("C:/Users/jonny/OneDrive/Desktop/Code/Research/Knockout/DiffExp/DE.Rdata")

save.image(file="C:/Users/jonny/OneDrive/Desktop/Code/Research/Knockout/DiffExp/DE.Rdata")

library('pheatmap')

sample_names <- c('WT 1','WT 2','WT 3','DSC1_KO 1', 'DSC1_KO 2', 'DSC1_KO 3')
# Create normalized count dataframe for only differentially expressed genes, relabel columns
diff_genes <- c(diff_genes, "LOC5577794")
sig_counts <- normalized_counts[diff_genes,]
colnames(sig_counts) <- sample_names

# Heatmap annotations
ann_df <- data.frame(Condition = rep(c('WT','DSC1 -/-'), c(3,3)))
row.names(ann_df) <- colnames(sig_counts)

ann_colors <- list(
  Condition = c("WT" = "#3A90E3",
                "DSC1 -/-" = "#FF8D42")
)

# Clustered heatmap of all 619 differentially expressed genes
heatmap <- pheatmap(
  sig_counts,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  main = 'Clustered Heatmap of All Differentially Expressed Genes',
  colorRampPalette(c(
    'blue',
    'lightyellow',
    'red'
  ))(25
  ),
  scale = 'row',
  angle_col = 0,
  annotation_col = ann_df,
  annotation_colors = ann_colors
)
heatmap

# Notable genes (from GO analysis)
#TODO: annotate rows of cluster heatmap with GO terms

gene_filter <- function(GO_type, GO_annotation) {
  if (GO_type == 'CC') {
    annot <- allAnnotations[allAnnotations$X %in% intersect(unlist(unname(genesInTerm(CC_GOdata, GO_annotation))), rownames(sig_counts)),]
  }
  else if (GO_type == 'BP') {
    #annot <- allAnnotations[allAnnotations$X %in% intersect(unlist(unname(genesInTerm(BP_GOdata, GO_annotation))), rownames(sig_counts)),]
    annot <- allAnnotations[allAnnotations$X %in% unlist(unname(genesInTerm(BP_GOdata, GO_annotation))),]
  }
  else if (GO_type == 'MF') {
    #annot <- Annotated_res_geneIDs[Annotated_res_geneIDs$X %in% intersect(unlist(unname(genesInTerm(MF_GOdata, GO_annotation))), rownames(sig_counts)),]
    #annot <- allAnnotations[allAnnotations$X %in% intersect(unlist(unname(genesInTerm(MF_GOdata, GO_annotation))), rownames(sig_counts)),]
    annot <- allAnnotations[allAnnotations$X %in% unlist(unname(genesInTerm(MF_GOdata, GO_annotation))),]
  }
  else {
    print("Invalid GO_type, must be CC/BP/MF")
    return -1
  }
  #return (annot[order(annot$log2FoldChange, decreasing=TRUE),]$X[1:3])  # return top 3 greatest absolute value log2FoldChange significant genes with the GO annotation
  #return (annot[order(annot$log2FoldChange, decreasing=TRUE),]$X)
  return (annot)
}

trans_genes <- gene_filter('BP','GO:0006412') # translation
proteolysis_genes <- gene_filter('BP','GO:0006508') # proteolysis
iron_bind_genes <- gene_filter('MF','GO:0005506') # iron ion binding

vg_potassium_channels <- gene_filter('MF','GO:0005249') # voltage-gated potassium channel
write.csv(vg_potassium_channels[,1:9], file='potassium_channels.csv')

ion_transport_genes <- gene_filter('MF','GO:0042625')
visual_genes <- gene_filter('BP','GO:0007601')
# RNA helicase GO didn't return any genes

#notable_genes <- c(trans_genes, proteolysis_genes, iron_bind_genes)

notable_genes <- c('LOC5571255','LOC5568798','LOC5573303','LOC5569124','LOC5576882','LOC5572301','LOC5571938','LOC5577794')

notable_counts <- sig_counts[notable_genes,]

# Notable annotations dataframes
#notable_ann_df <- data.frame(Function = rep(c('Translation','Proteolysis','Iron ion binding'), c(3,3,3)))
notable_ann_df <- data.frame(Function = rep(c('VGKC', 'Iron binding','Vision','TRP'), c(2,2,2,2)))
row.names(notable_ann_df) <- rownames(notable_counts)

notable_ann_colors <- list(
  Condition = c("WT" = "cyan3",
                "DSC1 -/-" = "lightpink2"),
  Function = c("VGKC" = "lightgreen",
               "Iron binding" = "slategray",
               "Vision" = 'skyblue',
               "TRP" = 'gold')
)

# notable_ann_colors <- list(
#   Condition = c("WT" = "cyan3",
#                 "DSC1 -/-" = "lightpink2"),
#   Function = c("Translation" = "lightgreen",
#                 "Proteolysis" = "lightpink3",
#                 "Iron ion binding" = "slategray")
# )

# try to average columns
# remove DSC1_KO 3, repeat DESeq2?
heatmap <- pheatmap(
  notable_counts,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  main = 'Clustered Heatmap of Selected Differentially Expressed Genes',
  colorRampPalette(c(
    'blue',
    'lightyellow',
    'red'
  ))(25
  ),
  scale = 'row',
  angle_col = 0,
  annotation_row = notable_ann_df,
  annotation_col = ann_df,
  annotation_colors = notable_ann_colors
)

