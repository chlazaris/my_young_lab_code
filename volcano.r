#!/bin/Rscript

# Load the required libraries
library(EnhancedVolcano)

# Read in the DESeq2 data
res <- read.table("input/all_deseq2_results.tsv", header=TRUE, stringsAsFactors=FALSE, check.names=F)
se_genes <- read.table("input/SE_uniq_genes.txt", header=FALSE)
colnames(se_genes) <- c("gene")

#lab_italics <- paste0("italic('", res$Gene, "')")
#selectLab_italics <- paste0("italic('", se_genes$gene, "')")

# Get the genes that change significantly
de_genes <- subset(res, abs(res$log2FoldChange)>1 & res$padj<0.05)
# Save the differentially expressed genes
write.table(de_genes, "diff_expressed_genes.tsv", row.names=F, col.names=T, sep="\t", quote=F)

# Isolate the differentially expressed, SE-bound genes
de_se_genes <- de_genes[de_genes$Gene %in% se_genes$gene,]
# Save the differentially expressed genes in a table
write.table(de_se_genes, "SE_regulated_diff_expressed_genes.tsv", row.names=F, col.names=T, sep="\t", quote=F)

# Isolate the differentially expressed, Ikaros-bound genes
ikaros_genes <- read.table("input/Ikaros_TSS_bound_protein_coding_genes.tsv", header=FALSE)
colnames(ikaros_genes) <- c("gene")
ikaros_genes
de_ikzf1_genes <- de_genes[de_genes$Gene %in% ikaros_genes$gene,]
de_ikzf1_genes

pdf("se_volcano.pdf", useDingbats=FALSE)
EnhancedVolcano(res,
  lab = res$Gene,
  selectLab = de_se_genes$Gene,
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-6, 6),
  title = '',
  #subtitle = paste0('p-value cutoff drawn ',
  #    'at equivalent of adjusted p=0.05'),
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,
  labSize = 3.0,
  colAlpha = 1,
  cutoffLineType = 'solid',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.5,
  col=c('gray', 'gray', 'gray', 'red3'),
  #hlineCol = c('black', 'black', 'black', 'black'),
  #hlineType = c('longdash', 'longdash', 'dotdash', 'dotdash'),
  #hlineWidth = c(0.4, 0.4, 0.8, 0.8),
  legendPosition = 'bottom',
  boxedLabels = FALSE,
  drawConnectors = TRUE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE)
dev.off()

pdf("ikzf1_target_volcano.pdf", useDingbats=FALSE)
EnhancedVolcano(res,
  lab = res$Gene,
  selectLab = de_ikzf1_genes$Gene,
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-6, 6),
  title = '',
  #subtitle = paste0('p-value cutoff drawn ',
  #    'at equivalent of adjusted p=0.05'),
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,
  labSize = 3.0,
  colAlpha = 1,
  cutoffLineType = 'solid',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.5,
  col=c('gray', 'gray', 'gray', 'red3'),
  #hlineCol = c('black', 'black', 'black', 'black'),
  #hlineType = c('longdash', 'longdash', 'dotdash', 'dotdash'),
  #hlineWidth = c(0.4, 0.4, 0.8, 0.8),
  legendPosition = 'bottom',
  boxedLabels = FALSE,
  drawConnectors = TRUE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE)
dev.off()

