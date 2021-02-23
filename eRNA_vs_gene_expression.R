#!/bin/Rscript 

# Load the required libraries
library(dplyr)
library(ggpubr)
library(gtools)

# Expression relationship between enhancers and genes

# Read in the input data
data <- read.table("R_input.txt", header=T, sep="\t")

# Replace all zeros to NA
data[data == 0] <- NA
head(data)

# Remove all the ones with NA
data <- na.omit(data)

# Add 1 to the data
# data[,2:3] <- data[,2:3] + 1

# Isolate the enhancer data
enhancer <- subset(data, subset = grepl("_te_|_se_", data$Feature_ID))

# Isolate the gene data
gene <- subset(data, subset = !grepl("_te_|_se_", data$Feature_ID))

# Replace everything after the underscore for enhancer name
enhancer$Feature_ID <- gsub("\\_.*","",enhancer$Feature_ID)

# Add a type to enhancer data
enhancer$type <- rep("enhancer", length(enhancer$Feature_ID))

# Add a type to gene data
gene$type <- rep("gene", length(gene$Feature_ID))

# Now use dplyr to join the gene data and the enhancer data based on the gene name
expression_data <- inner_join(enhancer, gene, by=c("Feature_ID"))

# Separate dox minus from dox plus
cols1 <- c(1,2,5)
dox_minus <- expression_data[,cols1]
#dox_minus[,2:3] <- log2(dox_minus[,2:3])
colnames(dox_minus) <- c("gene","enhancer_TPM","gene_TPM")
dox_minus <- mutate(dox_minus, quantile_rank = ntile(dox_minus$enhancer_TPM, 5)) %>% as.data.frame()
head(dox_minus)

# Make the plot
pdf("dox_minus.pdf")
ggscatter(dox_minus, x = "gene_TPM", y = "enhancer_TPM", color = "black", shape = 21, size = 3, add = "reg.line") +
  stat_cor(label.x = 1, label.y = 5) +
  stat_regline_equation(label.x = 1, label.y = 4) + xscale("log2", .format = TRUE) + yscale("log2", .format = TRUE)
dev.off()

write.table(dox_minus, "dox_minus_quantile_expression.tsv", row.names=F, col.names=T, sep="\t", quote=F)

cols2 <- c(1,3,6)
dox_plus <- expression_data[,cols2]
#dox_plus[,2:3] <- log2(dox_plus[,2:3])
colnames(dox_plus) <- c("gene","enhancer_TPM","gene_TPM")
dox_plus <- mutate(dox_plus, quantile_rank = ntile(dox_plus$enhancer_TPM, 5)) %>% as.data.frame()
head(dox_plus)

write.table(dox_plus, "dox_plus_quantile_expression.tsv", row.names=F, col.names=T, sep="\t", quote=F)

# Make the plot
pdf("dox_plus.pdf")
ggscatter(dox_plus, x = "gene_TPM", y = "enhancer_TPM", color = "black", shape = 21, size = 3, add = "reg.line") +
  stat_cor(label.x = 1, label.y = 5) +
  stat_regline_equation(label.x = 1, label.y = 4) + xscale("log2", .format = TRUE) + yscale("log2", .format = TRUE)
dev.off()

# Get the expression of genes based on the quantile of eRNA expression

pdf("quantile_dox_minus.pdf")
ggviolin(dox_minus, "quantile_rank", "gene_TPM",  color = "quantile_rank",
   add = "boxplot") + yscale("log2", .format = TRUE)
dev.off()

pdf("quantile_dox_plus.pdf")
ggviolin(dox_plus, "quantile_rank", "gene_TPM",  color = "quantile_rank",
   add = "boxplot") + yscale("log2", .format = TRUE)
dev.off()

# Get the common genes between dox_minus and dox_plus
dox_minus$treatment <- rep("dox_minus", nrow(dox_minus))
dox_plus$treatment <- rep("dox_plus", nrow(dox_plus))
dox <- rbind(dox_minus, dox_plus)

# Make the plot for dox
pdf("quantile_dox_unpaired.pdf")
ggviolin(dox, "quantile_rank", "gene_TPM",  color = "treatment",
   add = "boxplot") + yscale("log2", .format = TRUE) + stat_compare_means(aes(group = treatment), paired=FALSE)
dev.off()

pdf("quantile_dox_paired.pdf")
ggviolin(dox, "quantile_rank", "gene_TPM",  color = "treatment",
   add = "boxplot") + yscale("log2", .format = TRUE) + stat_compare_means(aes(group = treatment), paired=TRUE) 
dev.off()

# Make the plot for dox but now plotting eRNA TPM
pdf("quantile_eRNA_dox_unpaired.pdf")
ggviolin(dox, "quantile_rank", "enhancer_TPM",  color = "treatment",
   add = "boxplot") + yscale("log2", .format = TRUE) + stat_compare_means(aes(group = treatment), paired=FALSE)
dev.off()

pdf("quantile_eRNA_dox_paired.pdf")
ggviolin(dox, "quantile_rank", "enhancer_TPM",  color = "treatment",
   add = "boxplot") + yscale("log2", .format = TRUE) + stat_compare_means(aes(group = treatment), paired=TRUE)
dev.off()

# Isolate the unique genes in the fifth quantile for dox- and dox+
dox_plus_5th_quant <- subset(dox_plus, dox_plus$quantile_rank == 5)
dox_plus_5th_quant_genes <- unique(dox_plus_5th_quant$gene)
dox_plus_5th_quant_gene_df <- data.frame(dox_plus_5th_quant_genes)
colnames(dox_plus_5th_quant_gene_df) <- c("gene")
write.table(dox_plus_5th_quant_gene_df, "dox_plus_5th_quant_genes.tsv", col.names=T, row.names=F, sep="\t", quote=F)

dox_minus_5th_quant <- subset(dox_minus, dox_minus$quantile_rank == 5)
dox_minus_5th_quant_genes <- unique(dox_minus_5th_quant$gene)
dox_minus_5th_quant_gene_df <- data.frame(dox_minus_5th_quant_genes)
colnames(dox_minus_5th_quant_gene_df) <- c("gene")
write.table(dox_minus_5th_quant_gene_df, "dox_minus_5th_quant_genes.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# Now join the dox- and dox+ and plot eRNA and gene expression for dox+ based on dox- genes and quantiles
dox_minus_plus <- inner_join(dox_minus, dox_plus, by=c("gene"))
head(dox_minus_plus)
cols <- c(1,4,6,7)
dox_minus_plus <- dox_minus_plus[,cols]
colnames(dox_minus_plus) <- c("gene","quantile_rank","enhancer_TPM","gene_TPM")
dox_minus_plus$quantile_rank <- factor(dox_minus_plus$quantile_rank, levels=c(1,2,3,4,5))
head(dox_minus_plus)
str(dox_minus_plus)

# Get the expression of genes based on the quantile of eRNA expression
pdf("quantile_dox_minus_plus_eRNA_expression.pdf")
ggviolin(dox_minus_plus, "quantile_rank", "enhancer_TPM",  color = "quantile_rank",
   add = "boxplot") + yscale("log2", .format = TRUE)
dev.off()

# Get the expression of genes based on the quantile of eRNA expression
pdf("quantile_dox_minus_plus_gene_expression.pdf")
ggviolin(dox_minus_plus, "quantile_rank", "gene_TPM",  color = "quantile_rank",
   add = "boxplot") + yscale("log2", .format = TRUE)
dev.off()

# Get eRNA expression and gene expression both for dox- and dox+ based on the quantiles of eRNA from dox- condition
dox_minus_plus <- inner_join(dox_minus, dox_plus, by=c("gene"))
head(dox_minus_plus)
cols <- c(1,2,3,4,5)
dox_minus <- dox_minus_plus[,cols]
colnames(dox_minus) <- c("gene", "enhancer_TPM", "gene_TPM", "quantile_rank", "treatment")
cols <- c(1,6,7,4,9)
dox_plus <- dox_minus_plus[,cols]
colnames(dox_plus) <- c("gene", "enhancer_TPM", "gene_TPM", "quantile_rank", "treatment")
dox_total <- rbind(dox_minus, dox_plus)
tail(dox_total)
dox_total$quantile_rank <- factor(dox_total$quantile_rank, levels=c(1,2,3,4,5))
dox_total$treatment <- factor(dox_total$treatment)

# Make the plot for dox plus gene expression based on the genes on quantiles from dox minus eRNA
pdf("quantile_dox_total_unpaired.pdf")
ggviolin(dox_total, "quantile_rank", "gene_TPM",  color = "treatment",
   add = "boxplot") + yscale("log2", .format = TRUE) + stat_compare_means(aes(group = treatment), paired=FALSE)
dev.off()

# Make the plot for dox plus based on the genes on quantiles from dox minus
pdf("quantile_eRNA_dox_total_unpaired.pdf")
ggviolin(dox_total, "quantile_rank", "enhancer_TPM",  color = "treatment",
   add = "boxplot") + yscale("log2", .format = TRUE) + stat_compare_means(aes(group = treatment), paired=FALSE)
dev.off()

# Now import the data for log2fold change between dox+ and dox- (dox_minus/dox_plus)
fc <- read.table("input_for_foldchange.tsv", header=T, sep="\t")

# Remove the ones with not known fold-change
fc <- na.omit(fc)

# Isolate the enhancer data
fc_enhancer <- subset(fc, subset = grepl("_te_|_se_", fc$Gene))

# Isolate the gene data
fc_gene <- subset(fc, subset = !grepl("_te_|_se_", fc$Gene))

# Replace everything after the underscore for enhancer name
fc_enhancer$Gene <- gsub("\\_.*","",fc_enhancer$Gene)

# Add a type to enhancer data
fc_enhancer$type <- rep("enhancer", length(fc_enhancer$Gene))

# Add a type to gene data
fc_gene$type <- rep("gene", length(fc_gene$Gene))

# Now use dplyr to join the gene data and the enhancer data based on the gene name
fc_data <- inner_join(fc_enhancer, fc_gene, by=c("Gene"))

# Select the columns we care about from fc_data
cols3 <- c(1,2,4)
fc_data <- fc_data[,cols3]
colnames(fc_data) <- c("gene","enhancer_FC_expression","gene_FC_expression")

genes <- c('Nanog', 'Trim28', 'Pou5f1', 'miR290a', 'Sox2', 'Essrb')
select_data <- subset(fc_data, fc_data$gene %in% genes) 

# Create a scatterplot
pdf("fold_change_expression.pdf")
ggscatter(fc_data, x = "enhancer_FC_expression", y = "gene_FC_expression", color = "black", shape = 21, size = 3, add = "reg.line") +
  stat_cor(label.x = 1, label.y = 6) +
  stat_regline_equation(label.x = 1, label.y = 5) 
dev.off()


