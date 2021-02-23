#!/bin/Rscript

# module load r/3.3.0

# Check if the required packages are present
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggrepel")) install.packages("ggrepel")

# Load the required packages
library(ggplot2)
library(ggrepel)

# Remove the scientific notation
options(scipen=999)

# Read the allEnhancers file here (ENHANCER TO GENE)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {print ("USAGE: Rscript plot_superenhancergraph.R INPUT_FILE"); q(save="no")}

# Get the input
input <- args[1]

# Read in the data - they have header
data <- read.table(sprintf("%s", input), sep="\t", header=T)
# Get the cutoff to call the SEs
#cutoff <- as.numeric(as.character(args[2]))
#genes <- read.table(sprintf("%s", args[3]), sep="\t", header=F)
#genes <- as.vector(as.character(genes$V1))

# Name the data
colnames(data) <- c("REGION_ID","CHROM","START","STOP","NUM_LOCI",
	"CONSTITUENT_SIZE","VALUE","OVERLAP_GENES","PROXIMAL_GENES","CLOSEST_GENE","RANK","ISSUPER")
data$RANK2 <- c(1:nrow(data))

# Isolate only the SE data
se_data <- subset(data, data$ISSUPER==1)

# Mark the superenhancers of the 10 genes with highest SE rank
gene_data <- se_data[1:10,]
# Plot additional genes
additional_genes <- subset(se_data, se_data$CLOSEST_GENE %in% c("IKZF1","IRF4","MYC"))
gene_data <- rbind(gene_data, additional_genes)

# Isolate the SE with the highest rank for each one of the genes
gene_enh_data <- gene_data[!duplicated(gene_data$CLOSEST_GENE),]
# Get a vector with the info you need for the plotting
x <- as.vector(as.character(gene_enh_data$CLOSEST_GENE))
y <- as.vector(as.character(gene_enh_data$RANK))
enh_labels <- paste(x, sprintf("(%s)",y))

# Get the number of total enhancers and SEs
enhancer_no = dim(data)[1]
se_no=dim(se_data)[1]

# Plot the SE graph
out=strsplit(input,"[.]")[[1]][1] 
output <- paste(out,"SE_graph.pdf",sep="_")

# Get the maximum value for H3K27ac signal
max_value <- max(data$VALUE)

pdf(sprintf("%s", output), useDingbats=F)
ggplot(data, aes(x=RANK2,y=VALUE)) + geom_line(color='red') + geom_point(size=1) + geom_point(data=se_data, aes(x=RANK2,y=VALUE), color="red", size=1) + theme_bw() + xlab("Enhancers Ranked") + ylab("ChIP-seq signal") + geom_text_repel(data=gene_enh_data, aes(x=RANK,y=VALUE), label=enh_labels) + annotate("text", x = 11000, y = max_value + 20000, label = sprintf("Superenhancers detected: %s",se_no), color="red") + annotate("text", x = 11000, y = max_value, label = sprintf("Total enhancers detected: %s", enhancer_no), color="black") + scale_x_reverse()
dev.off()
