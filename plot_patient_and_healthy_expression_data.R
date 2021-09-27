#!/bin/Rscript

# Load the required libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(forcats)

# Check the number of command-line arguments
input <- commandArgs(trailingOnly=TRUE)

# Get the data
data <- read.table(sprintf("%s", input), header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)

# Isolate the columns that you need
df <- (data[,2:ncol(data)])
df2 <- melt(df, id.vars = c("gene_name"))
head(df2)

# Add the column names
colnames(df2) <- c("gene_name","sample","TPM")
df3 <- df2 %>% mutate(log2_TPM=log2(TPM+1))
df3$sample <- as.factor(df3$sample)
df3$sample <- fct_relevel(df3$sample, "Donor977","Donor978","Donor979","Donor981","Donor985","Donor988","Donor991","Donor994","Donor996","Donor999","Donor1000","Donor1001","Donor1002","Donor1004","Donor1005","Donor1007","Donor1011","Donor1012","Donor1014","Donor1016","Donor1017","Donor1022","Donor1026","Donor1028","Donor1033","Donor1041","Donor1044","Donor1045","NormalBM1","NormalBM2","NormalBM3","NormalBM4","NormalBM5","NormalBM6","NormalBM7","NormalBM8","NormalBM9","NormalBM10","NormalBM11","NormalBM12","NormalBM13")

# Isolate Ikaros for each one
df4 <- subset(df3, df3$gene_name=='IKZF1')
head(df4)

pdf("TPM_values.pdf")
ggplot(df3, aes(x=sample, y=log2_TPM)) + geom_boxplot(outlier.shape=NA) + geom_point(data=df4, aes(x=sample, y=log2_TPM), show.legend=c("IKZF1"), color="red") + scale_y_continuous(limits = c(0, 6.0)) + xlab("Sample") + ylab("log2(TPM + 1)") + theme(axis.text.x = element_text(angle = 90))
dev.off()
