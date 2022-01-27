#install the core bioconductor packages, if not already installed
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("hugene10stv1cdf", "hugene10stv1probe", "hugene10stprobeset.db", "hugene10sttranscriptcluster.db"))

# install additional bioconductor libraries, if not already installed
#biocLite("GEOquery")
#biocLite("affy")
#biocLite("gcrma")
#biocLite("hugene10stv1cdf")
#biocLite("hugene10stv1probe")
#biocLite("hugene10stprobeset.db")
#biocLite("hugene10sttranscriptcluster.db")

#Load the necessary libraries
library(GEOquery)
library(affy)
library(limma) # Used to plot densities and for differential expression
#library(gcrma)
library(hugene10stv1cdf)
library(hugene10stv1probe)
library(hugene10stprobeset.db)
library(hugene10sttranscriptcluster.db)

#Set working directory for download
#setwd("/lab/solexa_young/lazaris/raw_data/microarray")

#Download the CEL file package for this dataset (by GSE - Geo series id)
#getGEOSuppFiles("GSE31365")

#Unpack the CEL files
#setwd("/lab/solexa_young/lazaris/raw_data/microarray/GSE31365")
#untar("GSE31365_RAW.tar", exdir="data")
#cels = list.files("data/", pattern = "CEL")
#sapply(paste("data", cels, sep="/"), gunzip)
#cels = list.celfiles("/lab/solexa_young/lazaris/raw_data/microarray/GSE31365/data/")
#cels

#setwd("/lab/solexa_young/lazaris/raw_data/microarray/GSE31365/data")
raw.data=ReadAffy() #From bioconductor

#perform RMA normalization (I would normally use GCRMA but it did not work with this chip)
data.rma.not_norm=rma(raw.data, normalize=F, background=F)
data.rma.norm=rma(raw.data, normalize=T, background=T)

# Get rma for normalized and not normalized
rma.nonorm <- exprs(data.rma.not_norm) 
rma.norm <- exprs(data.rma.norm)

#Get the important stuff out of the data - the expression estimates for each array
pdf(file="DensityNoNorm.pdf", w=6, h=6)
plotDensities(rma.nonorm, main="Arrays Not Normalized")
dev.off()

pdf(file="DensityNorm.pdf", w=6, h=6)
plotDensities(rma.norm, main="Arrays Not Normalized")
dev.off()

# Boxplot of intensity values before normalization
pdf(file="BoxplotNoNorm.pdf", w=6, h=6)
boxplot(rma.nonorm, main="Arrays Not Normalized")
dev.off()

# Boxplot of intensity values after normalization
pdf(file="BoxplotNorm.pdf", w=6, h=6)
boxplot(rma.norm, main="Arrays Normalized")
dev.off()

#Format values to 5 decimal places
rma=format(rma.norm, digits=5)

#Map probe sets to gene symbols or other annotations
#To see all available mappings for this platform
#ls("package:hugene10stprobeset.db") #Annotations at the exon probeset level
#ls("package:hugene10sttranscriptcluster.db") #Annotations at the transcript-cluster level (more gene-centric view)

#Extract probe ids, entrez symbols, and entrez ids
probes=row.names(rma)
Symbols = unlist(mget(probes, hugene10sttranscriptclusterSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hugene10sttranscriptclusterENTREZID, ifnotfound=NA))

#Combine gene annotations with raw data
rma2=cbind(probes,Symbols,Entrez_IDs,rma)

#Write RMA-normalized, mapped data to file
write.table(rma2, file = "rma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Get the phenotypic data
treatment <- c("DMSO","DMSO","JQ1","JQ1")

eset <- rma(raw.data, normalize=T, background=T)
design <- model.matrix(~factor(treatment))
colnames(design) <- c("DMSO", "JQ1")

fit <- lmFit(eset, design)
fit <- eBayes(fit)
res <- topTable(fit, number=Inf, adjust.method="BH", coef=1)
write.table(res, 'diff_exp.txt', sep="\t", quote=F)
