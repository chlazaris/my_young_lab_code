# Goal:
# Take FASTA and produce barcode plots.

# Use:
# Rscript ./barcode_plot.r protein_sequences.fa barcode_plots.pdf

################################################################################
# Import libraries

# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
library(Biostrings)

################################################################################


################################################################################
# Define amino acid alphabet
# - Order of letters is what is displayed in output.

#alphabet <- c("N", "H", "C", "M", "L", "I", "Q", "T", "S", "P", "G", "Y", "W", "F", "D", "E", "V", "A", "R", "K")

#Define the alphabet having the negative AAs first, then positive, uncharged, hydrophobic and special cases
alphabet <- rev(c("G","E","P","D","R","H","K","C","S","T","N","Q","A","V","I","L","M","F","T","W"))


################################################################################
# Process command line arguments

args <- commandArgs(TRUE)

# Check for the number of arguments
if (length(args)!=2) {
	print("Please provide the correct number of arguments...")
	print("USAGE: Rscript barcode_plot.r [protein FASTA sequence] [barcode plot (PDF)]")
	q(save="no")
}

# Read the arguments
fasta <- args[1]
barcodePlots <- args[2]

################################################################################
# Import sequences

stringSet <- readAAStringSet(fasta)
print(stringSet)


################################################################################
# Make barcode plots

pdf(barcodePlots)
for (i in 1:length(stringSet)) {

	# Retrieve string
	string <- as.character(stringSet[[i]])
	name <- names(stringSet)[i]
	count <- nchar(string)


	# Make distribution matrix
	distMat <- matrix(0, length(alphabet), nchar(string))
	row.names(distMat) <- alphabet

	# Populate distribution matrix
	for (j in 1:length(alphabet)) {
		aa <- alphabet[j]
		for (k in 1:nchar(string)) {
			residue <- substring(string, k, k)
			if (residue == aa) {
				distMat[j,k] <- 1
			}
		}
	}

	# Make plot
	palette <- colorRampPalette(c("white","black"))(n=2)
	heatmap(distMat,
		main=name, 
		labCol=count,
		col=palette, Rowv=NA, Colv=NA, scale="none")
}
dev.off()




