#!/usr/bin/Rscript

# Load the required libraries
library(venneuler)

# Create the Venn diagram
pdf("venn.pdf", useDingbats=FALSE)
v <- venneuler(c(A=74000, B=63550, "A&B"=49098))
plot(v)
dev.off()
