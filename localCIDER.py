#!/bin/python

# Import the required modules, libraries etc.
import numpy
import scipy
import matplotlib
import localcider
from localcider.sequenceParameters import SequenceParameters
from Bio import SeqIO
import sys
import subprocess
import pandas as pd

# Function to calculate running average over a certain window
def running_mean(x, N):
    cumsum = numpy.cumsum(numpy.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)
    
################################################################################
# 1. Define inputs
prots = sys.argv[1]
window = int(sys.argv[2])

# Check number of arguments
if len(sys.argv)!=3:
    print("Provide the correct number of arguments")
    print("USAGE: python2 localCIDER.py [FASTA file] [w]")
    sys.exit(1)

# Read in the .fasta sequence
fasta_sequences = SeqIO.parse(open(prots),'fasta')

################################################################################
# 2. Process protein sequences

for fasta in fasta_sequences:

  name, sequence = fasta.id, str(fasta.seq)
  size = len(sequence)
  SeqOb = SequenceParameters(sequence)

  linFCR = SeqOb.get_linear_FCR(blobLen=window)
  linNCPR = SeqOb.get_linear_NCPR(blobLen=window)
  linSigma = SeqOb.get_linear_sigma(blobLen=window)
  linHydropathy = SeqOb.get_linear_hydropathy(blobLen=window)
  linOmegaList = list(SeqOb.get_Omega_sequence().replace("O","0").replace("X","1"))
  linOmegaBinary = numpy.asarray(linOmegaList, dtype=numpy.float32)
  linOmegaNoPad = running_mean(linOmegaBinary, window)
  linOmega = numpy.pad(linOmegaNoPad, ((window-1)/2, window-((window-1)/2)-1), "constant")

  combo = numpy.column_stack((linFCR[1], linNCPR[1], linSigma[1], linHydropathy[1], linOmega))
  output = pd.DataFrame(data = combo, columns=["FCR", "NCPR", "Sigma", "Hydropathy", "Omega"])
  outname = name + "_blob" + str(window) + ".csv"
  output.to_csv(outname, index=False)
#  filename = name + "_blob" + str(window) + ".txt"
#  numpy.savetxt(filename, output, delimiter="\t")

