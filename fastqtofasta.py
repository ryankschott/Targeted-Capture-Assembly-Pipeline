from Bio import SeqIO
import sys
import os

SeqIO.convert(sys.argv[1], "fastq", sys.argv[2], "fasta")
