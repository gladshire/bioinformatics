# Given a text file containing a string of DNA, count A C G T nucleotides and display in order

import sys

with open(sys.argv[1], "r") as dna_file:
    dna_strand = dna_file.read()

def nucleotide_count(dna_strand):
    nucleotide_count = {}
    nucleotides      = ['A', 'C', 'G', 'T']
    
    for base in nucleotides:
        nucleotide_count[base] = dna_strand.count(base)

    return nucleotide_count

print(nucleotide_count(dna_strand))
