# Given a text file containing a string of DNA, create a file containing its complement

import sys
sys.path.append('/Users/macbook/bioinformatics/')

import bioi

def generate_complement_strand(dna_strand):
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '\n': ''}
    return ''.join([pairs[base] for base in dna_strand[::-1]])

if len(sys.argv) == 2:
    with open(sys.argv[1], 'r') as dna_file_in:
        dna_data = dna_file_in.read()
    
    dna_parsed = bioi.fparse(dna_data)
    dna_data = list(dna_parsed.values())[0]

    dna_comp = generate_complement_strand(dna_data)

    with open("DNA_complementary.fa", "w") as dna_file_out: 
        dna_file_out.write(dna_comp)

