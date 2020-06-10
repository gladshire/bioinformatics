# Given an mRNA strand, translate and output single-letter protein string

import sys

def strans(mrna_strand):
    codon_protein_sym = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
                         'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
                         'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
                         'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W',
                         'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
                         'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                         'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                         'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                         'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
                         'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                         'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                         'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                         'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
                         'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                         'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                         'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
                         '\n': ''}
    
    rna_codons = []
    peps       = []
    curr_ind   = 0
    start_ind  = 0

    while curr_ind < len(mrna_strand):
        curr_ind = start_ind + 3
        rna_codons.append(mrna_strand[start_ind:curr_ind])
        start_ind = curr_ind
    
    protein_seq = ''.join([codon_protein_sym[codon] for codon in rna_codons])
    
    for prot_str in protein_seq.split('Stop'):
        if len(prot_str) != 0:
            peps.append(prot_str)

    return peps

with open(sys.argv[1], 'r') as rna_file:
    mrna = rna_file.read()
    mrna.strip('\n')

print(strans(mrna))
