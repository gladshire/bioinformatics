# Given DNA data, output locations and lengths of all palindromic sequences between lengths 4 and 12

import sys
sys.path.append('home/programming/python_/bioinformatics')

import bioi

with open(sys.argv[1]) as dna_file:
    dna_data = dna_file.read()

dna_parsed = bioi.fparse(dna_data)
dna_strand = list(dna_parsed.values())[0] 

def rsites(dna_strand, lower_limit, upper_limit):
    limit_list  = range(lower_limit, upper_limit + 1) 
    palindromes = []
    for i, base in enumerate(dna_strand):

        for seq_length in limit_list:
            
            if i + seq_length > len(dna_strand):
               break 

            dna_seq  = dna_strand[i:i + seq_length]
            comp_seq = bioi.rcomp(dna_seq)            
            
            if dna_seq == comp_seq:
                
                p_element  = i + 1
                p_length   = seq_length
                palindromes.append((p_element, p_length))  

    return palindromes

results = rsites(dna_strand, 4, 12)

print("\n".join([' '.join(map(str, r)) for r in results]))
