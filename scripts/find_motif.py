# Given a FASTA file containing a DNA dataset, prompt user to search for a motif
# Output indices of all locations of searched motif

import sys

with open(sys.argv[1]) as fasta_file:
    fasta_string = fasta_file.read()

def search_dna(motif, dna_strand):    
        
    motif_indices = []
    search_start  = 0

    while search_start < len(dna_strand):

        if dna_strand.find(motif, search_start) != -1:
        
            curr_found_ind = dna_strand.find(motif, search_start)
            motif_indices.append(curr_found_ind + 1) 
            search_start   = curr_found_ind + 1
        
        else:
            return motif_indices
    

print("Enter motif to search for:")
motif_searched  = input()

motif_locations = search_dna(motif_searched, fasta_string)

if len(motif_locations) != 0:
    print("Motif found at bp location(s):")
    print(*motif_locations)

else:
    print("Motif not found in given strand")
