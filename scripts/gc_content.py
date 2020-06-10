# Given a text file containing DNA segments in FASTA format, parse and calculate GC-content for each segment

import sys

with open(sys.argv[1], 'r') as fasta_file:
    fasta_raw_data = fasta_file.read()


def fasta_parse(raw_data):
    results = {}
    fasta_list = raw_data.strip().split('>')

    for i in fasta_list:

        if len(i) == 0:
            continue

        parts = i.split()
        label = parts[0]
        bases = ''.join(parts[1:])
    
        results[label] = bases
        
    return results
      
def gc_content(dna_strand):
    gc_count = dna_strand.count('G') + dna_strand.count('C')
    
    gc_percent = gc_count / len(dna_strand) * 100

    return gc_percent


fasta_dict = fasta_parse(fasta_raw_data)
for label in fasta_dict:
    gc_percent = gc_content(fasta_dict[label])
    print(label)
    print("{:.5f}".format(gc_percent))
    

