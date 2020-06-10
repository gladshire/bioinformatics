 # Given a text file in fasta format containing homologous seqs of DNA
 # Output consensus strand
nucleotides = ['A', 'C', 'G', 'T']
import sys
sys.path.append('/home/miles/programming/python_/bioinformatics')

import bioi

nucleotides = ['A', 'C', 'G', 'T']

def prof(seq_list):

    col = len(seq_list[0])
    row = len(seq_list)

    profile = {}

    for i, base in enumerate(nucleotides):
        base_list = []
        for j in range(col):
            curr_column = [seq_list[k][j] for k in range(row)]
            base_list.append(curr_column.count(base))
       
        profile[base] = base_list
    
    return profile
            
    
    
def cons(seq_list):

    ind_nt = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    lengt  = len(seq_list[0])
    
    seq_profile = prof(seq_list)
    cons_arr    = []

    for j in range(lengt):
        curr_column = [seq_profile[base][j] for base in nucleotides]
        cons_arr.append(ind_nt[curr_column.index(max(curr_column))])

    return ''.join([cons_arr[i] for i in range(lengt)])

with open(sys.argv[1], 'r') as fasta_file:
    fasta_data = bioi.fparse(fasta_file.read())
    seq_list   = list(fasta_data.values())

consensus = cons(seq_list)
profile   = prof(seq_list)
print(consensus)
for i, base in enumerate(nucleotides):
    print('{}: '.format(base))
    print(*profile[base])
    








