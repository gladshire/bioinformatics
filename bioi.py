 # bioi Bioinformatics library for python 3.7.3
 # By Miles Woodcock-Girard

#######################################################################

nucleotides = ['A', 'C', 'G', 'T']
pairs       = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

 # Generate a complementary strand to a given strand of DNA, in 5' - 3'
def comp(dna_strand):
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '\n': ''}

    return ''.join([pairs[base] for base in dna_strand])


 # Generate a complementary strand to a given strand of DNA, in 3' - 5'
def rcomp(dna_strand):
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '\n': ''}
    
    return ''.join([pairs[base] for base in dna_strand[::-1]])


 # Generate a strand of DNA of a given bp count
def gen_dna(dna_length):
    
    return ''.join([random.choice('ATCG') for base in range(n)])

 
 # Create a dictionary of DNA datasets from a fasta file
def fparse(dna_string):
    results = {}
    fasta_list = dna_string.strip().split('>')
    for i in fasta_list:
        if len(i) == 0:
            continue
    
        parts = i.split()
        label = parts[0]
        bases = ''.join(parts[1:])
        results[label] = bases        
    
    return results


 # Calculate the GC-content percentage of a given strand of DNA
def gcp(dna_strand):
    gc_count = dna_strand.count('G') + dna_strand.count('C')
    gc_percent = gc_count / len(dna_strand) * 100

    return gc_percent


 # Transcribe a given strand of DNA
def transcribe(dna_strand):
    rna_strand = dna_strand.replace('T', 'U')
    
    return rna_strand


 # Search for a motif in a given strand and return indices of all occurrences
def msearch(motif, dna_strand):    
    motif_indices = []
    search_start  = 0
    while search_start < len(dna_strand):
        if dna_strand.find(motif, search_start) != -1:
            curr_found_ind = dna_strand.find(motif, search_start)
            motif_indices.append(curr_found_ind + 1) 
            search_start   = curr_found_ind + 1    
        else:
            
            return motif_indices


 # Create a dictionary of base pair count in given strand of DNA
def ntcount(dna_strand):
    nucleotide_count = {}
    nucleotides      = ['A', 'C', 'G', 'T']
    for base in nucleotides:
        nucleotide_count[base] = dna_strand.count(base)

    return nucleotide_count


 # Find indices/lengths of all restriction sites within specified length range
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

 
 # Given an original and mutant strand of DNA return number of point mutations
def hdist(f_strand, s_strand):
    
    return sum(1 for a, b in zip(f_strand, s_strand) if a != b)


 # Translate a given strand of mRNA into single-letter protein string
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


 # Given a list of homologous DNA sequences, return profile of each nucleotide
def prof(seq_list):

    col = len(seq_list[0])
    row = len(seq_list)

    profile = {}

    for i, base in enumerate(nucleotides):
        base_list = []
        for j in range(col):
            curr_column = list([seq_list[k][j] for k in range(row)])
            base_list.append(curr_column.count(base))
       
        profile[base] = base_list
    
    return profile
            
    
 # Given a list of homologous DNA sequences, return the consensus sequence
def cons(seq_list):

    ind_nt = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    lengt  = len(seq_list[0])
    
    seq_profile = prof(seq_list)
    cons_arr    = []

    for j in range(lengt):
        curr_column = list([seq_profile[base][j] for base in nucleotides])
        cons_arr.append(ind_nt[curr_column.index(max(curr_column))])

    return ''.join([cons_arr[i] for i in range(lengt)])


 # Given a FASTA file, parse and create adjacency list of specified overlap
def alist(fasta_dict, overlap):

    adjac_list = []

    for seq_key in fasta_dict:
        curr_key = seq_key
        for seq_key in fasta_dict:
            if fasta_dict[curr_key] == fasta_dict[seq_key]:
                continue
            if fasta_dict[curr_key][-overlap:] == fasta_dict[seq_key][:overlap]:
                adjac_list.append("{} {}".format(curr_key, seq_key))
  
    return adjac_list
