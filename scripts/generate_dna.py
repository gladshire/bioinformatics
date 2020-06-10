

import random

def generate_dna(n):
    return ''.join([random.choice('ATCG') for base in range(n)])

print("Specify number of bases for strand to be generated")
num_bases  = int(input())

label = ">sequence_(Training)_(" + str(num_bases) + "bp)\n"

dna_string = label + generate_dna(num_bases) 

with open("dna_data.fa", "w") as dna_file:
    dna_file.write(dna_string)
