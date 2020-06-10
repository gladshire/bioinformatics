#Transcribe a string of DNA from a text file, and write a text file with the output

import sys

with open(sys.argv[1], "r") as DNAin
    DNA = DNAin.read()

RNA = DNA.replace('T', 'U')

with open("DNA_transcribed.txt", "w") as RNAout
    RNAout.write(RNA)
