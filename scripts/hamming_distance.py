# Given two strands of DNA, find number of point mutations

def hdist(f_strand, s_strand):
    return sum(1 for a, b in zip(f_strand, s_strand) if a != b)

print(hdist(s1, s2))

