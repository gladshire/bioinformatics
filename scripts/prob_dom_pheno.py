# Given a population with specified number of homozygous D, heterozygous, homozygous R
# Output probability of offspring of random parents harboring dominant allele

def calc_dom(homo_dom, hetero, homo_rec):
    
    pop_total  = homo_dom + hetero + homo_rec
    
    hdhd = (homo_dom / pop_total) * ((homo_dom - 1) / (pop_total - 1)) 
    
    hdhe = (homo_dom / pop_total) * (hetero / (pop_total - 1))
    hehd = (hetero / pop_total) * (homo_dom / (pop_total - 1))

    hehe = (hetero / pop_total) * ((hetero - 1) / (pop_total - 1)) * 0.75

    hdhr = (homo_dom / pop_total) * (homo_rec / (pop_total - 1))
    hrhd = (homo_rec / pop_total) * (homo_dom / (pop_total - 1))

    hehr = (hetero / pop_total) * (homo_rec / (pop_total - 1)) * 0.5
    hrhe = (homo_rec / pop_total) * (hetero / (pop_total - 1)) * 0.5

    return hdhd + hdhe + hehd + hehe + hdhr + hrhd + hehr + hrhe 


# print(calc_dom(27, 26, 22))

homo_dom = int(input('Number of homozygous dominant individuals:'))
hetero   = int(input('Number of heterozygous individuals:'))
homo_rec = int(input('Number of homozygous recessive individuals:'))

print(calc_dom(homo_dom, hetero, homo_rec))
