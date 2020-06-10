 # Given n (number of months) and k (number of rabbit pairs produced in a litter, every month)
 # Return the number of rabbit pairs present after n months

def recursive_rabbits(n, k):
    if n == 0:
        return 0
    if n == 1: 
        return 1
    else:
        return recursive_rabbits(n - 1, k) + k * recursive_rabbits(n - 2, k) 

print(recursive_rabbits(29, 4))
