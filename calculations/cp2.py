from math import comb

N = [1]

def calculate_term(d,k,l):
    return N[k-1]*N[l-1]*k**2*l*(l*comb(3*d-4,3*k-2) - k*comb(3*d-4,3*k-1))


def calculate_next():
    d = len(N)
    total = 0
    for k in range(1,d+1):
        l = d+1 - k
        total += calculate_term(d+1,k,l)
    return total



for _ in range(50):
    N.append(calculate_next())
    d = len(N)
    print(f"N({d}) = {N[d-1]}")