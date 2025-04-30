# %%
from math import comb
import numpy as np
import pandas as pd
from fractions import Fraction # GW invariants are rational, so this is more precise

def combb(n,r):
    # A binomial that returns zero whenever arguments are illegal
    if n<0 or r<0:
        return 0
    else:
        return comb(n,r)

def calc_GW(max_degree: int) -> np.array:
    # max_degree: maximum degree
    # Output: N[d] is N_0(d)
    N = np.zeros((max_degree+1),Fraction)
    
    # Initial values
    N[0] = np.nan
    N[1] = Fraction(1,1)
    
    for d in range(2,len(N)):
        # Calculating N_0(d)
        N[d] = 0
        for k in range(1,d):
            l = d-k
            N[d] += k**2*l*N[k]*N[l]*(l*combb(3*d-4,3*k-2)-k*combb(3*d-4,3*k-1))
    return N


def GW_table(max_degree: int) -> pd.DataFrame:
    # Turn GWs into table
    N = calc_GW(max_degree)
    table = pd.DataFrame(N[1:],
                         index = range(1,max_degree+1),
                         columns = ["N_0(d)"])
    table.columns.name = "d"
    return table

GW_table(8)
# %%