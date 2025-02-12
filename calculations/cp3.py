# %%
from math import comb
import numpy as np
import pandas as pd
from fractions import Fraction

np.set_printoptions(suppress=True)

def combb(n,r):
    if n<0 or r<0:
        return 0
    else:
        return comb(n,r)

def calc_GW(max_degree: int) -> np.array:
    # max_degree: maximum degree
    # Output: N[d][n] is N_0(n;d)
    N = np.zeros((max_degree+1,2*max_degree+1),Fraction)
    
    # Initial values
    N[0][0] = Fraction(1,1)
    N[1][2] = Fraction(1,1)
    
    for d in range(1,len(N)):
        # Calculating N_0(1;d)
        N[d][1] = d*N[d][2]
        for k in range(1,d):
            l = d-k
            N[d][1] += k*N[k][0]*N[l][1]*(l**2*combb(4*d-4,4*k-2)-k**2*combb(4*d-4,4*k))
        
        # Calculating N_0(0;d)
        N[d][0] = 2*d*N[d][1]
        for k in range(1,d):
            l = d-k
            N[d][0] += k**2*N[k][0]*N[l][0]*(l*combb(4*d-3,4*k-1)-k*combb(4*d-3,4*k))
        
        # Calculating all other N_0(n;d)
        for n in range(3,2*d+1):
            N[d][n] = N[d][n-1]/d
            for k in range(1,d):
                l = d-k
                for a in range(0,n-1):
                    b = n-2-a
                    N[d][n] -= combb(n-2,a)*k*N[k][a]*N[l][b+1]*(l**2*combb(4*d-2*n,4*k-2*a-2)-k**2*combb(4*d-2*n,4*k-2*a))/d
        
        # Calculating N_0(2;d+1)
        if d+1 <= max_degree:
            for k in range(1,d+1):
                l = d+1-k
                N[d+1][2] += 2*k**2*l*combb(4*d-1,4*k-2)*N[k][1]*N[l][1] - k**2*N[k][0]*N[l][2]*(k*combb(4*d-1,4*k)+l*combb(4*d-1,4*k-1))
    return N


def GW_table(max_degree: int) -> pd.DataFrame:
    N = calc_GW(max_degree)
    table = pd.DataFrame(N,
                         index = range(max_degree+1),
                         columns = range(2*max_degree+1))
    return table

GW_table(6)
# %%
