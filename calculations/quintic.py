from sympy import symbols, Symbol, bell
from sympy.polys.polytools import Poly
from math import factorial
from fractions import Fraction
from collections.abc import Callable

class TaylorSeries:
    # Creates a Taylor series, allowing for polynomial representations
    def __init__(self, nthterm: Callable[[int],Fraction]):
        self.nthterm = nthterm
    
    def as_polynomial(self, degree: int, variable: Symbol) -> Poly:
        polynomial = 0
        for n in range(degree+1):
            polynomial += self.nthterm(n)*variable**n
        return Poly(polynomial, variable)

def truncate_poly(f: Poly, degree: int) -> Poly:
    # Truncates a polynomial to some degree
    result = {}
    monomials = f.monoms()
    coeffs = f.coeffs()
    for i in range(len(monomials)):
        if all(d <= degree for d in monomials[i]):
            result[monomials[i]] = coeffs[i]
    return Poly.from_dict(result, *f.gens)

def division(f: TaylorSeries, g:TaylorSeries, degree: int, variable: Symbol) -> Poly:
    # Divides two Taylor series up to some degree
    f_poly = f.as_polynomial(degree=degree, variable=variable)
    g_poly = g.as_polynomial(degree=degree, variable=variable)
    g_const = g_poly.coeff_monomial(1)
    g_new = g_poly - g_const
    result = 0
    for n in range(0,degree+1):
        result += (-1)**n*Fraction(1,g_const**n)*g_new**n
    return truncate_poly(f_poly*result, degree)

def exponentiation(f: Poly, degree: int) -> Poly:
    # Exponentiates some polynomial up to some degree
    result = 0
    for n in range(0,degree+1):
        result += f**n*Fraction(1,factorial(n))
    return truncate_poly(result, degree)

def lagrange_bell_invert(f: Poly, old_variable: Symbol, new_variable: Symbol) -> Poly:
    # An implementation of the Lagrange inversion theorem using Bell polynomials
    deg = f.degree(gen=old_variable)
    f_list = [f.as_expr().coeff(old_variable,i)*factorial(i) for i in range(deg+1)]
    if deg >= 1 and f_list[1] != 0:
        if f_list[0] == 0:
            fhat_list = [Fraction(f_list[k],k*f_list[1]) for k in range(1,len(f_list))]
            g = new_variable*Fraction(1,f_list[1])
            for n in range(2,len(f_list)):
                res = 0
                for k in range(1,n):
                    funct_tup = tuple(fhat_list[i] for i in range(1,n-k+1))
                    res += (-1)**k*Fraction(factorial(n+k-1),factorial(n-1))*bell(n-1,k,funct_tup)
                g += res*Fraction(1,f_list[1]**n*factorial(n))*new_variable**n
            return g
        else:
            raise Exception("This function has a nonzero constant term, so the algorithm doesn't apply.")
    else:
        raise Exception("This function isn't invertible.")


if __name__ == "__main__":
    b, q = symbols('b q')
    
    omega0 = TaylorSeries(nthterm = lambda n : Fraction(factorial(5*n),(factorial(n))**5))

    def sum_of_recips(n: int) -> Fraction:
        result = 0
        for i in range(n+1,5*n+1):
            result += Fraction(1,i)
        return result
    psi = TaylorSeries(nthterm = lambda n : Fraction(5*factorial(5*n),(factorial(n))**5)*sum_of_recips(n))

    q_func = b*exponentiation(division(psi,omega0,5,b), 5)
    b_func = lagrange_bell_invert(q_func,b,q)

    print(b_func)