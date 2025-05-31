from sympy import diff, symbols, Function, latex, expand, collect, Derivative, Matrix, zeros, simplify, fraction
from os import system

# This scripts calculates the Picard-Fuchs equation of the quintic mirror

z, t = symbols('z t')
f = Function('f')

# theta = z d/dz. This applies it n times
def theta(f_sym, n: int=1):
    result = f_sym
    for _ in range(n):
        result = z*diff(result,z)
    return result

# Outputs (5*theta + i) applied to function
def fivethetaplus(f_sym, i: float):
    return 5*theta(f_sym)+i*f_sym

# Performs product of (5*theta + i) applied to function
def part_res(f_sym, n: int=4):
    result = f_sym
    for i in range(1,n+1):
        result = fivethetaplus(result,i)
    return result

# Collects result in terms of derivatives
def organise_result(expr, f_sym, deg: int=4):
    derivs = [Derivative(f_sym, z, i) for i in range(deg+1)]
    return collect(expand(expr), derivs)

# Outputs theta^4 - 5 z \prod_{i=1}^4 (5*theta + i) applied to function
def PF_operator(f_sym):
    return theta(f_sym,4) - 5*z*part_res(f_sym,4)

# Extracts the coefficients of the derivatives from expression
def extract_coeffs(expr, f_sym, deg: int=4, var = z):
    derivs = [Derivative(f_sym, var, i) for i in range(deg+1)]
    derivs_list = collect(expand(expr), derivs, evaluate=False)
    result = []
    for i in range(deg+1):
        result.append(derivs_list.get(derivs[i],0))
    return result

# Collects info in a matrix. Works as derivative on expressions that satisfy PF
def PF_matrix(f_sym, deg: int=4):
    deriv_coeffs = extract_coeffs(PF_operator(f_sym),f_sym,deg)
    result = zeros(deg,deg)
    for i in range(deg-1):
        result[i, i+1] = 1
    for j in range(deg):
        result[deg-1, j] = -deriv_coeffs[j]/deriv_coeffs[deg]
    return result

Omega = Function('Ω')
S = Function('S_0')
z_func = Function('z')
C = Function('C')
# Extracts matrix for the change of coordinates Ansatz
def coord_trans_matrix():
    result = []
    exprlist = [
        Omega(z_func(t))/S(t),
        diff(Omega(z_func(t))/S(t),t),
        C(z_func(t))**(-1)*diff(Omega(z_func(t))/S(t),t,2),
        diff(C(z_func(t))**(-1)*diff(Omega(z_func(t))/S(t),t,2),t)
    ]
    for expr in exprlist:
        result.append(extract_coeffs(expr, Omega(z_func(t)), var=z_func(t))[:-1])
    return Matrix(result)

# Transform the connection matrix under change of coords
def changed_PF_matrix(f_sym):
    trans_matrix = coord_trans_matrix()
    conn_matrix = PF_matrix(f_sym)
    result = (diff(trans_matrix,t) + diff(z_func(t),t)*trans_matrix*conn_matrix)*trans_matrix.inv()
    return simplify(result)

# Outputs list of differential equations
def get_PF_eqs(f_sym):
    output = []
    conn_matrix = changed_PF_matrix(f_sym)
    for i in range(4):
        for j in range(4):
            if (i,j) not in [(0,1),(1,2),(2,3)]:
                if conn_matrix[i,j] != 0:
                    output.append(fraction(conn_matrix[i,j])[0])
    return output



# Print result in LaTeX format
with open("picardfuchs.tex", "w") as f_out:
    f_out.write("""
\\documentclass{standalone}
\\usepackage{amsmath}
\\begin{document}
\\(\\begin{gathered}
""" +  ",\\\\\n".join(latex(eq) + "=0" for eq in get_PF_eqs(f(z))) + """.
\\end{gathered}\\)
\\end{document}
""")
system("lualatex picardfuchs.tex")