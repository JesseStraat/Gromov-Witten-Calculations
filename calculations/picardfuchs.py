from sympy import diff, symbols, Function, latex, expand, collect, Derivative

# This scripts calculates the Picard-Fuchs equation of the quintic mirror

z = symbols('z')
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
def organise_result(f_sym, deg: int=4):
    derivs = [Derivative(f(z), z, i) for i in range(deg+1)]
    return collect(expand(f_sym), derivs)

# Outputs theta^4 - 5 z \prod_{i=1}^4 (5*theta + i) applied to function
def result(f_sym):
    return theta(f_sym,4) - 5*z*part_res(f_sym,4)

# Print result in LaTeX format
print(latex(organise_result(result(f(z)))))