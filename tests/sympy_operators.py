from sympy import sin, cos, Plot, Basic, Symbol, sympify, zeronm, lambdify,\
     symbols
from sympy.abc import x, y, z, t

from numpy import arange

def default_space_vars( vars ):
    if vars is None:
        vars = [x, y, z]
    return vars

def grad(f, vars=None):
    vars = default_space_vars( vars )
    f = sympify( f )
    div = zeronm(len(vars), 1)
    for i in range(len(vars)):
        div[i, 0] = f.diff(vars[i])
    return div

def div(field, vars=None):
    vars = default_space_vars( vars )
    field = list(field)
    assert len(field) == len(vars)
    r = 0
    for A_i, x_i in zip(field, vars):
        r += sympify( A_i ).diff(x_i)
    return r

def laplace(f, vars=None):
    return div(grad(f, vars), vars)

def boundary(f, vars):
    lap = laplace(f, vars)
    l = lambdify(lap, vars)
    a = 5.
    for x in arange(-a, a, 0.1):
        y = -a
        print x, y, l(x, y)

if __name__ == '__main__':
    f = sin(x)*cos(y)
    boundary(f, [x, y])
    #Plot(f, [x, -10, 10], [y, -10, 10])
