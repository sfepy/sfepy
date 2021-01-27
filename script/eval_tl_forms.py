"""
Operators present in the FE discretization of hyperelastic terms in the total
Lagrangian formulation.
"""
from __future__ import print_function
from __future__ import absolute_import
import sympy as s

def main():
    u1, u2, u3 = s.symbols(['u1', 'u2', 'u3'], commutative=False)
    u = s.Matrix([[u1], [u2], [u3]])

    g1, g2, g3 = s.symbols(['g1', 'g2', 'g3'], commutative=False)
    gc = s.Matrix([[g1], [g2], [g3]])

    aux = s.symbols(['u11', 'u12', 'u13', 'u21', 'u22', 'u23', 'u31', 'u32', 'u33'],
                    commutative=False)
    u11, u12, u13, u21, u22, u23, u31, u32, u33 = aux

    aux = s.symbols(['f11', 'f12', 'f13', 'f21', 'f22', 'f23', 'f31', 'f32', 'f33'],
                    commutative=False)
    f11, f12, f13, f21, f22, f23, f31, f32, f33 = aux

    print(gc)

    ## z = s.zeros((3, 1))

    ## g = s.Matrix([[gc, z, z],
    ##               [z, gc, z],
    ##               [z, z, gc]])

    g = s.Matrix([[g1, 0, 0],
                  [g2, 0, 0],
                  [g3, 0, 0],
                  [0, g1, 0],
                  [0, g2, 0],
                  [0, g3, 0],
                  [0, 0, g1],
                  [0, 0, g2],
                  [0, 0, g3]])

    print(g)
    print(g * u)


    h = s.Matrix([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 1, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 1],
                  [0, 1, 0, 1, 0, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0, 1, 0, 0],
                  [0, 0, 0, 0, 0, 1, 0, 1, 0]])

    print(h)

    print('linear part:')
    print(h * g * u)
    print(h * g)

    a = s.Matrix([[u11, 0, 0, u21, 0, 0, u31, 0, 0],
                  [0, u12, 0, 0, u22, 0, 0, u32, 0],
                  [0, 0, u13, 0, 0, u23, 0, 0, u33],
                  [u12, u11, 0, u22, u21, 0, u32, u31, 0],
                  [u13, 0, u11, u23, 0, u21, u33, 0, u31],
                  [0, u13, u12, 0, u23, u22, 0, u33, u32]])

    print(a)

    print((h + a) * g * u)

    b = (h + a) * g

    s.pprint(b)

    a = s.Matrix([[u11+1, 0, 0, u21, 0, 0, u31, 0, 0],
                  [0, u12, 0, 0, u22+1, 0, 0, u32, 0],
                  [0, 0, u13, 0, 0, u23, 0, 0, u33+1],
                  [u12, u11+1, 0, u22+1, u21, 0, u32, u31, 0],
                  [u13, 0, u11+1, u23, 0, u21, u33+1, 0, u31],
                  [0, u13, u12, 0, u23, u22+1, 0, u33+1, u32]])
    print(a)

    print(a * g * u)

    b2 = a * g

    s.pprint(b2)

    print(b == b2)

    u11p, u22p, u33p = s.symbols(['u11p', 'u22p', 'u33p'], commutative=False)
    a = s.Matrix([[u11p, 0, 0, u21, 0, 0, u31, 0, 0],
                  [0, u12, 0, 0, u22p, 0, 0, u32, 0],
                  [0, 0, u13, 0, 0, u23, 0, 0, u33p],
                  [u12, u11p, 0, u22p, u21, 0, u32, u31, 0],
                  [u13, 0, u11p, u23, 0, u21, u33p, 0, u31],
                  [0, u13, u12, 0, u23, u22p, 0, u33p, u32]])
    print(a)

    print(a * g * u)

    b = a * g

    s.pprint(b)

if __name__ == '__main__':
    main()
