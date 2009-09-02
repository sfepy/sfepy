from sfepy.base.base import *

##
# c: 22.07.2008
def youngpoisson_to_lame( young, poisson, plane = 'stress' ):

    if plane == 'stress':
        lam = young*poisson/(1.0 - poisson*poisson)
        mu = young/(2.0*(1.0 + poisson))
    elif plane == 'strain':
        lam = young*poisson/((1.0 + poisson)*(1.0 - 2.0*poisson))
        mu = young/(2.0*(1.0 + poisson))

    return lam, mu

##
# c: 22.07.2008
def stiffness_tensor_lame( dim, lam, mu ):

    sym = (dim + 1) * dim / 2
    o = nm.array( [1.] * dim + [0.] * (sym - dim), dtype = nm.float64 )
    oot = nm.outer( o, o )

    return lam * oot + mu * nm.diag( o + 1.0 )
    
##
# c: 22.07.2008
def stiffness_tensor_youngpoisson( dim, young, poisson, plane = 'stress' ):

    lam, mu = youngpoisson_to_lame( young, poisson, plane )

    return stiffness_tensor_lame( dim, lam, mu )
 
##
# c: 10.08.2009
def stiffness_tensor_lame_mixed( dim, lam, mu ):

    sym = (dim + 1) * dim / 2
    o = nm.array( [1.] * dim + [0.] * (sym - dim), dtype = nm.float64 )
    oot = nm.outer( o, o )

    return 2.0/3.0*(lam-mu) * oot + mu * nm.diag( o + 1.0 )

##
# c: 10.08.2009
def stiffness_tensor_youngpoisson_mixed( dim, young, poisson, plane = 'stress' ):

    lam, mu = youngpoisson_to_lame( young, poisson, plane )

    return stiffness_tensor_lame_mixed( dim, lam, mu )

##
# c: 10.08.2009
def bulk_modulus_lame( lam, mu ):

    return 1.0/3.0 * (2*mu + lam)

##
# c: 10.08.2009
def bulk_modulus_youngpoisson( young, poisson, plane = 'stress' ):

    lam, mu = youngpoisson_to_lame( young, poisson, plane )

    return bulk_modulus_lame( lam, mu )

class TransformToPlane( Struct ):
    """Transformmations of constitutive law coefficients of 3D problems to 2D."""

    def __init__( self, iplane = None ):
        """`iplane` ... vector of indices denoting the plane, e.g.: [0, 1]"""
        if iplane is None:
            iplane = [0, 1]

        # Choose the "master" variables and the "slave" ones
        # ... for vectors
        i_m = nm.sort( iplane )
        i_s = nm.setdiff1d( nm.arange( 3 ), i_m )

        # ... for second order tensors (symmetric storage)
        i_ms = {(0, 1) : [0, 1, 3],
               (0, 2) : [0, 2, 4],
               (1, 2) : [1, 2, 5]}[tuple( i_m )]
        i_ss = nm.setdiff1d( nm.arange( 6 ), i_ms )
        
        Struct.__init__( self, iplane = iplane,
                         i_m = i_m, i_s = i_s,
                         i_ms = i_ms, i_ss = i_ss )

    def tensor_plane_stress( self, c3 = None, d3 = None, b3 = None ):
       """Transforms all coefficients of the piezoelectric constitutive law
       from 3D to plane stress problem in 2D: strain/stress ordering/ 11 22
       33 12 13 23. If `d3` is None, uses only the stiffness tensor `c3`.

       `c3` ... stiffness tensor
       `d3` ... dielectric tensor
       `b3` ... piezoelectric coupling tensor"""

       mg = nm.meshgrid

       cs = c3[mg(self.i_ss,self.i_ss)]
       cm = c3[mg(self.i_ss,self.i_ms)].T
       if d3 is None: # elasticity only.
           A = cs
           Feps = cm

           Ainv = nm.linalg.inv( A )
           c2 = c3[mg(self.i_ms,self.i_ms)] \
                - nm.dot( Feps.T, nm.dot( Ainv, Feps ) )
           return c2

       else:
           dm = d3[mg(self.i_s,self.i_m)].T
           ds = d3[mg(self.i_s,self.i_s)]

           ii = mg( self.i_s, self.i_ss )
           A = nm.r_[nm.c_[cs, b3[ii]],
                     nm.c_[b3[ii].T, -ds]] #=> sym !!!
           F = nm.r_[nm.c_[cm, b3[mg(self.i_m,self.i_ss)]],
                     nm.c_[b3[mg(self.i_s,self.i_ms)].T, -dm ]]

           Feps = F[:,:3]
           FE = F[:,3:]
           Ainv = nm.linalg.inv( A )
           c2 = c3[mg(self.i_ms,self.i_ms)] \
                - nm.dot( Feps.T, nm.dot( Ainv, Feps ) )
           d2 = d3[mg(self.i_m,self.i_m)] \
                - nm.dot( FE.T, nm.dot( Ainv, FE ) )
           b2 = b3[mg(self.i_m,self.i_ms)].T \
                - nm.dot( FE.T, nm.dot( Ainv, Feps ) )
           return c2, d2, b2
