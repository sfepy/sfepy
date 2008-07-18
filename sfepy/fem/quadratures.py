from sfepy.base.base import *
"""
This is a test file only! Orders of quadrature may be wrong!!!

Simplex quadratures for <0, 1> x dim
Tensor product quadratures for <-1, 1> x dim
... due to base functions used.
"""


##
# 13.11.2007, c
class Quadrature( Struct ):
    """Naming conventions:
abstract description: <family>_<max. order>_<dimension>
geometry description: <family>_<max. order>_<geometry>


    """
    def __call__( self ):
        return self.vals, self.weights

##
# 16.11.2007, c
class CustomQuadrature( Quadrature ):
    family = 'custom'
    name = 'custom'

    def from_conf( conf ):
        if conf.mode == 'custom':
            obj = CustomQuadrature( vals = nm.asarray( conf.vals, nm.float64 ),
                                    weights = nm.asarray( conf.weights,
                                                          nm.float64 ) )
        else:
            print 'unknown custom quadrature mode:', conf.mode
            raise ValueError

        return obj
    from_conf = staticmethod( from_conf )


##
# created:       03.12.2007
# last revision: 11.12.2007
class GaussSimplexO1G12( Quadrature ):
    family = 'gauss_o1_d1'
    name = 'gauss_s_o1_g1_2'

    def __init__( self ):
        self.vals = nm.array( [[0.5]], dtype = nm.float64 )
        self.weights = nm.array( [1.0], dtype = nm.float64 )
    
##
# created:       11.12.2007
# last revision: 11.12.2007
class GaussSimplexO2G12( Quadrature ):
    family = 'gauss_o2_d1'
    name = 'gauss_s_o2_g1_2'

    def __init__( self ):
        GAUSS3 = 0.7745966692414833770358531
        GAUSS3W1 = 0.5555555555555555555555556
        GAUSS3W2 = 0.8888888888888888888888889

        self.vals = (nm.array( [[-GAUSS3], [0.0], [GAUSS3]],
                               dtype = nm.float64 ) + 1.0) / 2.0
        self.weights = nm.array( [GAUSS3W1, GAUSS3W2, GAUSS3W1],
                                 dtype = nm.float64 ) / 2.0
##
# created:       11.12.2007
# last revision: 11.12.2007
class GaussTensorProductO2G12( Quadrature ):
    family = 'gauss_o2_d1'
    name = 'gauss_tp_o2_g1_2'

    def __init__( self ):
        a = nm.sqrt( 3.0 ) / 3.0
        self.vals = nm.array( [[-a], [a]], dtype = nm.float64 )
        self.weights = nm.array( [1.0, 1.0], dtype = nm.float64 )

##
# 13.11.2007, c
class GaussO2G23( Quadrature ):
    family = 'gauss_o2_d2'
    name = 'gauss_s_o2_g2_3'

    def __init__( self ):
        c = 1.0 / 6.0
        d = 2.0 / 3.0

        self.vals = nm.array( [[c, c], [d, c], [c, d]], dtype = nm.float64 )
        self.weights = nm.array( [c] * 3, dtype = nm.float64 )

##
# created:       11.12.2007
# last revision: 11.12.2007
class GaussO3G23( Quadrature ):
    family = 'gauss_o3_d2'
    name = 'gauss_s_o3_g2_3'

    def __init__( self ):
        d = 1.0 / 3.0
        e = 1.0 / 5.0
        f = 3.0 / 5.0
        w3 = -27.0 / 96.0
        w4 = 25.0 / 96.0

        self.vals = nm.array(  [[d, d],
                                [e, e],
                                [f, e],
                                [e, f]], dtype = nm.float64 )
        self.weights = nm.array( [w3, w4, w4, w4], dtype = nm.float64 )

##
# created:       13.11.2007
# last revision: 11.12.2007
class GaussO2G24( Quadrature ):
    family = 'gauss_o2_d2'
    name = 'gauss_tp_o2_g2_4'

    def __init__( self ):
        a = nm.sqrt( 3.0 ) / 3.0

        self.vals = nm.array( [[-a, -a], [a, -a], [a, a], [-a, a]],
                              dtype = nm.float64 )
        self.weights = nm.array( [1.0] * 4, dtype = nm.float64 )

##
# created: 03.12.2007
class GaussO1G38( Quadrature ):
    family = 'gauss_o1_d3'
    name = 'gauss_tp_o1_g3_8'

    def __init__( self ):
        a = nm.sqrt( 3.0 ) / 3.0

        self.vals = nm.array( [[-a, -a, -a],
                               [ a, -a, -a],
                               [ a,  a, -a],
                               [-a,  a, -a],
                               [-a, -a,  a],
                               [ a, -a,  a],
                               [ a,  a,  a],
                               [-a,  a,  a]],
                              dtype = nm.float64 )
        self.weights = nm.array( [1.0] * 8, dtype = nm.float64 )

##
# created: 06.12.2007
class GaussO2G38( Quadrature ):
    family = 'gauss_o2_d3'
    name = 'gauss_tp_o2_g3_8'

    def __init__( self ):
        a = nm.sqrt( 3.0 ) / 3.0

        self.vals = nm.array( [[-a, -a, -a],
                               [ a, -a, -a],
                               [ a,  a, -a],
                               [-a,  a, -a],
                               [-a, -a,  a],
                               [ a, -a,  a],
                               [ a,  a,  a],
                               [-a,  a,  a]],
                              dtype = nm.float64 )
        self.weights = nm.array( [1.0] * 8, dtype = nm.float64 )

##
# created: 03.12.2007
class GaussO1G34( Quadrature ):
    family = 'gauss_o1_d3'
    name = 'gauss_s_o1_g3_4'

    def __init__( self ):
        self.vals = nm.array( [[0.25, 0.25, 0.25]], dtype = nm.float64 )
        self.weights = nm.array( [1.0/6.0], dtype = nm.float64 )

##
# created: 06.12.2007
class GaussO2G34( Quadrature ):
    family = 'gauss_o2_d3'
    name = 'gauss_s_o2_g3_4'

    def __init__( self ):
        a = (5.0 - pow( 5.0, 0.5 )) / 20.0
        b = (5.0 + 3.0 * pow( 5.0, 0.5 )) / 20.0
        w = 1.0 / 24.0

        self.vals = nm.array( [[a, a, a],
                               [a, b, a],
                               [a, a, b],
                               [b, a, a]], dtype = nm.float64 )
        self.weights = nm.array( [w] * 4, dtype = nm.float64 )

##
# created: 06.12.2007
class GaussO3G34( Quadrature ):
    family = 'gauss_o3_d3'
    name = 'gauss_s_o3_g3_4'

    def __init__( self ):
        # 3D 5 point, order 3
        a = 0.25
        b = 0.5
        c = 1.0 / 6.0
        w1 = -2.0 / 15.0
        w2 = 3.0 / 40.0

        self.vals = nm.array( [[a, a, a],
                               [c, c, c],
                               [c, b, c],
                               [c, c, b],
                               [b, c, c]], dtype = nm.float64 )
        self.weights = nm.array( [w1, w2, w2, w2, w2], dtype = nm.float64 )

##
# 13.11.2007, c
def collect_quadratures():
    var_dict = globals().items()

    quad_table = {}
    for key, var in var_dict:
        try:
            if is_derived_class( var, Quadrature ):
                quad_table[var.name] = var
        except TypeError:
            pass
#    print quad_table

    family_table = {}
    for name, cls in quad_table.iteritems():
        family_table.setdefault( cls.family, {}  )[name] = cls

#    print family_table

    return family_table

quadratures = collect_quadratures()

if __name__ == "__main__":
    print collect_quadratures()
