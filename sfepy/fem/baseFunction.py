from sfepy.base.base import *

##
# 30.03.2005, c
class BaseFunction( Struct ):
    ##
    # 30.03.2005, c
    # 20.07.2005
    # 07.12.2005
    def __init__( self, bf, nodes, var_set = None ):
        self.nodes = nodes.copy()
        self.bf = bf
        self.var_set = var_set

        self.aos = nm.sum( nodes, 1 )
    ##
    # base_fun.value( qp_coor['v'], nodes )
    # base_fun.value( qp_coor['v'], nodes, (0,1,2) )
    #
    # 30.03.2005, c
    # 01.04.2005
    # 20.07.2005
    # 23.11.2005
    # 07.12.2005
    # 01.09.2007
    def value( self, coors, nodes, ivs = None, suppress_errors = False ):

##         print self
##         print self.bf
##         print coors
##         print nodes

        # Check nodes (simplified!).
        naos = nm.sum( nodes, 1 )
        if nm.any( naos - self.aos ):
            print "bad base function node in:"
            print nodes, naos
            raise AssertionError
        
        n_nod = len( nodes )
        if not( ivs ):
            val = nm.zeros( (coors.shape[0], 1, n_nod), nm.float64 )
            for ii, node in enumerate( nodes ):
                val[:,0,ii] = self.bf( coors, node,
                                       suppress_errors = suppress_errors )
        else:
            # Check iv.
            for iv in ivs:
                if not(iv in self.var_set):
                    print "bad base function variable:", iv
                    raise ValueError

            val = nm.zeros( (coors.shape[0], len( ivs ), n_nod), nm.float64 )
            for ir, iv in enumerate( ivs ):
                for ic, node in enumerate( nodes ):
                    val[:,ir,ic] = self.bf( coors, node, iv,
                                            suppress_errors = suppress_errors )

        return( val )
