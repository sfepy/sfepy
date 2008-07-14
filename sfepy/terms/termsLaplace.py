from terms import *
from utils import fixScalarConstant, fixScalarInEl

##
# 28.11.2005, c
class LaplaceTerm( Term ):
    r""":description: Laplace term with $c$ constant or constant per element.
    :definition: $c \int_{\Omega}\nabla s \cdot \nabla r$
    or $\sum_{K \in \Tcal_h}\int_{T_K} c_K\ \nabla s \cdot \nabla r$
    """
    name = 'dw_laplace'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    symbolic = {'expression': 'c * div( grad( u ) )',
                'map' : {'u' : 'state', 'c' : 'material'}}

    ##
    # c: 28.11.2005, r: 27.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

    ##
    # c: 27.03.2008, r: 27.03.2008
    def getShape( self, diffVar, chunkSize, apr, apc = None ):
        self.dataShape = apr.getVDataShape( self.integralName )
        nEl, nQP, dim, nEP = self.dataShape
        if diffVar is None:
            return (chunkSize, 1, nEP, 1), 0
        elif diffVar == self.getArgName( 'state' ):
            return (chunkSize, 1, nEP, nEP), 1
        else:
            raise StopIteration

    ##
    # c: 27.03.2008, r: 28.03.2008
    def buildCFunArgs( self, mat, state, ap, vg ):
        vec, indx = state()
        matArg = fixScalarConstant( mat, nm.float64 )
        if matArg is None:
            matArg = fixScalarInEl( mat, self.dataShape[0], nm.float64 )
            self.function = terms.dw_st_pspg_p
        else:
            self.function = terms.dw_laplace
        fargs = vec, indx.start, matArg, vg, ap.econn

        return fargs
        
    ##
    # c: 28.11.2005, r: 28.03.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        material, virtual, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )

        shape, mode = self.getShape( diffVar, chunkSize, ap )
        fargs = self.buildCFunArgs( material, state, ap, vg )
        
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, *fargs + (chunk, mode) )
            yield out, chunk, status

##
# 01.08.2006, c
class DiffusionTerm( LaplaceTerm ):
    r""":description: General diffusion term with permeability $K_{ij}$
    constant or given in mesh vertices.
    :definition: $\int_{\Omega} K_{ij} \nabla_i q  \nabla_j p$
    """
    name = 'dw_diffusion'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    useCaches = {'mat_in_qp' : [['material']]}
    symbolic = {'expression': 'div( K * grad( u ) )',
                'map' : {'u' : 'state', 'K' : 'material'}}

    ##
    # c: 03.08.2006?, r: 28.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_diffusion )
        
    ##
    # c: 28.03.2008, r: 28.03.2008
    def buildCFunArgs( self, mat, state, ap, vg ):
        vec, indx = state()

        nEl, nQP, dim, nEP = self.dataShape
        cache = self.getCache( 'mat_in_qp', 0 )
        matQP = cache( 'matqp', self.getCurrentGroup(), 0,
                       mat = nm.asarray( mat ), ap = ap,
                       assumedShapes = [(1, nQP, dim, dim),
                                        (nEl, nQP, dim, dim)],
                       modeIn = None )
        return 1.0, vec, indx.start, matQP, vg, ap.econn

##
# 12.03.2007, c
class DiffusionIntegratedTerm( Term ):
    r""":description: Integrated general diffusion term with permeability
    $K_{ij}$  constant or given in mesh vertices.
    :definition: $\int_{\Omega} K_{ij} \nabla_i \bar{p}  \nabla_j r$
    """
    name = 'd_diffusion'
    argTypes = ('material', 'parameter_1', 'parameter_2')
    geometry = [(Volume, 'parameter_1'), (Volume, 'parameter_2')]
    useCaches = {'grad_scalar' : [['parameter_1'], ['parameter_2']],
                 'mat_in_qp' : [['material']]}

    ##
    # c: 12.03.2007, r: 28.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_diffusion )
        
    ##
    # c: 12.03.2007, r: 28.03.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        mat, par1, par2 = self.getArgs( **kwargs )
        ap, vg = par1.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )
        shape = (chunkSize, 1, 1, 1)

        cache = self.getCache( 'grad_scalar', 0 )
        gp1 = cache( 'grad', self.getCurrentGroup(), 0, state = par1 )
        cache = self.getCache( 'grad_scalar', 1 )
        gp2 = cache( 'grad', self.getCurrentGroup(), 0, state = par2 )

        cache = self.getCache( 'mat_in_qp', 0 )
        matQP = cache( 'matqp', self.getCurrentGroup(), 0,
                       mat = mat, ap = ap,
                       assumedShapes = [(1, nQP, dim, dim),
                                        (nEl, nQP, dim, dim)],
                       modeIn = None )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, 1.0, gp1, gp2, matQP, vg, chunk )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

##
# 23.04.2007, c
class PermeabilityRTerm( Term ):
    r""":description: Special-purpose diffusion-like term with permeability
    $K_{ij}$ constant or given in mesh vertices (to use on a right-hand side).
    :definition: $\int_{\Omega} K_{ij} \nabla_j q$
    """
    name = 'dw_permeability_r'
    argTypes = ('material', 'virtual', 'index')
    geometry = [(Volume, 'virtual')]
    useCaches = {'mat_in_qp' : [['material']]}

    ##
    # c: 23.04.2007, r: 28.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_permeability_r )
        
    ##
    # c: 23.04.2007, r: 28.03.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        mat, virtual, index = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, nEP, 1)
        else:
            raise StopIteration

        cache = self.getCache( 'mat_in_qp', 0 )
        matQP = cache( 'matqp', self.getCurrentGroup(), 0,
                       mat = mat, ap = ap,
                       assumedShapes = [(1, nQP, dim, dim),
                                        (nEl, nQP, dim, dim)],
                       modeIn = None )
        matQP = nm.ascontiguousarray( matQP[...,index:index+1] )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, matQP, vg, ap.econn, chunk )
            yield out, chunk, status

##
# 07.09.2006, c
class DiffusionVelocityTerm( Term ):
    r""":description: Diffusion velocity averaged in elements.
    :definition: vector of $\forall K \in \Tcal_h: \int_{T_K} K_{ij} \nabla_j r
    / \int_{T_K} 1$
    """
    name = 'de_diffusion_velocity'
    argTypes = ('material','parameter')
    geometry = [(Volume, 'parameter')]
    useCaches = {'mat_in_qp' : [['material']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.de_diffusion_velocity )
        
    ##
    # c: 07.09.2006, r: 14.07.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        mat, parameter = self.getArgs( **kwargs )
        ap, vg = parameter.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, dim, 1)
        else:
            raise StopIteration

        vec, indx = parameter()
        cache = self.getCache( 'mat_in_qp', 0 )
        matQP = cache( 'matqp', self.getCurrentGroup(), 0,
                       mat = nm.asarray( mat ), ap = ap,
                       assumedShapes = [(1, nQP, dim, dim),
                                        (nEl, nQP, dim, dim)],
                       modeIn = None )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start,
                                    matQP, vg, ap.econn, chunk )
            out1 = out / vg.variable( 2 )[chunk]
            yield out1, chunk, status
