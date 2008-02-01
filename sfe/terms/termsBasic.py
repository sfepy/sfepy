from terms import *

##
# 12.04.2007, c
class IntegrateVolumeTerm( Term ):
    """:definition: $\int_\Omega y$"""
    name = 'd_volume_integrate'
    argTypes = ('parameter',)
    geometry = [(Volume, 'parameter')]
    useCaches = {'state_in_volume_qp' : [['parameter']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
    ##
    # created:       12.04.2007
    # last revision: 21.12.2007
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        par, = self.getArgs( **kwargs )
        ap, vg = par.getApproximation( self.getCurrentGroup(), 'Volume' )
        shape = (chunkSize, 1, 1, 1)

        cache = self.getCache( 'state_in_volume_qp', 0 )
        vec = cache( 'state', self.getCurrentGroup(), 0, state = par )

        for out, chunk in self.charFun( chunkSize, shape ):
            status = vg.integrateChunk( out, vec[chunk], chunk )
            out1 = nm.sum( out )
            yield out1, chunk, status

##
# 01.11.2007, c
class IntegrateVolumeOperatorTerm( Term ):
    """:definition: $\int_\Omega q$"""

    name = 'dw_volume_integrate'
    argTypes = ('virtual',)
    geometry = [(Volume, 'virtual')]

    ##
    # 01.11.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
    ##
    # created:       01.11.2007
    # last revision: 21.12.2007
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        virtual, = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, nEP, 1 )
            mode = 0
        else:
            raise StopIteration

        bf = ap.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            bfT = nm.tile( bf.transpose( (0, 2, 1) ), (chunkSize, 1, 1, 1) )
            status = vg.integrateChunk( out, bfT, chunk )
            yield out, chunk, 0

##
# 24.04.2007, c
class IntegrateSurfaceTerm( Term ):
    """:definition: $\int_\Gamma y$, for vectors: $\int_\Gamma \ul{y} \cdot
    \ul{n}$"""
    name = 'd_surface_integrate'
    argTypes = ('parameter',)
    geometry = [(Surface, 'parameter')]
    useCaches = {'state_in_surface_qp' : [['parameter']]}

    ##
    # 24.04.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        self.dofConnType = 'surface'

    ##
    # c: 24.04.2007, r: 15.01.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        """
        Integrates over surface.
        """
        par, = self.getArgs( **kwargs )
        ap, sg = par.getApproximation( self.getCurrentGroup(), 'Surface' )
        shape = (chunkSize, 1, 1, 1)

        sd = ap.surfaceData[self.region.name]

        cache = self.getCache( 'state_in_surface_qp', 0 )
        vec = cache( 'state', self.getCurrentGroup(), 0, state = par )
        
        for out, chunk in self.charFun( chunkSize, shape ):
            lchunk = self.charFun.getLocalChunk()
            status = sg.integrateChunk( out, vec[lchunk], lchunk )
            out1 = nm.sum( out )
            yield out1, chunk, status

##
# 26.09.2007, c
class DotProductVolumeTerm( Term ):
    """:description: Volume $L^2(\Omega)$ dot product for both scalar and
    vector fields.
    :definition: $\int_\Omega p r$, $\int_\Omega \ul{u} \cdot \ul{w}$"""
    name = 'd_volume_dot'
    argTypes = ('parameter_1', 'parameter_2')
    geometry = [(Volume, 'parameter_1'), (Volume, 'parameter_2')]
    useCaches = {'state_in_volume_qp' : [['parameter_1'], ['parameter_2']]}

    ##
    # 26.09.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
    ##
    # created:       26.09.2007
    # last revision: 13.12.2007
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        par1, par2 = self.getArgs( **kwargs )
        ap, vg = par1.getApproximation( self.getCurrentGroup(), 'Volume' )
        shape = (chunkSize, 1, 1, 1)

        cache = self.getCache( 'state_in_volume_qp', 0 )
        vec1 = cache( 'state', self.getCurrentGroup(), 0, state = par1 )
        cache = self.getCache( 'state_in_volume_qp', 1 )
        vec2 = cache( 'state', self.getCurrentGroup(), 0, state = par2 )

        for out, chunk in self.charFun( chunkSize, shape ):
            if vec1.shape[-1] > 1:
                vec = nm.sum( vec1[chunk] * vec2[chunk], axis = -1 )
            else:
                vec = vec1[chunk] * vec2[chunk]
            status = vg.integrateChunk( out, vec, chunk )
            out1 = nm.sum( out )
            yield out1, chunk, status

##
# 09.10.2007, c
class DotProductSurfaceTerm( Term ):
    """:description: Surface $L^2(\Gamma)$ dot product for both scalar and
    vector fields.
    :definition: $\int_\Gamma p r$, $\int_\Gamma \ul{u} \cdot \ul{w}$"""
    name = 'd_surface_dot'
    argTypes = ('parameter_1', 'parameter_2')
    geometry = [(Surface, 'parameter_1'), (Surface, 'parameter_2')]
    useCaches = {'state_in_surface_qp' : [['parameter_1'], ['parameter_2']]}

    ##
    # 09.10.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
    ##
    # c: 09.10.2007, r: 15.01.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        par1, par2 = self.getArgs( **kwargs )
        ap, sg = par1.getApproximation( self.getCurrentGroup(), 'Surface' )
        shape = (chunkSize, 1, 1, 1)

        cache = self.getCache( 'state_in_surface_qp', 0 )
        vec1 = cache( 'state', self.getCurrentGroup(), 0, state = par1 )
        cache = self.getCache( 'state_in_surface_qp', 1 )
        vec2 = cache( 'state', self.getCurrentGroup(), 0, state = par2 )

        for out, chunk in self.charFun( chunkSize, shape ):
            lchunk = self.charFun.getLocalChunk()
            if vec1.shape[-1] > 1:
                vec = nm.sum( vec1[lchunk] * vec2[lchunk], axis = -1 )
            else:
                vec = vec1[lchunk] * vec2[lchunk]
            status = sg.integrateChunk( out, vec, lchunk )
            out1 = nm.sum( out )
            yield out1, chunk, status

##
# 16.07.2007, c
class VolumeTerm( Term ):
    """:description: Volume of a domain. Uses approximation of the parameter
    variable.
    :definition: $\int_\Omega 1$"""
    name = 'd_volume'
    argTypes = ('parameter',)
    geometry = [(Volume, 'parameter')]
    useCaches = {'volume' : [['parameter']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
    ##
    # created:       16.07.2007
    # last revision: 13.12.2007
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        par, = self.getArgs( **kwargs )
        shape = (1, 1, 1, 1)

        cache = self.getCache( 'volume', 0 )
        volume = cache( 'volume', self.getCurrentGroup(), 0,
                        region = self.charFun.region, field = par.field )
        yield volume, 0, 0
