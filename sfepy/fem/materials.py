from sfepy.base.base import *

##
# 21.07.2006, c
class Materials( Container ):
    ##
    # 24.07.2006, c
    def from_conf( conf, wanted = None ):

        if wanted is None:
            wanted = conf.keys()

        objs = OneTypeList( Material )
        for key, val in conf.iteritems():
            if key in wanted:
                objs.append( Material( **val ) )
        obj = Materials( objs )
        return obj
    from_conf = staticmethod( from_conf )

    ##
    # 22.08.2006, c
    def setup_regions( self, regions ):
        for mat in self:
            mat.setup_regions( regions )

    ##
    # c: 01.08.2006, r: 20.02.2008
    def time_update( self, ts, funmod, domain, extra_mat_args = None ):
        extra_mat_args = get_default( extra_mat_args, {} )
        output( 'updating materials...' )
        tt = time.clock()
        for mat in self:
            output( ' ', mat.name )
            extra_args = extra_mat_args.setdefault( mat.name, {} )
            mat.time_update( ts, funmod, domain, **extra_args )
        output( '...done in %.2f s' % (time.clock() - tt) )

    ##
    # 22.08.2006, c
    def set_current_group( self, ig ):
        for mat in self:
            mat.set_current_group( ig )

##
# 21.07.2006, c
class Material( Struct ):
    """
    A class holding constitutive and other material parameters.

    Example input:

    material_2 = {
       'name' : 'm',
       'mode' : 'here',
       'region' : 'Omega',
       'E' : 1.0,
    }

    Material parameters are passed to terms using the dot notation,
    i.e. 'm.E' in our example case.

    >>> mat = Material( **material_2 )
    >>> print mat
Material:m
  E:
    1.0
  name:
    m
  region:
    None
  extra_args:
    {}
  mode:
    here
  region_name:
    Omega
    
    """
    ##
    # c: 22.08.2006, r: 02.07.2008
    def __init__( self, **kwargs ):
        kwargs.setdefault( 'extra_args', {} )
        Struct.__init__( self, **kwargs )

        self.region_name = self.region
        self.region = None
        self.datas = None
        self.kind = get_default_attr( self, 'kind', 'time-dependent' )

    ##
    # 22.08.2006, c
    def setup_regions( self, regions ):
        region = regions[self.region_name]
        self.igs = region.igs
        self.region = region 

    ##
    # c: 01.08.2006, r: 02.07.2008
    def time_update( self, ts, funmod, domain, **extra_args ):
        """coor is in region.vertices[ig] order (i.e. sorted by node number)"""
        if self.mode == 'function':
            self.data = None
            if (self.datas is not None) and (self.kind == 'stationary'):
                return
            
            self.datas = []

            if isinstance( self.function, str ):
                fun = getattr( funmod, self.function )
            else:
                fun = self.function

            kwargs = copy( self.extra_args )
            kwargs.update( extra_args )
            args = dict( ts = ts, region = self.region, **kwargs )

            for ig in self.igs:
                coor = domain.get_mesh_coors()[self.region.get_vertices( ig )]
                args.update( {'coor' : coor, 'ig' : ig} )
                self.datas.append( fun( **args ) )

    ##
    # 31.07.2007, c
    def set_data( self, datas ):
        self.mode = 'user'
        self.datas = datas
        self.data = None

    ##
    # 01.08.2007, c
    def set_function( self, function ):
        self.mode =  'function'
        self.function = function

    ##
    # 01.08.2007, c
    def set_extra_args( self, extra_args ):
        self.extra_args = extra_args
        
    ##
    # 22.08.2006, c
    # 22.02.2007
    # 31.07.2007
    def set_current_group( self, ig ):
        if (self.mode == 'function') or (self.mode == 'user'):
            try:
                ii = self.igs.index( ig )
                self.data = self.datas[ii]
            except:
                self.data = None

    ##
    # c: 02.08.2006, r: 02.05.2008
    def get_data( self, region_name, ig, name = None ):
        """Returns None in function mode if set_current_group() was not called."""
##         print 'getting', name

        if name is None:
            output( 'material arguments must use the dot notation!' )
            output( '(material: %s, region: %s)' % (self.name, region_name) )
            raise ValueError

        if self.mode == 'here':
            return getattr( self, name )
        else:
            ii = self.igs.index( ig )
            if isinstance( self.datas[ii], Struct ):
                return getattr( self.datas[ii], name )
            else:
                return self.datas[ii][name]

    ##
    # 01.08.2007, c
    def reduce_on_datas( self, reduce_fun, init = 0.0 ):
        out = {}.fromkeys( self.datas[0].keys(), init )

        for data in self.datas:
            for key, val in data.iteritems():
                out[key] = reduce_fun( out[key], val )

        return out
