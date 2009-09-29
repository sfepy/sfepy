from collections import deque

from sfepy.base.base import *

class DataCaches( Container ):
    def __init__( self, **kwargs ):
        Container.__init__( self, **kwargs )

        self.cache_info = {}

    def insert_cache( self, cache ):
        self.append( cache )

    def insert_term( self, cache_name, term_name, term_args ):
        self.cache_info.setdefault( cache_name, [] ).append( (term_name,
                                                              term_args) )

    def del_term( self, cache_name, term_name, term_args ):
        """Destroys a cache if its term count reaches zero."""
        ci = self.cache_info[cache_name]
        ci.remove( (term_name, term_args) )
        if len( ci ) == 0:
            del self.cache_info[cache_name]
            self.remove_name( cache_name )

##
# 24.11.2006, c
class DataCache( Struct ):
    """valid[key][ckey] ... for ih == 0,
       data[key][ckey][ih]
       mem_sizes[key], history_sizes[key] ... same for all ckey"""
    arg_types = ()
    
    ##
    # c: 19.10.2006, r: 02.04.2008
    def __init__( self, name, arg_names, keys, history_sizes = None,
                  function = None ):
        """history_sizes[key] = (history_size, max_mem_history_size)
        Both history_size and max_mem_history_size can be -1, i.e. growing as
        necessary"""
        self.name = name
        self.arg_names = arg_names
        self.function = function
#        self.initialized = dict.fromkeys( keys, False )
        self.valid = dict_from_keys_init( keys, dict )
        self.clear()


        self.data = dict_from_keys_init( keys, dict )
        self.override = False
        self.mem_sizes = {}.fromkeys( keys, 1 )
        self.history_sizes = {}.fromkeys( keys, 1 )
        self.mem_growing = {}.fromkeys( keys, False )
        self.history_growing = {}.fromkeys( keys, False )
        self.merge_history_sizes( history_sizes )
        self.step = 0

        self.dtype = nm.float64
        self.geometry = 'volume'
        self.region_name = None

    ##
    # 28.11.2006, c
    def keys( self ):
        return self.valid.keys()
    
    ##
    # 19.10.2006, c
    # 27.02.2007
    def clear( self ):
        for valids in self.valid.itervalues():
            for ckey in valids.iterkeys():
                valids[ckey] = False

    ##
    # c: 07.05.2008, r: 07.05.2008
    def reset( self ):
        """Complete reset: unlike clear(), reset() removes all ckeys,
           so that init_data() must be called."""
        self.valid = dict_from_keys_init( self.keys(), dict )

    ##
    # 02.03.2007, c
    def set_mode( self, override ):
        self.override = override

    ##
    # 27.02.2007, c
    def get_args( self, **kwargs ):
        """Caches can use an optional 'history' kwarg - it is up to an instance
        to obtain it from kwargs."""
        args = []
        for at in self.arg_types:
            args.append( kwargs[at] )
        return args

    ##
    # 08.06.2007, c
    # 11.06.2007
    # 12.06.2007
    def merge_history_sizes( self, history_sizes = None ):
        if history_sizes is not None:
            for key, hs in history_sizes.iteritems():
                if not key in self.keys():
                    msg = 'unknown cache key! (%s in %s)' % (key, self.keys())
                    raise ValueError( msg )
                if hs[0] >= 0 and (not self.history_growing[key]):
                    self.history_sizes[key] = max( hs[0],
                                                  self.history_sizes[key] )
                else:
                    self.history_sizes[key] = 1
                    self.history_growing[key] = True

                if hs[1] >= 0 and (not self.mem_growing[key]):
                    max_mem_history_size = max( hs[1], self.mem_sizes[key] )
                else:
                    self.mem_sizes[key] = max_mem_history_size = 1
                    self.mem_growing[key] = True

                    self.history_sizes[key] = 1
                    self.history_growing[key] = True

                self.mem_sizes[key] = min( max_mem_history_size,
                                          self.history_sizes[key] )

    ##
    # c: 02.04.2008, r: 02.04.2008
    def init_time( self, ts ):
        self.step = ts.step

    ##
    # c: 28.11.2006, r: 02.04.2008
    def init_data( self, key, ckey, shape ):

        self.valid[key][ckey] = False
        data = self.data[key]
        
        # Current time data.
        data[ckey] = deque()
        arr = nm.empty(shape, dtype=self.dtype)
        data[ckey].append( arr )
            
        # History data.
        for ih in xrange( self.mem_sizes[key] - 1 ):
            data[ckey].append( nm.empty_like( arr ) )

    def init_datas( self, ckey, shapes ):
        for key, shape in shapes.iteritems():
            DataCache.init_data( self, key, ckey, shape )
    ##
    # 30.11.2006, c
    # 08.06.2007
    # 11.06.2007
    # 12.06.2007
    def advance( self, step ):
        """History advancement."""
        # Should data be zeroed here?
        if step != self.step + 1:
            msg = 'bad step %d == %d' % (step, self.step + 1)
            raise ValueError(msg)

        self.step = step

        if hasattr(self, 'custom_advance'):
            for key in self.data.iterkeys():
                for ckey in self.valid[key].iterkeys():
                    self.custom_advance(ckey, step)

        else:
            for key, data in self.data.iteritems():
                mhs = self.mem_sizes[key]
                if self.mem_growing[key]:
                    for ckey in self.valid[key].iterkeys():
                        arr = nm.empty_like( self.data[key][ckey][0] )
                        self.data[key][ckey].appendleft( arr )
                        self.mem_sizes[key] += 1
                        self.history_sizes[key] += 1
                elif mhs == 1: continue
                else:
                    for ckey in self.valid[key].iterkeys():
                        # Rotate list.
                        self.data[key][ckey].rotate()

                        if mhs < self.history_sizes:
                            # Data falling to disk. To finish...
                            data = self.data[key][ckey][0]

    ##
    # c: 13.12.2007, r: 15.01.2008
    def g_to_c( self, group_indx ):
        if hasattr( self, 'region_matters' ) and self.region_matters:
            return group_indx
        else:
            return group_indx[0], group_indx[-1]
        

    def _call( self, key, group_indx, ih, **kwargs ):
        """group_indx : term.get_current_group() - term.region.name ignored
                  ih : history level index
                       0 .. current, 1, 2, 3...
        """
#        print key, group_indx, ih, kwargs
        if not key in self.valid.keys():
            err = 'invalid cache key: %s not in %s' % (key, self.keys())
            raise ValueError( err )

        ckey = self.g_to_c( group_indx )
        if not self.valid[key].has_key( ckey ):
            self.init_data( key, ckey, **kwargs )

        if (not self.valid[key][ckey] and (ih == 0)) or self.override:
            self.update( key, group_indx, ih, **kwargs )
            self.valid[key][ckey] = True
#            print self.name, key, ckey, ih, 'up'

#        print key, kwargs

        if ih < self.mem_sizes[key]:
#            print self.name, key, ckey, ih, 'get'
            return self.data[key][ckey][ih]
        else:
            print ih, self.mem_sizes[key]
            raise NotImplementedError

    def __call__( self, key, group_indx, ih, **kwargs ):
        """group_indx : term.get_current_group() - term.region.name ignored
                  ih : history level index
                       0 .. current, 1, 2, 3...
        """
        if isinstance( key, str ):
            out = self._call( key, group_indx, ih, **kwargs )

        else:
            assert_( isinstance( key, list ) )
            out = [self._call( k, group_indx, ih, **kwargs ) for k in key]

        return out
    
