from collections import deque

from sfe.base.base import *

##
# 24.11.2006, c
class DataCache( Struct ):
    """valid[key][ckey] ... for ih == 0,
       data[key][ckey][ih]
       memSizes[key], historySizes[key] ... same for all ckey"""
    argTypes = ()
    
    ##
    # c: 19.10.2006, r: 02.04.2008
    def __init__( self, name, argNames, keys, historySizes = None,
                  function = None ):
        """historySizes[key] = (historySize, maxMemHistorySize)
        Both historySize and maxMemHistorySize can be -1, i.e. growing as
        necessary"""
        self.name = name
        self.argNames = argNames
        self.function = function
#        self.initialized = dict.fromkeys( keys, False )
        self.valid = dictFromKeysInit( keys, dict )
        self.clear()


        self.data = dictFromKeysInit( keys, dict )
        self.override = False
        self.memSizes = {}.fromkeys( keys, 1 )
        self.historySizes = {}.fromkeys( keys, 1 )
        self.memGrowing = {}.fromkeys( keys, False )
        self.historyGrowing = {}.fromkeys( keys, False )
        self.mergeHistorySizes( historySizes )
        self.step = 0

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
           so that initData() must be called."""
        self.valid = dictFromKeysInit( self.keys(), dict )

    ##
    # 02.03.2007, c
    def setMode( self, override ):
        self.override = override

    ##
    # 27.02.2007, c
    def getArgs( self, **kwargs ):
        """Caches can use an optional 'history' kwarg - it is up to an instance
        to obtain it from kwargs."""
        args = []
        for at in self.argTypes:
            args.append( kwargs[at] )
        return args

    ##
    # 08.06.2007, c
    # 11.06.2007
    # 12.06.2007
    def mergeHistorySizes( self, historySizes = None ):
        if historySizes is not None:
            for key, hs in historySizes.iteritems():
                if not key in self.keys():
                    print 'unknown cache key! (%s in %s)' % (key, self.keys())
                    raise AssertionError
                if hs[0] >= 0 and (not self.historyGrowing[key]):
                    self.historySizes[key] = max( hs[0],
                                                  self.historySizes[key] )
                else:
                    self.historySizes[key] = 1
                    self.historyGrowing[key] = True

                if hs[1] >= 0 and (not self.memGrowing[key]):
                    maxMemHistorySize = max( hs[1], self.memSizes[key] )
                else:
                    self.memSizes[key] = maxMemHistorySize = 1
                    self.memGrowing[key] = True

                    self.historySizes[key] = 1
                    self.historyGrowing[key] = True

                self.memSizes[key] = min( maxMemHistorySize,
                                          self.historySizes[key] )

    ##
    # c: 02.04.2008, r: 02.04.2008
    def initTime( self, ts ):
        self.step = ts.step

    ##
    # c: 28.11.2006, r: 02.04.2008
    def initData( self, key, ckey, shape ):

        self.valid[key][ckey] = False
        data = self.data[key]
        
        # Current time data.
        data[ckey] = deque()
        arr = nm.empty( shape, dtype = nm.float64 )
        data[ckey].append( arr )
            
        # History data.
        for ih in xrange( self.memSizes[key] - 1 ):
            data[ckey].append( nm.empty_like( arr ) )

    ##
    # 30.11.2006, c
    # 08.06.2007
    # 11.06.2007
    # 12.06.2007
    def advance( self, step ):
        """History advancement."""
        # Should data be zeroed here?
        if step != self.step + 1:
            print 'bad step %d == %d' % (step, self.step + 1)
            raise RuntimeError

        self.step = step

        for key, data in self.data.iteritems():
            mhs = self.memSizes[key]
            if self.memGrowing[key]:
                for ckey in self.valid[key].iterkeys():
                    arr = nm.empty_like( self.data[key][ckey][0] )
                    self.data[key][ckey].appendleft( arr )
                    self.memSizes[key] += 1
                    self.historySizes[key] += 1
            elif mhs == 1: continue
            else:
                for ckey in self.valid[key].iterkeys():
                    # Rotate list.
                    self.data[key][ckey].rotate()

                    if mhs < self.historySizes:
                        # Data falling to disk. To finish...
                        data = self.data[key][ckey][0]

    ##
    # c: 13.12.2007, r: 15.01.2008
    def gToC( self, groupIndx ):
        if hasattr( self, 'regionMatters' ) and self.regionMatters:
            return groupIndx
        else:
            return groupIndx[0], groupIndx[-1]
        
    ##
    # c: 19.10.2006, r: 13.12.2007
    def __call__( self, key, groupIndx, ih, **kwargs ):
        """groupIndx : term.getCurrentGroup() - term.region.name ignored
                  ih : history level index
                       0 .. current, 1, 2, 3...
        """
        if not key in self.valid.keys():
            print 'invalid cache key: %s in %s' % (key, self.keys())
            raise ValueError

        ckey = self.gToC( groupIndx )
        if not self.valid[key].has_key( ckey ):
            self.initData( key, ckey, **kwargs )

        if (not self.valid[key][ckey] and (ih == 0)) or self.override:
            self.update( key, groupIndx, ih, **kwargs )
            self.valid[key][ckey] = True
#            print self.name, key, ckey, ih, 'up'

#        print key, kwargs

        if ih < self.memSizes[key]:
#            print self.name, key, ckey, ih, 'get'
            return self.data[key][ckey][ih]
        else:
            print ih, self.memSizes[key]
            raise NotImplementedError
