from sfe.base.base import *
from parseEq import createBNF
from materials import Materials
from sfe.terms import Terms, Term, termTable, cacheTable

"""
Note:
- create found materials, variables from configuration/input file data
  ... no - user should be able to create objects even if they are not
  used in equations
"""

##
# 03.01.2006. c
# 10.01.2006
# 21.03.2006
# 21.05.2006
# 21.07.2006
# 01.08.2006
# 11.08.2006
# 12.02.2007
# 20.03.2007
def parseTerms( regions, desc, itps ):
    # Parse.
    termDescs = []
    bnf = createBNF( termDescs, itps )
    try:
        bnf.parseString( desc )
    except:
        print 'cannot parse:\n', desc
        raise
    
    # Construct terms.
    terms = OneTypeList( Term )
    for td in termDescs:
##         print td
##         pause()
        try:
            constructor = termTable[td.name]
        except:
            print "term '%s' is not in %s" % (td.name,
                                              sorted( termTable.keys() ))
            raise ValueError
        region = regions[td.region]
        term = constructor( region, td.name, td.sign )
        term.argNames = td.args
        term.integralName = td.integral
        terms.append( term )

    return terms

##
# c: 27.11.2006, r: 18.01.2008
def setupTermArgs( terms, variables, materials, user = None ):
    """terms ... can be both Terms or Term class
       - checks term argument existence in variables, materials, user
       - checks equality of field and term subdomain lists (igs)"""
    terms.classifyArgs()
    for term in terms:
        igs = term.charFun.igs
        vns = term.getVariableNames()
        for name in vns:
            if name not in variables.names:
                output( 'variable "%s" not found' % name )
                raise IndexError
            field = variables[name].field
            if not set( igs ).issubset( set( field.aps.igs ) ):
                output( ('%s: incompatible regions: (term, field)'
                         + ' (%s(%s) in %s(%s)') %\
                        (term.name, igs, name, field.igs(), field.name) )
                raise ValueError

        mns = term.getMaterialNames()
        for name in mns:
            if name not in materials.names:
                output( 'material "%s" not found' % name )
                raise IndexError
            mat = materials[name]
            if not set( igs ).issubset( set( mat.igs ) ):
                output( ('%s: incompatible regions: (term, material)'
                         + ' (%s(%s) in %s(%s)') %\
                        (term.name, igs, name, mat.igs, mat.name) )
                raise ValueError

    if user is None:
        return

    uns = terms.getUserNames()
    uks = user.keys()
    for name in uns:
        if name not in uks:
            output( 'user data "%s" not found' % name )
            raise IndexError
#        print '********* ok'

##
# 24.07.2006, c
def buildArgs( term, variables, materials, **kwargs ):
    args = kwargs
    vns = term.getVariableNames()
    for vn in vns:
        args[vn] = variables[vn]
    mns = term.getMaterialNames()
    for mn in mns:
        args[mn] = materials[mn]
    return args

##
# 21.07.2006, c
class Equations( Container ):

    ##
    # c: 18.04.2006, r: 20.02.2008
    def fromConf( conf ):
        objs = OneTypeList( Equation )

        conf = copy( conf )
        tps = conf.pop( 'namespaces', {} )
        itps = invertDict( tps, True )

        ii = 0
        for name, desc in conf.iteritems():
            output( 'equation "%s":' %  name )
            output( desc )
            eq = Equation( name = name,
                           desc = desc,
                           itps = itps )
            objs.append( eq )
            ii += 1

        obj = Equations( objs, itps = itps )
        
        return obj
    fromConf = staticmethod( fromConf )

    ##
    # 21.07.2006, c
    # 27.11.2006
    # 29.11.2006
    # 12.02.2007
    # 13.02.2007
    # 27.02.2007
    # 02.03.2007
    def parseTerms( self, regions ):
        self.caches = {}
        for eq in self:
            eq.parseTerms( regions, self.caches )

    ##
    # 21.07.2006, c
    # 24.07.2006
    # 27.11.2006
    # 27.02.2007
    def setupTermArgs( self, variables, materials, user = None ):
        for eq in self:
            eq.setupTermArgs( variables, materials, user )

        self.materials = materials
        self.variables = variables

    ##
    # c: ??, r: 26.02.2008
    def describeGeometry( self, geometries, variables, integrals ):
        output( 'describing geometries...' )
        tt = time.clock()
        for eq in self:
            eq.describeGeometry( geometries, variables, integrals )
        output( '...done in %.2f s' % (time.clock() - tt) )
        

    ##
    # 24.08.2006, c
    # 24.04.2007
    def getTermGeometries( self ):
        tgs = set()
        for eq in self:
            tgs.update( eq.getTermGeometries() )
        return tgs

    ##
    # 16.11.2007, c
    def getTermIntegralNames( self ):
        iNames = set()
        for eq in self:
            iNames.update( eq.getTermIntegralNames() )
        return iNames

    ##
    # 27.02.2007, c
    def invalidateTermCaches( self ):
        for cache in self.caches.itervalues():
            cache.clear()

    ##
    # c: 07.05.2008, r: 07.05.2008
    def resetTermCaches( self ):
        for cache in self.caches.itervalues():
            cache.reset()

    ##
    # 02.03.2007, c
    def setCacheMode( self, cacheOverride ):
        for cache in self.caches.itervalues():
            cache.setMode( cacheOverride )

    ##
    # c: 02.04.2008, r: 02.04.2008
    def initTime( self, ts ):
        for cache in self.caches.itervalues():
            cache.initTime( ts )

    ##
    # 08.06.2007, c
    def advance( self, ts ):
        for cache in self.caches.itervalues():
            cache.advance( ts.step + 1 )

##
# 21.07.2006, c
class Equation( Struct ):

    ##
    # 25.07.2006, c
    # 28.08.2006
    # 12.02.2007
    def fromDesc( name, desc, termPrefixes = None ):
        if termPrefixes is None: termPrefixes = {}

        obj = Equation( name = name, desc = desc,
                        itps = invertDict( termPrefixes, True ) )
        return obj
    fromDesc = staticmethod( fromDesc )

    ##
    # 21.07.2006, c
    # 25.07.2006
    # 01.08.2006
    # 11.08.2006
    # 12.02.2007
    # 27.02.2007
    def parseTerms( self, regions, caches ):
        terms = parseTerms( regions, self.desc, self.itps )
        self.terms = Terms( terms )
        self.assignTermCaches( caches )

    ##
    # 29.11.2006, c
    # 27.02.2007
    # 08.06.2007
    # 11.06.2007
    def assignTermCaches( self, caches ):
        """
        History sizes for a particular cache instance are taken as maximum
        of historySizes requirements of all terms using the instance.
        """
        for term in self.terms:
            if not hasattr( term, 'useCaches' ): continue

#            print term.name
            for name, argLists in term.useCaches.iteritems():
##                print name, argLists
                for args in argLists:
                    # Order should be handled in terms...
                    args = copy( args )
                    if type( args[-1] ) == dict:
                        historySizes = args.pop()
                    else:
                        historySizes = None
                    ans = [term.getArgName( arg ) for arg in args]
                    cname = '_'.join( [name] + ans )
##                     print term.name, name, argLists, args, self.name, cname
##                     print historySizes
##                     debug()
                    if caches.has_key( cname ):
                        caches[cname].mergeHistorySizes( historySizes )
                        continue

#                    print 'new'
                    try:
                        constructor = cacheTable[name]
                    except:
                        raise RuntimeError, 'cache not found! %s in %s'\
                              % (name, sorted( cacheTable.keys() ))
                    caches[cname] = constructor( cname, ans, historySizes )
            term.caches = caches

    ##
    # 21.07.2006, c
    # 24.07.2006
    # 22.08.2006
    # 25.08.2006
    # 27.11.2006
    # 20.02.2007
    def setupTermArgs( self, variables, materials, user = None ):
        """- checks term argument existence in variables, materials, user
           - checks compatability of field and term subdomain lists (igs)"""
        setupTermArgs( self.terms, variables, materials, user )

    ##
    #
    def describeGeometry( self, geometries, variables, integrals ):
        for term in self.terms:
            term.describeGeometry( geometries, variables, integrals )

    ##
    # 24.04.2007, c
    def getTermGeometries( self ):
        tgs = set()
        for term in self.terms:
            for tg in term.getGeometry():
                tgs.add( tg )
        return tgs

    ##
    # 16.11.2007, c
    def getTermIntegralNames( self ):
        iNames = set()
        for term in self.terms:
            iNames.add( term.integralName )
        return iNames
