"""
== Region description ==

Regions serve for selection of certain parts of the computational domain (= selection of nodes and elements of a FE mesh). They are used to define the boundary conditions, the domains of terms and materials etc.
"""
from pyparsing import Literal, CaselessLiteral, Word, delimitedList,\
     Group, Optional, ZeroOrMore, nums, alphas, alphanums,\
     Combine, StringStart, StringEnd, Forward, oneOf, ParseException

opCodes = ['OA_SubN', 'OA_SubE', 'OA_AddN', 'OA_AddE',
           'OA_IntersectN', 'OA_IntersectE']
evalCodes = ['E_NIR', 'E_NOS', 'E_NBF', 'E_EOG', 'E_ONIR', 'E_NI']
kwCodes = ['KW_All', 'KW_Region']

##
# 11.05.2006, c
def toStack( stack ):
    def pushFirst( str, loc, toks ):
        if toks:
#            print stack
            stack.append( toks[0] )
#            print toks, '->', stack
#            print ''
        return toks
    return pushFirst

##
# 14.06.2006, c
# 15.06.2006
# 02.05.2007
def replace( what, keep = False ):
    def _replace( str, loc, toks ):
        ret = {'token' : what, 'orig' : []}
        if keep:
            ret['orig'] = list( toks[0] )
        return ret
    return _replace

##
# 02.05.2007, c
def replaceWithRegion( what, rIndex ):
    def _replace( str, loc, toks ):
        ret = {'token' : what, 'orig' : []}

        orig = toks[0]
        rOrig = orig[rIndex]
        if isinstance( rOrig, dict ) and (rOrig['token'] == 'KW_Region'):
            orig = list( orig[:rIndex] ) + rOrig['orig']
        ret['orig'] = orig
        return ret
    return _replace

##
# 14.06.2006, c
def joinTokens( str, loc, toks ):
#    print toks
    return [" ".join( toks[0] )]

##
# 14.06.2006, c
def visitStack( stack, opVisitor, leafVisitor ):

    def visit( stack, level ):
        op = stack.pop()

        token = op['token']
        if token in opCodes:
            res2 = visit( stack, level + 1 )
            res1 = visit( stack, level + 1 )
            return opVisitor( level, op, res1, res2 )

        elif token in evalCodes:
            return leafVisitor( level, op )

        elif token in kwCodes:
            return leafVisitor( level, op )

        else:
            raise ValueError, token
            
    return visit( stack, 0 )


##
# 14.06.2006, c
def printOp( level, op, item1, item2 ):
    print level * '  ' + (': %s' % op)

##
# 14.06.2006, c
def printLeaf( level, op ):
    print level * '  ' + ('< %s' % op)

##
# 14.06.2006, c
def printStack( stack ):
    visitStack( stack, printOp, printLeaf )

##
# 13.06.2006, c
# 14.06.2006
# 15.06.2006
# 20.06.2006
# 12.10.2006
# 19.02.2007
# 02.03.2007
# 10.04.2007
# 02.05.2007
def createBNF( stack ):
    point = Literal( "." )
    e = CaselessLiteral( "E" )
    inumber = Word( nums )
    fnumber = Combine( Word( "+-"+nums, nums ) + 
                       Optional( point + Optional( Word( nums ) ) ) +
                       Optional( e + Word( "+-"+nums, nums ) ) )
    _of = Literal( 'of' )
    _in = Literal( 'in' )
    _by = Literal( 'by' )
    _copy = Literal( 'copy' )

    _mn = Literal( '-n' ).setParseAction( replace( 'OA_SubN' ) )
    _me = Literal( '-e' ).setParseAction( replace( 'OA_SubE' ) )
    _pn = Literal( '+n' ).setParseAction( replace( 'OA_AddN' ) )
    _pe = Literal( '+e' ).setParseAction( replace( 'OA_AddE' ) )
    _inn = Literal( '*n' ).setParseAction( replace( 'OA_IntersectN' ) )
    _ine = Literal( '*e' ).setParseAction( replace( 'OA_IntersectE' ) )
    regop = (_mn | _me | _pn | _pe | _inn | _ine)

    lpar  = Literal( "(" ).suppress()
    rpar  = Literal( ")" ).suppress()

    _all = Literal( 'all' ).setParseAction( replace( 'KW_All' ) )
    node = Literal( 'node' )
    nodes = Literal( 'nodes' )
    elements = Literal( 'elements' )
    group = Literal( 'group' )
    surface = Literal( 'surface' )
    variable = Word( 'xyz', max = 1 )
    anyVar = Word( alphas + '_', alphanums + '_' ) | fnumber

    args = delimitedList( variable, combine = True )\
           + ZeroOrMore( Literal( ',' ) + anyVar )

    function = Word( alphas, alphanums + '_' )\
               + Literal( '(' ) + Optional( args ) + Literal( ')' )
    function = Group( function ).setParseAction( joinTokens )

    region = Combine( Literal( 'r.' ) + Word( alphas, '_' + alphas + nums ) )
    region = Group( Optional( _copy, default = 'nocopy' ) + region )
    region.setParseAction( replace( 'KW_Region', keep = True ) )

    coor = oneOf( 'x y z' )
    boolop = oneOf( '& |' )
    relop = oneOf( '< > <= >= != ==' )
    boolTerm = ZeroOrMore( '(' ) + (coor | fnumber ) + relop + (coor | fnumber)\
               + ZeroOrMore( ')' )
    relation = Forward()
    relation << ZeroOrMore( '(' )\
             + boolTerm + ZeroOrMore( boolop + relation )\
             + ZeroOrMore( ')' )
    relation = Group( relation ).setParseAction( joinTokens )

    nos = Group( nodes + _of + surface ).setParseAction( replace( 'E_NOS' ) )
    nir = Group( nodes + _in + relation ).setParseAction( \
        replace( 'E_NIR', keep = True ) )
    nbf = Group( nodes + _by + function ).setParseAction( \
        replace( 'E_NBF', keep = True ) )
    eog = Group( elements + _of + group + Word( nums ) ).setParseAction( \
        replace( 'E_EOG', keep = True ) )
    onir = Group( node + _in + region ).setParseAction( \
        replaceWithRegion( 'E_ONIR', 2 ) )
    ni = Group( node + inumber ).setParseAction( \
        replace( 'E_NI', keep = True ) )

    regionExpression = Forward()

    atom1 = (_all | region | ni | onir | nos | nir | nbf | eog)
    atom1.setParseAction( toStack( stack ) )
    atom2 = (lpar + regionExpression.suppress() + rpar)
    atom = (atom1 | atom2)

    aux = (regop + regionExpression)
    aux.setParseAction( toStack( stack ) )
    regionExpression << atom + ZeroOrMore( aux )
    regionExpression = StringStart() + regionExpression + StringEnd()

#    region.setDebug()
#    relation.setDebug()
#    regionExpression.setDebug()

    return regionExpression

_testStrs = ['nodes of surface -n r.egion_1',
            'r.egion_2 +n copy r.egion_1',
            'nodes in (y <= 0.00001) & (x < 0.11)',
            'nodes in ((y <= 0.00001) & (x < 0.11))',
            'nodes in (((y <= 0.00001) & (x < 0.11)))',
            'nodes in (((0.00001 < y) & (x < 0.11)))',
            'nodes in (y < 1.0)',
            'all -n nodes in (y == 0.00001)',
            'all -n nodes of surface',
            'all -e r.egion_100',
            'r.egion_1 -n nodes of surface *e r.egion_8 *n nodes in (y > 0)',
            'nodes of surface +n nodes by pokus( x, y, z )',
            'elements of group 6 +e nodes by fn2_3c( x )',
            """r.egion_1 *n (r.egion_2 +e (nodes in (y > 0) *n r.egion_32))
            -n nodes of surface -e r.egion_5""",
            'nodes by noargs()',
            'nodes by extraargs( x, y, z, abc,3 )',
            'node in r.region_3',
            'node 10']

if __name__ == "__main__":
    testStrs = _testStrs

    stack = []
    bnf = createBNF( stack )

    nFail = 0
    for testStr in testStrs:
        print testStr
        stack[:] = []

        try:
            out = bnf.parseString( testStr )
#            print out
#            print stack
        except:
            print '...failed!'
            nFail += 1
            continue

        printStack( stack )
    print 'failed: %d' % nFail
