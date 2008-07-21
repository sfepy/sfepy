from pyparsing import Combine, Literal, Word, delimitedList, Group, Optional,\
     ZeroOrMore, OneOrMore, nums, alphas, alphanums,\
     StringStart, StringEnd, CaselessLiteral


##
# 12.02.2007, c
class TermParse( object ):
    def __str__( self ):
        ss = "%s\n" % self.__class__
        for key, val in self.__dict__.iteritems():
            ss += "  %s:\n    %s\n" % (key, self.__dict__[key])
        return ss

##
# 11.05.2006, c
# 01.08.2006
# 09.08.2006
# 12.02.2007
# 10.04.2007
# 23.04.2007
# 13.11.2007
def collect_term( term_descs, lc, itps ):
    signs = {'+' : 1.0, '-': -1.0}
    def append( str, loc, toks ):
        sign = signs[toks.sign] * signs[lc[0]]
        term_prefix = itps.get( toks.term_desc.name, '' ) 
##        print toks, lc, term_prefix
##         print name, integral, region
        tp = TermParse()
        tp.integral = toks.term_desc.integral
        tp.region = toks.term_desc.region
        tp.flag = toks.term_desc.flag
        tp.sign = sign * float( toks.mul[0] )
        tp.name = term_prefix + toks.term_desc.name
        tp.args = toks.args
#        print tp
        term_descs.append( tp )
    return append

##
# 21.05.2006, c
def rhs( lc ):
    def aux( str, loc, toks ):
        if toks:
#            print toks
            lc[0] = '-'
    return aux
        
##
# c: 21.05.2006, r: 08.07.2008
def create_bnf( term_descs, itps ):
    """term_descs .. list of TermParse objects (sign, term_name, term_arg_names)"""

    lc = ['+'] # Linear combination context.
    equal = Literal( "=" ).setParseAction( rhs( lc ) )
    zero  = Literal( "0" ).suppress()

    point = Literal( "." )
    e = CaselessLiteral( "E" )
    inumber = Word( "+-"+nums, nums )
    fnumber = Combine( Word( "+-"+nums, nums ) + 
                       Optional( point + Optional( Word( nums ) ) ) +
                       Optional( e + Word( "+-"+nums, nums ) ) )
    lbracket = Literal( '[' ).suppress()
    rbracket = Literal( ']' ).suppress()


    ident = Word( alphas, alphanums + "_")
    history = Optional( lbracket + inumber + rbracket, default = 0 )( "history" )
    history.setParseAction( lambda str, loc, toks: int( toks[0] ) )
    variable = Group( Word( alphas, alphanums + '._' ) + history )
    flag = Literal( 'a' )

    term = Optional( Literal( '+' ) | Literal( '-' ), default = '+' )( "sign" )\
           + Optional( fnumber + Literal( '*' ).suppress(),
                       default = '1.0' )( "mul" ) \
           + Combine( ident( "name" )\
                      + Optional( "." + (ident( "integral" ) + "."
                                  + ident( "region" ) + "." + flag( "flag" ) |
                                  ident( "integral" ) + "." + ident( "region" ) |
                                  ident( "region" )
                                  )))( "term_desc" ) + "("\
                                  + Optional( delimitedList( variable ),
                                default = [] )( "args" ) + ")"
    term.setParseAction( collect_term( term_descs, lc, itps ) )

    rhs1 = equal + OneOrMore( term )
    rhs2 = equal + zero
    equation = StringStart() + OneOrMore( term )\
               + Optional( rhs1 | rhs2 ) + StringEnd()

#    term.set_debug()

    return equation
    
    
if __name__ == "__main__":

    test_str = """d_term1.Y( fluid, u, w, Nu, dcf, mode )
                 + 5.0 * d_term2.Omega( u, w, Nu, dcf, mode )
                 - d_another_term.Elsewhere( w, p[-1], Nu, dcf, mode )
                 = - dw_rhs.Y3.a( u, q, Nu, dcf, mode )"""
    
    term_descs = []
    bnf = create_bnf( term_descs, {} )
    out = bnf.parseString( test_str )

    print 'out:', out, '\n'
    for tp in term_descs:
        print tp
