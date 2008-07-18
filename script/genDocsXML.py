#!/usr/bin/env python
import sys
from optparse import OptionParser
import pyparsing as pp

sys.path.append( '.' )
from sfepy.base.base import *
from sfepy.terms import term_table, cache_table

##
# 09.11.2007, c
# 13.11.2007
def set_section( sec ):
    def action( str, loc, toks ):
        if toks:
            sec[0] = toks[0][1:-1]
        return toks
    return action

##
# 09.11.2007, c
def to_list( slist, sec ):
    def action( str, loc, toks ):
        if toks:
##             print toks
##             pause()
            slist.append( (sec[0], toks[0]) )
        return toks
    return action
    
##
# 09.11.2007, c
# 13.11.2007
def create_bnf( slist, current_section ):
    
    colon = pp.Literal( ':' )

    section = pp.Combine( colon
                          + pp.Word( pp.alphas, pp.alphanums + '_ ' )
                          + colon )
    section.set_parse_action( set_section( current_section ) )
    section.set_name( 'section' )
#    section.set_debug()

    text = pp.SkipTo( section | pp.StringEnd() )
    text.set_parse_action( to_list( slist, current_section ) )
    text.set_name( 'text' )
#    text.set_debug()

##     doc = pp.StringStart()\
##           + pp.ZeroOrMore( section + text )\
##           + pp.StringEnd()

    doc = pp.StringStart()\
           + pp.Optional( text ) + pp.ZeroOrMore( section + text )\
           + pp.StringEnd()

    return doc

begining = r"""
<article>
    <articleinfo>
        <title>SfePy Documentation</title>
        <mathinclude>
\setlength{\parindent}{0pt}
\def\dt{{\Delta t}}
\def\pdiff#1#2{\frac{\partial {#1}}{\partial {#2}}}
\def\tdiff#1#2{\frac{{\rm d} {#1}}{{\rm d} {#2}}}
\def\difd#1{\ {\rm d}#1}
\newcommand{\dvg}{\mathop{\rm div}}
\newcommand{\ul}[1]{\underline{#1}}
\newcommand{\uld}[1]{\dot{\underline{#1}}}
\newcommand{\ull}[1]{\underline{\underline{#1}}}
\def\Vcal{\mathcal{V}}
\def\Tcal{\mathcal{T}}
\def\Hcal{\mathcal{H}}
\def\Fcal{\mathcal{F}}
\def\Gcal{\mathcal{G}}
\def\figDir{../doc/tex/figures}
\newcommand{\sfepy}{SfePy}
        </mathinclude>
    </articleinfo>
"""

ending = r"""
</article>
"""

item_section = r"""
<em>%s</em>:
"""

term_syntax = r"""
<command>%s.%s( &lt;%s> )</command>
"""            

cache_syntax = r"""
<p>
<command>cache = term.get_cache( '%s', &lt;index> )</command>
</p><p>
<command>data = cache( &lt;data name>, &lt;ig>, &lt;ih>, %s )</command>
</p>
"""            

term_definition = """<center>%s</center>"""

##
# 29.10.2007, c
# 15.11.2007
def items_per_sections( table, sec_name_prefix, omit_list ):
    if omit_list is not None:
        omit_list = omit_list.split()

    ips = {}
    for name, cls in table.iteritems():
        module = cls.__module__.split( '.' )[-1]
        if module in omit_list: continue
        sec_name = '%s in %s' % (sec_name_prefix, module)
        ips.setdefault( sec_name, [] ).append( name )
    return ips

##
# c: 14.11.2007, r: 24.10.2008
def typeset_term_syntax( fd, cls, name ):
    fd.write( item_section % 'Syntax' )
    arg_types = '>, &lt;'.join( cls.arg_types )
    fd.write( term_syntax % (name, '&lt;i>.&lt;r>', arg_types) )

##
# 14.11.2007, c
def typeset_cache_syntax( fd, cls, name ):
    arg_types = ', '.join( cls.arg_types )
    fd.write( cache_syntax % (name, arg_types) )

##
# 15.11.2007, c
def typeset_arguments( fd, args_text ):
    fd.write( '<center><table><tgroup>\n' )
    aux = args_text.split( ',' )
    fd.write( '<tbody>\n' )
    for arg_desc in aux:
        print arg_desc
        args = arg_desc.split( r':' )
        fd.write( '<row><entry>%s</entry><entry>%s</entry></row>\n' % \
                (args[0], args[1].replace("$", "")) )

    fd.write( '</tbody></tgroup></table></center>\n' )

##
# c: 21.03.2008, r: 21.03.2008
def replace_dollars( text ):
    text = text.replace( '<', '&lt;' )
    while text.find("$") != -1:
        c0 = c1 = True
        text2 = text.replace("$", "<m>", 1)
        if text2 == text:
            c0 = False
        text = text2.replace("$", "</m>", 1)
        if text2 == text:
            c1 = False
        if c0 and not c1:
            print 'missing end $!'
            pause()
            raise ValueError
    return text

##
# c: 26.03.2008, r: 26.03.2008
def include( fd, file_name ):
    fd2 = open( file_name, 'r' )
    fd.write( fd2.read() )
    fd2.close()
    
##
# c: 21.03.2008, r: 07.05.2008
def typeset_item_table( fd, item_table ):
    sec_list = []
    current_section = [None]
    bnf = create_bnf( sec_list, current_section )

    fd.write( r'<section><title>_list of all terms</title>' )
    fd.write( '\n\n' )

    fd.write( '<center><table><tgroup>\n' )
    fd.write( '<tbody>\n' )

    fd.write( '<row><entry>section</entry><entry>name</entry>'
              '<entry format="p{5cm}">definition</entry></row>\n' )
    fd.write( '<row><entry></entry><entry></entry><entry></entry></row>\n' )

    row_format = '<row><entry><a ref="%s"/></entry><entry>%s</entry><entry format="p{5cm}">%s</entry></row>\n'

    keys = item_table.keys()
    keys.sort()
    for key in keys:
        item_class = item_table[key]
        doc = item_class.__doc__
        if doc is not None:
            sec_list[:] = []
            current_section[0] = None
            out = bnf.parse_string( doc )
            dd = [x[1] for x in sec_list if x[0] == 'definition']
            if len( dd ):
                dd = dd[0]
            else:
                dd = ''
            print dd
##             print replace_dollars( dd )
            fd.write( row_format % (item_class.name, item_class.name,
                                   replace_dollars( dd )) )

    fd.write( '</tbody></tgroup></table></center>\n' )
    fd.write( r'</section>' )
    

##
# c: 14.11.2007, r: 07.05.2008
def typeset( fd, items_per_section, item_table, typeset_syntax ):
    sec_list = []
    current_section = [None]
    bnf = create_bnf( sec_list, current_section )

    for sec_name, item_names in items_per_section.iteritems():
        fd.write( r'<section><title>%s</title>' % sec_name )
        fd.write( '\n\n' )

        item_names.sort()
        for item_name in item_names:
            item_class = item_table[item_name]
            name = item_class.name#.replace( '_', '\_' )
            fd.write( '\n<section id="%s">\n    <title>%s</title>\n'\
                      % (name, name) )
            fd.write( '<p>\n' )
            fd.write( r'<em>_class</em>: %s' % item_class.__name__ )
            fd.write( '</p>\n' )

            doc = item_class.__doc__
            if doc is not None:
                print doc
                sec_list[:] = []
                current_section[0] = None
                out = bnf.parse_string( doc )
##                 print sec_list
##                 pause()

                for sec_name, sec_text in sec_list:
#                    print sec_name
#                    print sec_text
                    if sec_name is None and not sec_text: continue
                    if sec_name is not None:
                        fd.write("\n<p>")
                        fd.write( item_section % sec_name.capitalize() )
                        if sec_name == 'definition':
                            sec_text = replace_dollars( sec_text )
                            fd.write( term_definition % sec_text )
                        elif sec_name == 'arguments':
                            sec_text = replace_dollars( sec_text )
                            fd.write("\n<p>")
                            typeset_arguments( fd, sec_text )
                            fd.write("\n</p>")
                        else:
                            while sec_text.find("$") != -1:
                                sec_text = sec_text.replace("$", "<m>", 1)
                                sec_text = sec_text.replace("$", "</m>", 1)
                            fd.write( sec_text )
                        fd.write("\n</p>")

            fd.write("\n<p>")
            typeset_syntax( fd, item_class, item_name )
            fd.write("\n</p>")
            fd.write( r'</section>' )
        fd.write( r'</section>' )

usage = """%prog [options]"""
help = {
    'omit' :
    'omit listed sections',
    'output_file_name' :
    'output file name',
}

# ./gen_docs.py --omit="terms_adjoint_navier_stokes terms_hdpm  caches_hdpm  caches_basic"

##
# c: 29.10.2007, r: 26.03.2008
def main():
    parser = OptionParser( usage = usage, version = "%prog" )
    parser.add_option( "", "--omit", metavar = 'list of sections',
                       action = "store", dest = "omit_list",
                       default = "", help = help['omit'] )
    parser.add_option( "-o", "--output", metavar = 'output_file_name',
                       action = "store", dest = "output_file_name",
                       default = "terms.xml", help = help['output_file_name'] )
    (options, args) = parser.parse_args()

    tps = items_per_sections( term_table, 'Terms', options.omit_list )
    cps = items_per_sections( cache_table, 'Term caches', options.omit_list )

    fd = open( options.output_file_name, 'w' )
    fd.write( begining )

    include( fd, 'doc/pages/notation.xml' )
    typeset_item_table( fd, term_table )
    include( fd, 'doc/pages/intro.xml' )
    typeset( fd, tps, term_table, typeset_term_syntax )
    typeset( fd, cps, cache_table, typeset_cache_syntax )
            
    fd.write( ending )
    fd.close()

if __name__ == '__main__':
    main()
