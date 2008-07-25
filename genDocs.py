#!/usr/bin/env python
import os
import os.path as op
from optparse import OptionParser
import pyparsing as pp

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

latex = 'pdflatex'
latex_options = r'-interaction=\nonstopmode'
latex_runs = 3
output_dir = './tmp'

defines = r"""
\def\dt{{\Delta t}}
\def\pdiff#1#2{\frac{\partial {#1}}{\partial {#2}}}
\def\difd#1{\ {\rm d}#1}
\newcommand{\dvg}{\mathop{\rm div}}
\newcommand{\ul}[1]{\underline{#1}}
\newcommand{\uld}[1]{\dot{\underline{#1}}}
\newcommand{\ull}[1]{\underline{\underline{#1}}}
\def\Vcal{\mathcal{V}}
\def\Tcal{\mathcal{T}}
\def\figDir{../doc/tex/figures}
\newcommand{\sfepy{SFE}
"""

header = r"""
\documentclass[10pt,a4paper]{article}

\usepackage{bm}
\usepackage{a4wide}

%s
"""  % defines

begining = r"""
\setlength{\parindent}{0pt}

\newif\ifPDF
\ifx\pdftexversion\undefined \PDFfalse
\else \PDFtrue
\fi

\ifPDF
   \usepackage[pdftex]{graphicx}
   \DeclareGraphicsExtensions{.pdf,.png}
   \pdfinfo{
    /Title      (Overview of terms in \sfepy{})
    /Author     (R. Cimrman)
  }
\else
  \usepackage{graphicx}
  \DeclareGraphicsExtensions{.eps}
\fi

\begin{document}

%s

%%\setcounter{page}{4}

\tableofcontents

\newpage

%s
"""

ending = r"""
\end{document}
"""

item_section = r"""
\textbf{%s}:
"""

term_syntax = r"""{
\small
\verb|%s.%s( <%s> )|
}"""            

cache_syntax = r"""{
\small
\begin{verbatim}
cache = term.get_cache( '%s', <index> )
data = cache( <data name>, <ig>, <ih>, %s )
\end{verbatim}
}"""            

term_definition = r"""
\begin{center}
  %s
\end{center}
"""

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
    arg_types = '>, <'.join( cls.arg_types )
    fd.write( term_syntax % (name, '<i>.<r>', arg_types) )

##
# 14.11.2007, c
def typeset_cache_syntax( fd, cls, name ):
    arg_types = ', '.join( cls.arg_types )
    fd.write( cache_syntax % (name, arg_types) )

##
# 15.11.2007, c
def typeset_arguments( fd, args_text ):
    fd.write( r'\begin{center}' )
    fd.write( '\n' )
    aux = args_text.split( ',' )
    fd.write( r'\begin{tabular}{rcl}' )
    fd.write( '\n' )
    for arg_desc in aux:
        print arg_desc
        args = arg_desc.split( r':' )
        fd.write( r'%s & \dots & %s \\' % (args[0], args[1]) )
        fd.write( '\n' )

    fd.write( r'\end{tabular}' )
    fd.write( '\n' )
    fd.write( r'\end{center}' )
    fd.write( '\n' )

##
# 14.11.2007, c
def typeset( fd, items_per_section, item_table, typeset_syntax ):
    sec_list = []
    current_section = [None]
    bnf = create_bnf( sec_list, current_section )

    for sec_name, item_names in items_per_section.iteritems():
        fd.write( r'\section{%s}' % sec_name )
        fd.write( '\n\n' )

        item_names.sort()
        for item_name in item_names:
            item_class = item_table[item_name]
            name = item_class.name.replace( '_', '\_' )
            fd.write( r'\subsection{%s}' % name )
            fd.write( '\n' )
            fd.write( r'\textbf{_class}: %s' % item_class.__name__ )
            fd.write( '\n\n' )

            doc = item_class.__doc__
            if doc is not None:
                print doc
                sec_list[:] = []
                current_section[0] = None
                out = bnf.parseString( doc )
##                 print sec_list
##                 pause()

                for sec_name, sec_text in sec_list:
#                    print sec_name
#                    print sec_text
                    if sec_name is None and not sec_text: continue
                    if sec_name is not None:
                        fd.write( item_section % sec_name.capitalize() )
                        if sec_name == 'definition':
                            fd.write( term_definition % sec_text )
                            fd.write( '\n\n' )
                        elif sec_name == 'arguments':
                            typeset_arguments( fd, sec_text )
                        else:
                            fd.write( sec_text )
                            fd.write( '\n\n' )

            typeset_syntax( fd, item_class, item_name )

usage = """%prog [options]"""
help = {
    'omit' :
    'omit listed sections',
}

# ./gen_docs.py --omit="terms_adjoint_navier_stokes terms_hdpm  caches_hdpm  caches_basic"

##
# c: 29.10.2007, r: 24.01.2008
def main():
    parser = OptionParser( usage = usage, version = "%prog" )
    parser.add_option( "", "--omit", metavar = 'list of sections',
                       action = "store", dest = "omit_list",
                       default = "", help = help['omit'] )
    (options, args) = parser.parse_args()

    tps = items_per_sections( term_table, 'Terms', options.omit_list )
    cps = items_per_sections( cache_table, 'Term caches', options.omit_list )

    latex_filename = 'terms.tex'
    latex_filename_complete = op.join( output_dir, latex_filename )
    print latex_filename_complete

    fd = open( 'doc/pages/title_sfepy.tex', 'r' )
    title_src = fd.read()
    fd.close()

    fd = open( 'doc/pages/intro.tex', 'r' )
    intro_src = fd.read()
    fd.close()

    
    fd = open( latex_filename_complete, 'w' )
    fd.write( header )
    fd.write( begining % (title_src, intro_src) )

    typeset( fd, tps, term_table, typeset_term_syntax )
    typeset( fd, cps, cache_table, typeset_cache_syntax )
            
    fd.write( ending )
    fd.close()

    cd_in = 'cd %s;' % output_dir
    cd_out = 'cd %s;' % os.curdir
    
    cmd = ' '.join( (cd_in,
                     latex, latex_options, latex_filename, ';',
                     cd_out) )
    for ii in xrange( latex_runs ):
        os.system( cmd )

if __name__ == '__main__':
    main()
