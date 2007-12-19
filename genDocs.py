#!/usr/bin/env python
import os
import os.path as op
from optparse import OptionParser
import pyparsing as pp

from sfe.base.base import *
from sfe.terms import termTable, cacheTable

##
# 09.11.2007, c
# 13.11.2007
def setSection( sec ):
    def action( str, loc, toks ):
        if toks:
            sec[0] = toks[0][1:-1]
        return toks
    return action

##
# 09.11.2007, c
def toList( slist, sec ):
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
def createBNF( slist, currentSection ):
    
    colon = pp.Literal( ':' )

    section = pp.Combine( colon
                          + pp.Word( pp.alphas, pp.alphanums + '_ ' )
                          + colon )
    section.setParseAction( setSection( currentSection ) )
    section.setName( 'section' )
#    section.setDebug()

    text = pp.SkipTo( section | pp.StringEnd() )
    text.setParseAction( toList( slist, currentSection ) )
    text.setName( 'text' )
#    text.setDebug()

##     doc = pp.StringStart()\
##           + pp.ZeroOrMore( section + text )\
##           + pp.StringEnd()

    doc = pp.StringStart()\
           + pp.Optional( text ) + pp.ZeroOrMore( section + text )\
           + pp.StringEnd()

    return doc

latex = 'pdflatex'
latexOptions = r'-interaction=\nonstopmode'
latexRuns = 3
outputDir = './tmp'

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
"""

header = r"""
\documentclass[10pt,a4paper]{article}

\usepackage{bm}
\usepackage{a4wide}

%s

\setlength{\parindent}{0pt}

\newif\ifPDF
\ifx\pdftexversion\undefined \PDFfalse
\else \PDFtrue
\fi

\ifPDF
   \usepackage[pdftex]{graphicx}
   \DeclareGraphicsExtensions{.pdf,.png}
   \pdfinfo{
    /Title      (Overview of terms in SFE)
    /Author     (R. Cimrman)
  }
\else
  \usepackage{graphicx}
  \DeclareGraphicsExtensions{.eps}
\fi

\begin{document}

\begin{titlepage}
  \begin{center}
    \ \\
    [40mm]
    \Huge Overview of terms in SFE \\
    [10mm]
    \huge Robert Cimrman \& genDocs.py \\
    [20mm]
    \Large
    Department of Mechanics \& New Technologies Research Centre \\
    University of West Bohemia \\
    Plze\v{n}, Czech Republic\\
    [10mm]
%%    \includegraphics[height=1.5cm,viewport=96 122 504 527,clip]{\figDir/zcu}
    \includegraphics[height=0.1\linewidth]{\figDir/znak_zcu}
    \hspace{5mm}
    \includegraphics[height=0.1\linewidth]{\figDir/logo_fav_en}
    \hspace{5mm}
    \includegraphics[height=0.1\linewidth]{\figDir/logoen_bw}
    \vfill
    \today
  \end{center}
\end{titlepage}

\tableofcontents

\newpage

\section{Introduction}


\textbf{Notation}:
\begin{center}
  \begin{tabular}{rl}
    $\Omega$ & volume (sub)domain \\
    $\Gamma$ & surface (sub)domain \\
    $t$ & time \\
    $y$ & any function \\
    $\ul{y}$ & any vector function \\
    $\ul{n}$ & unit outward normal \\
    $q$, $s$ & scalar test function \\
    $p$, $r$ & scalar unknown or parameter function \\
    $\bar{p}$ & scalar parameter function \\
    $\ul{v}$ & vector test function \\
    $\ul{w}$, $\ul{u}$ & vector unknown or parameter function \\
    $\ul{b}$ & vector parameter function \\
    $\ull{e}(\ul{u})$ & Cauchy strain tensor ($\frac{1}{2}((\nabla u) + (\nabla
    u)^T)$) \\
    $\ul{f}$ & vector volume forces \\
    $\rho$ & density \\
    $\nu$ & kinematic viscosity \\
    $c$ & any constant \\
    $\delta_{ij}$, $\ull{I}$ & Kronecker delta, identity matrix \\
  \end{tabular}
\end{center}

The suffix "$_0$" denotes a quatity related to a previous time step.

General syntax of a term call:
\begin{verbatim}
<term_name>.<region>( <arg1>, <arg2>, ... )
\end{verbatim}

In the following, \verb|<virtual>| corresponds to a test function,
\verb|<state>| to a unknown function and \verb|<parameter>| to a known function
arguments.

""" % defines

ending = r"""
\end{document}
"""

itemSection = r"""
\textbf{%s}:
"""

termSyntax = r"""{
\small
\verb|%s.%s( <%s> )|
}"""            

cacheSyntax = r"""{
\small
\begin{verbatim}
cache = term.getCache( '%s', <index> )
data = cache( <data name>, <ig>, <ih>, %s )
\end{verbatim}
}"""            

termDefinition = r"""
\begin{center}
  %s
\end{center}
"""

##
# 29.10.2007, c
# 15.11.2007
def itemsPerSections( table, secNamePrefix, omitList ):
    if omitList is not None:
        omitList = omitList.split()

    ips = {}
    for name, cls in table.iteritems():
        module = cls.__module__.split( '.' )[-1]
        if module in omitList: continue
        secName = '%s in %s' % (secNamePrefix, module)
        ips.setdefault( secName, [] ).append( name )
    return ips

##
# 14.11.2007, c
def typesetTermSyntax( fd, cls, name ):
    fd.write( itemSection % 'Syntax' )
    argTypes = '>, <'.join( cls.argTypes )
    fd.write( termSyntax % (name, '<region>', argTypes) )

##
# 14.11.2007, c
def typesetCacheSyntax( fd, cls, name ):
    argTypes = ', '.join( cls.argTypes )
    fd.write( cacheSyntax % (name, argTypes) )

##
# 15.11.2007, c
def typesetArguments( fd, argsText ):
    fd.write( r'\begin{center}' )
    fd.write( '\n' )
    aux = argsText.split( ',' )
    fd.write( r'\begin{tabular}{rcl}' )
    fd.write( '\n' )
    for argDesc in aux:
        print argDesc
        args = argDesc.split( r':' )
        fd.write( r'%s & \dots & %s \\' % (args[0], args[1]) )
        fd.write( '\n' )

    fd.write( r'\end{tabular}' )
    fd.write( '\n' )
    fd.write( r'\end{center}' )
    fd.write( '\n' )

##
# 14.11.2007, c
def typeset( fd, itemsPerSection, itemTable, typesetSyntax ):
    secList = []
    currentSection = [None]
    bnf = createBNF( secList, currentSection )

    for secName, itemNames in itemsPerSection.iteritems():
        fd.write( r'\section{%s}' % secName )
        fd.write( '\n\n' )

        itemNames.sort()
        for itemName in itemNames:
            itemClass = itemTable[itemName]
            name = itemClass.name.replace( '_', '\_' )
            fd.write( r'\subsection{%s}' % name )
            fd.write( '\n' )
            fd.write( r'\textbf{Class}: %s' % itemClass.__name__ )
            fd.write( '\n\n' )

            doc = itemClass.__doc__
            if doc is not None:
                print doc
                secList[:] = []
                currentSection[0] = None
                out = bnf.parseString( doc )
##                 print secList
##                 pause()

                for secName, secText in secList:
#                    print secName
#                    print secText
                    if secName is None and not secText: continue
                    if secName is not None:
                        fd.write( itemSection % secName.capitalize() )
                        if secName == 'definition':
                            fd.write( termDefinition % secText )
                            fd.write( '\n\n' )
                        elif secName == 'arguments':
                            typesetArguments( fd, secText )
                        else:
                            fd.write( secText )
                            fd.write( '\n\n' )

            typesetSyntax( fd, itemClass, itemName )

usage = """%prog [options]"""
help = {
    'omit' :
    'omit listed sections',
}

# ./genDocs.py --omit="termsAdjointNavierStokes termsHDPM  cachesHDPM  cachesBasic"

##
# 29.10.2007, c
# 04.11.2007
# 13.11.2007
# 14.11.2007
# 15.11.2007
def main():
    parser = OptionParser( usage = usage, version = "%prog" )
    parser.add_option( "", "--omit", metavar = 'list of sections',
                       action = "store", dest = "omitList",
                       default = "", help = help['omit'] )
    (options, args) = parser.parse_args()

    tps = itemsPerSections( termTable, 'Terms', options.omitList )
    cps = itemsPerSections( cacheTable, 'Term caches', options.omitList )

    latexFileName = 'terms.tex'
    latexFileNameComplete = op.join( outputDir, latexFileName )
    print latexFileNameComplete
    
    fd = open( latexFileNameComplete, 'w' )
    fd.write( header )

    typeset( fd, tps, termTable, typesetTermSyntax )
    typeset( fd, cps, cacheTable, typesetCacheSyntax )
            
    fd.write( ending )
    fd.close()

    cdIn = 'cd %s;' % outputDir
    cdOut = 'cd %s;' % os.curdir
    
    cmd = ' '.join( (cdIn,
                     latex, latexOptions, latexFileName, ';',
                     cdOut) )
    for ii in xrange( latexRuns ):
        os.system( cmd )

if __name__ == '__main__':
    main()
