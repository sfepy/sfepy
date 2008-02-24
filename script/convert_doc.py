#! /usr/bin/python

"""
This module uses lxml heavily. Read the excellent documentation at

http://codespeak.net/lxml

that explains everything needed to understand how this file works.
"""

import os
import sys
import tempfile
from optparse import OptionParser

import pexpect
from lxml.etree import parse, SubElement, Element, ElementTree, Comment, dump
from lxml.builder import E

sys.path.append( '.' )
from sfe.base.progressbar import progressbar
#from style import style_string

def replace(old, new):
    p = old.getparent()
    i = p.index(old)
    p.remove(old)
    p.insert(i, neweq)

def convert_xml_docbook(infile, outfile):
    """Converts our "simplified DocBook" to the regular DocBook.

    Converted things:

    <a ref="s"/>    ->    <xref linkend="s"/>
    <p>ss</p> -> <para>ss</para>
    <m>x^2</m> -> <inlineequation><alt>x^2</alt></inlineequation>
    <e>x^2</e> -> <informalequation><alt>x^2</alt></informalequation>
    <e id="ss">x^2</e> -> <equation id="ss"><alt>x^2</alt></equation>

    All whitespace is normalized like this:
    <a> ss </a>   -> <s>ss</a>
    <a> ss <b> x </b>   </a>   -> <a>ss <b>x</b></a>
    <a> ss <b> x </b> y </a>   -> <a>ss <b>x</b> y</a>

    """
    tree = parse(infile)
    root = tree.getroot()

    for e in root.getiterator("m"):
        e.tag = "inlineequation"
        SubElement(e, "alt").text = e.text
        e.text = ""

    for e in root.getiterator("e"):
        if "id" in e.attrib:
            e.tag = "equation"
            SubElement(e, "alt").text = e.text
        else:
            e.tag = "informalequation"
            SubElement(e, "alt").text = e.text
        e.text = ""
        if e.tail:
            e.tail = e.tail.lstrip()

    for e in root.getiterator("p"):
        e.tag = "para"

    for e in root.getiterator("a"):
        e.tag = "xref"
        e.attrib["linkend"] = e.attrib["ref"]
        del e.attrib["ref"]


    def remove_whitespace(text):
        text = text.strip()
        text = text.replace("\n", " ")
        while text.find("  ") != -1:
            text = text.replace("  ", " ")
        return text
    for e in root.getiterator():
        if e.text:
            preserve_space_end = (e.text[-1] in [" ", "\n"]) and len(e) > 0
            e.text = remove_whitespace(e.text)
            if preserve_space_end and len(e.text) > 0:
                e.text = e.text + " "
        if e.tail:
            preserve_space_beginning = e.tail[0] in [" ", "\n"]
            preserve_space_end = (e.tail[-1] in [" ", "\n"]) and \
                    (e.getnext() is not None)
            e.tail = remove_whitespace(e.tail)
            if len(e.tail) > 0:
                if preserve_space_beginning: 
                    e.tail = " " + e.tail
                if preserve_space_end: 
                    e.tail = e.tail + " "

    tree.write(outfile)

class Converter(object):

    def __init__(self, root):
        self.root = root

class LaTeXConverter(Converter):
    """
    Converts DocBook to LaTeX. 

    It doesn't touch whitespace (if you want your whitespace to be normalized,
    use the convert_docbook converter).
    """

    def __init__(self, root):
        Converter.__init__(self, root)
        self.data = {}
        self.section_level = 0
        self.align = False

    def convert(self):
        return self.convert_node(self.root)

    def escape(self, text):
        replacements = {
                "&": r"\&",
                ">": r"\hbox{$>$}",
                "<": r"\hbox{$<$}",
                "_": r"\_",
                }
        for old, new in replacements.iteritems():
            text = text.replace(old, new)
        return text

    def convert_node(self, node):
        """Converts any type of node."""
        if node.tag == "article":
            return self.convert_article(node)
        elif node.tag == "articleinfo":
            return self.convert_articleinfo(node)
        elif node.tag == "title":
            return self.convert_articleinfo(node, onlytitle=True)
        elif node.tag == "section":
            return self.convert_section(node)
        elif node.tag == "sect1":
            return self.convert_section(node, 1)
        elif node.tag == "sect2":
            return self.convert_section(node, 2)
        elif node.tag == "sect3":
            return self.convert_section(node, 3)
        elif node.tag == Comment:
            return self.convert_comment(node)
        elif node.tag == "para":
            return self.convert_para(node)
        elif node.tag == "note":
            return self.convert_note(node)
        elif node.tag == "tip":
            return self.convert_tip(node)
        elif node.tag == "programlisting":
            return self.convert_programlisting(node)
        elif node.tag == "em":
            return self.convert_em(node)
        elif node.tag == "command":
            return self.convert_command(node)
        elif node.tag == "inlineequation":
            return self.convert_inlineequation(node)
        elif node.tag == "informalequation":
            return self.convert_equation(node, numbered=False)
        elif node.tag == "equation":
            return self.convert_equation(node, numbered=True)
        elif node.tag == "xref":
            return self.convert_xref(node)
        elif node.tag == "align":
            return self.convert_align(node)
        elif node.tag == "indexterm":
            return self.convert_skip(node)
        elif node.tag in ["email", "ulink", "application", "table",
                "itemizedlist", "guibutton", "guilabel", "guimenu", 
                "menuchoice", "variablelist", "keycombo", "orderedlist",
                "figure", "guimenuitem"]:
            return self.convert_text_tail(node)
        return self.convert_default(node)

    def convert_articleinfo_node(self, node):
        """Converts any type of node."""
        if node.tag == "title":
            return self.convert_articleinfo_title(node)
        elif node.tag == Comment:
            return self.convert_comment(node)
        elif node.tag == "abstract":
            return self.convert_articleinfo_abstract(node)
        elif node.tag in ["releaseinfo", "revhistory", "authorgroup",
                "publisher", "copyright"]:
            return self.convert_skip(node)
        return self.convert_default(node)

    def convert_default(self, node):
        print "Unimplemented element %s. Skipped." % (node.tag)
        dump(node)
        return self.convert_text_tail(node)

    def convert_article(self, node):
        assert node.tag == "article"
        self.check_zero_tail(node.tail)
        r = r"""\documentclass[12pt]{article}
\usepackage{amsmath}
"""
        for x in node:
            r += self.convert_node(x)
        r += """\n\\vfill\n\end{document}\n"""

        return r

    def convert_articleinfo(self, node, onlytitle=False):
        if onlytitle:
            assert node.tag == "title"
            self.check_zero_tail(node.tail)
            r = "\\title{%s}\n" % self.escape(node.text)
            r += "\\begin{document}\n"
            r += "\\maketitle\n"
        else:
            assert node.tag == "articleinfo"
            self.check_zero_tail(node.tail)
            for x in node:
                # this doesn't produce any output, only fills the self.data dict
                self.convert_articleinfo_node(x)
            r = ""
            if "title" in self.data:
                r += "\\title{%s}\n" % self.data["title"]

            r += "\\begin{document}\n"
            if "title" in self.data:
                r += "\\maketitle\n"
            if "abstract" in self.data:
                r += "\\begin{abstract}\n%s\n\\end{abstract}\n" % \
                        self.data["abstract"]
        r += "\\tableofcontents\n"
        return r

    def default_label(self, node):
        if "id" in node.attrib:
            label = "\label{%s}" % node.attrib["id"]
        else:
            label = ""
        return label

    def convert_articleinfo_title(self, node):
        assert node.tag == "title"
        self.check_zero_tail(node.tail)
        self.data["title"] = node.text

    def convert_articleinfo_abstract(self, node):
        assert node.tag == "abstract"
        self.check_zero_tail(node.tail)
        self.data["abstract"] = node.text
        r = ""
        for x in node: 
            r += self.convert_node(x)
        self.data["abstract"] = r

    def convert_section(self, node, level = None):
        self.section_level += 1
        if level is None:
            assert node.tag == "section"
        else:
            assert level == self.section_level
            assert node.tag == "sect%d" % self.section_level
        self.check_zero_tail(node.tail)
        title = node.find("title")
        label = self.default_label(node)
        if self.section_level == 1:
            r = """\n\section{%s}%s\n""" % (self.escape(title.text), label)
        elif self.section_level == 2:
            r = """\n\subsection{%s}%s\n""" % (self.escape(title.text), label)
        elif self.section_level == 3:
            r = """\n\subsubsection{%s}%s\n""" % (self.escape(title.text), label)
        else:
            r = ""
        node.remove(title)
        for x in node: 
            r += self.convert_node(x)
        self.section_level -= 1
        return r

    def convert_comment(self, node):
        assert node.tag == Comment
        self.check_zero_tail(node.tail)
        r = "%% %s\n" % node.text.replace("\n", "\n% ")
        return r

    def convert_para(self, node):
        assert node.tag == "para"
        self.check_zero_tail(node.tail)
        r = "\n"
        if node.text:
            r += self.escape(node.text)
        for x in node: 
            r += self.convert_node(x)
        r += "\n"
        return r

    def convert_align(self, node):
        """
        Rules:

        <align sep="ss"> </align> 

        The "ss" is a separator, before which the & should be inserted (first
        match). If not found, & is put at the and of the equation. Default
        value for "sep" is "=".

        TODO: maybe find a better name for "sep".
        """
        assert node.tag == "align"
        r = ""
        if node.text:
            r += self.escape(node.text)
        if "sep" in node.attrib:
            self.align_sep = node.attrib["sep"]
        else:
            self.align_sep = "="
        r += "\n\\begin{eqnarray*}"
        self.align = True
        for x in node: 
            r += self.convert_node(x)
        self.align = False
        r += "\n\\end{eqnarray*}"
        if node.tail:
            r += "\n" + self.escape(node.tail)
        return r

    def convert_note(self, node):
        assert node.tag == "note"
        self.check_zero_tail(node.tail)
        r = "\n"
        title = node.find("title")
        if title is not None:
            r += """{\\bf %s}\n""" % self.escape(title.text)
            node.remove(title)
        if node.text:
            r += node.text
        for x in node: 
            r += self.convert_node(x)
        return r

    def convert_em(self, node):
        assert node.tag == "em"
        r = r"\textbf{"
        if node.text is not None:
            r += self.escape(node.text)
        r += "}"
        r += self.default_label(node)
        if node.tail is not None:
            r += self.escape(node.tail)
        return r

    def convert_command(self, node):
        assert node.tag == "command"
        r = r"\texttt{"
        if node.text is not None:
            r += self.escape(node.text)
        r += "}"
        r += self.default_label(node)
        if node.tail is not None:
            r += self.escape(node.tail)
        return r

    def convert_tip(self, node):
        assert node.tag == "tip"
        self.check_zero_tail(node.tail)
        r = "\n"
        if node.text:
            r += node.text
        for x in node: 
            r += self.convert_node(x)
        return r

    def convert_programlisting(self, node):
        assert node.tag == "programlisting"
        r = "\n"
        if node.text:
            r += self.escape(node.text)
        if node.tail:
            r += self.escape(node.tail)
        return r

    def convert_inlineequation(self, node):
        assert node.tag == "inlineequation"
        r = "$%s$" % node.find("alt").text
        if node.tail:
            r += self.escape(node.tail)
        return r

    def convert_equation(self, node, numbered=True):
        if self.align:
            r = "\n%s" % (node.find("alt").text)
            align_sep = self.align_sep
            if "sep" in node.attrib:
                align_sep = node.attrib["sep"]
            i = r.find(align_sep)
            if i != -1:
                r = r[:i] + "&" + r[i:i+len(align_sep)] + "&" + \
                        r[i+len(align_sep):]
            else:
                r = r + "&"
            r = r + r" \\"
        else:
            if numbered:
                assert node.tag == "equation"
                label = self.default_label(node)
                r = "\n\\begin{equation}\n  %s  %s\n\\end{equation}" \
                        % (node.find("alt").text, label)
            else:
                assert node.tag == "informalequation"
                r = """\n\\begin{equation*}\n  %s\n\\end{equation*}""" % \
                        node.find("alt").text
        if node.tail:
            r += "\n"+self.escape(node.tail)
        return r

    def convert_xref(self, node):
        assert node.tag == "xref"
        linkend = node.attrib["linkend"]
        r = "\\ref{%s}" % linkend
        if node.tail:
            r += self.escape(node.tail)
        return r

    def convert_skip(self, node):
        return ""

    def convert_text_tail(self, node):
        r = ""
        if node.text is not None:
            r += self.escape(node.text)
        r += self.default_label(node)
        if node.tail is not None:
            r += self.escape(node.tail)
        return r

    def check_zero_tail(self, tail):
        if tail is None:
            return
        if tail.strip() != "":
            print "Unhandled tail: '%s'" % tail

class SfePyDocConverter(LaTeXConverter):

    def convert_article(self, node):
        assert node.tag == "article"
        self.check_zero_tail(node.tail)
        r = r"""\documentclass[10pt]{article}
\usepackage{amsmath}
\setlength{\parindent}{0pt}
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
\newcommand{\sfe}{SFE}
"""
        for x in node:
            r += self.convert_node(x)
        r += """\n\\vfill\n\end{document}\n"""

        return r

    def convert_xref(self, node):
        assert node.tag == "xref"
        linkend = node.attrib["linkend"]
        r = "(\\ref{%s})" % linkend
        if node.tail:
            r += self.escape(node.tail)
        return r

    def convert_command(self, node):
        assert node.tag == "command"
        r = r"{\small\verb|"
        if node.text is not None:
            r += node.text
        r += "|}"
        r += self.default_label(node)
        if node.tail is not None:
            r += self.escape(node.tail)
        return r

class XHTMLConverter(Converter):
    """
    Converts DocBook to HTML. 

    It doesn't touch whitespace (if you want your whitespace to be normalized,
    use the convert_docbook converter).
    """

    def __init__(self, root):
        Converter.__init__(self, root)
        self.data = {}
        self.section_level = 0
        self.align = False
        self.eqlist = []

    def convert(self):
        return self.convert_node(self.root)

    def escape(self, text):
        replacements = {
                "&": r"\&",
                ">": r"\hbox{$>$}",
                }
        for old, new in replacements.iteritems():
            text = text.replace(old, new)
        return text

    def convert_node(self, node):
        """Converts any type of node."""
        if node.tag == "article":
            return self.convert_article(node)
        elif node.tag == "articleinfo":
            return self.convert_articleinfo(node)
        elif node.tag == "title":
            return self.convert_articleinfo(node, onlytitle=True)
        elif node.tag == "section":
            return self.convert_section(node)
        elif node.tag == "sect1":
            return self.convert_section(node, 1)
        elif node.tag == "sect2":
            return self.convert_section(node, 2)
        elif node.tag == "sect3":
            return self.convert_section(node, 3)
        elif node.tag == Comment:
            return self.convert_comment(node)
        elif node.tag == "para":
            return self.convert_para(node)
        elif node.tag == "note":
            return self.convert_note(node)
        elif node.tag == "tip":
            return self.convert_tip(node)
        elif node.tag == "programlisting":
            return self.convert_programlisting(node)
        elif node.tag == "inlineequation":
            return self.convert_inlineequation(node)
        elif node.tag == "informalequation":
            return self.convert_equation(node, numbered=False)
        elif node.tag == "equation":
            return self.convert_equation(node, numbered=True)
        elif node.tag == "xref":
            return self.convert_xref(node)
        elif node.tag == "align":
            return self.convert_align(node)
        elif node.tag == "indexterm":
            return self.convert_skip(node)
        elif node.tag in ["email", "ulink", "application", "table",
                "itemizedlist", "guibutton", "guilabel", "guimenu", 
                "menuchoice", "variablelist", "keycombo", "orderedlist",
                "figure", "command", "guimenuitem"]:
            return self.convert_text_tail(node)
        return self.convert_default(node)

    def convert_articleinfo_node(self, node):
        """Converts any type of node."""
        if node.tag == "title":
            return self.convert_articleinfo_title(node)
        elif node.tag == Comment:
            return self.convert_comment(node)
        elif node.tag == "abstract":
            return self.convert_articleinfo_abstract(node)
        elif node.tag in ["releaseinfo", "revhistory", "authorgroup",
                "publisher", "copyright"]:
            return self.convert_skip(node)
        return self.convert_default(node)

    def convert_default(self, node):
        print "Unimplemented element %s. Skipped." % (node.tag)
        dump(node)
        return self.convert_text_tail(node)

    def convert_article(self, node):
        assert node.tag == "article"
        self.check_zero_tail(node.tail)
        r = r"""<html>"""
        for x in node:
            r += self.convert_node(x)
        r += """</html>"""

        return r

    def convert_articleinfo(self, node, onlytitle=False):
        if onlytitle:
            assert node.tag == "title"
            self.check_zero_tail(node.tail)
            r = "<head><title>%s</title></head>" % self.escape(node.text)
        else:
            assert node.tag == "articleinfo"
            self.check_zero_tail(node.tail)
            for x in node:
                # this doesn't produce any output, only fills the self.data dict
                self.convert_articleinfo_node(x)
            r = ""
            if "title" in self.data:
                r += "\\title{%s}\n" % self.data["title"]

            r += "\\begin{document}\n"
            if "title" in self.data:
                r += "\\maketitle\n"
            if "abstract" in self.data:
                r += "\\begin{abstract}\n%s\n\\end{abstract}\n" % \
                        self.data["abstract"]
        r += "\n"
        return r

    def default_label(self, node):
        if "id" in node.attrib:
            label = ' id="%s"' % node.attrib["id"]
        else:
            label = ""
        return label

    def convert_articleinfo_title(self, node):
        assert node.tag == "title"
        self.check_zero_tail(node.tail)
        self.data["title"] = node.text

    def convert_articleinfo_abstract(self, node):
        assert node.tag == "abstract"
        self.check_zero_tail(node.tail)
        self.data["abstract"] = node.text
        r = ""
        for x in node: 
            r += self.convert_node(x)
        self.data["abstract"] = r

    def convert_section(self, node, level = None):
        self.section_level += 1
        if level is None:
            assert node.tag == "section"
        else:
            assert level == self.section_level
            assert node.tag == "sect%d" % self.section_level
        self.check_zero_tail(node.tail)
        title = node.find("title")
        label = self.default_label(node)
        r = """<h%d%s>%s</h%d>\n""" % (self.section_level, label, 
                self.escape(title.text), self.section_level)
        node.remove(title)
        for x in node: 
            r += self.convert_node(x)
        self.section_level -= 1
        return r

    def convert_comment(self, node):
        assert node.tag == Comment
        self.check_zero_tail(node.tail)
        r = "%% %s\n" % node.text.replace("\n", "\n% ")
        return r

    def convert_para(self, node):
        assert node.tag == "para"
        self.check_zero_tail(node.tail)
        r = "\n<p>\n"
        if node.text:
            r += self.escape(node.text)
        for x in node: 
            r += self.convert_node(x)
        r += "\n</p>"
        return r

    def convert_align(self, node):
        """
        Rules:

        <align sep="ss"> </align> 

        The "ss" is a separator, before which the & should be inserted (first
        match). If not found, & is put at the and of the equation. Default
        value for "sep" is "=".

        TODO: maybe find a better name for "sep".
        """
        assert node.tag == "align"
        r = ""
        if node.text:
            r += self.escape(node.text)
        if "sep" in node.attrib:
            self.align_sep = node.attrib["sep"]
        else:
            self.align_sep = "="
        #r += "\n\\begin{eqnarray*}"
        self.align = True
        for x in node: 
            r += self.convert_node(x)
        self.align = False
        #r += "\n\\end{eqnarray*}"
        if node.tail:
            r += "\n" + self.escape(node.tail)
        return r

    def convert_note(self, node):
        assert node.tag == "note"
        self.check_zero_tail(node.tail)
        r = "\n"
        title = node.find("title")
        if title is not None:
            r += """{\\bf %s}\n""" % self.escape(title.text)
            node.remove(title)
        if node.text:
            r += node.text
        for x in node: 
            r += self.convert_node(x)
        return r

    def convert_tip(self, node):
        assert node.tag == "tip"
        self.check_zero_tail(node.tail)
        r = "\n"
        if node.text:
            r += node.text
        for x in node: 
            r += self.convert_node(x)
        return r

    def convert_programlisting(self, node):
        assert node.tag == "programlisting"
        r = "\n"
        if node.text:
            r += self.escape(node.text)
        if node.tail:
            r += self.escape(node.tail)
        return r

    def convert_inlineequation(self, node):
        assert node.tag == "inlineequation"
        eq = node.find("alt").text
        eqfile = "figures/eq%d.png" % len(self.eqlist)
        self.eqlist.append((eq, eqfile, True))
        r = """<img src="%s"/>""" % eqfile
        if node.tail:
            r += self.escape(node.tail)
        return r

    def convert_equation(self, node, numbered=True):
        if self.align:
            eq = node.find("alt").text
            eqfile = "figures/eq%d.png" % len(self.eqlist)
            self.eqlist.append((eq, eqfile, False))
            r = """<br/><img src="%s"/><br/>""" % eqfile
        else:
            if numbered:
                assert node.tag == "equation"
                label = self.default_label(node)
                eq = node.find("alt").text
                eqfile = "figures/eq%d.png" % len(self.eqlist)
                self.eqlist.append((eq, eqfile, False))
                r = """<br/><img src="%s"%s/><br/>""" % (eqfile, label)
            else:
                assert node.tag == "informalequation"
                eq = node.find("alt").text
                eqfile = "figures/eq%d.png" % len(self.eqlist)
                self.eqlist.append((eq, eqfile, False))
                r = """<br/><img src="%s"/><br/>""" % eqfile
        if node.tail:
            r += "\n"+self.escape(node.tail)
        return r

    def convert_xref(self, node):
        assert node.tag == "xref"
        linkend = node.attrib["linkend"]
        r = '<a href="#%s">LINK</a>' % linkend
        if node.tail:
            r += self.escape(node.tail)
        return r

    def convert_skip(self, node):
        return ""

    def convert_text_tail(self, node):
        r = ""
        if node.text is not None:
            r += self.escape(node.text)
        r += self.default_label(node)
        if node.tail is not None:
            r += self.escape(node.tail)
        return r

    def check_zero_tail(self, tail):
        if tail is None:
            return
        if tail.strip() != "":
            print "Unhandled tail: '%s'" % tail


def convert_docbook_latex(infile, outfile):
    print "Converting text..."
    tree = parse(infile)
    root = tree.getroot()

    f = open(outfile, "w")
    #f.write(LaTeXConverter(root).convert())
    f.write(SfePyDocConverter(root).convert())

def create_image(filename, eq, inline=False):
    r"""Runs "eq" through TeX and saves the result into the "filename" as a png.

    Example:

    create_image("figures/eq1.png", "x^2 + y^2 = z^2")

    With inline=True, it is transormed to "x^2 + y^2 = z^2", otherwise to
    "\begin{equation*}x^2 + y^2 = z^2\end{equation*}".
    """

    filename, ext = os.path.splitext(filename)
    assert ext == ".png"
    filename = os.path.abspath(filename)


    euler = False

    if not euler:
        format = r"""\documentclass[12pt]{article}
                     \usepackage{amsmath}
                     \begin{document}
                     \pagestyle{empty}
                     %s
                     \vfill
                     \end{document}
                 """
    else:
        format = r"""\documentclass[12pt]{article}
                     \usepackage{amsmath}
                     \usepackage{eulervm}
                     \begin{document}
                     \pagestyle{empty}
                     %s
                     \vfill
                     \end{document}
                 """

    def latex(e, inline=True):
        if inline:
            return "$" + e + "$"
        else:
            return r"\begin{equation*}%s\end{equation*}" % e

    tmp = tempfile.mktemp()

    tex = open(tmp + ".tex", "w")
    tex.write(format % latex(eq, inline=inline))
    tex.close()

    cwd = os.getcwd()
    try:
        os.chdir(tempfile.gettempdir())

        def execute(cmd):
            return pexpect.run(cmd, withexitstatus=True)[1]

        if execute("latex -halt-on-error %s.tex" % tmp) != 0:
            raise SystemError("Failed to generate DVI output.")

        os.remove(tmp + ".tex")
        os.remove(tmp + ".aux")
        os.remove(tmp + ".log")

        command = "dvipng -T tight -z 9 --truecolor -o %s.png %s.dvi"

        output = "png"
        pardir = filename.rstrip(os.path.basename(filename))
        if not os.path.exists(pardir):
            os.mkdir(pardir)
        try:
            if execute(command % (filename, tmp)) != 0:
                raise SystemError("Failed to generate output.")
            else:
                os.remove(tmp + ".dvi")
        except KeyError:
            raise SystemError("Invalid output format: %s" % output)
    finally:
        os.chdir(cwd)

def convert_docbook_xhtml(infile, outfile):
    print "Converting text..."
    tree = parse(infile)
    root = tree.getroot()

    f = open(outfile, "w")
    c = XHTMLConverter(root)
    f.write(c.convert())
    p = progressbar("Creating images", maxval = len(c.eqlist))
    i = 0
    for eq, eq_file, inline in c.eqlist:
        try:
            create_image(eq_file, eq, inline)
        except SystemError:
            print "Warning: file %s, equation '%s' failed to generate." % \
                    (eq_file, eq)
            pass
        i += 1
        p.update(i)
    p.finish()

def convert_docbook_xhtml_old(infile, outfile):
    print "Converting text..."
    tree = parse(infile)
    root = tree.getroot()
    eqlist = []

    # convert nonrecursive things first:

    for e in root.getiterator("m"):
        e.tag = "div"
        e.attrib["class"] = "inline-equation"
        imgfile = "figures/eq%d.png" % (len(eqlist) + 1)
        eqlist.append((e.text, imgfile))
        img = Element("img")
        img.attrib["src"] = imgfile
        img.attrib["alt"] = e.text
        e.text = ""
        e.append(img)

    for e in root.getiterator("eq"):
        e.tag = "div"
        if "id" in e.attrib:
            e.attrib["class"] = "equation"
        else:
            e.attrib["class"] = "informal-equation"
        imgfile = "figures/eq%d.png" % (len(eqlist) + 1)
        eqlist.append((e.text, imgfile))
        img = Element("img")
        img.attrib["src"] = imgfile
        img.attrib["alt"] = e.text
        e.text = ""
        e.append(img)

    # convert everything else:
    def convert(x, level, eqlist):
        if x.tag == "article":
            title = x.find("title")
            html = Element("html")
            head = SubElement(html, "head")
            head.append(title)
            style = SubElement(head, "style")
            style.attrib["type"] = "text/css"
            style.text = style_string
            body = Element("body")
            for y in x:
                body.extend(convert(y, level, eqlist))
            html.append(body) 
            html.tail = "\n"
            return html
        elif x.tag == "section":
            title = x.find("title")
            h = Element("h%d" % level)
            h.text = title.text
            h.tail = title.tail
            x.remove(title)
            r = [h]
            for y in x: r.extend(convert(y, level + 1, eqlist))
            return r
        elif x.tag == "p":
            p = Element("p")
            p.text = x.text
            p.tail = x.tail
            for y in x:
                p.extend(convert(y, level + 1, eqlist))
            return [p]
        elif x.tag in ["a", "div", "alt", "graphic", "img"]:
            return [x]
        raise Exception("unknown tag %r" % x)

    ElementTree(convert(root, 1, eqlist)).write(outfile)

    p = progressbar("Creating images", maxval = len(eqlist))
    i = 0
    for eq, eq_file in eqlist:
        create_image(eq_file, eq, False)
        i += 1
        p.update(i)
    p.finish()

def convert_tex_xml(infile, outfile):
    """
    Converts some easy things from tex to our simplified DocBook.

    <       ->  &lt;
    &       ->  &amp;
    $x^2$   ->  <m>x^2</m>
    $$x^2$$ ->  <e>x^2</e>
    empty line -> </p><p>
    
    """
    s = "".join(open(infile).readlines())
    s = s.replace("<", "&lt;")
    s = s.replace("&", "&amp;")
    while s.find("$$") != -1:
        s = s.replace("$$", "<e>", 1) 
        s = s.replace("$$", "</e>", 1) 
    while s.find("$") != -1:
        s = s.replace("$", "<m>", 1) 
        s = s.replace("$", "</m>", 1) 
    s = s.replace("\n\n", "\n</p>\n<p>\n")
    open(outfile, "w").write(s)

def main():
    usage = """%prog [options] filein fileout

Converts "filename" DocBook to other formats.

    Examples:

    Convert to TeX:

    $ %prog m.xml m.tex

    Convert to LaTeX:

    $ %prog --latex m.xml m.tex

    Convert to xhtml:

    $ %prog m.xml m.html

    Converts our simpler format to DocBook:

    $ %prog m.xml master.xml
    """

    parser = OptionParser(usage = usage)
    parser.add_option( "-l", "--latex",
                       action = "store_true", dest = "latex",
                       default = False, help = "Converts to LaTeX" )
    options, args = parser.parse_args()
    if len(args) == 2:
        name, ext = os.path.splitext(args[1])
        if ext == ".xml":
            name, ext = os.path.splitext(args[0])
            if ext == ".xml":
                convert_xml_docbook(args[0], args[1])
                return
            elif ext == ".tex":
                convert_tex_xml(args[0], args[1])
                return
        if ext == ".tex":
            convert_docbook_latex(args[0], args[1])
            return
        if ext == ".html":
            convert_docbook_xhtml(args[0], args[1])
            return
    parser.print_help()

if __name__ == "__main__":
    main()
