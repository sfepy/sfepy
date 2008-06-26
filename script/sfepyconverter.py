from convert_doc import LaTeXConverter

class SfePyDocConverter(LaTeXConverter):

    def convert_article(self, node):
        assert node.tag == "article"
        self.check_zero_tail(node.tail)
        r = r"""\documentclass[10pt]{article}
\usepackage{bm}
\usepackage{a4wide}
\usepackage{longtable}
\usepackage{colortbl}
\usepackage{amsmath}
\usepackage[utf8]{inputenc}
\usepackage{hyperref}
\def\mybackslash{$\backslash$}
"""
        for x in node:
            r += self.convert_node(x)
        r += """\n\\vfill\n\end{document}\n"""

        return r

    def convert_command(self, node):
        assert node.tag in ["command", "citetitle"]
        r = r"{\small\verb|"
        if node.text is not None:
            r += node.text
        r += "|}"
        r += self.default_label(node)
        if node.tail is not None:
            r += self.escape(node.tail)
        return r

    def convert_table(self, node):
        assert node.tag == "table"
        self.check_zero_tail(node.tail)
        r = "\n"
        body = node.find("tgroup").find("tbody")
        cols = "|"
        ncolumns = len(body.find("row"))
        for i in range(ncolumns-1):
            cols += " l | "
        cols += "p{5cm} |"
        r += r"\renewcommand{\arraystretch}{1.4}" + '\n'
        r += "\\begin{longtable}{%s}\n" % cols
        r += r"""
  \hline
  \hline
  \endfirsthead
  \hline
  \multicolumn{%d}{|c|}{\dots{}\emph{continued}} \\
  \hline
  \endhead
  \hline
  \multicolumn{%d}{|c|}{\emph{continued}\dots} \\
  \hline
  \endfoot
  \hline
  \hline
  \endlastfoot
""" % (ncolumns, ncolumns)
        for ii, row in enumerate( body ):
            assert row.tag == "row"
#            r += "\hline\n"
            r += r"\rowcolor[gray]{%2.1f}[.9\tabcolsep]" % {0 : 1.0,
                                                            1: 0.9}[ii % 2]
            r += "\n"
            for entry in row:
                assert row.tag == "row"
                self.check_zero_tail(row.tail)
                if entry.text:
                    r += self.escape(entry.text)
                for x in entry: 
                    r += self.convert_node(x)
                r += " & "
            if r[-2] == "&":
                r = r[:-3]
            r += "\\\\ \n"
#        r += "\hline\n"
        r += "\end{longtable}\n"
        r += r"\renewcommand{\arraystretch}{1.0}" + '\n'
        return r
