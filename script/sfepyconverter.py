from convert_doc import LaTeXConverter

class SfePyDocConverter(LaTeXConverter):

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
        r += r"\renewcommand{\arraystretch}{1.3}" + '\n'
        r += "\\begin{tabular}{%s}\n" % cols
        for row in body:
            assert row.tag == "row"
            r += "\hline\n"
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
        r += "\hline\n"
        r += "\end{tabular}\n"
        r += r"\renewcommand{\arraystretch}{1.0}" + '\n'
        return r
