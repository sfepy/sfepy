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

