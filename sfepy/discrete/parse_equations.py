from pyparsing import Combine, Literal, Word, delimitedList, Group, Optional,\
     ZeroOrMore, OneOrMore, nums, alphas, alphanums,\
     StringStart, StringEnd, CaselessLiteral, Forward, oneOf

class TermParse:
    def __str__(self):
        ss = "%s\n" % self.__class__
        for key, val in self.__dict__.items():
            ss += "  %s:\n    %s\n" % (key, self.__dict__[key])
        return ss

def collect_term(term_descs, lc):
    signs = {'+': 1.0, '-': -1.0}

    def append(str, loc, toks):
        tp = TermParse()
        if len(toks.mul):
            tp.sign = eval(''.join(toks.mul)) * signs[lc[0]]
        else:
            tp.sign = signs[toks.sign] * signs[lc[0]]
        tp.integral = toks.term_desc.integral
        if not tp.integral:
            tp.integral = 'a'
        tp.region = toks.term_desc.region
        tp.flag = toks.term_desc.flag
        tp.name = toks.term_desc.name
        tp.args = ', '.join(toks.args[0])
        term_descs.append(tp)
    return append

def rhs(lc):
    def aux(str, loc, toks):
        if toks:
            lc[0] = '-'
    return aux

def create_bnf(term_descs):
    """term_descs .. list of TermParse objects
    (sign, term_name, term_arg_names), where sign can be real or complex
    multiplier"""

    lc = ['+']  # Linear combination context.
    equal = Literal("=").setParseAction(rhs(lc))
    zero = Literal("0").suppress()

    point = Literal(".")
    e = CaselessLiteral("E")
    inumber = Word("+-" + nums, nums)
    fnumber = Combine(Word("+-" + nums, nums) +
                      Optional(point + Optional(Word(nums))) +
                      Optional(e + Word("+-" + nums, nums)))
    number = fnumber + Optional(Literal('j'), default='')
    add_op = oneOf('+ -')
    number_expr = Forward()
    number_expr << Optional(add_op) + ZeroOrMore('(') + number \
                + ZeroOrMore(add_op + number_expr) \
                + ZeroOrMore(')')

    ident = Word(alphas, alphanums + "_")

    integral = Combine((Literal('i') + Word(alphanums)) | Literal('i')
                       | Literal('a') | Word(nums))("integral")

    history = Optional('[' + inumber + ']', default='')("history")

    variable = Combine(Word(alphas, alphanums + '._') + history)

    derivative = Combine(Literal('d') + variable
                         + Literal('/') + Literal('dt'))
    # Not needed here, just for correspondence with create_arg_parser().
    derivative2 = Combine(OneOrMore(Literal('d')) + variable)

    trace = Combine(Literal('tr') + '('
                    + Optional(ident + Literal(','))
                    + variable
                    + ')', adjacent=False)

    generalized_var = derivative | derivative2 | trace | variable
    args = Group(delimitedList(generalized_var))

    flag = Literal('a')

    term = ((Optional(Literal('+') | Literal('-'), default='+')("sign")
             ^ Optional(number_expr + Literal('*').suppress(),
                        default=['1.0', ''])("mul"))
            + Combine(
                ident("name") + Optional(
                    "."
                    + (integral + "." + ident("region") + "." + flag("flag")
                       | integral + "." + ident("region")
                       | ident("region"))
                ))("term_desc") + "("
            + Optional(args, default=[''])("args") + ")")
    term.setParseAction(collect_term(term_descs, lc))

    rhs1 = equal + OneOrMore(term)
    rhs2 = equal + zero
    equation = StringStart() + OneOrMore(term) \
               + Optional(rhs1 | rhs2) + StringEnd()
    ## term.setDebug()
    return equation

if __name__ == "__main__":

    test_str = """d_term1.Y(fluid, u, w, Nu, dcf, mode)
                  + 5.0 * d_term2.Omega(u, w, Nu, dcf, mode)
                  - d_another_term.Elsewhere(w, p[-1], Nu, dcf, mode)
                  = - dw_rhs.a.Y3(u, q, Nu, dcf, mode)"""

    term_descs = []
    bnf = create_bnf(term_descs)
    out = bnf.parseString(test_str)

    print('out:', out, '\n')
    for tp in term_descs:
        print(tp)
