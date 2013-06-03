"""
Grammar for selecting regions of a domain.

Regions serve for selection of certain parts of the computational domain (=
selection of nodes and elements of a FE mesh). They are used to define the
boundary conditions, the domains of terms and materials etc.

Notes
-----
History: pre-git versions already from from 13.06.2006.
"""
from pyparsing import Literal, CaselessLiteral, Word, delimitedList,\
     Group, Optional, ZeroOrMore, nums, alphas, alphanums,\
     Combine, StringStart, StringEnd, Forward, oneOf, ParseException

ParseException # Needed for importing elsewhere.

op_codes = ['OA_SubV', 'OA_SubE', 'OA_SubF', 'OA_SubC', 'OA_SubS',
            'OA_AddV', 'OA_AddE', 'OA_AddF', 'OA_AddC', 'OA_AddS',
            'OA_IntersectV', 'OA_IntersectE', 'OA_IntersectF',
            'OA_IntersectC', 'OA_IntersectS']
eval_codes = ['E_NIR', 'E_NOS', 'E_NBF', 'E_NOG', 'E_ONIR', 'E_NI', 'E_NOSET',
              'E_EBF', 'E_EOG', 'E_EI1', 'E_EI2', 'E_EOSET']
kw_codes = ['KW_All', 'KW_Region']

def to_stack(stack):
    def push_first(str, loc, toks):
        if toks:
            stack.append(toks[0])
        return toks
    return push_first

def replace(what, keep=False):
    def _replace(str, loc, toks):
        ret = {'token' : what, 'orig' : []}
        if keep:
            ret['orig'] = list(toks[0])
        return ret
    return _replace

def replace_with_region(what, r_index):
    def _replace(str, loc, toks):
        ret = {'token' : what, 'orig' : []}

        orig = toks[0]
        r_orig = orig[r_index]
        if isinstance(r_orig, dict) and (r_orig['token'] == 'KW_Region'):
            orig = list(orig[:r_index]) + r_orig['orig']
        ret['orig'] = orig
        return ret
    return _replace

def join_tokens(str, loc, toks):
    return [" ".join(toks[0])]

def visit_stack(stack, op_visitor, leaf_visitor):

    def visit(stack, level):
        op = stack.pop()

        token = op['token']
        if token in op_codes:
            res2 = visit(stack, level + 1)
            res1 = visit(stack, level + 1)
            return op_visitor(level, op, res1, res2)

        elif token in eval_codes:
            return leaf_visitor(level, op)

        elif token in kw_codes:
            return leaf_visitor(level, op)

        else:
            raise ValueError, token

    return visit(stack, 0)

def print_op(level, op, item1, item2):
    print level * '  ' + (': %s' % op)

def print_leaf(level, op):
    print level * '  ' + ('< %s' % op)

def print_stack(stack):
    visit_stack(stack, print_op, print_leaf)

def create_bnf(stack):
    point = Literal(".")
    comma = Literal(",")
    e = CaselessLiteral("E")
    inumber = Word(nums)
    fnumber = Combine(Word("+-"+nums, nums) +
                       Optional(point + Optional(Word(nums))) +
                       Optional(e + Word("+-"+nums, nums)))
    _of = Literal('of')
    _in = Literal('in')
    _by = Literal('by')
    _copy = Literal('copy')

    _mv = Literal('-v').setParseAction(replace('OA_SubV'))
    _me = Literal('-e').setParseAction(replace('OA_SubE'))
    _mf = Literal('-f').setParseAction(replace('OA_SubF'))
    _mc = Literal('-c').setParseAction(replace('OA_SubC'))
    _ms = Literal('-s').setParseAction(replace('OA_SubS'))
    _pv = Literal('+v').setParseAction(replace('OA_AddV'))
    _pe = Literal('+e').setParseAction(replace('OA_AddE'))
    _pf = Literal('+f').setParseAction(replace('OA_AddF'))
    _pc = Literal('+c').setParseAction(replace('OA_AddC'))
    _ps = Literal('+s').setParseAction(replace('OA_AddS'))
    _inv = Literal('*v').setParseAction(replace('OA_IntersectV'))
    _ine = Literal('*e').setParseAction(replace('OA_IntersectE'))
    _inf = Literal('*f').setParseAction(replace('OA_IntersectF'))
    _inc = Literal('*c').setParseAction(replace('OA_IntersectC'))
    _ins = Literal('*s').setParseAction(replace('OA_IntersectS'))
    regop = (_mv | _me | _mf | _mc | _ms |
             _pv | _pe | _pf | _pc | _ps |
             _inv | _ine | _inf | _inc | _ins)

    lpar  = Literal("(").suppress()
    rpar  = Literal(")").suppress()

    _all = Literal('all').setParseAction(replace('KW_All'))
    node = Literal('node')
    nodes = Literal('nodes')
    element = Literal('element')
    elements = Literal('elements')
    group = Literal('group')
    _set = Literal('set')
    surface = Literal('surface')

    ident = Word(alphas + '_.', alphanums + '_.')
    set_name = Word(nums) | ident

    function = Word(alphas + '_', alphanums + '_')
    function = Group(function).setParseAction(join_tokens)

    region = Combine(Literal('r.') + Word(alphas + '_',
                                          '_' + alphas + nums + '.'))
    region = Group(Optional(_copy, default='nocopy') + region)
    region.setParseAction(replace('KW_Region', keep=True))

    coor = oneOf('x y z')
    boolop = oneOf('& |')
    relop = oneOf('< > <= >= != ==')
    bool_term = (ZeroOrMore('(') + (coor | fnumber) + relop + (coor | fnumber)
                 + ZeroOrMore(')'))
    relation = Forward()
    relation << (ZeroOrMore('(')
                 + bool_term + ZeroOrMore(boolop + relation)
                 + ZeroOrMore(')'))
    relation = Group(relation).setParseAction(join_tokens)

    nos = Group(nodes + _of + surface).setParseAction(replace('E_NOS'))
    nir = Group(nodes + _in + relation).setParseAction(
        replace('E_NIR', keep=True))
    nbf = Group(nodes + _by + function).setParseAction(
        replace('E_NBF', keep=True))
    ebf = Group(elements + _by + function).setParseAction(
        replace('E_EBF', keep=True))
    eog = Group(elements + _of + group + Word(nums)).setParseAction(
        replace('E_EOG', keep=True))
    nog = Group(nodes + _of + group + Word(nums)).setParseAction(
        replace('E_NOG', keep=True))
    onir = Group(node + _in + region).setParseAction(
        replace_with_region('E_ONIR', 2))
    ni = Group(node + delimitedList(inumber)).setParseAction(
        replace('E_NI', keep=True))
    ei1 = Group(element + delimitedList(inumber)).setParseAction(
        replace('E_EI1', keep=True))
    etuple = (lpar.suppress() + inumber + comma.suppress()
              + inumber + rpar.suppress())
    ei2 = Group(element + delimitedList(etuple)).setParseAction(
        replace('E_EI2', keep=True))
    noset = Group(nodes + _of + _set + set_name).setParseAction(
        replace('E_NOSET', keep=True))
    eoset = Group(elements + _of + _set + set_name).setParseAction(
        replace('E_EOSET', keep=True))

    region_expression = Forward()

    atom1 = (_all | region | ni | onir | nos | nir | nbf
             | ei1 | ei2 | ebf | eog | nog | noset | eoset)
    atom1.setParseAction(to_stack(stack))
    atom2 = (lpar + region_expression.suppress() + rpar)
    atom = (atom1 | atom2)

    aux = (regop + region_expression)
    aux.setParseAction(to_stack(stack))
    region_expression << atom + ZeroOrMore(aux)
    region_expression = StringStart() + region_expression + StringEnd()

    return region_expression

_test_strs = ['nodes of surface -n r.egion_1',
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
             'nodes of surface +n nodes by pokus',
             'elements of group 6 +e nodes by fn2_3c',
             """r.egion_1 *n (r.egion_2 +e (nodes in (y > 0) *n r.egion_32))
             -n nodes of surface -e r.egion_5""",
             'node in r.region_3',
             'node 10',
             'elements by afun']

if __name__ == "__main__":
    test_strs = _test_strs

    stack = []
    bnf = create_bnf(stack)

    n_fail = 0
    for test_str in test_strs:
        print test_str
        stack[:] = []

        try:
            out = bnf.parseString(test_str)
        except:
            print '...failed!'
            n_fail += 1
            continue

        print_stack(stack)
    print 'failed: %d' % n_fail
