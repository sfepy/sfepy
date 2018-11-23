"""
Create pyparsing grammar for problem configuration and options.
"""
from pyparsing import (Word, Group, Suppress, Combine, Optional,
                       Forward, Empty, quotedString, oneOf, removeQuotes,
                       delimitedList, nums, alphas, alphanums,
                       Keyword, CaselessLiteral)

(lparen, rparen, lbrack, rbrack,
 lbrace, rbrace, colon, equal_sign) = map(Suppress, '()[]{}:=')

word_free_wb = Word(alphas + '@_-/.+**' + alphanums)
word_free = Forward()
word_free = word_free_wb + Optional(lbrace + word_free + rbrace + word_free)
word_free.setParseAction(lambda toks: ''.join(toks[0]))

word_strict = Word(alphas, alphas + alphanums + '_' )

integer = Combine(Optional(oneOf('+ -')) + Word(nums)).setName('integer')
cvt_int = lambda toks: int(toks[0])
integer.setParseAction(cvt_int)

boolean_true = Keyword('True', caseless=True)
boolean_true.setParseAction(lambda x: True)
boolean_false = Keyword('False', caseless=True)
boolean_false.setParseAction(lambda x: False)

boolean = boolean_true | boolean_false

none = Keyword('None', caseless=True)

cvt_none = lambda toks: [None]
none.setParseAction(cvt_none)

e = CaselessLiteral("e")
real = (Combine(Optional(oneOf('+ -')) + Word(nums) +
               '.' + Optional(Word(nums)) +
               Optional(e + Optional(oneOf('+ -')) + Word(nums))
       )| Combine(Optional(
            oneOf('+ -')) + Word(nums) +
            Optional('.') + Optional(Word(nums)) +
            e + Optional(oneOf('+ -')) + Word(nums))
       ).setName('real')
cvt_real = lambda toks: float(toks[0])
real.setParseAction(cvt_real)

cmplx = real + CaselessLiteral('j')
cvt_cmplx = lambda toks: complex(toks[0])
cmplx.setParseAction(cvt_cmplx)

array_index = integer + Optional(colon + integer +
                                 Optional(colon + integer))
cvt_array_index = lambda toks: int(toks[0]) if len(toks) == 1 \
                  else slice(*toks)
array_index.setParseAction(cvt_array_index)
array_braces = lbrack + array_index + rbrack

def create_bnf(allow_tuple=False, free_word=False):
    word = word_free if free_word else word_strict
    defs = get_standard_type_defs(word)
    empty = Empty()
    empty.setParseAction( lambda toks: [{}])
    if allow_tuple:
        return defs['dict'].inner | defs['tuple'].inner | empty
    else:
        return defs['dict'].inner | empty

def list_of(element, *elements):
    """
    Return lexical element that parses a list of items. The items can be a one
    or several lexical elements. For example, result of ``list_of(real,
    integer)`` parses list of real or integer numbers.
    """
    for e in elements:
        element ^= e
    lst = delimitedList(element)
    return lst + Optional(Suppress(','))

def get_standard_type_defs(word=word_free):
    """
    Return dict of the pyparsing base lexical elements.

    The compound types (tuple, list, dict) can contain compound types or simple
    types such as integers, floats and words.

    Parameters
    ----------
    word : lexical element
        A custom lexical element for word.

    Returns
    -------
    defs : dict
        The dictionary with the following items:

        - tuple: (..., ..., ...)
        - list: [..., ...., ...]
        - dict: {...:..., ...:..., ....} or {...=..., ...=..., ....}
        - list_item: any of preceding compound types or simple types
    """
    tuple_str = Forward()
    list_str = Forward()
    dict_str = Forward()
    cvt_tuple = lambda toks : [ tuple(toks.asList()) ]
    cvt_dict = lambda toks: [ dict(toks.asList())]

    list_item = (none ^ boolean ^ cmplx ^ real ^ integer ^
                 list_str ^ tuple_str ^ dict_str ^
                 quotedString.setParseAction(removeQuotes) ^
                 word)
    list_item2 = list_item | Empty().setParseAction(lambda: [None])

    tuple_str.inner = Empty() ^ list_of(list_item)
    list_str.inner = tuple_str.inner.copy()
    tuple_str.inner.setParseAction(cvt_tuple)
    tuple_str << (lparen + tuple_str.inner + rparen)

    list_str.inner.setParseAction(lambda toks: [toks.asList()])
    list_str << (lbrack + list_str.inner + rbrack)

    dict_entry = Group(list_item + (colon | equal_sign) + list_item2)
    dict_str.inner = Empty() ^ list_of(dict_entry)
    dict_str.inner.setParseAction(cvt_dict)
    dict_str << (lbrace + (dict_str.inner |
                           Empty().setParseAction( lambda x: [{}] ) ) + rbrace)

    defs = {'tuple' : tuple_str,
            'list' : list_str,
            'dict' : dict_str,
            'list_item' : list_item}

    return defs

def list_dict(word=word_free):
    """
    Return the pyparsing lexical element, that parses a string either as a list
    or as a dictionary.

    Parameters
    ----------
    word : lexical element
        A custom lexical element for word.

    Returns
    -------
    ld : lexical element
        The returned lexical element parses a string in the form
        ``..., ..., ...`` or ``key1:..., key2=..., key3: ...``
        where ``...`` is a ``list_item`` from :func:`get_standard_type_defs()`
        and interprets it as a list or a dictionary.
    """
    defs = get_standard_type_defs(word)
    i = defs['list_item']
    arg = i.copy()
    arg.setParseAction(lambda t: (t[0],))
    narg = word_strict + (colon | equal_sign) + i
    narg.setParseAction(lambda t: (t[0], t[1]))

    ld = Group(list_of(narg ^ arg) | Empty() )
    ld.setParseAction(lambda t:  ([x[0] for x in t[0] if len(x) == 1],
                                 dict([x for x in t[0] if len(x) > 1]))
                      )
    return ld
