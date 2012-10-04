"""
Create pyparsing grammar for problem configuration and options.
"""
from pyparsing import (Word, Group, Suppress, Combine, Optional,
                       Forward, Empty, quotedString, oneOf, removeQuotes,
                       delimitedList, nums, alphas, alphas8bit, alphanums,
                       Keyword)

def create_bnf(allow_tuple=False, free_word=False):
    cvt_int = lambda toks: int(toks[0])
    cvt_real = lambda toks: float(toks[0])
    cvt_bool =  lambda toks: toks[0].lower == 'true'
    cvt_none =  lambda toks: [None]
    cvt_tuple = lambda toks : tuple(toks.asList())
    cvt_dict = lambda toks: dict(toks.asList())


    # define punctuation as suppressed literals
    (lparen, rparen, lbrack, rbrack,
     lbrace, rbrace, colon) = map(Suppress,"()[]{}:")

    integer = Combine(Optional(oneOf("+ -")) + Word(nums)).setName("integer")
    integer.setParseAction(cvt_int)

    boolean = Keyword("False", caseless = True)
    boolean.setParseAction(cvt_bool)

    none = Keyword("None", caseless = True)
    none.setParseAction(cvt_none)

    real = Combine(Optional(oneOf("+ -"))+ Word(nums)
                   + "." + Optional(Word(nums))
                   + Optional("e" + Optional(oneOf("+ -"))
                              + Word(nums))).setName("real")
    real.setParseAction(cvt_real)

    tuple_str = Forward()
    list_str = Forward()
    dict_str = Forward()

    if free_word:
        string = Word(alphas8bit + "_-/.+**" + alphanums)

    else:
        string = Word(alphas8bit + alphas, alphas8bit + alphanums + "_" )

    list_item = (none | boolean | real | integer | list_str | tuple_str
                 | dict_str
                 | quotedString.setParseAction(removeQuotes)
                 | string )
    list_item2 = list_item | Empty().setParseAction(lambda: [None])

    tuple_inner = Optional(delimitedList(list_item)) + Optional(Suppress(","))
    tuple_inner.setParseAction(cvt_tuple)
    tuple_str << (Suppress("(") + tuple_inner  + Suppress(")"))

    list_inner = Optional(delimitedList(list_item) + Optional(Suppress(",")))
    list_inner.setParseAction(lambda toks: list(toks))
    list_str << (lbrack + list_inner + rbrack)

    dict_entry = Group(list_item + colon + list_item2)
    dict_inner = delimitedList(dict_entry) + Optional(Suppress(","))
    dict_inner.setParseAction(cvt_dict)
    dict_str << (lbrace + Optional(dict_inner) + rbrace)

    dict_or_tuple =  dict_inner | tuple_inner

    if allow_tuple:
        return dict_or_tuple

    else:
        return dict_inner
