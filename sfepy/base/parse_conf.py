"""
Create pyparsing grammar for problem configuration and options.
"""
from pyparsing import (Word, Group, Suppress, Combine, Optional,
                       Forward, Empty, quotedString, oneOf, removeQuotes,
                       delimitedList, nums, alphas, alphas8bit)

def create_bnf():
    cvt_int = lambda toks: int(toks[0])
    cvt_real = lambda toks: float(toks[0])
    cvt_tuple = lambda toks : tuple(toks.asList())
    cvt_dict = lambda toks: dict(toks.asList())

    # define punctuation as suppressed literals
    (lparen, rparen, lbrack, rbrack,
     lbrace, rbrace, colon) = map(Suppress,"()[]{}:")

    integer = Combine(Optional(oneOf("+ -")) + Word(nums)).setName("integer")
    integer.setParseAction(cvt_int)

    real = Combine(Optional(oneOf("+ -"))
                   + Word(nums) + "." + Optional(Word(nums))).setName("real")
    real.setParseAction(cvt_real)

    tuple_str = Forward()
    list_str = Forward()
    dict_str = Forward()

    list_item = (real | integer | Group(list_str) | tuple_str | dict_str
                 | quotedString.setParseAction(removeQuotes)
                 | Word(alphas8bit + alphas + "_"))
    list_item2 = list_item | Empty().setParseAction(lambda: [None])

    tuple_str << (Suppress("(") + Optional(delimitedList(list_item)) +
                  Optional(Suppress(",")) + Suppress(")"))
    tuple_str.setParseAction(cvt_tuple)

    list_str << (lbrack + Optional(delimitedList(list_item) +
                                   Optional(Suppress(","))) + rbrack)

    dict_entry = Group(list_item + colon + list_item2)
    dict_inner = delimitedList(dict_entry) + Optional(Suppress(","))
    dict_inner.setParseAction(cvt_dict)

    dict_str << (lbrace + Optional(dict_inner) + rbrace)

    return dict_inner
