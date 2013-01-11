"""
Create pyparsing grammar for problem configuration and options.
"""
from pyparsing import (Word, Group, Suppress, Combine, Optional,
                       Forward, Empty, quotedString, oneOf, removeQuotes,
                       delimitedList, nums, alphas, alphas8bit, alphanums,
                       Keyword)


(lparen, rparen, lbrack, rbrack,
lbrace, rbrace, colon) = map(Suppress,"()[]{}:")

integer = Combine(Optional(oneOf("+ -")) + Word(nums)).setName("integer")
cvt_int = lambda toks: int(toks[0])
integer.setParseAction(cvt_int)

boolean = Keyword("False", caseless = True)
cvt_bool =  lambda toks: toks[0].lower == 'true'
boolean.setParseAction(cvt_bool)

none = Keyword("None", caseless = True)
cvt_none =  lambda toks: [None]
none.setParseAction(cvt_none)

real = Combine(Optional(oneOf("+ -"))+ Word(nums)
	   + "." + Optional(Word(nums))
	   + Optional("e" + Optional(oneOf("+ -"))
		      + Word(nums))).setName("real")
cvt_real = lambda toks: float(toks[0])
real.setParseAction(cvt_real)

array_index = Word(nums) + Optional(colon + integer + Optional(colon + integer))
cvt_array_index = lambda toks: int(toks[0]) if len(toks)==1 else slice(*toks)
array_index.setParseAction(cvt_array_index)

word_free = Word(alphas8bit + "_-/.+**" + alphanums)
word_strict = Word(alphas8bit + alphas, alphas8bit + alphanums + "_" )


def create_bnf(allow_tuple=False, free_word=False):
    word = word_free if free_word else word_strict
    tuple, list, dict, 
    defs = standard_structs(word)

    if allow_tuple:
       return defs['dict'].inner | defs['tuple'].inner
    else:
       return defs['dict'].inner

def list_of(element, *elements):
    for e in elements: 
        element|=e 
    lst = delimitedList(element)
    return lst + Optional(Suppress(","))

def standard_structs(word = None):
    if word is None:
       word = strict_word

    tuple_str = Forward()
    list_str = Forward()
    dict_str = Forward()
    cvt_tuple = lambda toks : tuple(toks.asList())
    cvt_dict = lambda toks: dict(toks.asList())

    list_item = (none | boolean | real | integer | list_str | tuple_str
                 | dict_str
                 | quotedString.setParseAction(removeQuotes)
                 | word )
    list_item2 = list_item | Empty().setParseAction(lambda: [None])
   

    tuple_str.inner = list_of(list_item)
    tuple_str.inner.setParseAction(cvt_tuple)
    tuple_str << (Suppress("(") + tuple_str.inner  + Suppress(")"))

    list_str.inner = tuple_str.inner.copy()
    list_str.inner.setParseAction(lambda toks: list(toks))
    list_str << (lbrack + list_str.inner + rbrack)

    dict_entry = Group(list_item + colon + list_item2)
    dict_str.inner = list_of(dict_entry) 
    dict_str.inner.setParseAction(cvt_dict)
    dict_str << (lbrace + Optional(dict_str.inner) + rbrace)

    return {'tuple' : tuple_str, 'list' : list_str, 'dict' : dict_str, 'lsit_item' : list_item }
