from __future__ import absolute_import
import sfepy
from . import terms
from . import extmods
from .terms import Terms, Term
from .terms_th import THTerm, ETHTerm
from sfepy.base.base import load_classes

term_files = sfepy.get_paths('sfepy/terms/terms*.py')
term_table = load_classes(term_files, [Term], ignore_errors=True)

del sfepy

def register_term(cls):
    """
    Register a custom term.
    """
    term_table[cls.name] = cls
