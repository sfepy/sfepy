import re, os
from glob import glob

import terms
import extmods
from terms import Terms, Term, CharacteristicFunction, vector_chunk_generator
from cache import DataCache, DataCaches

def load_classes( filenames, is_class ):
    table = {}
    for filename in filenames:
        name = os.path.splitext( filename )[0]
#        print filename, name
        mod = __import__( name )
#        print mod
        for key, var in mod.__dict__.iteritems():
            if is_class( key ):
                table[var.name] = var
    return table

term_files = glob( os.path.join( 'sfepy', 'terms', 'terms*.py' ) )
is_term = re.compile( '[a-zA-Z_0-9]+Term$' ).match
term_table = load_classes( term_files, is_term )

cache_files = glob( os.path.join( 'sfepy', 'terms', 'caches*.py' ) )
is_cache = re.compile( '[a-zA-Z_0-9]+DataCache$' ).match
cache_table = load_classes( cache_files, is_cache )
