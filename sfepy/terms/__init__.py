import re, os
from glob import glob

import sfepy
import terms
import extmods
from terms import Terms, Term, CharacteristicFunction, vector_chunk_generator
from cache import DataCache, DataCaches

def load_classes( filenames, is_class ):
    table = {}
    for filename in filenames:
        name = os.path.splitext( filename )[0]
        parts = name.split( os.path.sep )
        mod, name = '.'.join( parts[-3:] ), parts[-1:]
        mod = __import__( mod, globals(), locals(), name )
        for key, var in mod.__dict__.iteritems():
            if is_class( key ):
                table[var.name] = var
    return table

def get_paths(pattern):
    if not sfepy.in_source_tree:
        pattern = '../' + pattern

    files = glob(os.path.normpath(os.path.join(sfepy.top_dir, pattern)))
    return files

term_files = get_paths('sfepy/terms/terms*.py')
is_term = re.compile( '[a-zA-Z_0-9]+Term$' ).match
term_table = load_classes( term_files, is_term )

cache_files = get_paths('sfepy/terms/caches*.py')
is_cache = re.compile( '[a-zA-Z_0-9]+DataCache$' ).match
cache_table = load_classes( cache_files, is_cache )
