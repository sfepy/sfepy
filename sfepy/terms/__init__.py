import os
from glob import glob

import sfepy
import terms
import extmods
from terms import Terms, Term, CharacteristicFunction, vector_chunk_generator
from cache import DataCache, DataCaches

def load_classes(filenames, cls):
    from sfepy.base.base import import_file, find_subclasses
    table = {}
    for filename in filenames:
        mod = import_file(filename)
        table.update(find_subclasses(vars(mod), [cls], omit_unnamed=True))

    return table

def get_paths(pattern):
    if not sfepy.in_source_tree:
        pattern = '../' + pattern

    files = glob(os.path.normpath(os.path.join(sfepy.top_dir, pattern)))
    return files

term_files = get_paths('sfepy/terms/terms*.py')
term_table = load_classes(term_files, Term)

cache_files = get_paths('sfepy/terms/caches*.py')
cache_table = load_classes(cache_files, DataCache)
