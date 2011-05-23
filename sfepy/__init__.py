import os, glob

from version import __version__, in_source_tree, top_dir

if in_source_tree:
    data_dir = top_dir
else:
    op = os.path
    data_dir = op.normpath(op.join(top_dir, '../../../../share/sfepy/'))

base_dir = os.path.dirname(__file__)

def get_paths(pattern):
    """
    Get files/paths matching the given pattern in the sfepy source tree.
    """
    if not in_source_tree:
        pattern = '../' + pattern

    files = glob.glob(os.path.normpath(os.path.join(top_dir, pattern)))
    return files
