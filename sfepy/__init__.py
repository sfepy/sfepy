from __future__ import absolute_import
import os, glob

from .config import Config
from .version import __version__, in_source_tree, top_dir

data_dir = os.path.realpath(top_dir)
base_dir = os.path.dirname(os.path.normpath(os.path.realpath(__file__)))

def get_paths(pattern):
    """
    Get files/paths matching the given pattern in the sfepy source tree.
    """
    if not in_source_tree:
        pattern = '../' + pattern

    files = glob.glob(os.path.normpath(os.path.join(top_dir, pattern)))
    return files
