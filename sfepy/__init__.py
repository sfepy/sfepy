from sfepy.base.base import *
from version import version as __version__
from version import in_source_tree, top_dir

if in_source_tree:
    data_dir = '.'
else:
    import os.path as op
    data_dir = op.normpath(op.join(top_dir, '../../../../share/sfepy/'))

base_dir = os.path.dirname(__file__)
