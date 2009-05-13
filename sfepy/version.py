import os.path as op

# If installed, up_dir is '.', otherwise (in (git) source directory) '..'.
for up_dir in ['..', '.']:
    top_dir = op.normpath(op.join(op.dirname(__file__), up_dir))
    version_file = op.join(top_dir, 'VERSION')
    if op.isfile(version_file):
        break
else:
    raise RuntimeError('cannot find VERSION file!')

fd = open(version_file, 'r')
version = fd.readline().strip()
fd.close()

master = op.join(top_dir, '.git/refs/heads/master')
if op.isfile(master):
    fd = open(master, 'r')
    version += ' (%s)' % fd.readline().strip()
    fd.close()
else:
    version += ''

in_source_tree = up_dir == '..'

del fd, up_dir
