import os.path as op

top_dir = op.normpath(op.join(op.dirname(__file__), '..'))

fd = open(op.join(top_dir, 'VERSION'), 'r')
version = fd.readline().strip()
fd.close()

master = op.join(top_dir, '.git/refs/heads/master')
if op.isfile(master):
    fd = open(master, 'r')
    version += ' (%s)' % fd.readline().strip()
    fd.close()
else:
    version = ''

del fd
