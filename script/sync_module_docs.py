#!/usr/bin/env python
"""
Synchronize the documentation files in a given directory ``doc_dir`` with the
actual state of the SfePy sources in ``top_dir``. Missing files are created,
files with no corresponding source file are removed, other files are left
untouched.

Notes
-----
The developer guide needs to be edited manually to reflect the changes.
"""
import os
import fnmatch
from optparse import OptionParser

from sfepy.base.base import output
from sfepy.base.ioutils import locate_files, edit_filename, ensure_path

usage = '%prog [options] doc_dir top_dir\n' + __doc__.rstrip()

omits = [
    '__init__.py',
    '__config__.py',
    'debug.py',
    'friction_slip.py',
    'gitwash_dumper.py',
    'setup.py',
    'site_cfg.py',
    'site_cfg_template.py',
]

doc_template = """%s
%s

.. automodule:: %s
   :members:
   :undoc-members:
"""

help = {
    'dry_run' :
    'only show what changes would be made',
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-n', '--dry-run',
                      action='store_true', dest='dry_run',
                      default=False, help=help['dry_run'])
    options, args = parser.parse_args()

    if len(args) == 2:
        doc_dir, top_dir = [os.path.realpath(ii) for ii in args]
    else:
        parser.print_help(),
        return

    docs = set(ii for ii in locate_files('*.rst', root_dir=doc_dir))

    sources = set(ii for ii in
                  locate_files('*.py',
                               root_dir=os.path.join(top_dir, 'sfepy'))
                  if os.path.basename(ii) not in omits)
    sources.update(ii for ii in
                   locate_files('*.pyx',
                                root_dir=os.path.join(top_dir, 'sfepy')))
    scripts = set(ii for ii in
                  locate_files('*.py',
                               root_dir=os.path.join(top_dir, 'script'))
                  if os.path.basename(ii) not in omits)
    top_scripts = set(os.path.realpath(ii)
                      for ii in fnmatch.filter(os.listdir(top_dir), '*.py')
                      if os.path.basename(ii) not in omits)

    all_sources = set()
    all_sources.update(sources, scripts, top_scripts)

    cwd = os.path.realpath(os.path.curdir) + os.path.sep

    output.prefix = 'smd:'
    output('removing unneeded rst files in "%s"...' % doc_dir)
    for doc in sorted(docs):
        aux = edit_filename(doc, new_ext='.py')
        src1 = os.path.normpath(aux.replace(doc_dir, top_dir))

        aux = edit_filename(doc, new_ext='.pyx')
        src2 = os.path.normpath(aux.replace(doc_dir, top_dir))

        if (src1 not in all_sources) and (src2 not in all_sources):
            output('remove: %s' % doc.replace(cwd, ''))
            if not options.dry_run:
                os.remove(doc)
    output('...done')

    output('creating missing rst files in "%s"...' % doc_dir)
    for src in sorted(all_sources):
        aux = edit_filename(src, new_ext='.rst')
        doc = os.path.normpath(aux.replace(top_dir, doc_dir))

        if doc not in docs:
            output('create: %s' % doc.replace(cwd, ''))
            if not options.dry_run:
                mod_filename = src.replace(top_dir + os.path.sep, '')
                mod_name = mod_filename.replace(os.path.sep, '.')
                mod_name = edit_filename(mod_name, new_ext='')
                if mod_name.startswith('sfepy'): # Module.
                    title = mod_name + ' module'

                else: # Script.
                    title = mod_filename + ' script'
                    mod_name = mod_name.split('.')[-1]

                underlines = '=' * len(title)

                contents = doc_template % (title, underlines, mod_name)

                ensure_path(doc)
                fd = open(doc, 'w')
                fd.write(contents)
                fd.close()

    output('...done')

if __name__ == '__main__':
    main()
