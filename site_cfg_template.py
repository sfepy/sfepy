##
# Template file for site configuration - copy it to site_cfg.py:
# $ cp site_cfg_template.py site_cfg.py

# Force python version (e.g. '2.4'), or left 'auto' for autodetection.
python_version = 'auto'

# Operating system - one of 'posix', 'windows' or None for automatic
# determination.
system = None

# Extra flags added to the flags supplied by distutils to compile C
# extension modules.
compile_flags = '-g -O2'

# Extra flags added to the flags supplied by distutils to link C
# extension modules.
link_flags = ''

# Can be '' or one or several from '-DDEBUG_FMF', '-DDEBUG_MESH'. For
# developers internal use only.
debug_flags = ''

# Sphinx documentation uses numpydoc extension. Set the path here in case it is
# not installed in a standard location.
numpydoc_path = None

# True for a release, False otherwise. If False, current git commit hash
# is appended to version string, if the sources are in a repository.
is_release = False

# Tetgen executable path.
tetgen_path = '/usr/bin/tetgen'
