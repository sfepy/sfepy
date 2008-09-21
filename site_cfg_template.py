##
# Template file for site configuration - copy it to site_cfg.py:
# $ cp site_cfg_template.py site_cfg.py

# Tetgen executable path.
tetgen_path = '/usr/bin/tetgen'

# 'lib' or 'lib64' depending on your architecture (32bit or 64bit)
archlib = 'lib'

# Force python version (e.g. '2.4'), or left 'auto' for autodetection.
python_version = 'auto'

# Installation prefix. Ignore if you do not install numpy/scipy on local account.
# Numpy headers should be in
# <numpy_prefix>/usr/<archlib>/python<python_version>/site-packages/numpy/core/include
numpy_include = None

# Flags passed to C compiler to compile C extension modules.
# Below is a Linux gcc default.
# For an Intel Mac, use
# '-g -O2 -fPIC -DPIC -fno-strict-aliasing -fno-common -dynamic'
opt_flags = '-g -O2 -fPIC -DPIC'

# Flags passed to C compiler to link C extension modules.
# Below is a Linux gcc default.
# For an Intel Mac, use '-dynamiclib -undefined dynamic_lookup -fPIC -DPIC'
link_flags = '-shared -fPIC -DPIC'

# Can be a combination of '-DDEBUG_FMF' '-DDEBUG_MESH'.
debug_flags = ''
