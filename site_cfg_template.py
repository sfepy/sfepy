##
# Template file for site configuration - copy it to site_cfg.py:
# $ cp site_cfg_template.py site_cfg.py

# Tetgen executable path.
tetgen_path = '/home/share/software/bin/tetgen'

# 'lib' or 'lib64' depending on your architecture (32bit or 64bit)
archlib = 'lib'

# Force python version (e.g. '2.4'), or left 'auto' for autodetection.
python_version = 'auto'

# Installation prefix. Ignore if you do not install numpy/scipy on local account.
# Numpy headers should be in
# <numpy_prefix>/usr/<archlib>/python<python_version>/site-packages/numpy/core/include
numpy_prefix = '/home/share/software'
