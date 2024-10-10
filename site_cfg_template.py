##
# Template file for site configuration - copy it to site_cfg.py:
# $ cp site_cfg_template.py site_cfg.py

# Set python version (e.g. '2.7', '3.7'}), or leave '3.*' for autodetection.
python_version = '3.*'

# Operating system - one of 'posix', 'windows' or None for automatic
# determination.
system = None

# Extra flags added to the flags supplied by scikit-build to compile C
# extension modules.
compile_flags = ['-g', '-O2']

# Can be '' or one or several from 'DEBUG_FMF', 'DEBUG_MESH'. For
# developers internal use only.
debug_flags = []

# Sphinx documentation uses numpydoc extension. Set the path here in case it is
# not installed in a standard location.
numpydoc_path = None

# True for a release, False otherwise. If False, current git commit hash
# is appended to version string, if the sources are in a repository.
is_release = False

# Tetgen executable path.
tetgen_path = '/usr/bin/tetgen'

# Reference mapping memory allocation safety factor (float), requires psutil
# installed. Use None for skipping the memory check.
refmap_memory_factor = None

# If True, the large deformation terms print cells with a negative deformation
# gradient determinant, save the negative volume indicator to
# warped_cells.vtk in the working directory and raise RuntimeError.
debug_warped_cells = False
