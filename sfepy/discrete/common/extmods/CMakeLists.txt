set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(common_src fmfield.c refmaps.c geommech.c common_python.c)
add_library(sfepy_common SHARED ${common_src})
python_extension_module(sfepy_common)
target_include_directories(sfepy_common PRIVATE ${NumPy_INCLUDE_DIRS})

add_cython_target(_fmfield _fmfield.pyx)
python_extension_module(_fmfield)
add_library(_fmfield MODULE ${_fmfield})
target_include_directories(_fmfield PRIVATE ${NumPy_INCLUDE_DIRS})
target_link_libraries(_fmfield sfepy_common)

install(TARGETS _fmfield LIBRARY DESTINATION _fmfield)
