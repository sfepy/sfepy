set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(common_path ../../discrete/common/extmods)

file(GLOB terms_src "*.c")
add_library(sfepy_terms SHARED ${terms_src})
python_extension_module(sfepy_terms)
target_link_libraries(sfepy_terms sfepy_common)
target_include_directories(sfepy_terms PRIVATE ${NumPy_INCLUDE_DIRS} ${common_path})

#add_cython_target(_fmfield _fmfield.pyx)
#python_extension_module(_fmfield)
#add_library(_fmfield MODULE ${_fmfield})
#target_include_directories(_fmfield PRIVATE ${NumPy_INCLUDE_DIRS})
#target_link_libraries(_fmfield sfepy_common)

install(TARGETS sfepy_terms LIBRARY DESTINATION sfepy)
