
set(common_path ../../discrete/common/extmods)

file(GLOB terms_src "*.c")
add_library(sfepy_terms STATIC ${terms_src})
python_extension_module(sfepy_terms)
target_link_libraries(sfepy_terms sfepy_common)
target_include_directories(sfepy_terms PRIVATE ${NumPy_INCLUDE_DIRS} ${common_path})

add_cython_target(terms terms.pyx)
add_library(terms MODULE ${terms})
python_extension_module(terms)
target_include_directories(terms PRIVATE ${NumPy_INCLUDE_DIRS} ${common_path})
target_link_libraries(terms sfepy_common sfepy_terms)


install(TARGETS terms LIBRARY DESTINATION sfepy/terms/extmods)
