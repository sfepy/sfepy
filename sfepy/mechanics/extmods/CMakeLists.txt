
set(common_path ../../discrete/common/extmods)

add_cython_target(ccontres ccontres.pyx)
add_library(ccontres MODULE ${ccontres} contres.c)
python_extension_module(ccontres)
target_include_directories(ccontres PRIVATE ${NumPy_INCLUDE_DIRS} ${common_path})
target_link_libraries(ccontres sfepy_common)


install(TARGETS ccontres LIBRARY DESTINATION sfepy/mechanics/extmods)
