set(common_path ../../common/extmods)

add_cython_target(igac igac.pyx)
add_library(igac MODULE ${igac} nurbs.c)
python_extension_module(igac)
target_include_directories(igac PRIVATE ${NumPy_INCLUDE_DIRS} ${common_path})
target_link_libraries(igac sfepy_common)

install(
    TARGETS igac
    LIBRARY DESTINATION sfepy/discrete/iga/extmods
)
