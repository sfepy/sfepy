set(common_path ../../common/extmods)

add_cython_target(bases bases.pyx)
add_library(bases MODULE ${bases} lagrange.c)
python_extension_module(bases)
target_include_directories(bases PRIVATE ${NumPy_INCLUDE_DIRS} ${common_path})
target_link_libraries(bases sfepy_common)

add_cython_target(lobatto_bases lobatto_bases.pyx)
add_library(lobatto_bases MODULE ${lobatto_bases} lobatto.c lobatto1d.c)
python_extension_module(lobatto_bases)
target_include_directories(lobatto_bases PRIVATE ${NumPy_INCLUDE_DIRS} ${common_path})
target_link_libraries(lobatto_bases sfepy_common)


install(
    TARGETS bases lobatto_bases
    LIBRARY DESTINATION sfepy/discrete/fem/extmods
)
