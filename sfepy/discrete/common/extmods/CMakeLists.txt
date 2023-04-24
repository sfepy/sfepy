if(DEBUG_MESH)
    add_compile_definitions(DEBUG_MESH)
endif()

set(common_src fmfield.c geommech.c common_python.c)
add_library(sfepy_common STATIC ${common_src})
python_extension_module(sfepy_common)
target_include_directories(sfepy_common PRIVATE ${NumPy_INCLUDE_DIRS})

add_cython_target(_fmfield _fmfield.pyx)
add_library(_fmfield MODULE ${_fmfield})
python_extension_module(_fmfield)
target_include_directories(_fmfield PRIVATE ${NumPy_INCLUDE_DIRS})
target_link_libraries(_fmfield sfepy_common)

add_cython_target(cmapping cmapping.pyx)
add_library(cmapping MODULE ${cmapping})
python_extension_module(cmapping)
target_include_directories(cmapping PRIVATE ${NumPy_INCLUDE_DIRS})
target_link_libraries(cmapping sfepy_common)

add_cython_target(assemble assemble.pyx)
add_library(assemble MODULE ${assemble})
python_extension_module(assemble)
target_include_directories(assemble PRIVATE ${NumPy_INCLUDE_DIRS})

add_cython_target(cmesh cmesh.pyx)
add_library(cmesh MODULE ${cmesh} geomtrans.c mesh.c meshutils.c sort.c common_python.c)
python_extension_module(cmesh)
target_include_directories(cmesh PRIVATE ${NumPy_INCLUDE_DIRS})

add_cython_target(crefcoors crefcoors.pyx)
add_library(crefcoors MODULE ${crefcoors} refcoors.c geomtrans.c mesh.c)
python_extension_module(crefcoors)
target_include_directories(crefcoors PRIVATE ${NumPy_INCLUDE_DIRS})
target_link_libraries(crefcoors sfepy_common)

add_cython_target(_geommech _geommech.pyx)
add_library(_geommech MODULE ${_geommech})
python_extension_module(_geommech)
target_include_directories(_geommech PRIVATE ${NumPy_INCLUDE_DIRS})
target_link_libraries(_geommech sfepy_common)

install(
    TARGETS _fmfield cmapping assemble cmesh crefcoors _geommech
    LIBRARY DESTINATION sfepy/discrete/common/extmods
)
