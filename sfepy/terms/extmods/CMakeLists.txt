set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2")

set(common_path ../../discrete/common/extmods)

file(GLOB terms_src "*.c")
add_library(sfepy_terms SHARED ${terms_src})
python_extension_module(sfepy_terms)
target_link_libraries(sfepy_terms sfepy_common)
target_include_directories(sfepy_terms PRIVATE ${NumPy_INCLUDE_DIRS} ${common_path})

install(TARGETS sfepy_terms LIBRARY DESTINATION sfepy/terms/extmods)
