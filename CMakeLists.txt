cmake_minimum_required(VERSION 3.4...3.22)
project(sfepy)

set(CMAKE_VERBOSE_MAKEFILE:BOOL ON)

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# this needs discussion, see PR description
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2")

if(MSVC)
    add_compile_definitions(SFEPY_PLATFORM=1 inline=__inline)
else()
   add_compile_definitions(SFEPY_PLATFORM=0 inline=inline)
endif()

add_subdirectory(sfepy/discrete/common/extmods)
add_subdirectory(sfepy/discrete/fem/extmods)
add_subdirectory(sfepy/discrete/iga/extmods)
add_subdirectory(sfepy/mechanics/extmods)
add_subdirectory(sfepy/terms/extmods)
