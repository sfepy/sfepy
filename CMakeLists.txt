cmake_minimum_required(VERSION 3.4...3.22)
project(sfepy)

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)

add_subdirectory(sfepy/discrete/common/extmods)
add_subdirectory(sfepy/discrete/fem/extmods)
add_subdirectory(sfepy/discrete/iga/extmods)
add_subdirectory(sfepy/mechanics/extmods)
add_subdirectory(sfepy/terms/extmods)
