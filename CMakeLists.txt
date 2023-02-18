cmake_minimum_required(VERSION 3.4...3.22)
project(sfepy)

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)

add_subdirectory(sfepy/discrete/common/extmods)