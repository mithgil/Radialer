cmake_minimum_required(VERSION 3.10)
project(cnpytest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Find cnpy package
find_path(CNPY_INCLUDE_DIR cnpy.h)
find_library(CNPY_LIBRARY NAMES cnpy)

include_directories(${CNPY_INCLUDE_DIR})

# Add your executable
add_executable(cnpytest npytest.cpp)

# Link against cnpy
target_link_libraries(cnpytest ${CNPY_LIBRARY})
