cmake_minimum_required(VERSION 3.10)
project(cnpytest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Find ZLIB package
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

# Find cnpy package
find_path(CNPY_INCLUDE_DIR cnpy.h)
find_library(CNPY_LIBRARY NAMES cnpy)

if(NOT CNPY_INCLUDE_DIR OR NOT CNPY_LIBRARY)
    message(FATAL_ERROR "Could not find cnpy. Please set CNPY_INCLUDE_DIR and CNPY_LIBRARY.")
endif()

include_directories(${CNPY_INCLUDE_DIR})

# Add your executable
add_executable(cnpytest npytest.cpp)

# Link against cnpy and zlib
target_link_libraries(cnpytest ${CNPY_LIBRARY} ${ZLIB_LIBRARIES})

