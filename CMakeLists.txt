cmake_minimum_required(VERSION 3.10)

project(RadialProject)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find HDF5 package (both C and C++ component
find_package(HDF5 COMPONENTS C CXX REQUIRED)

# Find cnpy package
find_path(CNPY_INCLUDE_DIR cnpy.h)
find_library(CNPY_LIBRARY NAMES cnpy)

include_directories(${CNPY_INCLUDE_DIR})

# Specify the directories containing header files
include_directories(include)
# Include HDF5 header directories
include_directories(${HDF5_INCLUDE_DIRS})

# Add the library for your main source
add_library(radial SHARED src/radial.cpp)

# Add the executable for your tests
add_executable(runradial tests/main.cpp)

target_link_libraries(runradial radial ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${CNPY_LIBRARY})
