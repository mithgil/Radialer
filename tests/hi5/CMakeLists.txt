cmake_minimum_required(VERSION 3.10)
project(HDF5Example)

# Find HDF5 package (both C and C++ components)
find_package(HDF5 COMPONENTS C CXX REQUIRED)

# Include HDF5 header directories
include_directories(${HDF5_INCLUDE_DIRS})

# Create an executable from main.cpp
add_executable(h5complex h5complex.cpp)

# Link HDF5 libraries
target_link_libraries(h5complex ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
