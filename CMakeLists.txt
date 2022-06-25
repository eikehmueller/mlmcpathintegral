cmake_minimum_required (VERSION 3.12)

project(MLMCMCPathIntegral)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define user options
option(USE_MPI "Use MPI" ON)

# Find dependencies
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)

if(USE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Initialize global list of cc files
set(SOURCES "" CACHE INTERNAL "sources")

include_directories(src)

add_subdirectory(src)