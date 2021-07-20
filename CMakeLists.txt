cmake_minimum_required (VERSION 3.10)

project(MLMCMCPathIntegral)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)

include_directories(src)

set(SOURCES "" CACHE INTERNAL "sources")

add_subdirectory(src)
