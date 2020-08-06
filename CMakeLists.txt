cmake_minimum_required(VERSION 2.6)

project(Soliton)
include(FindPkgConfig)


set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -ansi -Wwrite-strings -Wuninitialized -Wunused-result -Wunreachable-code -Wmissing-braces -Wmissing-field-initializers -m64 -O3 -g -std=c++14" )

message(STATUS "COMPILE OPTIONS : ${CMAKE_CXX_FLAGS}\n")

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release")
endif()

message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}")
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DDEBUG)
ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug")


option(USE_VERBOSE "VERBOSE MODE" OFF)
message(STATUS "VRBOSE MODE : ${USE_VERBOSE}")

if (USE_VERBOSE)
  add_definitions(-DVERBOSE)
endif()

set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/libs)

#set(CMAKE_BUILD_TYPE Release)
set(TODAY "24-Jul-2020")

#set(LIB_DIR ${CMAKE_SOURCE_DIR}/libs)

#include(${LIB_DIR}/lib.cmake)

find_package(Eigen3 REQUIRED)

if (NOT EIGEN3_FOUND)
    message(STATUS "find_package : NOT found EIGEN library")
else()
    message(STATUS "find_package : found EIGEN library")
    include(${EIGEN3_USE_FILE})
    include_directories(${EIGEN3_INCLUDE_DIR})
    include_directories("${EIGEN3_INCLUDE_DIR}/Eigen/unsupported")
endif()

set(SOLITON_DIR ${CMAKE_SOURCE_DIR}/soliton_src)
add_subdirectory (${SOLITON_DIR})
add_subdirectory (${CMAKE_SOURCE_DIR}/poisson_src)

