cmake_minimum_required(VERSION 2.6)

project(Soliton)
include(FindPkgConfig)


set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -ansi -Wwrite-strings -Wuninitialized -Wunused-result -Wunreachable-code -Wmissing-braces -Wmissing-field-initializers -ftrack-macro-expansion=0 -fno-diagnostics-show-caret -m64 -O3 -g -std=c++14" )

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


option(USE_VERBOSE "VERBOSE MODE" ON)
message(STATUS "VERBOSE MODE : ${USE_VERBOSE}")

if (USE_VERBOSE)
  add_definitions(-DVERBOSE)
endif()

option(PRINT_HASHMAP "Write Hashmap in file (can take a long time)" OFF)
message(STATUS "PRINT_HASHMAP : ${PRINT_HASHMAP}")

if (PRINT_HASHMAP)
  add_definitions(-DPRINTHASHMAP)
endif()

set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/libs)

set(TODAY "24-Jul-2020")

find_package(Boost COMPONENTS filesystem system unit_test_framework REQUIRED)

# Assign the include directories
if (NOT Boost_FOUND)
    message(STATUS "find_package : NOT found Boost library")
else()
    message(STATUS "find_package : found Boost library")
    include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(Eigen3 REQUIRED)

if (NOT EIGEN3_FOUND)
    message(STATUS "find_package : NOT found EIGEN library")
else()
    message(STATUS "find_package : found EIGEN library")
    include(${EIGEN3_USE_FILE})
    include_directories(${EIGEN3_INCLUDE_DIR})
    include_directories("${EIGEN3_INCLUDE_DIR}/Eigen/unsupported")
endif()

# Shared library Soliton
set(SOLITON_DIR ${CMAKE_SOURCE_DIR}/soliton_src)
add_subdirectory (${SOLITON_DIR})

# Unit tests
option(USE_UNIT_TEST "Build unit tests " ON)
if (USE_UNIT_TEST)
    message ("-- Build unit tests")
    set(SOLITON_TEST_UNIT_DIR ${CMAKE_SOURCE_DIR}/soliton_test_unit)
    add_subdirectory (${SOLITON_TEST_UNIT_DIR})
endif()


# Executable
add_subdirectory (${CMAKE_SOURCE_DIR}/genexe_src)
add_subdirectory (${CMAKE_SOURCE_DIR}/meshconvergence_src)


