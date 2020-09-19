cmake_minimum_required(VERSION 3.5)

# Project settings
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ".")
set(PROJECT_UNIT_TESTS_NAME soliton_unit_tests)

# Dependencies
find_package(Boost COMPONENTS filesystem system unit_test_framework REQUIRED)

# Gather the header and source files
file(GLOB UNIT_TESTS_INC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file(GLOB UNIT_TESTS_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Include paths
set(UNIT_TESTS_INCLUDES ${SOLITON_DIR} ${Boost_INCLUDE_DIR})

# Assign the include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${UNIT_TESTS_INCLUDES})

# Build unit tests
add_definitions(-DBOOST_TEST_DYN_LINK)
#add_executable(${PROJECT_UNIT_TESTS_NAME} ${UNIT_TESTS_SRC_FILES} $<TARGET_OBJECTS:demo_objects> )
add_executable(${PROJECT_UNIT_TESTS_NAME} ${UNIT_TESTS_SRC_FILES} ${SOLITON_SOURCE_FILES})
target_link_libraries(${PROJECT_UNIT_TESTS_NAME} ${Boost_LIBRARIES})
target_include_directories (${PROJECT_UNIT_TESTS_NAME} PUBLIC ${SOLITON_DIR})
target_link_libraries(${PROJECT_UNIT_TESTS_NAME} Soliton)
target_link_libraries(${PROJECT_UNIT_TESTS_NAME} ${EIGEN3_LIBRARIES})

enable_testing()
add_test(SolitonTest ${PROJECT_UNIT_TESTS_NAME})
add_custom_command(
     TARGET ${PROJECT_UNIT_TESTS_NAME}
     COMMENT "Run tests"
     POST_BUILD
     COMMAND ${PROJECT_UNIT_TESTS_NAME} -l message
)
