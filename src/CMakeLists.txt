file(GLOB_RECURSE SOLITON_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE SOLITON_HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)

add_library(
    Soliton
    SHARED
    ${SOLITON_SOURCE_FILES}
    ${SOLITON_HEADER_FILES}
    )

set_target_properties(Soliton PROPERTIES OUTPUT_NAME "soliton" COMPILE_FLAGS ${CMAKE_CXX_FLAGS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

option(BUILD_DOC "Build documentation" ON)

if(BUILD_DOC)
    message("-- Build documentation ON")
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/docs/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
#        message("Doxygen build started")

        add_custom_target(DocDoxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
        message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
else(BUILD_DOC)
    message("-- Build documentation OFF")
endif(BUILD_DOC)
