
option(DFUNSTRUCTURED_EXEC "Build mesh convergence executable" ON)
message("-- dfunstructured will be built : ${DFUNSTRUCTURED_EXEC}")

if(DFUNSTRUCTURED_EXEC)
    file(GLOB_RECURSE DFUNSTRUCTURED_EXEC_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
    file(GLOB_RECURSE DFUNSTRUCTURED_EXEC_HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)

    add_executable (dfunstructured ${DFUNSTRUCTURED_EXEC_SOURCE_FILES} ${DFUNSTRUCTURED_EXEC_HEADER_FILES})

    include_directories (${CMAKE_SOURCE_DIR}/src)
    add_library (libsoliton SHARED IMPORTED)
    set_target_properties(libsoliton PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/libs/libsoliton.so )
    target_link_libraries (dfunstructured libsoliton)
endif(DFUNSTRUCTURED_EXEC)

