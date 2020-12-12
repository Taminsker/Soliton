
option(CONV4MESH_EXEC "Build mesh convergence executable" ON)
message("-- conv4mesh will be built : ${CONV4MESH_EXEC}")

if(CONV4MESH_EXEC)
    file(GLOB_RECURSE CONV4MESH_EXEC_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
    file(GLOB_RECURSE CONV4MESH_EXEC_HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)

    add_executable (conv4mesh ${CONV4MESH_EXEC_SOURCE_FILES} ${CONV4MESH_EXEC_HEADER_FILES})

    include_directories (${CMAKE_SOURCE_DIR}/src)
    add_library (libsoliton SHARED IMPORTED)
    set_target_properties(libsoliton PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/libs/libsoliton.so )
    target_link_libraries (conv4mesh libsoliton)
endif(CONV4MESH_EXEC)

