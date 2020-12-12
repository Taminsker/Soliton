option(POISSON_EXEC "Build Poisson executable" ON)
message("-- poisson will be built : ${POISSON_EXEC}")

if(POISSON_EXEC)
    file(GLOB_RECURSE POISSON_EXEC_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
    file(GLOB_RECURSE POISSON_EXEC_HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)

    add_executable (poisson ${POISSON_EXEC_SOURCE_FILES} ${POISSON_EXEC_HEADER_FILES})

    include_directories (${CMAKE_SOURCE_DIR}/src)
    add_library (libsoliton SHARED IMPORTED)
    set_target_properties(libsoliton PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/libs/libsoliton.so )
    target_link_libraries (poisson libsoliton)
endif(POISSON_EXEC)

