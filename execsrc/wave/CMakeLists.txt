
option(WAVE_EXEC "Build Wave executable" ON)
message("-- wave will be built : ${WAVE_EXEC}")

if(WAVE_EXEC)
    file(GLOB_RECURSE WAVE_EXEC_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
    file(GLOB_RECURSE WAVE_EXEC_HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

    add_executable (wave ${WAVE_EXEC_SOURCE_FILES} ${WAVEE_EXEC_HEADER_FILES})

    include_directories (${CMAKE_SOURCE_DIR}/soliton_src)
    add_library (libsoliton SHARED IMPORTED)
    set_target_properties(libsoliton PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/libs/libsoliton.so )
    target_link_libraries (wave libsoliton)
endif(WAVE_EXEC)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

