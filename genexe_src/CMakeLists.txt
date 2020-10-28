
option(GENERIC_EXECUTABLE "Build generic executable (work with input.dat)" ON)

message("-- generic executable will be built : ${GENERIC_EXECUTABLE}")
if(GENERIC_EXECUTABLE)
    add_executable (GenExe main.cpp)
    target_include_directories (GenExe PUBLIC ${SOLITON_DIR})
    target_link_libraries (GenExe Soliton)
    target_link_libraries (GenExe ${EIGEN3_LIBRARIES})
endif(GENERIC_EXECUTABLE)

