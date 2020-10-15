
option(MESH_CONVERGENCE "Build mesh convergence" ON)

message("-- mesh convergence will be built : ${MESH_CONVERGENCE}")

if(MESH_CONVERGENCE)
    add_executable(MeshConvergence main.cpp functions.cpp)
    target_include_directories (MeshConvergence PUBLIC ${SOLITON_DIR})
    target_link_libraries(MeshConvergence Soliton)
    target_link_libraries(MeshConvergence ${EIGEN3_LIBRARIES})
endif(MESH_CONVERGENCE)

