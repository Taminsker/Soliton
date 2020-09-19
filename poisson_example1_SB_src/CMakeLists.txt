
option(EXECUTABLE_POISSON_EXAMPLE_1_SB "Build example 1 for Poisson equation (with SB)" ON)

if(EXECUTABLE_POISSON_EXAMPLE_1_SB)
    message("-- example 1 (with SB) for Poisson equation will be built : ON")
    add_executable(PoissonExample1_sb main.cpp functions.cpp)
    target_include_directories (PoissonExample1_sb PUBLIC ${SOLITON_DIR})
    target_link_libraries(PoissonExample1_sb Soliton)
    target_link_libraries(PoissonExample1_sb ${EIGEN3_LIBRARIES})
else(EXECUTABLE_POISSON_EXAMPLE_1_SB)
    message("-- example 1 (with SB) for Poisson equation will be built : OFF")
endif(EXECUTABLE_POISSON_EXAMPLE_1_SB)

