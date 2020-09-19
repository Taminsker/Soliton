
option(EXECUTABLE_POISSON_EXAMPLE_1_NSB "Build example 1 for Poisson equation (not SB)" ON)

if(EXECUTABLE_POISSON_EXAMPLE_1_NSB)
    message("-- example 1 (not SB) for Poisson equation will be built : ON")
    add_executable(PoissonExample1_nsb main.cpp functions.cpp)
    target_include_directories (PoissonExample1_nsb PUBLIC ${SOLITON_DIR})
    target_link_libraries(PoissonExample1_nsb Soliton)
    target_link_libraries(PoissonExample1_nsb ${EIGEN3_LIBRARIES})
else(EXECUTABLE_POISSON_EXAMPLE_1_NSB)
    message("-- example 1 (not SB) for Poisson equation will be built : OFF")
endif(EXECUTABLE_POISSON_EXAMPLE_1_NSB)

