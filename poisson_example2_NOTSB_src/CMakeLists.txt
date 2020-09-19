
option(EXECUTABLE_POISSON_EXAMPLE_2_NSB "Build example 2 for Poisson equation (not SB)" ON)

if(EXECUTABLE_POISSON_EXAMPLE_2_NSB)
    message("-- example 2 (not SB) for Poisson equation will be built : ON")
    add_executable(PoissonExample2_nsb main.cpp functions.cpp)
    target_include_directories (PoissonExample2_nsb PUBLIC ${SOLITON_DIR})
    target_link_libraries(PoissonExample2_nsb Soliton)
    target_link_libraries(PoissonExample2_nsb ${EIGEN3_LIBRARIES})
else(EXECUTABLE_POISSON_EXAMPLE_2_NSB)
    message("-- example 2 (not SB) for Poisson equation will be built : OFF")
endif(EXECUTABLE_POISSON_EXAMPLE_2_NSB)

