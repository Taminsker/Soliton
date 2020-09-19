
option(EXECUTABLE_POISSON_CONVERGENCE_NSB "Build example convergence for Poisson equation (not SB)" ON)

if(EXECUTABLE_POISSON_CONVERGENCE_NSB)
    message("-- example convergence (not SB) for Poisson equation will be built : ON")
    add_executable(PoissonExampleConv_nsb main.cpp functions.cpp)
    target_include_directories (PoissonExampleConv_nsb PUBLIC ${SOLITON_DIR})
    target_link_libraries(PoissonExampleConv_nsb Soliton)
    target_link_libraries(PoissonExampleConv_nsb ${EIGEN3_LIBRARIES})
else(EXECUTABLE_POISSON_CONVERGENCE_NSB)
    message("-- example convergence (not SB) for Poisson equation will be built : OFF")
endif(EXECUTABLE_POISSON_CONVERGENCE_NSB)

