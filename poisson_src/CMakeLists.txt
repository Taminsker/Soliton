
option(EXECUTABLE_POISSON "Build executable for Poisson equation" ON)

if(EXECUTABLE_POISSON)
    message("-- the executable for Poisson equation will be built : ON")
    add_executable(
        Poisson
        main.cpp
        )

    target_include_directories (Poisson PUBLIC ${SOLITON_DIR})
    target_link_libraries(Poisson Soliton)
    target_link_libraries(Poisson ${EIGEN3_LIBRARIES})
else(EXECUTABLE_POISSON)
    message("-- the executable for Poisson equation will be built : OFF")
endif(EXECUTABLE_POISSON)

