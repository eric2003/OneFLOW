# src/CMakeLists.txt
find_package(Threads REQUIRED)

set(RECONSTRUCTOR_SOURCES
    reconstructor/basic/FirstOrderReconstructor.cpp
    reconstructor/basic/SecondOrderReconstructor.cpp
    reconstructor/basic/MUSCLReconstructor.cpp
    reconstructor/advanced/EnoReconstructor.cpp
    reconstructor/advanced/WenoReconstructor.cpp
    reconstructor/advanced/MPWenoReconstructor.cpp
    reconstructor/ReconstructorFactory.cpp
)

set(FLUX_SOURCES
    flux/basic/RusanovFlux.cpp
    flux/basic/EngquistOsherFlux.cpp
    flux/basic/UpwindFlux.cpp
    flux/basic/CentralFlux.cpp
    flux/basic/LaxWendroffFlux.cpp
    flux/advanced/RoeFlux.cpp
    flux/advanced/HLLFlux.cpp
    flux/FluxFactory.cpp
)

set(TIME_INTEGRATOR_SOURCES
    time_integrator/TimeIntegrator.cpp
    time_integrator/rk1.cpp
    time_integrator/rk2.cpp
    time_integrator/rk3.cpp
    time_integrator/TimeIntegratorFactory.cpp
)

set(INITIAL_CONDITION_SOURCES
    initial_condition/StepInitialCondition.cpp
    initial_condition/SinusoidalInitialCondition.cpp
    initial_condition/GaussianInitialCondition.cpp
    initial_condition/CustomInitialCondition.cpp
    initial_condition/InitialConditionFactory.cpp
	initial_condition/InitialCondition.cpp
)

set(SOURCES
    config.cpp
    mesh.cpp
    domain.cpp
    solution.cpp
    ${FLUX_SOURCES}
    ${RECONSTRUCTOR_SOURCES}
    ${INITIAL_CONDITION_SOURCES}
    boundary.cpp
    residual.cpp
    ${TIME_INTEGRATOR_SOURCES}
    cfd_solver.cpp
    analysis.cpp
    main.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/initial_condition
    ${CMAKE_SOURCE_DIR}/include/time_integrator
    ${CMAKE_SOURCE_DIR}/include/flux
    ${CMAKE_SOURCE_DIR}/include/flux/basic
    ${CMAKE_SOURCE_DIR}/include/flux/advanced
    ${CMAKE_SOURCE_DIR}/include/reconstructor
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced
)

if(MSVC)
	target_compile_definitions( ${PROJECT_NAME} PRIVATE _USE_MATH_DEFINES )
    target_compile_options(${PROJECT_NAME} PRIVATE 
        "$<$<C_COMPILER_ID:MSVC>:/utf-8>"
        "$<$<CXX_COMPILER_ID:MSVC>:/utf-8>"
    )
endif()


