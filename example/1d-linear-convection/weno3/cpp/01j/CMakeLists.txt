# CMakeLists.txt
cmake_minimum_required(VERSION 4.2.1)

project(OneFLOW_CFD LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 启用数学库
find_package(Threads REQUIRED)

set(RECONSTRUCTOR_SOURCES
    src/reconstructor/basic/FirstOrderReconstructor.cpp
    src/reconstructor/basic/SecondOrderReconstructor.cpp
    src/reconstructor/basic/MUSCLReconstructor.cpp
    src/reconstructor/advanced/EnoReconstructor.cpp
    src/reconstructor/advanced/WenoReconstructor.cpp
    src/reconstructor/advanced/MPWenoReconstructor.cpp
    src/reconstructor/ReconstructorFactory.cpp
)

set(FLUX_SOURCES
    src/flux/basic/RusanovFlux.cpp
    src/flux/basic/EngquistOsherFlux.cpp
    src/flux/basic/UpwindFlux.cpp
    src/flux/basic/CentralFlux.cpp
    src/flux/basic/LaxWendroffFlux.cpp
    src/flux/advanced/RoeFlux.cpp
    src/flux/advanced/HLLFlux.cpp
    src/flux/FluxFactory.cpp
)

set(TIME_INTEGRATOR_SOURCES
    src/time_integrator/TimeIntegrator.cpp
    src/time_integrator/rk1.cpp
    src/time_integrator/rk2.cpp
    src/time_integrator/rk3.cpp
    src/time_integrator/TimeIntegratorFactory.cpp
)

set(INITIAL_CONDITION_SOURCES
    src/initial_condition/StepInitialCondition.cpp
    src/initial_condition/SinusoidalInitialCondition.cpp
    src/initial_condition/GaussianInitialCondition.cpp
    src/initial_condition/CustomInitialCondition.cpp
    src/initial_condition/InitialConditionFactory.cpp
	src/initial_condition/InitialCondition.cpp
)

# 更新总源文件列表
set(SOURCES
    src/config.cpp
    src/mesh.cpp
    src/domain.cpp
    src/solution.cpp
    ${FLUX_SOURCES}
    ${RECONSTRUCTOR_SOURCES}
    ${INITIAL_CONDITION_SOURCES}
    src/boundary.cpp
    src/residual.cpp
    ${TIME_INTEGRATOR_SOURCES}
    src/cfd_solver.cpp
    src/analysis.cpp
    src/main.cpp
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
	set_property(TARGET ${PROJECT_NAME} 
				 PROPERTY VS_STARTUP_PROJECT 
				 ${PROJECT_NAME}
	)
endif()

