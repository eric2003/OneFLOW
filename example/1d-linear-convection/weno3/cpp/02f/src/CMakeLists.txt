# src/CMakeLists.txt
find_package(Threads REQUIRED)

set(CORE_SOURCES
    analysis.cpp
    boundary.cpp
    cfd_solver.cpp
    config.cpp
    domain.cpp
    mesh.cpp
    residual.cpp
    solution.cpp
)

set(FLUX_SOURCES
    flux/FluxFactory.cpp
    flux/advanced/HLLFlux.cpp
    flux/advanced/RoeFlux.cpp
    flux/basic/CentralFlux.cpp
    flux/basic/EngquistOsherFlux.cpp
    flux/basic/LaxWendroffFlux.cpp
    flux/basic/RusanovFlux.cpp
    flux/basic/UpwindFlux.cpp
)

set(INITIAL_CONDITION_SOURCES
    initial_condition/ComplexWaveInitialCondition.cpp
    initial_condition/CustomInitialCondition.cpp
    initial_condition/GaussianInitialCondition.cpp
    initial_condition/InitialCondition.cpp
    initial_condition/InitialConditionFactory.cpp
    initial_condition/SinusoidalInitialCondition.cpp
    initial_condition/StepInitialCondition.cpp
)

set(RECONSTRUCTOR_SOURCES
    reconstructor/ReconstructorFactory.cpp
    reconstructor/advanced/EnoReconstructor.cpp
    reconstructor/advanced/MPWenoReconstructor.cpp
    reconstructor/advanced/WenoReconstructor.cpp
    reconstructor/basic/FirstOrderReconstructor.cpp
    reconstructor/basic/MUSCLReconstructor.cpp
    reconstructor/basic/SecondOrderReconstructor.cpp
)

set(TIME_INTEGRATOR_SOURCES
    time_integrator/TimeIntegrator.cpp
    time_integrator/TimeIntegratorFactory.cpp
    time_integrator/rk1.cpp
    time_integrator/rk2.cpp
    time_integrator/rk3.cpp
)

set(CORE_HEADERS
    ${CMAKE_SOURCE_DIR}/include/analysis.hpp
    ${CMAKE_SOURCE_DIR}/include/boundary.hpp
    ${CMAKE_SOURCE_DIR}/include/cfd_solver.hpp
    ${CMAKE_SOURCE_DIR}/include/config.hpp
    ${CMAKE_SOURCE_DIR}/include/domain.hpp
    ${CMAKE_SOURCE_DIR}/include/mesh.hpp
    ${CMAKE_SOURCE_DIR}/include/residual.hpp
    ${CMAKE_SOURCE_DIR}/include/solution.hpp
)

set(FLUX_HEADERS
    ${CMAKE_SOURCE_DIR}/include/flux/FluxCalculator.h
    ${CMAKE_SOURCE_DIR}/include/flux/FluxFactory.h
    ${CMAKE_SOURCE_DIR}/include/flux/advanced/HLLFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/advanced/RoeFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/CentralFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/EngquistOsherFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/LaxWendroffFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/RusanovFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/UpwindFlux.h
)

set(INITIAL_CONDITION_HEADERS
    ${CMAKE_SOURCE_DIR}/include/initial_condition/ComplexWaveInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/CustomInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/GaussianInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/InitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/InitialConditionFactory.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/SinusoidalInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/StepInitialCondition.h
)

set(RECONSTRUCTOR_HEADERS
    ${CMAKE_SOURCE_DIR}/include/reconstructor/Reconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/ReconstructorFactory.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/EnoReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/MPWenoReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/WenoReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/FirstOrderReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/MUSCLReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/SecondOrderReconstructor.h
)

set(TIME_INTEGRATOR_HEADERS
    ${CMAKE_SOURCE_DIR}/include/time_integrator/TimeIntegrator.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/TimeIntegratorFactory.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk1.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk2.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk3.h
)

# Unified source list
set(CORE_SOURCES
    ${CORE_SOURCES}
    ${FLUX_SOURCES}
    ${INITIAL_CONDITION_SOURCES}
    ${RECONSTRUCTOR_SOURCES}
    ${TIME_INTEGRATOR_SOURCES}
)

add_library(cfd_core STATIC ${CORE_SOURCES})

target_sources(cfd_core PRIVATE
    ${CMAKE_SOURCE_DIR}/include/analysis.hpp
    ${CMAKE_SOURCE_DIR}/include/boundary.hpp
    ${CMAKE_SOURCE_DIR}/include/cfd_solver.hpp
    ${CMAKE_SOURCE_DIR}/include/config.hpp
    ${CMAKE_SOURCE_DIR}/include/domain.hpp
    ${CMAKE_SOURCE_DIR}/include/mesh.hpp
    ${CMAKE_SOURCE_DIR}/include/residual.hpp
    ${CMAKE_SOURCE_DIR}/include/solution.hpp
    ${CMAKE_SOURCE_DIR}/include/flux/FluxCalculator.h
    ${CMAKE_SOURCE_DIR}/include/flux/FluxFactory.h
    ${CMAKE_SOURCE_DIR}/include/flux/advanced/HLLFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/advanced/RoeFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/CentralFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/EngquistOsherFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/LaxWendroffFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/RusanovFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/UpwindFlux.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/ComplexWaveInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/CustomInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/GaussianInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/InitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/InitialConditionFactory.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/SinusoidalInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/StepInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/Reconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/ReconstructorFactory.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/EnoReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/MPWenoReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/WenoReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/FirstOrderReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/MUSCLReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/SecondOrderReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/TimeIntegrator.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/TimeIntegratorFactory.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk1.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk2.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk3.h
)

# === VS Folder View ===
source_group("include" FILES
    ${CMAKE_SOURCE_DIR}/include/analysis.hpp
    ${CMAKE_SOURCE_DIR}/include/boundary.hpp
    ${CMAKE_SOURCE_DIR}/include/cfd_solver.hpp
    ${CMAKE_SOURCE_DIR}/include/config.hpp
    ${CMAKE_SOURCE_DIR}/include/domain.hpp
    ${CMAKE_SOURCE_DIR}/include/mesh.hpp
    ${CMAKE_SOURCE_DIR}/include/residual.hpp
    ${CMAKE_SOURCE_DIR}/include/solution.hpp
)
source_group("include/flux" FILES
    ${CMAKE_SOURCE_DIR}/include/flux/FluxCalculator.h
    ${CMAKE_SOURCE_DIR}/include/flux/FluxFactory.h
)
source_group("include/flux/advanced" FILES
    ${CMAKE_SOURCE_DIR}/include/flux/advanced/HLLFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/advanced/RoeFlux.h
)
source_group("include/flux/basic" FILES
    ${CMAKE_SOURCE_DIR}/include/flux/basic/CentralFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/EngquistOsherFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/LaxWendroffFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/RusanovFlux.h
    ${CMAKE_SOURCE_DIR}/include/flux/basic/UpwindFlux.h
)
source_group("include/initial_condition" FILES
    ${CMAKE_SOURCE_DIR}/include/initial_condition/ComplexWaveInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/CustomInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/GaussianInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/InitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/InitialConditionFactory.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/SinusoidalInitialCondition.h
    ${CMAKE_SOURCE_DIR}/include/initial_condition/StepInitialCondition.h
)
source_group("include/reconstructor" FILES
    ${CMAKE_SOURCE_DIR}/include/reconstructor/Reconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/ReconstructorFactory.h
)
source_group("include/reconstructor/advanced" FILES
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/EnoReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/MPWenoReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/WenoReconstructor.h
)
source_group("include/reconstructor/basic" FILES
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/FirstOrderReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/MUSCLReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/SecondOrderReconstructor.h
)
source_group("include/time_integrator" FILES
    ${CMAKE_SOURCE_DIR}/include/time_integrator/TimeIntegrator.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/TimeIntegratorFactory.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk1.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk2.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk3.h
)
source_group("src" FILES
    analysis.cpp
    boundary.cpp
    cfd_solver.cpp
    config.cpp
    domain.cpp
    mesh.cpp
    residual.cpp
    solution.cpp
)
source_group("src/flux" FILES
    flux/FluxFactory.cpp
)
source_group("src/flux/advanced" FILES
    flux/advanced/HLLFlux.cpp
    flux/advanced/RoeFlux.cpp
)
source_group("src/flux/basic" FILES
    flux/basic/CentralFlux.cpp
    flux/basic/EngquistOsherFlux.cpp
    flux/basic/LaxWendroffFlux.cpp
    flux/basic/RusanovFlux.cpp
    flux/basic/UpwindFlux.cpp
)
source_group("src/initial_condition" FILES
    initial_condition/ComplexWaveInitialCondition.cpp
    initial_condition/CustomInitialCondition.cpp
    initial_condition/GaussianInitialCondition.cpp
    initial_condition/InitialCondition.cpp
    initial_condition/InitialConditionFactory.cpp
    initial_condition/SinusoidalInitialCondition.cpp
    initial_condition/StepInitialCondition.cpp
)
source_group("src/reconstructor" FILES
    reconstructor/ReconstructorFactory.cpp
)
source_group("src/reconstructor/advanced" FILES
    reconstructor/advanced/EnoReconstructor.cpp
    reconstructor/advanced/MPWenoReconstructor.cpp
    reconstructor/advanced/WenoReconstructor.cpp
)
source_group("src/reconstructor/basic" FILES
    reconstructor/basic/FirstOrderReconstructor.cpp
    reconstructor/basic/MUSCLReconstructor.cpp
    reconstructor/basic/SecondOrderReconstructor.cpp
)
source_group("src/time_integrator" FILES
    time_integrator/TimeIntegrator.cpp
    time_integrator/TimeIntegratorFactory.cpp
    time_integrator/rk1.cpp
    time_integrator/rk2.cpp
    time_integrator/rk3.cpp
)

target_include_directories(cfd_core PRIVATE 
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

# MSVC 通用设置（放这里，示例会自动继承）
if(MSVC)
    target_compile_definitions(cfd_core PRIVATE _USE_MATH_DEFINES)
    target_compile_options(cfd_core PRIVATE /utf-8)  # 彻底解决 C4819 警告
endif()

# 可选：把源文件按文件夹组织（VS 中更好看）
set_property(TARGET cfd_core PROPERTY FOLDER "Core")


