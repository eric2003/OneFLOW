# src/CMakeLists.txt
find_package(Threads REQUIRED)

# ===============================
# Auto-generated CMake snippet
# ===============================

# === Source directory variables ===
set(CORE_SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(FLUX_SRC_DIR ${CMAKE_SOURCE_DIR}/src/flux)
set(INITIAL_CONDITION_SRC_DIR ${CMAKE_SOURCE_DIR}/src/initial_condition)
set(RECONSTRUCTOR_SRC_DIR ${CMAKE_SOURCE_DIR}/src/reconstructor)
set(TIME_INTEGRATOR_SRC_DIR ${CMAKE_SOURCE_DIR}/src/time_integrator)

# === Header directory variables ===
set(CORE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(FLUX_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/flux)
set(INITIAL_CONDITION_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/initial_condition)
set(RECONSTRUCTOR_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/reconstructor)
set(TIME_INTEGRATOR_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/time_integrator)

set(CORE_SOURCES
    ${CORE_SRC_DIR}/analysis.cpp
    ${CORE_SRC_DIR}/boundary.cpp
    ${CORE_SRC_DIR}/cfd_solver.cpp
    ${CORE_SRC_DIR}/config.cpp
    ${CORE_SRC_DIR}/domain.cpp
    ${CORE_SRC_DIR}/mesh.cpp
    ${CORE_SRC_DIR}/residual.cpp
    ${CORE_SRC_DIR}/solution.cpp
)

set(FLUX_SOURCES
    ${FLUX_SRC_DIR}/FluxFactory.cpp
    ${FLUX_SRC_DIR}/advanced/HLLFlux.cpp
    ${FLUX_SRC_DIR}/advanced/RoeFlux.cpp
    ${FLUX_SRC_DIR}/basic/CentralFlux.cpp
    ${FLUX_SRC_DIR}/basic/EngquistOsherFlux.cpp
    ${FLUX_SRC_DIR}/basic/LaxWendroffFlux.cpp
    ${FLUX_SRC_DIR}/basic/RusanovFlux.cpp
    ${FLUX_SRC_DIR}/basic/UpwindFlux.cpp
)

set(INITIAL_CONDITION_SOURCES
    ${INITIAL_CONDITION_SRC_DIR}/ComplexWaveInitialCondition.cpp
    ${INITIAL_CONDITION_SRC_DIR}/CustomInitialCondition.cpp
    ${INITIAL_CONDITION_SRC_DIR}/GaussianInitialCondition.cpp
    ${INITIAL_CONDITION_SRC_DIR}/InitialCondition.cpp
    ${INITIAL_CONDITION_SRC_DIR}/InitialConditionFactory.cpp
    ${INITIAL_CONDITION_SRC_DIR}/SinusoidalInitialCondition.cpp
    ${INITIAL_CONDITION_SRC_DIR}/StepInitialCondition.cpp
)

set(RECONSTRUCTOR_SOURCES
    ${RECONSTRUCTOR_SRC_DIR}/ReconstructorFactory.cpp
    ${RECONSTRUCTOR_SRC_DIR}/advanced/EnoReconstructor.cpp
    ${RECONSTRUCTOR_SRC_DIR}/advanced/MPWenoReconstructor.cpp
    ${RECONSTRUCTOR_SRC_DIR}/advanced/WenoReconstructor.cpp
    ${RECONSTRUCTOR_SRC_DIR}/basic/FirstOrderReconstructor.cpp
    ${RECONSTRUCTOR_SRC_DIR}/basic/MUSCLReconstructor.cpp
    ${RECONSTRUCTOR_SRC_DIR}/basic/SecondOrderReconstructor.cpp
)

set(TIME_INTEGRATOR_SOURCES
    ${TIME_INTEGRATOR_SRC_DIR}/TimeIntegrator.cpp
    ${TIME_INTEGRATOR_SRC_DIR}/TimeIntegratorFactory.cpp
    ${TIME_INTEGRATOR_SRC_DIR}/rk1.cpp
    ${TIME_INTEGRATOR_SRC_DIR}/rk2.cpp
    ${TIME_INTEGRATOR_SRC_DIR}/rk3.cpp
)

set(CORE_HEADERS
    ${CORE_INCLUDE_DIR}/analysis.hpp
    ${CORE_INCLUDE_DIR}/boundary.hpp
    ${CORE_INCLUDE_DIR}/cfd_solver.hpp
    ${CORE_INCLUDE_DIR}/config.hpp
    ${CORE_INCLUDE_DIR}/domain.hpp
    ${CORE_INCLUDE_DIR}/mesh.hpp
    ${CORE_INCLUDE_DIR}/residual.hpp
    ${CORE_INCLUDE_DIR}/solution.hpp
)

set(FLUX_HEADERS
    ${FLUX_INCLUDE_DIR}/FluxCalculator.h
    ${FLUX_INCLUDE_DIR}/FluxFactory.h
    ${FLUX_INCLUDE_DIR}/advanced/HLLFlux.h
    ${FLUX_INCLUDE_DIR}/advanced/RoeFlux.h
    ${FLUX_INCLUDE_DIR}/basic/CentralFlux.h
    ${FLUX_INCLUDE_DIR}/basic/EngquistOsherFlux.h
    ${FLUX_INCLUDE_DIR}/basic/LaxWendroffFlux.h
    ${FLUX_INCLUDE_DIR}/basic/RusanovFlux.h
    ${FLUX_INCLUDE_DIR}/basic/UpwindFlux.h
)

set(INITIAL_CONDITION_HEADERS
    ${INITIAL_CONDITION_INCLUDE_DIR}/ComplexWaveInitialCondition.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/CustomInitialCondition.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/GaussianInitialCondition.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/InitialCondition.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/InitialConditionFactory.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/SinusoidalInitialCondition.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/StepInitialCondition.h
)

set(RECONSTRUCTOR_HEADERS
    ${RECONSTRUCTOR_INCLUDE_DIR}/Reconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/ReconstructorFactory.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/advanced/EnoReconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/advanced/MPWenoReconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/advanced/WenoReconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/basic/FirstOrderReconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/basic/MUSCLReconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/basic/SecondOrderReconstructor.h
)

set(TIME_INTEGRATOR_HEADERS
    ${TIME_INTEGRATOR_INCLUDE_DIR}/TimeIntegrator.h
    ${TIME_INTEGRATOR_INCLUDE_DIR}/TimeIntegratorFactory.h
    ${TIME_INTEGRATOR_INCLUDE_DIR}/rk1.h
    ${TIME_INTEGRATOR_INCLUDE_DIR}/rk2.h
    ${TIME_INTEGRATOR_INCLUDE_DIR}/rk3.h
)

# Unified lists
set(ALL_SOURCES
    ${CORE_SOURCES}
    ${FLUX_SOURCES}
    ${INITIAL_CONDITION_SOURCES}
    ${RECONSTRUCTOR_SOURCES}
    ${TIME_INTEGRATOR_SOURCES}
)

set(ALL_HEADERS
    ${CORE_HEADERS}
    ${FLUX_HEADERS}
    ${INITIAL_CONDITION_HEADERS}
    ${RECONSTRUCTOR_HEADERS}
    ${TIME_INTEGRATOR_HEADERS}
)

add_library(cfd_core STATIC)

target_sources(cfd_core PRIVATE ${ALL_SOURCES} ${ALL_HEADERS})

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


