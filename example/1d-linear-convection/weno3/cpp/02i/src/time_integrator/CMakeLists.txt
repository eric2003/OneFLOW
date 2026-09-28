# src/reconstructor/CMakeLists.txt

set(TIME_INTEGRATOR_SRC_DIR ${CMAKE_SOURCE_DIR}/src/time_integrator)
set(TIME_INTEGRATOR_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/time_integrator)

set(TIME_INTEGRATOR_HEADERS
    ${TIME_INTEGRATOR_INCLUDE_DIR}/TimeIntegrator.h
    ${TIME_INTEGRATOR_INCLUDE_DIR}/TimeIntegratorFactory.h
    ${TIME_INTEGRATOR_INCLUDE_DIR}/rk1.h
    ${TIME_INTEGRATOR_INCLUDE_DIR}/rk2.h
    ${TIME_INTEGRATOR_INCLUDE_DIR}/rk3.h
)

set(TIME_INTEGRATOR_SOURCES
    ${TIME_INTEGRATOR_SRC_DIR}/TimeIntegrator.cpp
    ${TIME_INTEGRATOR_SRC_DIR}/TimeIntegratorFactory.cpp
    ${TIME_INTEGRATOR_SRC_DIR}/rk1.cpp
    ${TIME_INTEGRATOR_SRC_DIR}/rk2.cpp
    ${TIME_INTEGRATOR_SRC_DIR}/rk3.cpp
)

add_library(cfd_time_integrator STATIC)
target_sources(cfd_time_integrator PRIVATE ${TIME_INTEGRATOR_HEADERS} ${TIME_INTEGRATOR_SOURCES})

target_include_directories(cfd_time_integrator PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)


source_group("include/time_integrator" FILES
    ${CMAKE_SOURCE_DIR}/include/time_integrator/TimeIntegrator.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/TimeIntegratorFactory.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk1.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk2.h
    ${CMAKE_SOURCE_DIR}/include/time_integrator/rk3.h
)

source_group("src/time_integrator" FILES
    time_integrator/TimeIntegrator.cpp
    time_integrator/TimeIntegratorFactory.cpp
    time_integrator/rk1.cpp
    time_integrator/rk2.cpp
    time_integrator/rk3.cpp
)
