# src/initial_condition/CMakeLists.txt

set(INITIAL_CONDITION_SRC_DIR ${CMAKE_SOURCE_DIR}/src/initial_condition)
set(INITIAL_CONDITION_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/initial_condition)

set(INITIAL_CONDITION_HEADERS
    ${INITIAL_CONDITION_INCLUDE_DIR}/ComplexWaveInitialCondition.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/CustomInitialCondition.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/GaussianInitialCondition.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/InitialCondition.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/InitialConditionFactory.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/SinusoidalInitialCondition.h
    ${INITIAL_CONDITION_INCLUDE_DIR}/StepInitialCondition.h
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


add_library(cfd_initial_condition STATIC)
target_sources(cfd_initial_condition PRIVATE ${INITIAL_CONDITION_HEADERS} ${INITIAL_CONDITION_SOURCES})

target_include_directories(cfd_initial_condition PRIVATE
    ${CMAKE_SOURCE_DIR}/include
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

source_group("src/initial_condition" FILES
    initial_condition/ComplexWaveInitialCondition.cpp
    initial_condition/CustomInitialCondition.cpp
    initial_condition/GaussianInitialCondition.cpp
    initial_condition/InitialCondition.cpp
    initial_condition/InitialConditionFactory.cpp
    initial_condition/SinusoidalInitialCondition.cpp
    initial_condition/StepInitialCondition.cpp
)