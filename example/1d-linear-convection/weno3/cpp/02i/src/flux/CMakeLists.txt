# src/flux/CMakeLists.txt
# 添加基础和高阶通量子目录

set(FLUX_SRC_DIR ${CMAKE_SOURCE_DIR}/src/flux)
set(FLUX_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/flux)

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

# 创建通量库（包含公共部分）
add_library(cfd_flux STATIC)
target_sources(cfd_flux PRIVATE ${FLUX_HEADERS} ${FLUX_SOURCES})

target_include_directories(cfd_flux PRIVATE
    ${CMAKE_SOURCE_DIR}/include
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
