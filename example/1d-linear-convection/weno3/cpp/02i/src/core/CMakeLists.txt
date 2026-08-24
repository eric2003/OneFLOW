# src/core/CMakeLists.txt
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

# 创建核心库
add_library(cfd_core_impl STATIC)

target_sources(cfd_core_impl PRIVATE ${CORE_SOURCES} ${CORE_HEADERS})


# 包含目录
target_include_directories(cfd_core_impl PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

# VS源文件分组
source_group("src" FILES ${CORE_SOURCES})
source_group("include" FILES ${CORE_HEADERS})