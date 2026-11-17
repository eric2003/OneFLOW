# src/infrastructure/CMakeLists.txt
message(STATUS "配置基础设施模块...")

add_library(infrastructure
    config.f90
    mesh.f90
)

target_link_libraries(infrastructure PUBLIC core)

target_include_directories(infrastructure PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

install(TARGETS infrastructure
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

message(STATUS "基础设施模块配置完成")