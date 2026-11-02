# src/numerics/reconstructor/CMakeLists.txt
message(STATUS "配置重构器模块...")

add_library(reconstructor
    base.f90
    eno.f90
    weno3.f90
    weno5.f90
    factory.f90
)

target_link_libraries(reconstructor PUBLIC core infrastructure)

target_include_directories(reconstructor PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

install(TARGETS reconstructor
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

message(STATUS "重构器模块配置完成")