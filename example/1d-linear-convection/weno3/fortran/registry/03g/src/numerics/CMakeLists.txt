# src/numerics/CMakeLists.txt
message(STATUS "配置数值方法模块...")

add_subdirectory(reconstructor)
add_subdirectory(flux)
add_subdirectory(time_integration)  # 新增

message(STATUS "数值方法模块配置完成")