# src/CMakeLists.txt
# 源代码目录的CMake配置

message(STATUS "配置源代码目录...")

# 添加核心模块子目录
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/core AND 
   EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/core/CMakeLists.txt)
    add_subdirectory(core)
else()
    message(WARNING "core目录或CMakeLists.txt不存在")
endif()

# 可选：添加其他模块子目录
# if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/infrastructure)
#     add_subdirectory(infrastructure)
# endif()

message(STATUS "源代码目录配置完成")