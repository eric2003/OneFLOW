# CMakeLists.txt
cmake_minimum_required(VERSION 4.2.1)

project(OneFLOW_CFD LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

if(MSVC)
    foreach(CONFIG Debug Release RelWithDebInfo MinSizeRel)
        string(TOUPPER ${CONFIG} CONFIG_UPPER)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG_UPPER} ${CMAKE_BINARY_DIR}/bin)
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG_UPPER} ${CMAKE_BINARY_DIR}/lib)
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG_UPPER} ${CMAKE_BINARY_DIR}/lib)
    endforeach()
endif()

# 定义项目版本
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")

# 打印项目信息
message(STATUS "==========================================")
message(STATUS "  OneFLOW-CFD Project")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "==========================================")

# ==================== 其他 ====================
# 启用文件夹组织（Visual Studio）
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# 示例专用输出目录（子目录，便于区分）
set(EXAMPLES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/bin/examples)

# 添加子目录 - 先构建库，再构建示例
add_subdirectory(src)        # 先构建核心库
add_subdirectory(examples)   # 再构建示例程序

# ==================== 启动项目（放松：不强制设置） ====================
# 注释掉或删除下面的代码，让用户手动选择想调试的示例
# if(MSVC)
#     set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT eno_weno_example)
# endif()

message(STATUS "==========================================")
message(STATUS "OneFLOW-CFD Project - Multi-Example Build")
message(STATUS "Version: ${PROJECT_VERSION}(示例数量将由 examples/CMakeLists.txt 自动检测)")
message(STATUS "==========================================")