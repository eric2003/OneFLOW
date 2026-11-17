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

# 添加子目录 - 先构建库，再构建示例
add_subdirectory(src)          # 构建CFD核心库
#add_subdirectory(examples)     # 构建示例程序

# 启用文件夹组织（Visual Studio）
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(MSVC)
    # 正确方式：设置根目录的 VS_STARTUP_PROJECT 属性
    set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})

endif()
