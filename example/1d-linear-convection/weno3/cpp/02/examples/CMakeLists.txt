# examples/CMakeLists.txt

add_executable(OneFLOW_CFD_Example main.cpp)

# 链接核心库
target_link_libraries(OneFLOW_CFD_Example PRIVATE OneFLOW_CFD_Core)

# 包含头文件（通常继承自库的 PUBLIC include）
target_include_directories(OneFLOW_CFD_Example PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    # ... 如果需要额外路径
)

if(MSVC)
    target_compile_options(OneFLOW_CFD_Example PRIVATE /utf-8)
endif()

# 可选：把这个 target 放进 Examples 文件夹
set_property(TARGET OneFLOW_CFD_Example PROPERTY FOLDER "Examples")