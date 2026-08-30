# src/CMakeLists.txt
find_package(Threads REQUIRED)


# 创建主库，链接所有子模块
add_library(cfd_core INTERFACE)

# 链接所有子模块
target_link_libraries(cfd_core INTERFACE
    cfd_core_impl
    cfd_flux
    cfd_reconstructor
    cfd_time_integrator
    cfd_initial_condition
)

target_include_directories(cfd_core INTERFACE
	${CMAKE_SOURCE_DIR}/../include
)

# MSVC 通用设置（放这里，示例会自动继承）
if(MSVC)
    #target_compile_definitions(cfd_core INTERFACE _USE_MATH_DEFINES)
    #target_compile_options(cfd_core INTERFACE /utf-8)  # 彻底解决 C4819 警告
	add_compile_definitions(_USE_MATH_DEFINES)
	add_compile_options(/utf-8)
endif()

# 可选：把源文件按文件夹组织（VS 中更好看）
set_property(TARGET cfd_core PROPERTY FOLDER "Core")

# 添加子目录
add_subdirectory(core)
add_subdirectory(flux)
add_subdirectory(reconstructor)
add_subdirectory(time_integrator)
add_subdirectory(initial_condition)



