# examples/CMakeLists.txt

message(STATUS "Building example programs...")

# 函数：添加一个示例程序
function(add_cfd_example example_name source_file)
    add_executable(${example_name} ${source_file})
    
    target_link_libraries(${example_name} PRIVATE cfd_core)
    
    target_include_directories(${example_name} PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        # 如果示例需要额外头文件，可在这里加
    )
    
    # 输出到 bin/examples/ 子目录（便于区分主程序和示例）
    set_target_properties(${example_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${EXAMPLES_OUTPUT_DIR}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG ${EXAMPLES_OUTPUT_DIR}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE ${EXAMPLES_OUTPUT_DIR}
        FOLDER "Examples"
    )
    
    # MSVC 设置（示例专属，可从库继承）
    if(MSVC)
        target_compile_options(${example_name} PRIVATE /utf-8)
    endif()
    
    message(STATUS "  Added example: ${example_name} (${source_file})")
endfunction()

# 自动发现子目录中的 main.cpp
file(GLOB EXAMPLE_SUBDIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} LIST_DIRECTORIES true "*")
set(AVAILABLE_EXAMPLES "")

foreach(subdir ${EXAMPLE_SUBDIRS})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${subdir})
        set(main_file ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/main.cpp)
        if(EXISTS ${main_file})
            set(target_name "${subdir}_example")  # 如 eno_weno_example
            add_cfd_example(${target_name} ${main_file})
            list(APPEND AVAILABLE_EXAMPLES ${target_name})
        endif()
    endif()
endforeach()

# 导出变量供其他地方使用（可选）
set(AVAILABLE_EXAMPLES ${AVAILABLE_EXAMPLES} CACHE INTERNAL "List of available CFD examples")

message(STATUS "Available examples: ${AVAILABLE_EXAMPLES}")