# src/numerics/flux/CMakeLists.txt
message(STATUS "Configuring flux calculator module...")

# Start with just the base module
add_library(flux STATIC
    base.f90
	rusanov.f90
)

#target_link_libraries(flux PRIVATE core infrastructure)

target_include_directories(flux PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# 明确设置不使用 /Qm64 选项
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM|Intel")
    set_target_properties(flux PROPERTIES
        Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /Qm64"  # 明确设置，避免继承问题
    )
endif()

install(TARGETS flux
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

message(STATUS "Flux calculator module configured")