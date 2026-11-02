# src/numerics/flux/CMakeLists.txt
message(STATUS "Configuring flux calculator module...")

add_library(flux
    base.f90
    rusanov.f90
)

target_link_libraries(flux PUBLIC core infrastructure)

target_include_directories(flux PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

install(TARGETS flux
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

message(STATUS "Flux calculator module configured")