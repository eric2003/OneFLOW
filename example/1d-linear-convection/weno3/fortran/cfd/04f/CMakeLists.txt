# CMakeLists.txt
cmake_minimum_required(VERSION 4.2.1)

project(OneFLOW_CFD LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)  # 共享库(.so/.dylib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(SOURCES
    src/kinds.f90
	src/string_utils.f90
	src/io_utils.f90
	src/config.f90               # 配置模块
    src/mesh.f90                 # 网格模块
    src/domain.f90               # 计算域模块
    src/solution.f90             # 解变量模块
    src/flux/rusanov_flux.f90
	src/flux/engquist_osher_flux.f90
	src/flux/inviscid_flux.f90
	src/flux/flux.f90                 # 通量模块
    src/boundary.f90             # 边界条件模块
    src/initial_condition/sine.f90
	src/initial_condition/step.f90
	src/initial_condition/initial_condition.f90    # 初始条件和解析解模块
	src/reconstruct/recon_coefficients.f90
	src/reconstruct/eno.f90
	src/reconstruct/weno3.f90
	src/reconstruct/reconstruct.f90
	src/residual.f90
	src/time_integration/rk1.f90
	src/time_integration/rk2.f90
	src/time_integration/time_integration.f90
    src/cfd_solver.f90           # 主求解器
	src/run_eno_weno.f90
)



add_executable(oneflow_cfd ${SOURCES})

target_include_directories(oneflow_cfd PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
