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
    kinds.f90
	string_utils.f90
	io_utils.f90
	config.f90               # 配置模块
    mesh.f90                 # 网格模块
    domain.f90               # 计算域模块
    solution.f90             # 解变量模块
    flux/rusanov_flux.f90
	flux/engquist_osher_flux.f90
	flux/inviscid_flux.f90
	flux/flux.f90                 # 通量模块
    boundary.f90             # 边界条件模块
    initial_condition/step.f90
	initial_condition/initial_condition.f90    # 初始条件和解析解模块
	reconstruct/recon_coefficients.f90
	reconstruct/eno.f90
	reconstruct/weno3.f90
	reconstruct/reconstruct.f90
	residual.f90
	time_integration/rk1.f90
	time_integration/rk2.f90
	time_integration/time_integration.f90
    cfd_solver.f90           # 主求解器
	run_eno_weno.f90
)



add_executable(oneflow_cfd ${SOURCES})

target_include_directories(oneflow_cfd PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
