set(CMAKE_FOLDER ${CMAKE_FOLDER}/euler)

set(ONEFLOW_EULER_BACKEND_DIR
  ${CMAKE_SOURCE_DIR}/ports/kunshan/oneflow_1d_hip
)

add_executable(oneflow_euler_backend_cpu_test
  euler_backend_contract_test.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEuler.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEulerBackend.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDWeno5.cpp
)

target_include_directories(oneflow_euler_backend_cpu_test
  PRIVATE
    ${ONEFLOW_EULER_BACKEND_DIR}
)

target_compile_options(oneflow_euler_backend_cpu_test
  PRIVATE
    -ffp-contract=off
)

target_link_libraries(oneflow_euler_backend_cpu_test
  PRIVATE
    GTest::gtest_main
)

gtest_discover_tests(oneflow_euler_backend_cpu_test
  DISCOVERY_MODE PRE_TEST
)

add_executable(oneflow_flux_backend_bridge_test
  flux_backend_bridge_test.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEuler.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEulerBackend.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDWeno5.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/CpuFluxBackend.cpp
)

target_include_directories(oneflow_flux_backend_bridge_test
  PRIVATE
    ${ONEFLOW_EULER_BACKEND_DIR}
    ${CMAKE_SOURCE_DIR}/codes/accel/include
    ${CMAKE_SOURCE_DIR}/codes/project/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_compile_options(oneflow_flux_backend_bridge_test
  PRIVATE
    -ffp-contract=off
)

target_link_libraries(oneflow_flux_backend_bridge_test
  PRIVATE
)

# Bridge test is a standalone main(), not a gtest
# Run manually: ./bin/oneflow_flux_backend_bridge_test


add_executable(oneflow_euler_domain_contract_test
  ${CMAKE_SOURCE_DIR}/tests/euler_domain_contract_test.cpp)
target_include_directories(oneflow_euler_domain_contract_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_domain_contract_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_domain_contract_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_domain_backend_adapter_test
  euler_domain_backend_adapter_test.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEuler.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEulerBackend.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEulerDomainAdapter.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDWeno5.cpp)
target_include_directories(oneflow_euler_domain_backend_adapter_test PRIVATE
  ${ONEFLOW_EULER_BACKEND_DIR}
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_compile_options(oneflow_euler_domain_backend_adapter_test PRIVATE
  -ffp-contract=off)
target_link_libraries(oneflow_euler_domain_backend_adapter_test PRIVATE
  GTest::gtest_main)
gtest_discover_tests(oneflow_euler_domain_backend_adapter_test
  DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_domain_state_registry_test
  ${CMAKE_SOURCE_DIR}/tests/euler_domain_state_registry_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateRegistry.cpp)
target_include_directories(oneflow_euler_domain_state_registry_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_domain_state_registry_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_domain_state_registry_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_domain_state_lifecycle_test
  ${CMAKE_SOURCE_DIR}/tests/euler_domain_state_lifecycle_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateRegistry.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateLifecycle.cpp)
target_include_directories(oneflow_euler_domain_state_lifecycle_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_domain_state_lifecycle_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_domain_state_lifecycle_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_cpu_domain_backend_test
  ${CMAKE_SOURCE_DIR}/tests/euler_cpu_domain_backend_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/CpuEulerDomainBackend.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateRegistry.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateLifecycle.cpp)
target_include_directories(oneflow_euler_cpu_domain_backend_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_cpu_domain_backend_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_cpu_domain_backend_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_domain_mrfield_adapter_test
  ${CMAKE_SOURCE_DIR}/tests/euler_domain_mrfield_adapter_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/EulerDomainMrFieldAdapter.cpp)
target_include_directories(oneflow_euler_domain_mrfield_adapter_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/main/include
  ${CMAKE_SOURCE_DIR}/codes/array/include
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_domain_mrfield_adapter_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_domain_mrfield_adapter_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_cpu_adapter_test
  ${CMAKE_SOURCE_DIR}/tests/euler_cpu_adapter_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerCpuAdapter.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/CpuFluxBackend.cpp)
target_include_directories(oneflow_euler_cpu_adapter_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_cpu_adapter_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_cpu_adapter_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_inv_flux_capability_test
  ${CMAKE_SOURCE_DIR}/tests/euler_inv_flux_capability_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerInvFluxCapability.cpp)
target_include_directories(oneflow_euler_inv_flux_capability_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_inv_flux_capability_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_inv_flux_capability_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_runge_kutta_capability_test
  ${CMAKE_SOURCE_DIR}/tests/euler_runge_kutta_capability_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/solver/src/EulerRungeKuttaCapability.cpp)
target_include_directories(oneflow_euler_runge_kutta_capability_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/solver/include
  ${CMAKE_SOURCE_DIR}/codes/global/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_runge_kutta_capability_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_runge_kutta_capability_test DISCOVERY_MODE PRE_TEST)

# Hardware-dependent HIP/DCU contract. This is opt-in so ordinary CPU builds
# keep the root CTest suite reproducible on machines without a DCU. On a target
# node, configure with CMAKE_HIP_COMPILER and CMAKE_HIP_ARCHITECTURES, then the
# same target is discovered by the root CTest tree with the `HIP.` prefix and
# `hip`/`dcu` labels.
# The resident-state generation/ownership contract is backend-neutral and
# runs without a DCU.  It gates the HIP production path's device/host ghost
# epoch bookkeeping, so it stays in the CPU test set.
add_executable(oneflow_hip_resident_state_contract_test
  ${CMAKE_SOURCE_DIR}/tests/hip_resident_state_contract_test.cpp)
target_include_directories(oneflow_hip_resident_state_contract_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_hip_resident_state_contract_test PRIVATE
  GTest::gtest_main)
gtest_discover_tests(oneflow_hip_resident_state_contract_test
  DISCOVERY_MODE PRE_TEST)

if (ONEFLOW_ENABLE_HIP_TESTS)
  if (NOT CMAKE_HIP_ARCHITECTURES)
    if (ONEFLOW_HIP_ARCHITECTURES)
      set(CMAKE_HIP_ARCHITECTURES "${ONEFLOW_HIP_ARCHITECTURES}"
          CACHE STRING "HIP/DCU architecture list" FORCE)
    else()
      message(FATAL_ERROR
        "ONEFLOW_ENABLE_HIP_TESTS=ON requires CMAKE_HIP_ARCHITECTURES "
        "or ONEFLOW_HIP_ARCHITECTURES")
    endif()
  endif()

  enable_language(HIP)
  find_package(hip CONFIG REQUIRED)
  include(${CMAKE_SOURCE_DIR}/cmake/OneFLOWEulerContract.cmake)
  oneflow_add_euler_hip_contract_test(
    oneflow_euler_backend_hip_test
    PORT_DIR ${ONEFLOW_EULER_BACKEND_DIR}
    ACCEL_DIR ${CMAKE_SOURCE_DIR}/codes/accel
    PROJECT_INC ${CMAKE_SOURCE_DIR}/codes/project/include
    BASIC_INC ${CMAKE_SOURCE_DIR}/codes/basic/include
    TEST_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/euler_backend_contract_test.cpp
    ADAPTER_TEST_SOURCE
      ${CMAKE_CURRENT_SOURCE_DIR}/euler_domain_backend_adapter_test.cpp)

  # Single-face Roe contract: the device flux must equal the host flux that
  # NsInvFluxPolicy produces for the same state.  Built HIP sources come from
  # the production backend so the test exercises the shipped kernel.
  set_source_files_properties(
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipFluxBackend.hip
    PROPERTIES LANGUAGE HIP)
  add_executable(oneflow_device_roe_contract_test
    ${CMAKE_CURRENT_SOURCE_DIR}/device_roe_contract_test.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipFluxBackend.hip
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipBackend.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/AccelBackend.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/AccelRuntime.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/CpuBackend.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/StageProfiler.cpp
    ${CMAKE_SOURCE_DIR}/codes/parallel/src/BoundaryExchangeContract.cpp
    ${CMAKE_SOURCE_DIR}/codes/ns/src/NsInvFluxPolicy.cpp)
  target_include_directories(oneflow_device_roe_contract_test PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/accel/include
    ${CMAKE_SOURCE_DIR}/codes/parallel/include
    ${CMAKE_SOURCE_DIR}/codes/ns/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/project/include)
  target_link_libraries(oneflow_device_roe_contract_test PRIVATE
    GTest::gtest_main hip::host)
  gtest_discover_tests(oneflow_device_roe_contract_test
    DISCOVERY_MODE PRE_TEST)

  add_executable(oneflow_device_limiter_contract_test
    ${CMAKE_CURRENT_SOURCE_DIR}/device_limiter_contract_test.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipFluxBackend.hip
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipBackend.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/AccelBackend.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/AccelRuntime.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/CpuBackend.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/StageProfiler.cpp
    ${CMAKE_SOURCE_DIR}/codes/parallel/src/BoundaryExchangeContract.cpp
    ${CMAKE_SOURCE_DIR}/codes/usolver/src/ULimiterPolicy.cpp)
  target_include_directories(oneflow_device_limiter_contract_test PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/accel/include
    ${CMAKE_SOURCE_DIR}/codes/parallel/include
    ${CMAKE_SOURCE_DIR}/codes/usolver/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/project/include)
  target_link_libraries(oneflow_device_limiter_contract_test PRIVATE
    GTest::gtest_main hip::host)
  gtest_discover_tests(oneflow_device_limiter_contract_test
    DISCOVERY_MODE PRE_TEST)

  add_test(
    NAME HIP.MainSolverFluxOneCall
    COMMAND $<TARGET_FILE:OneFLOWHipSmoke>)
  set_tests_properties(
    HIP.MainSolverFluxOneCall
    PROPERTIES LABELS "hardware;hip;dcu")

  # Derived temperature field for the viscous/heat-conduction port.  The device
  # must reproduce NsCalcGamaT's relation, including the nscom.statecoef form of
  # gamma, rather than only allocate a buffer.
  add_executable(oneflow_device_temperature_contract_test
    ${CMAKE_CURRENT_SOURCE_DIR}/device_temperature_contract_test.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipFluxBackend.hip
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipBackend.cpp
    # HipBackend.cpp's Initialize() calls RunHipBackendSelfTest(), defined in
    # HipKernel.hip.  Without ONEFLOW_ENABLE_HIP that call site is compiled out
    # and the symbol is never needed; this target defines it below, so the
    # definition has to come along or the link fails.
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipKernel.hip
    ${CMAKE_SOURCE_DIR}/codes/accel/src/AccelBackend.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/AccelRuntime.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/CpuBackend.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/StageProfiler.cpp
    ${CMAKE_SOURCE_DIR}/codes/parallel/src/BoundaryExchangeContract.cpp)
  target_include_directories(oneflow_device_temperature_contract_test PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/accel/include
    ${CMAKE_SOURCE_DIR}/codes/parallel/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/project/include)
  # HipBackend.cpp only registers the HIP backend inside ONEFLOW_ENABLE_HIP.
  # Without this definition the registry contains the CPU backend alone, the
  # temperature contract cannot bind a device state, and every case is skipped
  # instead of failing loudly.  The test allocates through DeviceBuffer, so it
  # must run on the real execution space rather than a host stand-in.
  target_compile_definitions(oneflow_device_temperature_contract_test PRIVATE
    ONEFLOW_ENABLE_HIP
    ONEFLOW_DEFAULT_ACCEL_BACKEND="HIP")
  target_link_libraries(oneflow_device_temperature_contract_test PRIVATE
    GTest::gtest_main hip::host)
  gtest_discover_tests(oneflow_device_temperature_contract_test
    DISCOVERY_MODE PRE_TEST)

  # Spalart-Allmaras one-equation contract: single-cell source and
  # single-face scalar fluxes, compared against the production host functions
  # (TurbCom::CalcSrcSa, SaComputeEddyViscosity, the TurbInvFluxPolicy and
  # TurbVisFluxPolicy helpers).
  #
  # Linked against OneFLOW_RuntimeObjects rather than a hand-picked source
  # list.  The host reference for the source term is the real
  # TurbCom::CalcSrcSa driven through the `turbcom` singleton, and that
  # translation unit drags in nscom, turb_trans, vis_model, chem and the
  # parallel/zone globals; enumerating them here would be a fragile list that
  # breaks whenever the host chain gains a dependency.  RuntimeObjects already
  # carries HipTurbulenceBackend.hip when HIP is enabled, so the device functions
  # under test are the shipped ones and are not compiled twice.
  add_executable(oneflow_device_sa_contract_test
    ${CMAKE_CURRENT_SOURCE_DIR}/device_sa_contract_test.cpp)
  # OneFLOW_RuntimeObjects publishes its include dirs PUBLIC, but a target that
  # only links the OBJECT library does not inherit them unless it asks; the
  # host reference compiles TurbCom.cpp's own include graph here, and that
  # reaches NsCom.h -> HXArray.h.  Naming array/ and global/ explicitly is what
  # keeps the test's compilation independent of the object target's private
  # dirs.
  target_include_directories(oneflow_device_sa_contract_test PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/accel/include
    ${CMAKE_SOURCE_DIR}/codes/turb/include
    ${CMAKE_SOURCE_DIR}/codes/ns/include
    ${CMAKE_SOURCE_DIR}/codes/array/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/project/include)
  target_link_libraries(oneflow_device_sa_contract_test PRIVATE
    OneFLOW_RuntimeObjects
    GTest::gtest_main
    hip::host)
  gtest_discover_tests(oneflow_device_sa_contract_test
    DISCOVERY_MODE PRE_TEST)

  # Single-face viscous-flux contract: the device kernel must equal the host
  # UNsVisFlux chain for the same face, both sides populated, to 1e-14 relative.
  # The host reference links NsVisFluxPolicy, UGradPolicy and Sutherland, i.e.
  # the functions the solver itself calls, so neither side can drift alone.
  # ViscousFluxArithmetic.h is the shared inline decomposition the kernel
  # compiles; this target does not depend on the SA track's files.
  add_executable(oneflow_device_viscous_flux_contract_test
    ${CMAKE_CURRENT_SOURCE_DIR}/device_viscous_flux_contract_test.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipFluxBackend.hip
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipBackend.cpp
    # HipBackend.cpp's Initialize() calls RunHipBackendSelfTest(), defined in
    # HipKernel.hip.  Without ONEFLOW_ENABLE_HIP that call site is compiled out
    # and the symbol is never needed; this target defines it below, so the
    # definition has to come along or the link fails.
    ${CMAKE_SOURCE_DIR}/codes/accel/src/HipKernel.hip
    ${CMAKE_SOURCE_DIR}/codes/accel/src/AccelBackend.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/AccelRuntime.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/CpuBackend.cpp
    ${CMAKE_SOURCE_DIR}/codes/accel/src/StageProfiler.cpp
    ${CMAKE_SOURCE_DIR}/codes/parallel/src/BoundaryExchangeContract.cpp
    ${CMAKE_SOURCE_DIR}/codes/ns/src/NsVisFluxPolicy.cpp
    ${CMAKE_SOURCE_DIR}/codes/usolver/src/UGradPolicy.cpp
    ${CMAKE_SOURCE_DIR}/codes/model/src/Sutherland.cpp)
  target_include_directories(oneflow_device_viscous_flux_contract_test PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/accel/include
    ${CMAKE_SOURCE_DIR}/codes/parallel/include
    ${CMAKE_SOURCE_DIR}/codes/ns/include
    ${CMAKE_SOURCE_DIR}/codes/usolver/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
    ${CMAKE_SOURCE_DIR}/codes/model/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/project/include)
  # HipBackend.cpp only registers the HIP backend inside ONEFLOW_ENABLE_HIP.
  # Without these definitions the registry holds the CPU backend alone,
  # AccelRuntime selects it, and every device case compared against a CPU
  # look-alike - which is a false pass rather than a failure.  Observed on this
  # target: "OneFLOW compute backend: CPU" and device flux 0 against host
  # 2.28e-12, i.e. the kernel never ran.
  target_compile_definitions(oneflow_device_viscous_flux_contract_test PRIVATE
    ONEFLOW_ENABLE_HIP
    ONEFLOW_DEFAULT_ACCEL_BACKEND="HIP")
  target_link_libraries(oneflow_device_viscous_flux_contract_test PRIVATE
    GTest::gtest_main hip::host)
  gtest_discover_tests(oneflow_device_viscous_flux_contract_test
    DISCOVERY_MODE PRE_TEST)
endif()
