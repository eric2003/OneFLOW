# RFC: 共享的生产运行时目标（`oneflow_case_runtime`）

Proposer module: runtime-io
Date: 2026-11-01
Status: proposed (only a human changes this)

> 文件名说明：`.agentlab/scope/runtime-io.txt` 只允许 `.agentlab/RFC/runtime-io-*`，
> 因此日期放在主题之后（`runtime-io-2026-11-01-...`）；日期与主题仍按
> `.agentlab/RFC/<日期>-<module>-<主题>.md` 的约定记录在文件头。
> 若人工批准后希望统一改名，请在批准时同时更新 `scope/runtime-io.txt`。

## Background and problem

审计基线 `7328dc504342f6f8027065e9b8ac2894234e5878`。

`docs/AI_CONTEXT.md` 的 R1 路线与 `docs/plans/setup-case-environment-integration-test.md`
的最后一节把剩余缺口列为：生产 `SimuContext::TeardownCase()`、完整的
`SimuImp::RunCase()` 异常边界、进程运行时保持。审计确认这三项都无法用现有测试目标
表达，原因是**链接边界**，而不是断言缺失：

1. `SimuContext::TeardownCase()` 定义在 `codes/main/src/SimuContextEnv.cpp:100`。
   该翻译单元同时定义了 `SetupProcessEnvironment`/`SetupEnvironment`/
   `FinalizeEnvironment`/`TeardownEnvironment`/`ResolveTaskFromControl`，
   并引用下列模块的符号（`#include` 见 `SimuContextEnv.cpp:25-56`）：
   `Parallel`、`AccelRuntime`、`SolverMap`、`SolverNameList`、`Zone`、`FieldManager`、
   `DataBase`/`DataBaseType`、`HeatFlux`、`NsCom`、`NsSolver`、`INsSolver`、
   `TurbSolver`、`TurbCom`、`Tolerance`、`LogFile`、`ActionState`、`Iteration`、
   `Prj`、`ParaFile`、`SimuCaseLifecycle`（头文件模板）、`SimuTask`。
   只要链接这个 TU，以上符号都必须可解析，即使测试只调用 `TeardownCase()`。
2. 现有 `control_configuration_test`（`tests/io/CMakeLists.txt:59-94`）刻意只列
   `SimuContextBootstrap.cpp`、`SimuContextIO.cpp`、`SimuContextCaseEnvironment.cpp`，
   **不列** `SimuContextEnv.cpp`。因此 `TeardownCase` 既未链接也无覆盖；
   集成测试只能手工"镜像"三步清理
   （`tests/io/control_configuration_ingestion_test.cpp:340-344`）。
3. 生产可执行目标 `OneFLOW` 用 `ConstructSolutionDirTree` 递归收集 `codes/` 全部源码
   （`codes/CMakeLists.txt:347-414`）再加全局编译定义，**没有可复用的库目标**；
   测试目标只能手写源清单。手写清单已经产生两次构建事故（参见
   `docs/AI_CONTEXT.md`：`control_configuration_test` 源文件名错误、
   `OneFLOW_GeometryCore` 漏列 `PIOFileAccess.cpp`）。

结论：缺少一个"生产与测试共用、源清单只写一次"的构建边界。审计不修改 CMake，
本 RFC 只提出目标定义与迁移步骤，供人工决定。

## Proposal (specific to interface signatures)

**P1｜新增 CMake source-group 模块，作为唯一源清单来源**

新增 `cmake/OneFlowRuntimeSources.cmake`（人工所有），只定义变量，不建目标：

```cmake
# cmake/OneFlowRuntimeSources.cmake  (proposal, ASCII)
set(ONEFLOW_CASE_CONFIG_SOURCES
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlManifestReader.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlScriptReader.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlScriptParser.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlFileReader.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlDatabaseIO.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ConfigLoader.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ConfigDatabaseAdapter.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlParameterResolver.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/LegacyParameterSyntax.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/TextFileParser.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/CommentLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContextBootstrap.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContextIO.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContextCaseEnvironment.cpp
)

set(ONEFLOW_CASE_LIFECYCLE_SOURCES
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContextEnv.cpp
  # + SimuContextEnv.cpp 的链接闭包，需在迁移时用链接错误迭代确定，
  #   首次提案不臆造清单
)

set(ONEFLOW_PATH_FILE_SOURCES
  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/FileUtils.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/FileO.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PIO.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PIOFileAccess.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/LogFile.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/SimulationCommandLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PrjCmdLineArgs.cpp
)
```

**P2｜新增静态库目标（仅当 P1 被接受后启用）**

在 `codes/CMakeLists.txt` 的 `PRJ_*` 变量就绪之后、`add_executable(${PROJECT_NAME} ...)`
之前加入（草案）：

```cmake
add_library(oneflow_case_runtime STATIC
  ${ONEFLOW_CASE_CONFIG_SOURCES}
  ${ONEFLOW_CASE_LIFECYCLE_SOURCES}
  ${ONEFLOW_PATH_FILE_SOURCES}
)
target_include_directories(oneflow_case_runtime PUBLIC ${PRJ_INCLUDE_DIRS})
target_compile_features(oneflow_case_runtime PUBLIC ${PRJ_COMPILE_FEATURES})
target_compile_definitions(oneflow_case_runtime PUBLIC ${PRJ_COMPILE_DEF})
target_link_libraries(oneflow_case_runtime PUBLIC ${PRJ_LIBRARIES})
```

要点（这些是接口契约，不是实现细节）：

- `PUBLIC` 传播 include/compile definitions：`SimuContextEnv.cpp` 有
  `#ifdef ONEFLOW_ENABLE_HIP`（`SimuContextEnv.cpp:34,107`）与
  `HX_PARALLEL`/`ENABLE_MPI` 等条件编译，测试目标必须看到与生产一致的宏，
  否则出现"同一 TU 两种编译结果"的 ABI/行为分歧。
- 目标**不**参与数值内核编译；`FieldSimu`、`GridGeneration`、`ns`/`ins`/`turb`
  的实现仍在可执行目标里。因此该库只对"不执行任务体"的测试有意义。
- 库不得导出新的公共头文件或改名任何函数；`codes/**` 的对外符号保持不变。
- 若 P1 只落地 `ONEFLOW_CASE_CONFIG_SOURCES`，则先建
  `oneflow_case_config`（仅摄取/装配），让 `control_configuration_test` 变成
  `add_executable(... tests... ) + target_link_libraries(... oneflow_case_config)`；
  `ONEFLOW_CASE_LIFECYCLE_SOURCES` 的第二阶段再决定是否值得为清理契约付出闭包。

**P3｜可选的最小化前提（需另行批准）**

若人工认为 `SimuContextEnv.cpp` 的闭包太大，可先做一次结构变更：把
`TeardownCase` 的清理步骤列表提升为可枚举的声明
（例如 `codes/main/include/SimuContextCaseCleanup.h` 中的
`constexpr`/函数指针表 `CaseCleanupSteps()`），`TeardownCase` 只负责遍历执行。
这样清理**顺序与项数**可以在不链接 solver 的前提下断言。
评估：这会让测试驱动一次生产结构变更，且"清理步骤表"可能被误读为可插拔策略；
按 `AGENTS.md`"不要仅为引入模式而重构"的要求，建议**先走 P1/P2**，
只有在 P2 的闭包被证明不可维护时才考虑 P3。

## Affected modules and existing code

| 变更点 | 现状 | 影响 |
|---|---|---|
| `codes/CMakeLists.txt:389` | 单一 `add_executable(${PROJECT_NAME} ...)`，递归源收集 | 增加库目标 + 让可执行链接该库；源文件仍只编译一次，但 ARM/Windows 中文路径与 `_CRT_SECURE_NO_WARNINGS` 等定义必须保持一致 |
| `tests/io/CMakeLists.txt:59-94` | 手写 30 个源 | 可缩减为测试源 + 链接库；`tests/main/CMakeLists.txt:174-192`（`simu_context_parameter_test`）同理 |
| `tests/geometry/CMakeLists.txt` | 手写 `OneFLOW_GeometryCore` 源清单，曾漏 `PIOFileAccess.cpp` | 改为复用 `ONEFLOW_PATH_FILE_SOURCES` 后该事故类型消失 |
| `codes/main/src/SimuContextEnv.cpp` | 生产唯一清理路径 | 不修改；只是新增可链接性 |
| 平台脚本 | 无 | 不改 |

风险与回滚：库目标完全是加法；若失败，删除 `add_library` 与 `target_link_libraries`
即回到现状，production 行为不变（同一批源仍由可执行目标编译）。

## Alternatives

1. **继续复制源清单到每个测试目标**（现状）。成本：三次事故已证明不可控
   （另外 `tests/task/CMakeLists.txt` 也手写 40+ 源）。
2. **让测试直接执行真实可执行**（端到端策略 B）。优点：不改构建。缺点：
   `test/grid/rae2822convert` 的 `simutask=Grid` 会真的生成网格；
   `test/cfd/rae2822_roe_sa` 的 `simutask=Solve` 会真的求解，不能用于 setup-only 断言；
   失败定位差，且无法表达"清理失败后仍保留原异常"这类白盒契约。
3. **P3：清理步骤表化**。见上文评估，建议延后。
4. **给 `SimuContextEnv.cpp` 做编译期拆分**（把 `TeardownCase` 单独放进一个 TU）。
   成本低，但 `TeardownCase` 引用的符号一个都不会减少，因此**不解决**链接闭包问题；
   只能减少编译单元，属无关变更。

## Impact on verification cases

- 不新增数值 verification case；`docs/AI_CONTEXT.md` 的记录显示
  `.agentlab/verification/run_all.sh` 仍是 `UNVERIFIED` 占位，本 RFC 不改变这一点。
- 生效后新增的契约测试（见报告 §6 的 M1/M2）：
  `CaseTeardownContract.ClearsCaseStateButKeepsProcessBinding`、
  `RunCaseFailureStillTearsDownCaseAndKeepsProcessRuntime` 等，
  与 `ControlReaderIntegration` 共用临时 case 夹具，不依赖 `test/` 下的算例，
  因此不引入新的运行时间或数值基线。
- 多 rank 契约（报告 §6 M3）不受本 RFC 影响，仍是独立的
  `mpirun -np 2` 证据（本地允许，`AGENTS.domain.md` 规则 1）。
- 与 R-1（`io_paths.txt` 无写者过滤）的关系：该缺陷的修复会改
  `codes/main/src/SimuContextIO.cpp`，属 P1 (`case-config-ingestion`) 范围，
  但方案需要人来定（是改用 `PIO::OpenCaseFile` 语义，还是由 rank 判断包一层），
  建议在本 RFC 批准后单独开一条 RFC，不要与构建变更混在一个提交里。
