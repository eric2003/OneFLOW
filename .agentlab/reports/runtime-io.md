# Report: runtime-io

status: need_decision
updated: 2026-11-01
commit: ff5a34c6e766c62c369e2d3342e4d152bdd6e4f0 (本报告所在提交；后续仅追加提交标识的元数据提交)
branch: agent/runtime-io
model: scnetx/deepseek-v4.1-flash-261017
pi_version: 1.4.0
gate_mode: audit
gate_status: PASS

## 0. 审计范围与证据边界

本次为 **只读审计**（`AGENTS.domain.md`、`.agentlab/AGENTS.target.md` 规则 3）。

- 本地检出的审计基线：`agent/runtime-io` HEAD `7328dc504342f6f8027065e9b8ac2894234e5878`
  （`dev` 为 `df8676fc24a35439d355eff20a6fdca71b342313`，即 HEAD 的后续提交；审计对象为 HEAD 内容）。
  未切分支、未 rebase、未 push、未 merge。
- 阅读并交叉核对的文档：仓库根 `AGENTS.md`、`docs/AI_CONTEXT.md`、
  `docs/source/runtime_architecture.rst`、`docs/plans/setup-case-environment-integration-test.md`、
  `docs/source/database_architecture.rst`（存在性）、`docs/README.md` 索引约定。
- 阅读的代码：`codes/main/**`、`codes/io/**`、`codes/basic/src/Fatal.cpp`、
  `codes/parallel/src/Parallel.cpp`、`codes/database/{include,src}` 的公开接口、
  `codes/CMakeLists.txt`、`tests/CMakeLists.txt`、`tests/main/CMakeLists.txt`、`tests/io/CMakeLists.txt`。
- 证据等级：**全部结论来自源码与测试源码的静态阅读**。本次 **没有** 配置、编译、运行
  CTest、运行 `mpirun`，也没有查看任何 CI 运行结论（环境中未做网络访问）。
  因此本报告不使用 "verified by test run" 一类表述；凡"现有测试覆盖 X"均指
  "源码中存在断言 X 的测试用例"，不是"X 已在本机通过"。
- 依据 `AGENTS.md` "source + tests 决定实现状态"：本报告不采信 `docs/AI_CONTEXT.md`
  中可能滞后的状态描述，凡与代码不一致处显式标出（见 §8）。

## 1. 入口与调用者（生产调用链）

`main()` 只有一条生产入口，两条调度分支：

`codes/main/src/OneFlow.cpp:26` `main` → `Simulation(argc, argv)` →
`Simulation::Run()`（`codes/main/src/Simulation.cpp:169`，捕获异常并返回
`EXIT_FAILURE`）→ `Simulation::RunImpl()`（`Simulation.cpp:119`）。

| 阶段 | 入口符号 | 位置 | 直接调用者 | 被调用的 IO/运行时操作 |
|---|---|---|---|---|
| 分支判定 | `ClassifySimulationInvocation` | `codes/io/src/SimulationCommandLine.cpp` | `Simulation::RunImpl` | 纯函数，argv 个数决定 Lightweight / 参数错误 / Full |
| 轻量测试路径 | `Simulation::MakeDefaultSimulation` | `Simulation.cpp:66-117` | `RunImpl` | `TestRegistry::Instance()`（`codes/main/src/SimuBase.cpp:26`）、环境变量 `ONEFLOW_DEFAULT_TEST` |
| 全仿真参数 | `ParseSimulationCommandLineOptions` | `codes/io/src/SimulationCommandLine.cpp` | `RunImpl:145`、`SimuContext::ProcessCommandLine:63` | 无 IO |
| 进程路径状态 | `Prj::Init()` | `codes/io/src/Prj.cpp:179` | `RunImpl:149` **唯一调用者** | 填 `execute_dir/current_dir/system_root`（`FindSystemRoot` 向上查找 `system/action/actionFileList.txt`） |
| 进程运行时 | `SimuContext::SetupProcessEnvironment` | `codes/main/src/SimuContextEnv.cpp:75` | `RunImpl:157`（进程 context）、`SimuImp::Run:53` | `SetUpParallelEnvironment()`（`codes/parallel/src/Parallel.cpp:168`）、`InitializeAccelRuntime` |
| 多 case 顺序执行 | `RunCasesWithProcessLifecycle` | `codes/main/include/SimuCaseLifecycle.h:129` | `RunImpl:156` | 头文件模板，无 IO |
| 单 case 对象 | `SimuImp(caseDir, debug)` | `codes/main/src/SimuImp.cpp:41` | `RunImpl:161` | 构造 `SimuContext(caseDir, debug)`（`codes/main/src/SimuContextCaseEnvironment.cpp:21`）→ `Prj::ResolveCaseDir` |
| case 生命周期 | `SimuImp::RunCase` | `SimuImp.cpp:58` | `RunImpl:162`、`SimuImp::Run:54` | `ExecuteWithCleanup(PreProcess+MainProcess, PostProcess)` |
| case 配置 | `SimuContext::SetupCaseEnvironment` | `SimuContextCaseEnvironment.cpp:28` | `SimuImp::InitSimu:94` | `SimuContextDetail::InitializeCaseConfiguration(caseDir_)` |
| 配置装配 | `InitializeCaseConfiguration` | `codes/main/src/SimuContextBootstrap.cpp:27` | `SetupCaseEnvironment:34` | `Prj::SetPrjBaseDir` → `logFile.SetCaseDir` → `ReadControlInfo(caseDir)` → `ReadCaseIOSettings` → `ApplyCaseIOSettings` |
| 控制文件读取 | `ReadControlInfo(const std::string&)` | `codes/io/src/ControlFileReader.cpp:41` | `SimuContextBootstrap.cpp:32` **唯一生产调用者** | `[server] ReadPrjScript` → `TestSayHelloFromEveryProcess` → `BroadcastControlParameterToAllProcessors` → `DumpDataBase(caseDir)` |
| 脚本摄取 | `ReadPrjScript(caseDir)` | `codes/io/src/ControlScriptReader.cpp:49` | `ControlFileReader.cpp:45` | `ReadScriptFileNameList(caseDir, …)` + `ReadMultiScriptFiles(…)` |
| 清单发现 | `ReadScriptFileNameList(caseDir, …)` | `codes/io/src/ControlManifestReader.cpp:53` | `ControlScriptReader.cpp:52` | `TextFileParser::OpenCaseFile(caseDir, "script/control.txt")`，按序 push |
| 多文件解析 | `ReadMultiScriptFiles` | `codes/io/src/ControlScriptParser.cpp:46` | `ControlScriptReader.cpp:39/53` | 逐文件 `ReadOneFLOWScriptFile(fileName)` |
| 单文件提交 | `ReadOneFLOWScriptFile` + `ConfigDatabaseAdapter::Commit` | `ControlScriptParser.cpp:39`、`codes/io/src/ConfigDatabaseAdapter.cpp:29` | 同上 | `ConfigLoader(GetParameterArraySize)` → `ProcessData` 写入全局 `DataBase` |
| 数组长度解析 | `GetParameterArraySize` | `codes/io/src/ControlParameterResolver.cpp:30` | `ConfigLoader` 注入点 | 数字字面量或 `GetDataValue<int>(word)` |
| 参数广播 | `BroadcastControlParameterToAllProcessors` | `codes/io/src/ControlDatabaseIO.cpp:55` | `ControlFileReader.cpp:37/49` | `HXBcast(CompressData, DecompressData, serverid)`；`nProc<=1` 直接 return（`Parallel.cpp:196`） |
| 参数落盘 | `DumpDataBase(caseDir)` | `ControlDatabaseIO.cpp:46` | `ControlFileReader.cpp:50` | `PIO::OpenCaseFile(… "log/database.log")` → `DataPara::DumpData` → `PIO::CloseFile` |
| 路径应用 | `ReadCaseIOSettings` / `ApplyCaseIOSettings` | `codes/main/include/SimuContextParameter.h:94`、`codes/main/src/SimuContextIO.cpp:34` | `SimuContextBootstrap.cpp:35,37` | `Prj::SetProjectBaseDir/SetIOBaseDirs/SetDomainIOBaseDirs` + `Prj::WriteIOPathInfo`（mesh、flow 输出根，`SimuContextIO.cpp:63,70`） |
| 任务解析 | `SimuContext::ResolveTaskFromControl` | `codes/main/src/SimuContextEnv.cpp:181` | `SimuImp::RunSimu:99` | `simu_state.Init()`（全局，`codes/global/src/SimuDef.cpp:29`）→ `TaskEnum` |
| 任务分派 | `TaskRegistry::Create(taskName)` → `ISimuTask::Execute(ctx)` | `codes/main/src/SimuTask.cpp:37`、`codes/main/src/SimuTaskReg.cpp:43-123` | `SimuImp::RunSimu:110,121` | `Solve`→`FieldSimuRunPipeline(ctx)`；`Grid/Partition`→`GenerateGrid(ctx.CaseDir())`；`WallDist/Theory/ToyModel/PostTask/FunctionTest`→各自的全局自由函数 |
| case 清理 | `SimuContext::TeardownCase` | `SimuContextEnv.cpp:100` | `SimuImp::PostProcess:84` | 19 个 `ExecuteAllCleanupSteps` 步骤，含 `logFile.ClearCaseDir`、`Prj::ClearPrjBaseDir`、`DataPara/DataField.Clear`、`envReady_=false` |
| 进程收尾 | `SimuContext::FinalizeEnvironment` | `SimuContextEnv.cpp:138` | `RunImpl:164`、`SimuImp::FinalizeEnvironment:89`，异常路径经 `TeardownEnvironment` | `FinalizeAccelRuntime()` → `HXFinalize()`（`codes/parallel/src/BasicParallel.cpp:30`） |

无生产调用者的兼容 API（静态阅读确认，`grep` 全域无调用点）：

- `ReadControlInfo()`（`ControlFileReader.cpp:29`）、`DumpDataBase()`（`ControlDatabaseIO.cpp:36`）、
  `ReadPrjScript()`（`ControlScriptReader.cpp:35`）、`Prj::ProcessCmdLineArgs`（`Prj.cpp:103`）、
  `SimuImp(std::vector<std::string>&)` + `SimuImp::Run()`（`SimuImp.cpp:34,50`）、
  `SimuContext::SetupEnvironment/TeardownEnvironment`（`SimuContextEnv.cpp:92,174`）。

## 2. 状态 / 资源所有权

| 资源 | 存储位置 | 生命周期所有者 | 跨 case 是否复用 | 释放点 |
|---|---|---|---|---|
| 进程路径锚点 `execute_dir/current_dir/system_root` | `Prj` 静态字段（`codes/io/src/Prj.cpp:89-101`） | 进程 | 是（`Prj::Init` 只调用一次，`Prj::ClearPrjBaseDir` 不清这些） | 进程退出 |
| 当前 case 根与 IO 根 | `Prj::prjBaseDir/inputBaseDir/outputBaseDir/projectBaseDir/mesh*BaseDir/flow*BaseDir` | **全局单份**，被当作 case 状态 | 否，靠 `TeardownCase` 复位 | `SimuContextEnv.cpp:126` `Prj::ClearPrjBaseDir()` |
| 控制参数（`DataPara`）与字段注册（`DataField`） | 全局单例 `RequireGlobalDataBase()`（`codes/database/...`） | 进程静态、内容为 case 级 | 否，靠 teardown 清空 | `SimuContextEnv.cpp:123-124` |
| 日志流与 case 绑定 | 全局 `LogFile logFile`（`codes/io/src/LogFile.cpp:30`）的 `my_fstream/caseDir_/ifReWrite_` | 全局对象、case 级绑定 | 否，`SetCaseDir` 先 `Close()` 再重置 | `SimuContextEnv.cpp:125` |
| 并行状态 `pid/nProc/serverid/mode` | `Parallel` 静态字段（`codes/parallel/src/Parallel.cpp:38-42`） | 进程 | 是 | `HXFinalize()` |
| 任务态 `simu_state` / `TaskState` / `ActionState` / `ZoneState` | 全局 | `simu_state` 是进程静态、内容为 case 级；其余多数在 `TeardownCase` 清除 | 部分 | `SimuContextEnv.cpp:119-122` |
| 加速器状态 | `SimuContext::accelStates_`（`codes/main/include/SimuContext.h:137`），context 私有 | **case/上下文**（真正按实例拥有） | 否 | `ClearAccelStates()` + `FinalizeAccelRuntime()` |
| `expandedSolverNames_`、`task_`、`taskName_`、`envReady_`、`caseDir_` | `SimuContext` 成员 | case/上下文 | 否 | `TeardownCase` 末步复位；`caseDir_` 在 teardown 中**不复位**（有意：teardown 后仍可查询） |
| 文件句柄 | 各调用点栈上 `std::fstream`；`logFile.my_fstream`（public）；`ActionState::file` | 调用点 / 全局 | — | `Prj::CloseFile`、`PIO::CloseFile`、`LogFile::Close` |
| 注册表（`TaskRegistry`/`TestRegistry`） | 函数内 static 单例（`SimuTask.cpp:28`、`SimuBase.cpp:28`） | 进程 | 是 | 进程退出 |

结论性观察（关系到 §10 的 scope 切分）：

1. **唯一真正按 case 实例化拥有的是 `SimuContext` 的少数成员**（`accelStates_`、
   `expandedSolverNames_`、`task_`、`envReady_`、`caseDir_`）；其余 case 语义状态都在
   进程级静态里，正确性依赖 `TeardownCase` 的 19 步清理顺序。
2. `TeardownCase` 与 `Prj` 复位是 `logFile.ClearCaseDir()` → `Prj::ClearPrjBaseDir()` 的
   先后关系，顺序被代码固定（`SimuContextEnv.cpp:125-126`），但没有任何测试固定该顺序。
3. `DumpDataBase(caseDir)` 的文件所有权由 `PIO` 按 `Parallel::pid != Parallel::GetFid()`
   过滤（`codes/io/src/PIOFileAccess.cpp:44`）；`ApplyCaseIOSettings` 写 `io_paths.txt`
   走的是 `Prj::WriteIOPathInfo` → `Prj::OpenPrjFile`（**未经 PIO 过滤**，
   `codes/main/src/SimuContextIO.cpp:63`），见 §7 R-1。

## 3. 从代码推导的边界（不使用 starter 的 spatial/temporal 划分）

本仓库的 `codes/` 目录是"按物理/算法模块 + 按基础设施能力"混合切分的，
`spatial`/`temporal` 不是它的边界（例如 `codes/time`、`codes/residual`、`codes/grad`
都是独立目录，且共享 `codes/solver`、`codes/lhs` 的全局状态）。从 runtime-io 视角，
代码与测试支持的边界是下面三条，且它们**互不相同**：

- **B1 配置摄取（纯数据变换 + 全局参数写入）**
  `TextFileParser`/`CommentLine`/`ConfigLoader`/`LegacyParameterSyntax` → `ConfigDocument`
  → `ConfigDatabaseAdapter` → 全局 `DataBase`。不触碰文件系统根、不触碰并行。
  证据：`tests/io/ConfigLoader_test.cpp` 与 `tests/io/control_configuration_ingestion_test.cpp`
  可以在不初始化并行环境的情况下断言类型与覆盖顺序（例如
  `ControlConfigurationIngestion.ReadsManifestScriptsAndCommitsTypedValuesInOrder` 只调
  `ReadPrjScript`，该函数仅做清单+解析，见 `ControlScriptReader.cpp:49`）。
- **B2 路径与文件访问（进程级、全局可变状态、被大量模块依赖）**
  `Prj` + `FileUtils` + `PIO`/`PIOFileAccess` + `LogFile`。这是"共享基础"：
  41 个源文件引用 `Prj::`，分布在 17 个模块目录（见 §4）。它同时承载"解析路径规则"
  与"打开/截断文件"两种副作用。
- **B3 运行时生命周期与任务分派（跨模块的编排）**
  `Simulation`/`SimuImp`/`SimuContext*`/`SimuTask*`。`SimuTaskReg.cpp` 把字符串任务名
  绑定到 `FieldSimuRunPipeline`、`GenerateGrid`、`WalldistSimu` 等 **numerics 入口**，
  因此 B3 的改动天然跨 workstream（见 §4、§10 的规则）。

`codes/main` 与 `codes/io` 之间**不是**分层关系而是双向依赖：`codes/io` 的
`Prj`/`LogFile` 被 `codes/main` 调用，而 `codes/main` 的
`SimuContextParameter.h`/`SimuContextIO.cpp` 又承载了"控制参数 → 路径"的解释。
`tests/main/simu_context_parameter_test` 直接链接 `codes/io/src/Prj.cpp`
（`tests/main/CMakeLists.txt:190`），这是该双向依赖的测试侧证据。

## 4. 共享文件：不能被单个 worker 独占

统计方式：`grep -rl` 覆盖 `codes/**/*.cpp|*.h`。

| 文件 | 外部依赖面 | 为什么不能归一个 worker |
|---|---|---|
| `codes/io/src/Prj.cpp` + `codes/io/include/Prj.h` | 41 个文件 / 17 个目录（io 8、geometry 6、scalar 5、main 5、special 4、global 3、zone/visual/usolver/theory/test/solver/register/physics/force/cgns 各 1） | 路径规则一改，全部 mesh/solver/post 的输入输出位置同时变化；`geometry/src/CalcGrid.cpp:172,422`、`scalar/src/ScalarGrid.cpp:825,1815` 直接决定网格文件落盘位置 |
| `codes/io/include/LogFile.h` + `src/LogFile.cpp` | io 3、main 2、zone/parallel/geometry 各 1 | 全局 `logFile` 的 case 绑定是 `SetupCaseEnvironment` 与 `TeardownCase` 双方的契约 |
| `codes/io/src/PIOFileAccess.cpp` | 生产可执行 + `tests/geometry/CMakeLists.txt` 的 `OneFLOW_GeometryCore` | 曾因源清单未同步导致 `GeometryCore.lib(Zone.obj)` 链接失败（`docs/AI_CONTEXT.md` 已记录该历史），属跨 workstream 的构建边界 |
| `codes/io/src/{ControlFileReader,ControlScriptReader,ControlScriptParser,ControlManifestReader,ControlDatabaseIO}.cpp` + `codes/io/include/ParaFile.h` | 一个头文件承载 9 个函数的声明，生产者是 io、主要消费者是 `codes/main/src/SimuContextBootstrap.cpp` | 任何签名变动同时影响 main 与 tests/io 的 4 个测试目标 |
| `codes/main/src/SimuTaskReg.cpp` | 同时包含 io 语义（`ctx.CaseDir()`）与 numerics 入口（`FieldSimu.h`、`GridGeneration.h`、`Theory.h`） | 任务名表是 runtime 与 numerics 之间唯一的字符串契约；两侧都会想改 |
| `codes/main/include/SimuContext.h` | 被 `codes/main/*`、`codes/main/src/SimuTaskReg.cpp`、`codes/time`? 否；被 `SimuTask.h`/`SimuTaskRequire.h` 共享 | `SetupCaseEnvironment`（配置侧）、`ResolveTaskFromControl`（生命周期侧）、`accelStates_`（数值侧）三类消费者 |
| `codes/database/include/{DataBase,DataPara,DataObject}.h` | main 2、io 1、database 2 + 全部数值模块通过 `GetDataValue<T>(name)` 读参数 | 控制文件 schema 的"类型与默认值"实际写在数值模块里，契约跨层 |
| `codes/CMakeLists.txt`、`tests/*/CMakeLists.txt` | 全部 | 单一 `OneFLOW` 可执行目标 + 递归源收集 + 全局编译定义（`codes/CMakeLists.txt:347-414`），测试目标靠手写源清单，已出现两次清单漂移事故 |

## 5. 现有测试（源码中存在这些用例）

合计 `tests/main` 122 个 `TEST`、`tests/io` 26 个 `TEST`（`grep -c "^TEST("`）。

| 层 | 目标 / 文件 | 覆盖的契约 | 明确**不**覆盖 |
|---|---|---|---|
| 纯解析 | `ConfigLoader_test`（5）、`TextFileParser_test`（7） | 标量类型/顺序、跨行数组、注入式命名长度、缺省 resolver 报错、解析新文件替换旧 `ConfigDocument`；注释行、`funcMap/fileMap` 块、重复计数简写 | 真实 case 目录、全局 `DataBase` 提交 |
| 组件 | `TextFileParser_test` 内 `ControlManifestReader.PreservesListedScriptOrderAndCaseRoot` | 清单**追加**语义（预置元素保留）+ 顺序 + case 根解析 | 缺清单文件、空清单（`maxCases=0` 的 `control.txt`）的失败语义 |
| 组件 | `control_configuration_test` 内 `ControlParameterResolverTest`（3）、`ConfigDatabaseAdapterTest`（2） | 字面量/命名数组长度、缺名报错；类型化提交、后声明覆盖先前值 | `DataBase` 未初始化、重复类型冲突之外的提交失败 |
| 组件 | `simu_context_parameter_test`（20） | 缺失/合法/类型错误的 8 个根键 + `sourceGridCaseDir` 回退矩阵；`ApplyCaseIOSettings` 写到 `results/io_paths.txt`；显式 `ioPathInfoFileName` 写 mesh/flow 两个输出根 | `ApplyCaseIOSettings` 的 mesh/flow 同路径去重分支（`SimuContextIO.cpp:65-72`）；失败清理 |
| 集成（摄取） | `control_configuration_ingestion_test` 的 4 个 `ControlConfigurationIngestion` 用例（`:172,191,218,250`） | 缺脚本保留前序文件提交；解析失败不提交当前文件前缀；**adapter 提交非原子**可留下当前文件的更早条目；真实清单→解析→全局 `DataBase` 类型化提交 | 整清单回滚（测试刻意固定"无回滚"现状） |
| 集成（读取编排） | 同上 `ControlReaderIntegration.ReadsCaseScriptsAndWritesDatabaseDump`（`:406`） | 真实 `ReadControlInfo(caseDir)`：server 侧读取、`log/database.log` 内容、`log/log0.log` 生成 | 多 rank 广播与单写者（serial 下 `HXBcast` 在 `nProc<=1` 直接 return，`Parallel.cpp:196`） |
| 集成（case 装配） | `ControlReaderIntegration.CaseBootstrapAppliesSettingsAndRebindsNextCase`（`:293`） | `InitializeCaseConfiguration` 应用 7 个根 + `io_paths.txt` 落在 mesh/flow 输出根；第二个 case 不继承前一个 case 的根 | **未调用 `TeardownCase`**：用 `logFile.ClearCaseDir(); Prj::ClearPrjBaseDir(); ClearRuntimeParameters();` 手工"镜像"了 19 步清理里的 3 步（注释 `:340-341` 自述 "Mirror"） |
| 集成（就绪位） | `ContextSetupPublishesReadinessAfterBootstrap`（`:363`）、`FailedContextSetupDoesNotPublishReadiness`（`:388`） | `SimuContext::SetupCaseEnvironment` 成功后才 `IsEnvironmentReady()`；缺脚本时抛异常且仍为 false | `SimuImp::RunCase` 的 `ExecuteWithCleanup` 边界、`TeardownCase` 是否真的执行 |
| 生命周期辅助 | `simu_case_lifecycle_test`（17） | 头文件模板 `ExecuteWithCleanup/ExecuteAllCleanupSteps/ExecuteWithRollback/RunCasesWithProcessLifecycle` 的顺序、异常保留、rollback 只在 setup 失败时发生、重复路径重复执行、stop-on-failure | 任何真实 IO 副作用（纯 lambda） |
| 上下文 | `simu_context_test`（16） | 并行参数校验、任务名映射、`accelStates_` 拥有与清除、move 语义、`EnsureExpandedSolverNames` | `SetupCaseEnvironment`（该目标不链接 `SimuContextCaseEnvironment.cpp`） |
| CLI / 路径 | `prj_cmdline_test`（43） | argv 解析、`ResolveCaseDir` 相对/绝对、`OpenPrjFile/OpenCaseFile` 的 mode 相关根、`PrjIOPaths` 的继承/独立/诊断/`DomainRootsResetWhenBindingNextCase`（`:880`）、相对 `projectRoot` 取自 `current_dir`（`:952`） | `Prj::Init()`/`FindSystemRoot()`；`Prj::OpenFile` 失败路径（`Fatal` → `runtime_error`）在 case 上下文中的传播 |
| 任务分派 | `task_registry_test`（9）、`simu_base_interface_test`（3） | 注册/创建/未知返回 nullptr、`TaskEnum`↔字符串、`SolveFieldTask` 的前置条件（`RequireSolveFieldContext`） | 真实任务体、`SimuTaskReg.cpp` 的字符串表（测试用 `FakeTask`，不校验 `SimuTaskReg.cpp` 的表与 `TaskFilter` 一致） |

## 6. 缺失的契约测试

按"能否用现有测试目标链接"排序，M1/M2 需要新的链接边界（见 §11 RFC）。

- **M1 `TeardownCase` 生产清理契约**。现状：`SimuContextEnv.cpp` **不被任何测试目标编译**
  （`tests/io/CMakeLists.txt:59-94` 只列 `SimuContextCaseEnvironment.cpp`、
  `SimuContextBootstrap.cpp`、`SimuContextIO.cpp`），因此 `TeardownCase` 既未链接也未覆盖。
  应固定：成功装配后调用 `TeardownCase()` 得到
  `Prj::prjBaseDir/GetMeshOutputBaseDir()…` 为空、`DataPara/DataField` 为空、
  `logFile` 不再持句柄、`envReady_==false`、`task_==SOLVE_FIELD`、`TaskName()=="Solve"`；
  且 `Prj::current_dir`/`system_root`（进程级）**不被**清空。
- **M2 `SimuImp::RunCase` 异常边界**。应固定：`SetupCaseEnvironment` 抛异常时
  `PostProcess()` 仍执行（即 `TeardownCase` 被调用），原异常类型/消息保留，
  且进程级并行运行时**未** finalize（以便下一个 case 继续）。现有
  `simu_case_lifecycle_test.SetupFailureRunsRollbackAndPreservesFailure` 只固定了模板行为。
- **M3 多 rank 文件所有权与广播**。`mpirun -np 2`（本地允许）应固定：
  非 server rank 从广播得到参数；`log/database.log` 只由 server 打开/截断；
  `results/io_paths.txt` **单写者**——该断言在当前实现下预计会失败（§7 R-1），
  因此应先按"期望契约"写测试并标记，再由人工决定修复方式。
- **M4 `ioPathInfoFileName` 去重分支**。`SimuContextIO.cpp:65-72` 的
  mesh/flow 输出根相同时只写一次，当前无断言（现有两个 `SimuContextIOTest`
  用的都是不同根）。
- **M5 进程入口路径契约**。`Prj::Init()` 只被 `Simulation::RunImpl` 调用，
  无测试；`FindSystemRoot()` 找不到 `system/` 时 `Fatal`（`Prj.cpp:189-193`），
  而 `Prj::ResolveCaseDir` 在 `current_dir` 为空时 `Fatal`（`Prj.cpp:205-212`）。
  现有测试靠手工 `Prj::current_dir = std::filesystem::current_path().string();`
  绕开（`prj_cmdline_test.cpp:891`、`control_configuration_ingestion_test.cpp:325`），
  该绕开本身是"生产初始化缺失"的证据。
- **M6 兼容 API 的存废契约**。`ReadControlInfo()`/`DumpDataBase()`/`ReadPrjScript()`/
  `SimuContext::SetupEnvironment/TeardownEnvironment`/`SimuImp::Run()` 均无调用者。
  要么删除（需 human 批准，属跨模块），要么补一条最小契约测试；不能长期保持
  "公开但无验证"的状态。
- **M7 清单/脚本失败语义的边界补全**：缺 `script/control.txt` 本身（现有用例缺的是
  清单里**列出的**脚本）、空清单、同一脚本重复列出（重复提交语义）。

## 7. 风险

- **R-1（高）`io_paths.txt` 在多 rank 下无写者过滤。** `ApplyCaseIOSettings` 在
  **所有 rank** 上执行（`SetupCaseEnvironment` 无 rank 判断，
  `SimuContextCaseEnvironment.cpp:28`），却经 `Prj::WriteIOPathInfo` →
  `Prj::OpenPrjFile`（`Prj.cpp:437`）以 `out|trunc` 打开同一路径
  （`SimuContextIO.cpp:63,70`），没有 `PIO::OpenPrjFile` 的
  `pid != GetFid()` 过滤（`PIOFileAccess.cpp:39-47`）。这与 `docs/AI_CONTEXT.md`
  记录的 `database.log` 修复是**同一类缺陷**，`database.log` 已改用
  `PIO::OpenCaseFile`（`ControlDatabaseIO.cpp:50`），`io_paths.txt` 尚未。
  并发 `trunc` + 并发写相同内容，实际后果是内容交错/半写；且
  `Prj::CloseOutputFile` 在 `WriteIOPathInfo` 失败时**抛出**
  （`Prj.cpp:482-495`），单 rank 失败会让该 rank 的 case 装配失败。
  目录并发创建本身是安全的（`FileUtils.cpp` 的 `HX_CreateDirectory` 已处理竞态）。
- **R-2（高）`TeardownCase` 的 19 步清理无测试。** 集成测试用手工 3 步"镜像"
  代替（`:340-344`），生产清理一旦漏项或多清，测试不会发现；而
  "第二个 case 不继承第一个 case"这一结论正依赖这些步骤。
- **R-3（中）测试与生产的源清单重复。** `control_configuration_test`（`tests/io/CMakeLists.txt:59-94`）
  手写 30 个源文件，等于把生产依赖闭包复制了一份；历史上已两次漂移
  （`control_configuration_test` 源文件名错误、`OneFLOW_GeometryCore` 缺
  `PIOFileAccess.cpp`）。新增 `TeardownCase` 覆盖会再次需要 `SimuContextEnv.cpp`
  的完整闭包（Parallel、AccelRuntime、SolverMap、heat_flux、nscom、NsSolver、TurbSolver、
  Zone、FieldManager…），复制清单的做法不可持续。
- **R-4（中）全局 `DataBase` 的测试串扰。** `control_configuration_test` 把
  resolver/adapter 用例（`control_parameter_resolver_test.cpp:39` 用
  `SetDataInt`、`:53` 写 `dimension`/`global_dt`/`simutask`）与摄取/集成用例编进
  同一个可执行文件，只有摄取用例显式 `ClearRuntimeParameters()`；
  adapter 用例留下的键会跨用例存活，断言依赖"唯一名字"这一约定而非隔离。
- **R-5（中）`LogFile` 的自由函数路径有进程级 static `ifReWrite`**
  （`LogFile.cpp:34`），只有 `logFile.caseDir_` 为空时才走；该 static 在 case 之间
  不复位，进程级日志从第二次写开始变为 `app`。当前无测试覆盖该分支。
- **R-6（中）`DumpDataBase()`（无参）写 `"/log/database.log"`
  （`ControlDatabaseIO.cpp:40`）。** 该字符串在 POSIX 上
  `std::filesystem::path::is_absolute()` 为真，`ResolvePathUnderRoot`
  （`Prj.cpp:71-86`）因此**绕过配置的输出根**，直接落到文件系统 `/log/`；
  在 Windows 上则不视为绝对路径。该函数只被无调用者的 `ReadControlInfo()` 触发，
  属"死 API + 平台相关路径语义"，但仍是可被后续代码踩到的陷阱。
- **R-7（中）`SetupCaseEnvironment` 无回滚**（文档已声明，`runtime_architecture.rst`
  §Case lifecycle）。直接调用它的测试只能断言"未发布就绪"，实际清理依赖
  `SimuImp::RunCase`/`SetupEnvironment`；M2 未覆盖前，"失败后状态干净"没有证据。
- **R-8（低）临时目录固定路径。** `simu_context_parameter_test` 用固定名
  （`OneFLOW_CaseIOSettings_Default` 等），`control_manifest_reader_test:14` 亦然；
  与 `TemporaryControlCase` 的"时间戳+序号"策略不一致，并行/残留时会互相污染。
- **R-9（低）`ReadScriptFileNameList` 把清单中的空行/无词行视为结束**
  （`ControlManifestReader.cpp:66-68` 用 `if (!ReadNextNonEmptyLine()) break;`），
  中间出现空行之后的条目会被静默丢弃；现有测试的空行恰好在末尾之后/中间且能继续读取，
  未固定"空行之后的条目是否必须被读取"。

## 8. 当前 R1 gap（与文档的差异）

`docs/AI_CONTEXT.md` 与 `docs/plans/setup-case-environment-integration-test.md` 的
最新记录（"Implementation update: exercise the production context setup method"）
已经把 R1 分解为"读 → 装配 → 就绪位"，并指出还剩三项：生产 `TeardownCase()`、
`SimuImp::RunCase()` 异常边界、进程运行时保持。

对照代码，R1 的现状应更新为：

- **已被测试固定**（源码级证据，未在本机运行）：真实 `ReadControlInfo(caseDir)` 的
  serial 编排与 `log/database.log`；真实 manifest→解析→`DataBase` 提交及覆盖顺序；
  三种摄取失败语义；`InitializeCaseConfiguration` 的根应用与跨 case 重绑定；
  `SetupCaseEnvironment` 成功/失败的就绪位；`Prj` 根在重绑定时复位。
- **仍然缺失**：`TeardownCase` 生产清理契约（M1）、`SimuImp::RunCase` 异常边界（M2）、
  多 rank 广播与单写者（M3，并带出 R-1 缺陷）、`io_paths.txt` 去重分支（M4）、
  进程入口 `Prj::Init`（M5）。
- `docs/AI_CONTEXT.md` 的 "Remaining R1 gap" 段落（"A real integration test still needs
  to exercise … effective `Prj` root application, `io_paths.txt`/database diagnostics,
  sequential-case reset"）**已部分过时**：其中的 effective root application、
  `io_paths.txt` 诊断、跨 case 复位现在都有断言
  （`ControlReaderIntegration.CaseBootstrapAppliesSettingsAndRebindsNextCase`）。
  但该文件"Next single action"仍停留在实现这些测试上，与最新 plan 记录不一致。
  本报告不修改 `docs/AI_CONTEXT.md`（不在 audit scope 内），仅指出漂移，供人工在
  下一个工作切片中一并修正。

## 9. 可并行的未来 scope 提案

约束：`Prj` 有 41 个依赖文件，`codes/CMakeLists.txt` 是单一可执行目标，测试目标是手写
源清单；因此"两个 worker 各改一半 `codes/io`"必然冲突。下面三个切分满足
"文件集合互斥"，但 P2 建议**冻结**（见 §10 的问题）。

| 提案 scope | 独占文件（建议写入 `scope/<name>.txt`） | 依赖的冻结契约 | 首个可验证切片 |
|---|---|---|---|
| **P1 `case-config-ingestion`** | `codes/io/src/{ControlManifestReader,ControlScriptReader,ControlScriptParser,ControlFileReader,ControlDatabaseIO,ConfigLoader,ConfigDatabaseAdapter,ControlParameterResolver,LegacyParameterSyntax,TextFileParser,CommentLine}.cpp`、对应 `codes/io/include/{ParaFile,ConfigLoader,ConfigDocument,ConfigDatabaseAdapter,TextFileParser,CommentLine,LegacyParameterSyntax}.h`、`codes/main/src/{SimuContextBootstrap,SimuContextIO,SimuContextCaseEnvironment}.cpp`、`codes/main/include/{SimuContextBootstrap,SimuContextIO,SimuContextParameter}.h`、`tests/io/*` | 只读使用 `Prj` 的 `Set*BaseDir`/`WriteIOPathInfo`/`GetCaseFileName`，不改签名；只读使用 `DataBase` | M4（去重分支）+ M7（清单失败语义）；R-1 的修复方案需 RFC |
| **P2 `path-file-io`**（建议冻结，不作为 worker scope） | `codes/io/src/{Prj,FileUtils,FileO,PIO,PIOFileAccess,LogFile,StringUtils}.cpp`、`codes/io/include/{Prj,FileUtils,FileO,PIO,LogFile,StringUtils}.h`、`tests/main/prj_cmdline_test.cpp` | 本身是冻结契约 | 仅人工/共享基础维护；任何签名变更走 RFC |
| **P3 `runtime-lifecycle-dispatch`** | `codes/main/src/{Simulation,SimuImp,SimuContext,SimuContextEnv,SimuTask,SimuTaskReg,SimuBase}.cpp`、`codes/main/include/{Simulation,SimuImp,SimuContext,SimuContextEnv? 无,SimuTask,SimuTaskRequire,SimuCaseLifecycle,SimuProcessLifecycle,SimuBase}.h`、`tests/main/{simu_case_lifecycle,simu_context,simu_base_interface,simulation_command_line,task_registry}_test.cpp` | `SimuContext.h` 的构造/移动/就绪位语义；`SimuTaskReg.cpp` 的任务名表 | M1（`TeardownCase`）+ M2（`SimuImp::RunCase` 异常边界） |
| **P4 `result-restart-io`**（推迟） | `codes/restart/**`、`codes/unsteady`/`ins`/`ns`/`uturb` 的 restart 文件、`codes/visual/src/HeatFluxTask.cpp` | 需要先有结果/重启动文件的路径与所有权契约 | 不在 v0.1/v0.2 起步范围 |

`SimuContext.h` 是被 P1 与 P3 共享的头文件，必须在 scope 文本里写死归属：
建议 **P1 只拥有 `SimuContextCaseEnvironment.cpp` 与 `SimuContextParameter.h`**，
`SimuContext.h`/`SimuContextEnv.cpp` 归 P3；P1 若需改 `SimuContext` 成员，走 RFC。

## 10. 需要人确认的一个具体问题

**问题**：缺失的 M1（`TeardownCase` 生产清理）与 M2（`SimuImp::RunCase` 异常边界）
应采用哪种验证策略——
**(A)** 新建一个由人工维护的共享生产运行时库目标（暂名 `oneflow_case_runtime`），
让测试直接链接并调用 `SimuContext::TeardownCase` / `SimuImp::RunCase`；
还是 **(B)** 只用真实 `OneFLOW` 可执行 + `test/` 算例做端到端断言（不改构建系统）？

- 建议：**选 (A)**，并把库范围限定为 `codes/io` 的摄取源 + `codes/main` 的
  `SimuContextBootstrap`/`SimuContextIO`/`SimuContextCaseEnvironment`/`SimuContextEnv`；
  由人工（共享基础所有者）改 `codes/CMakeLists.txt` 与 `tests/io/CMakeLists.txt`，
  worker 只写测试与契约文本。
- 理由：`TeardownCase` 位于 `codes/main/src/SimuContextEnv.cpp:100`，其对象文件引用
  `Parallel`、`AccelRuntime`、`SolverMap`、`heat_flux`、`nscom`、`NsSolver`、`INsSolver`、
  `TurbSolver`、`turbcom`、`Tolerance`、`Zone`、`ActionState`、`Iteration`、`FieldManager`
  等十余个模块的符号，**整个 TU 必须可链接**才能调用其中任一函数。现状是这份闭包
  被排除在所有测试目标之外（`tests/io/CMakeLists.txt:59-94` 不列 `SimuContextEnv.cpp`），
  所以清理路径无法测试；若靠复制源清单补上（R-3），已经两次造成构建事故。
  (B) 不改构建，但必须跑真实任务（`test/grid/rae2822convert` 的 `simutask=Grid`
  会真的生成网格，`test/cfd/rae2822_roe_sa` 是 `Solve`，不能用于 setup-only 断言），
  失去快速、确定性的失败定位，且无法单独验证"清理后再失败"的异常边界。
- (B) 仍应保留给 M3（多 rank 广播/单写者），因为那类契约必须有真实进程。
- 具体目标定义、源分组与替代方案见
  `.agentlab/RFC/runtime-io-2026-11-01-shared-case-runtime-target.md`
  （文件名遵循 `.agentlab/scope/runtime-io.txt` 的 `runtime-io-*` 前缀，否则 `scope_check` 会判越界；
  RFC 内部记录了日期与主题）。

## 11. 与本次审计无关但需人工记录的事项

- `docs/AI_CONTEXT.md` 的 R1 状态段落已滞后于 `tests/io/control_configuration_ingestion_test.cpp`
  的 `ControlReaderIntegration` 四个用例（见 §8）。修正该文档不在 audit scope 内。
- 本报告 §5 的“现有测试”是源码中存在性清单；**未在本机运行**，不得引用为本机通过证据。
  用户侧最近一次报告的本地基线（`docs/AI_CONTEXT.md` 记录为 987/987）属**用户报告**，
  与本审计基线 `7328dc504` 不是同一 revision。

## 12. 交付记录

- 本报告所在提交：`ff5a34c6e766c62c369e2d3342e4d152bdd6e4f0`（branch `agent/runtime-io`，未 amend、未 push）。
  同一提交包含本报告、RFC 与 gate 生成的 `.agentlab/repro/runtime-io-20261101T110035Z.json`。
- Gate：`GATE_MODE=audit …/gate.sh runtime-io` → 结果见下方 gate 输出。
- 未 push、未 merge、未提交集群作业。

## Gate result (tail of gate.sh)

```
$ GATE_MODE=audit /Users/luql/Downloads/agent/repo/open-research-agents/platform/scripts/gate.sh runtime-io
== scope
scope: runtime-io OK (2 changed files)
repro record: .agentlab/repro/runtime-io-20261101T110035Z.json
----
GATE PASS (mode: audit, read-only)
```

（audit 模式只做 `scope_check`，不执行 build/ctest/verification。gate 退出码 0。）

## Files changed

- `.agentlab/reports/runtime-io.md`（本文件）
- `.agentlab/RFC/runtime-io-2026-11-01-shared-case-runtime-target.md`
- `.agentlab/repro/runtime-io-*.json`（gate 生成的快照）

## Repro record (required for every accuracy/performance conclusion)

本审计不产生精度/性能结论，因此不写 `*.case.json`。gate 生成的
`.agentlab/repro/runtime-io-*.json` 只记录 revision、环境与 gate 结果：

- commit / dirty: 审计基线 `7328dc504342f6f8027065e9b8ac2894234e5878`；gate 快照记录时工作树 dirty（本报告 + RFC 未提交）
- machine and hardware: macOS 本地检出（无 MPI/METIS/CGNS/HIP 配置记录）
- compiler and flags: 未编译（audit 模式禁止 build/ctest）
- commands and environment: 仅只读 `git`/`grep`/`sed` 与 `GATE_MODE=audit gate.sh runtime-io`
- input data/version: 审计基线即上表 revision；未使用任何网格或算例
- result files under `.agentlab/repro/`: `runtime-io-<UTC>.json`
- domain fields（`.agentlab/repro_fields.yaml`）：不适用（无数值运行）

## Known limitations

- 未编译任何目标，未运行 GoogleTest/CTest，未运行 `mpirun`，未查看 CI 结论。
  §5 的覆盖判断均来自测试源码阅读，可能与"实际能否编译/通过"存在差异。
- 未审计 `codes/ns` 等数值模块内部；`spatial`/`temporal` 的职责划分属 `numerics`
  workstream（A2）。§3 只给出 runtime-io 视角下可用的三条边界。
- `docs/plans/setup-case-environment-integration-test.md` 的完整 1183 行只按章节标题与
  末尾进展段落阅读，未逐段核对历史细节。
