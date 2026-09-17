# OneFLOW hexin 主路径重构 — Handoff

| 项 | 内容 |
|---|---|
| Branch | `hexin` |
| Date | 2026-09-15 |
| Scope | 主程序编排 / 任务分发 / solver 命名 / Cmx 任务名入口 |
| Principles | 不改初始化与求解数值；可测纯逻辑优先；Linux（GCC）与 Windows（MSVC）均需通过 |

---

## 1. 已完成阶段

### 阶段 0 — `SimuImp : SimuBase`

- 仿真入口继承清晰基类。
- `Run` = Pre / Main / Post。

### 阶段 1 — `ISimuTask` + `TaskRegistry`

- 控制文件 `simutask` 字符串 → 注册表创建任务。
- 生产注册：`codes/main/src/SimuTaskReg.cpp`（`Solve` / `Grid` / `WallDist` / …）。
- 单测：`tests/main/task_registry_test.cpp`。

### 阶段 2 — `SimuContext`（逻辑与环境拆分）

- 上下文持有：`args`、rank/size、task 名与枚举、env 就绪标志。
- `SetupEnvironment` / `TeardownEnvironment` / `ResolveTaskFromControl`。
- 测试可注入：`SetTaskByName`、`MarkEnvironmentReady` 等。
- 单测：`tests/main/simu_context_test.cpp`。

### 阶段 2.1 — `Execute(const SimuContext&)`

- `SimuImp::RunSimu`：`Create(taskName)` → 可选 `ConstructSystemMap` → `Execute(*ctx_)`。
- `SolveFieldTask`：先 `RequireSolveFieldContext(ctx)`，再跑场求解管线。

---

## 2. 本轮追加（Solver / Cmx / Field 管线）

### 2.1 Solver 命名策略

| 项 | 说明 |
|---|---|
| Header | `codes/solver/include/SolverNamePolicy.h` |
| API | 纯逻辑：`MakeUnstructuredSolverName` / `MakeStructuredSolverName`、`Expand*`、`FillExpandedSolverNames` |
| Production | `SolverNameClass::ReadSolverNames()` → `FillExpandedSolverNames` |
| Parse loop | `ReadNextMeaningfulLine`（与 MessageMap 对齐：跳过注释、避免 EOF 空串） |
| Tests | `tests/solver/solver_name_policy_test.cpp`（含 DualLists / Empty） |
| GCC note | 模板内用到的 `Make*` 必须在模板定义**之前**可见（两阶段查找；MSVC 更宽松） |

### 2.2 Cmx 任务名

| 项 | 说明 |
|---|---|
| Header | `codes/task/include/CmxTaskNames.h`（唯一来源） |
| Examples | `kInitFlowFieldTaskName`、`kPostProcessTaskName`、MG 相关名等 |
| FieldSimu | Init → `MultiSolverMultiGridTask(kInitFlowFieldTaskName)` |
| Multigrid | 原字面量改为同一头文件常量 |
| Contract tests | `MessageMapImp` 对 `INIT_FLOWFIELD` / `POST_PROCESS` 做 name↔id 往返 |

### 2.3 Field 管线编排

```text
FieldSimuRunPipeline()
  SetupGlobals → LoadGrid → PrepareWallDist
  → CreateSolvers → InitFlowField → Run

SolveFieldTask::Execute  → Require* + FieldSimuRunPipeline()
FieldSimu()              → FieldSimuRunPipeline()
FieldSimuInitFlowField() → MultiSolverMultiGridTask(kInitFlowFieldTaskName)
InitializeSolver()       → FieldSimuInitFlowField()   // compatibility alias
```

### 2.4 任务名如何进入 CmxTask（语义未改）

```text
MultiSolverMultiGridTask(name)
  → for each solver × grid level:
       SingleSolverSingleGridTask(name)
         → operationId = MessageMap::GetMsgId(name)
         → GenerateCmdList(operationId)
         → CMD::ExecuteCmd()
```

字符串仅做 name↔id 映射，再按当前 `SolverState::solverType` 查注册实现；**不改变初始化数值**。

---

## 3. 关键约定（后续不要破坏）

1. **数值冻结**：Init / MG / residual 相关逻辑本轮未改算法，只改入口与命名。
2. **纯策略可单测**：`SolverNamePolicy`、`CmxTaskNames`、`MessageMapImp` 不依赖 MPI/网格。
3. **GCC 两阶段查找**：模板内用到的自由函数必须在模板定义前可见。
4. **双平台**：改头文件/模板后必须 Linux + Windows 都编过。
5. **回归**：接口改动后跑全量 gtest；数值内核另走 `test/` 与 kunshan 套件（见仓库根目录 `AGENTS.md`）。

---

## 4. 刻意未做（Backlog）

| 项 | 说明 |
|---|---|
| `solver.txt` 的 Prj 路径单测 | 生产用 `OpenPrjFile`；可选抽 `OpenFile` 版再测，当前略过 |
| `CreateSolvers` 内部再拆纯函数 | S1–S4 仅注释分层，未抽可注入列表的工厂 |
| `SimuContext` 携带 solver 名列表 | 仍读全局 `SolverNameClass` |
| Multigrid 以外 TU 的 Cmx 字面量 | 仅 FieldSimu + Multigrid 已收编 |
| 删除 `InitializeSolver` | 仍作兼容别名保留 |

---

## 5. 建议的下一刀

1. **轻**：全库搜索仍残留的 `"INIT_FLOWFIELD"` / `"POST_PROCESS"` 等，统一走 `CmxTaskNames.h`。
2. **中**：`CreateSolvers`：名字列表 → clone 循环的可测边界（mock `SafeClone` 或只测 index map）。
3. **重**：solver 列表进入 `SimuContext`，`CreateSolvers` 读 context 而非隐式全局。

---

## 6. 关键文件索引

| 路径 | 角色 |
|---|---|
| `codes/main/src/SimuImp.cpp` | 注册表分发 + context |
| `codes/main/src/SimuTaskReg.cpp` | 生产 `ISimuTask` |
| `codes/main/include/SimuContext.h` | 运行上下文 |
| `codes/global/src/FieldSimu.cpp` | 六阶段 + Init→Cmx |
| `codes/global/include/FieldSimu.h` | 阶段 API / `FieldSimuRunPipeline` |
| `codes/solver/include/SolverNamePolicy.h` | U/S 前缀策略 |
| `codes/solver/src/SolverMap.cpp` | 读表 + `CreateSolvers` |
| `codes/task/include/CmxTaskNames.h` | Cmx 操作名 |
| `codes/task/src/CmxTask.cpp` | Multi/Single Solver×Grid 入口 |
| `codes/multigrid/src/Multigrid.cpp` | MG 内 Cmx 调用 |
| `tests/main/task_registry_test.cpp` | 任务注册契约 |
| `tests/main/simu_context_test.cpp` | 上下文契约 |
| `tests/solver/solver_name_policy_test.cpp` | solver 命名策略 |
| `tests/register/MessageMapTest.cpp` | MessageMap + Cmx 名契约 |

---

## 7. 一句话状态

主路径已从「巨型 switch + 散落字符串」收成：**Context 解析任务 → Registry 执行 → Field 六阶段管线 → Policy 展开 solver 名 → CmxTaskNames 进入 MessageMap/Cmx**；数值行为保持不变，gtest 全绿。

---

## 8. 放置说明

建议提交路径（与 `AGENTS.md` 文档地图一致）：

```text
doc/reports/architecture/oneflow-hexin-main-path-handoff-20260915.md
```

若仓库已有 `doc/reports/README.md` 的索引表，请在同一变更中追加本文件条目。

## 9. 续：2026-09-16

### 测试边界
- `SolverNameList` 从 `SolverMap` 物理拆出；`solver_name_policy_test` 仅依赖 Policy 头。
- `solver_name_list_test` 测 `LoadFromBaseNames` / `Reset`（轻依赖）。

### SolverMap
- `SolverBucket` / `BuildSolversInBucket`。
- `CreateSolvers(gridType, const StringField*)`：非空=已展开注册名；`nullptr`=默认 `GetSolverNames`。

### SimuContext
- `SetExpandedSolverNames` / `HasExpandedSolverNames` / `ClearExpandedSolverNames`。
- `FieldSimuRunPipeline(ctx)` / `FieldSimuCreateSolvers(ctx)`；`SolveFieldTask` 走带 ctx 管线。
- 生产默认仍不注入 → 行为与改前一致。

### CmxTaskNames
- 覆盖 FieldSimu、Multigrid、TimeIntegral、SolverState、MultiBlock。
- 生产 `Single/MultiSolver*Task("...")` 字面量已清扫。

## 10. 续：2026-09-18 — AddCmdToList 字面量收编

### 目标
将仍散落在 Restart / SolverImp / Ns / INs / Turb 中的 `AddCmdToList("...")` 字面量统一到 `CmxTaskNames.h`，与既有 Multi/SingleSolver 路径同一来源。

### 变更
| 文件 | 说明 |
|---|---|
| `codes/task/include/CmxTaskNames.h` | 新增 Restart 初始化、Interface 交换、Dump/Visual/Unsteady 常量 |
| `codes/restart/src/RestartTaskReg.cpp` | `InitFlowField` 使用 `kInitFirst*` / `kReadRestart*` 等 |
| `codes/solver/src/SolverImp.cpp` | `CommInterfaceData` 使用 upload/update/download 常量 |
| `codes/ns/src/NsSolverImp.cpp` | `NsPostprocess` / `NsFinalPostprocess` |
| `codes/ins/src/INsSolverImp.cpp` | `INsPostprocess` / `INsFinalPostprocess` |
| `codes/turb/src/TurbSolverImp.cpp` | `TurbPostprocess` / `TurbFinalPostprocess` |

### 新增常量（节选）
- Restart: `kInitFirstTaskName`, `kInitRestartTaskName`, `kReadRestartTaskName`, `kInitInsRestartTaskName`, `kReadInsRestartTaskName`, `kInitFinalTaskName`
- Interface: `kUploadInterfaceDataTaskName`, `kUpdateInterfaceDataTaskName`, `kDownloadInterfaceDataTaskName`
- Dump/Visual: `kDumpResidualTaskName`, `kDumpAerodynamicTaskName`, `kDumpPressureCoeffTaskName`, `kDumpHeatfluxCoeffTaskName`, `kDumpRestartTaskName`, `kDumpLaminarPlateTaskName`, `kDumpTurbPlateTaskName`, `kVisualizationTaskName`, `kUpdateUnsteadyFlowTaskName`

### 语义
纯字符串替换 → 注册表 / MessageMap 行为不变；无数值路径改动。

### 建议后续
1. 全库再扫一次 `AddCmdToList("` 确认无遗漏（当前生产路径已清）。
2. MessageMap 契约测试可按需追加新名 ↔ id 往返（与现有 INIT_FLOWFIELD 同模式）。
3. 中/重项仍按原 backlog：CreateSolvers 可测边界深化、solver 列表正式进 SimuContext 生产路径。
