# Report: numerics

status: ready
updated: 2026-11-01
commit: e936b8191（审计交付提交；§12 记录回填与复核提交）
branch: agent/numerics
model: audit worker (agent-lab / pi)
pi_version: n/a
gate_mode: audit
gate_status: PASS

任务：`.agentlab/TASKS.md` A2 —— OneFLOW numerics 边界只读审计。
基线提交：`df8676fc2`（`agent/numerics` == `dev`，`git diff --name-only dev...HEAD` 为空）。

---

## 0. 结论摘要（三句话）

1. OneFLOW 的 numerics **不是 spatial/temporal 二分**，而是四层：**(a) 纯函数 policy 层**（`*Policy.h/.cpp`，唯一有契约测试的层）、**(b) 逐单元/逐面 kernel 类**（`NsInvFlux`、`NsLusgs`、`NsUpdate`…，读写全局 `inv`/`gcom`/`nscom`）、**(c) 驱动/装配层**（`UNsInvFlux`、`UNsLusgs`、`UNsUpdate`、`UGrad` 的循环 + 设备分派）、**(d) 任务名调度层**（`TimeIntegral` + `TaskRegister`/`MessageMap`）。真正的"空间算子"和"时间算子"是 (c) 里的两个**不同模块**、却共用一个任务序列。
2. 全部 5 个 CPU 回归算例都走 **LU-SGS 隐式（`time_integral=2`）**，**显式多阶段 RK 在 `test/` 里没有任何残差基线**；`tests/` 里**没有任何测试实例化** `UNsInvFlux`/`UNsGrad`/`UNsUpdate`/`UNsLusgs` 中任一实际循环。因此目前"有测试"只等于"policy 纯函数有测试"。
3. 时间推进的**迭代驱动程序不在 numerics 目录里**：外层/内层循环由 `codes/multigrid/src/Multigrid.cpp` 的 `MG::Run()` 拥有，`codes/time/` 只是墙钟工具（`TimeSpan`/`TimeTest`）。据此提出 5 个文件级互斥的未来 scope（§10），并把 OneFLOW 专用契约写成 RFC（`.agentlab/RFC/numerics-2026-11-01-boundary-contracts.md`）。

---

## 1. 方法与证据边界

- 本审计为 **read-only**：没有编辑求解器代码/CMake/contracts/verification，**没有构建、没有运行任何 gtest/ctest/算例**。
- 因此本报告**不对任何测试是否通过作出声明**。文中出现的 `cell`/`line` 引用均来自源码阅读。
- 未复用 starter 的 `spatial`/`temporal` 划分；边界全部由实际调用链推导（§3、§4）。
- 用户/文档记录的基线仅作为"改前基线"引用，不作为本轮证据：`docs/AI_CONTEXT.md` 记录 982/982（用户报告）；`docs/plans/oneflow-development-todo.md` §0 记录 CPU 五算例 normal 5/5 + strict 5/5（上一轮证据）。
- 审计范围以 `.agentlab/TASKS.md` A2 列出的目录为主；实际调用链把它强制扩到 `codes/uns`、`codes/usolver`、`codes/multigrid`、`codes/interface`、`codes/accel`，因为没有这些目录就无法描述 numerics 的真实入口。这一点本身就是主要发现之一。

---

## 2. 入口与调用者

### 2.1 生产调用链（稳态 LU-SGS，即当前 5 个回归算例的路径）

```
FieldPipeline::Run(ctx)                     codes/global/src/FieldSimu.cpp
├─ FieldSimuSetupGlobals()                  → InitFlowSimuGlobal()
│    ├─ vis_model.Init(); ctrl.Init();      codes/global/src/Ctrl.cpp
│    ├─ Iteration::Init()
│    └─ timeIntegration.Init()              codes/unsteady/src/TimeIntegration.cpp
├─ FieldSimuLoadGrid(); FieldSimuPrepareWallDist()
├─ FieldSimuCreateSolvers(ctx)              → SolverCatalog::CreateDefault
├─ FieldSimuInitFlowField()                 → "INIT_FLOWFIELD" 任务
└─ FieldSimuRun(ctx)                        → MultigridSolve(ctx)
     └─ MG::MultigridSolve                    codes/multigrid/src/Multigrid.cpp:104
        ├─ MG::Allocate()
        ├─ MG::Run()                          ← 外层/内层迭代循环在这里，不在 numerics 目录
        │   while (SimuIterState::Running())          // 外层 = 物理/伪时间步
        │     while (!SolverState::Converge())        // 内层；稳态下恒为 1 次（见 §5.4）
        │       MG::SolveInnerIter()
        │         ├─ MG::PreprocessMultigridFlowField()      // "STORE_RHS"
        │         ├─ MG::PreRelaxationCycle(gl)
        │         │    └─ TimeIntegral::Relaxation(MG::nPre[,ctx])   codes/solver/src/TimeIntegral.cpp:243
        │         ├─ MG::MultiGridCycle(gl)                  // 粗网格
        │         └─ MG::PostRelaxationCycle(gl)
        │              └─ TimeIntegral::Relaxation(MG::nPost[,ctx])
        │                   └─ TimeIntegral::timeIntegral()          // 函数指针
        │                        ├─ RungeKutta() / RungeKutta(ctx)   // MULTI_STAGE == 1
        │                        ├─ Lusgs()                          // 其它（含 2）
        │                        └─ Simple()                         // SIMPLE == 3
        │                           └─ SingleSolverSingleGridTask("<TASK_NAME>", ctx)
        └─ MG::Deallocate()
```

关键点：`TimeIntegral::Relaxation()` 的**唯一生产调用者**是 `codes/multigrid/src/Multigrid.cpp:294/298/308/312`。即使不做多重网格，迭代驱动也由 multigrid 模块承担。

### 2.2 任务名 → 数值实现（真实"接口"）

`TimeIntegral::Lusgs()` 不发号施令，而是按字符串任务名派发（`codes/task/include/CmxTaskNames.h`）：

| 任务名 | 注册处 | 落地实现 |
|---|---|---|
| `LOAD_Q` | `codes/solver/src/FieldTaskReg.cpp` | `LoadQ`（受 `DeviceOwnsExplicitRkBookkeeping()` 门控） |
| `CALC_TIME_STEP` | `codes/ns/src/NsSolverImp.cpp` → `REGISTER_DATA_CLASS(NsCalcTimeStep)` | `UTimeStep::CalcTimeStep`（`codes/uns/src/UTimeStep.cpp`） |
| `UPDATE_RESIDUALS` | 同上 → `NsUpdateResiduals` | `NsRhs::UpdateResiduals` → `NsCalcRHS`（`codes/ns/src/NsRhs.cpp`） |
| `LOAD_RESIDUALS` | `FieldTaskReg.cpp` → `LoadResiduals` | `NegField(res, FIELD_RHS)`，被设备链短路 |
| `ZERO_DQ_FIELD` / `INIT_LUSGS` / `LUSGS_LOWER_SWEEP` / `EXCHANGE_INTERFACE_DQ` / `LUSGS_UPPER_SWEEP` / `UPDATE_FLOWFIELD_LUSGS` | `codes/implicit/src/ImplicitTaskReg.cpp`、`codes/interface/src/InterfaceTaskReg.cpp` | `LusgsState::GetLusgsSolver()` → `UNsLusgs` / `NsLusgs`；`EXCHANGE_INTERFACE_DQ` 由 `codes/interface` 拥有 |
| `CALC_LHS` | `codes/lhs/src/LhsTaskReg.cpp` | `CreateLhs()` → `ULhs::CalcLHS`（`codes/usolver/src/ULhs.cpp`） |
| `UPDATE_FLOWFIELD` | `codes/update/src/UpdateTaskReg.cpp` | `CreateUpdate()` → `UNsUpdate::UpdateFlowField` |
| `CALC_BOUNDARY` | `codes/ns/src/NsSolverImp.cpp` → `NsCalcBoundary` | `NsCalcGamaT` + `NsCalcBc`(`UNsBcSolver`) 或 `CalcHipBoundary()` |
| `POST_PROCESS` | `codes/global` → `NsPostprocess` | `CommInterfaceData()` + `AddCmdToList(DUMP_RESIDUAL/AERO/VISUAL/RESTART)` |

也就是说：**numerics 的稳定接口是"任务名字符串 + 全局 DataBase 字段名"**，不是 C++ 抽象类。这解释了为什么 `.agentlab/contracts/interfaces.md` 的 `ISpatialOperator`/`ITimeIntegrator` 在这里找不到对应物。

### 2.3 空间算子实际入口（面通量）

```
NsCalcInvFlux()                                  codes/ns/src/NsRhs.cpp:186
└─ static UNsInvFlux uNsInvFlux;                 // 函数级 static，跨 RK 阶段复用
   uNsInvFlux.CalcFlux()                         codes/uns/src/UNsInvFlux.cpp:458
      ├─ inv.Init(); ug.Init(); unsf.Init()      // 全局量重绑定
      ├─ Alloc(); SetPointer(nscom.ischeme)      → &NsInvFlux::Roe / Ausmp / Slau2 / ...
      ├─ CalcInvFace()                           // gradient → limiter → reconstruction → BC 修正
      ├─ CalcInvFlux()|CalcInvFluxCpuBatch()|CalcAndAddInvFluxHipBatch()
      └─ AddInvFlux() → AddF2CField(res, invflux) // 有向沉积
```

逐面契约（隐式，靠全局量传递）在 `UNsInvFlux::PrepareFaceValue()`：

- 输入：`ug.fId`、`ug.lc`、`ug.rc`；`gcom.xfn/yfn/zfn/vfn/farea`；`inv.prim1[iEqu]`、`inv.prim2[iEqu]`（来自 `limiter->GetLeftField()/GetRightField()`，面索引）；`inv.gama1/2`（来自 `unsf.gama` 的单元值）。
- 输出：`inv.flux[0..nEqu)` = **未乘面积**的面通量（`NsInv::Init()` 只 `resize`，不 assign，等价零初始化）。
- 面积乘与符号由 `UNsInvFlux::UpdateFaceInvFlux()`（`= gcom.farea * inv.flux`）与 `AddF2CField`（`DepositOrientedFaceFlux`：left `-=`、right `+=`，boundary 面只 left）完成。

同一件事在 `FluxBackend` 路径上语义不同：`CpuFluxBackend::CalcInvFlux` **内部已乘 `state.faceArea`**（`codes/accel/src/CpuFluxBackend.cpp:304`）。即"面通量"这个量在两条实现里量纲不同（有/无面积）。

---

## 3. 责任划分：spatial / temporal / shared（由代码推导）

### 3.1 Spatial（离散在空间上，与时间推进无关）

| 责任 | 实现文件 | 备注 |
|---|---|---|
| 梯度（Green-Gauss + 反距离修正） | `codes/usolver/src/UGrad.cpp`、`codes/grad/src/Grad.cpp`、`codes/uns/src/UNsGrad.cpp` | 纯策略在 `UGradPolicy.h` |
| 限制器/重构（barth / vencat / 0 阶） | `codes/usolver/src/ULimiter.cpp`、`ULimiterRange.h`、`codes/uns/src/UNsLimiter.cpp` | |
| 无粘面通量（12 个 scheme） | `codes/ns/src/NsInvFlux.cpp`（1291 行，`Roe`/`Vanleer`/`Steger`/`Hlle`/`LaxFriedrichs`/`Ausmp`/`AusmpUp`/`Ausmdv`/`Ausmw`/`Ausmpw`/`Slau2`） | 纯函数拆在 `NsInvFluxPolicy.cpp` |
| 粘性面通量 | `codes/ns/src/NsVisFlux.cpp`、`codes/uns/src/UNsVisFlux.cpp` | 策略 `NsVisFluxPolicy.h` |
| 边界条件（固壁/远场/进口/出口/周期/对称…） | `codes/ns/src/NsBcSolver.cpp`、`codes/uns/src/UNsBcSolver.cpp`、`codes/ns/src/NsBcPolicy.cpp` | ghost cell 写入 |
| 面法向/面积/网格速度几何 | mesh 层（`ug.xfn/yfn/zfn/vfn/farea`） | numerics 只消费 |
| 残差有向沉积 | `codes/usolver/src/UResidual.cpp` + `UResidualPolicy.h` | 全日志器路径共享 |
| 谱半径（面/单元） | `codes/ns/src/TimeStep.cpp`、`codes/uns/src/UNsSpectrum.cpp`、`codes/ns/src/NsSpectrum.cpp` | 见 §3.3：既被显式 dt 用，也被隐式谱用 |

### 3.2 Temporal

| 责任 | 实现文件 | 备注 |
|---|---|---|
| 阶段序列/循环 | `codes/solver/src/TimeIntegral.cpp`（`RungeKutta`/`Lusgs`/`Simple`） | 由任务名驱动，无状态 |
| RK 系数 | `codes/global/src/Ctrl.cpp:85-90`（`rk_coef[i] = 1/(i+1)`） | 另在 `codes/ns/src/NsCtrl.cpp:62`、`codes/ins/src/INsCtrl.cpp`、`codes/turb/src/TurbCtrl.cpp` **各复制一遍** |
| 局部时间步 | `codes/ns/src/TimeStep.cpp` + `TimeStepPolicy.cpp`、`codes/uns/src/UTimeStep.cpp` | |
| 状态更新 | `codes/ns/src/NsUpdate.cpp`（`CalcFlowField` = `q0 + dq`，非物理回退 `WeekSolutionFix`/`SolutionFix`）、`codes/uns/src/UNsUpdate.cpp` | |
| 显式 LHS（实为残差缩放） | `codes/usolver/src/ULhs.cpp`（`res *= dt/cvol*lhscoef`） | 名为 LHS，MULTI_STAGE 下**没有矩阵** |
| 隐式 LU-SGS 前/后扫 | `codes/ns/src/NsLusgs.cpp`、`codes/uns/src/UNsLusgs.cpp`、`codes/implicit/src/*` | `UNsLusgs::LowerSweep` 每次调用 `Init()` → `CellTopo::CalcC2f()` |
| 双时间源项/线性两步系数 | `codes/unsteady/src/TimeIntegration.cpp`、`codes/usolver/src/UUnsteady.cpp` | 见 §8 风险 R4/R5 |
| 迭代收敛判据 | `codes/solver/src/SolverState.cpp:117`（`SolverState::Converge`）、`codes/unsteady/src/UnsteadyConvergence.cpp` | 见 §5.4 |

### 3.3 真正"shared"的东西（既不属于 spatial 也不属于 temporal）

| 共享对象 | 定义处 | 被谁读写 |
|---|---|---|
| `Ctrl ctrl` | `codes/global/include/Ctrl.h`、`Ctrl.cpp` | 115 处 `ctrl.` 引用；`ilim`/`ischeme`/`rk_coef`/`idualtime` 同时被空间与时间读 |
| `NsCom nscom`（**818 处 `nscom.` 引用**） | `codes/ns/include/NsCom.h` | 混装：配置（`ischeme`/`visSRModel`/`nTModel`）、物性（`gama_ref`/`reynolds`/`prl`）、**逐面临时标量**（`q1`/`q2`/`gama1`/`gama2`/`invsr`/`vissr`/`timestep`）、守恒态（`q`/`q0`/`dq`/`prim`） |
| `NsInv inv` | `codes/ns/include/NsInvFlux.h`、`NsInvFlux.cpp:36` | 空间核的工作区单例 |
| `UCom ug` / `UNsCom unsf` | `codes/usolver/include/UCom.h`、`codes/uns/include/UNsCom.h` | 网格视图与字段指针缓存 |
| `gcom` | 面几何全局（`PrepareFaceValue` 每面刷新） | 空间（通量）、时间（`NsLusgs::CalcFaceEigenValue`）、谱半径共用 |
| `Iteration::*` | `codes/global/include/Iteration.h` | CFL、内外步计数、输出节奏，被 numerics 与 IO 同时读 |
| `timeIntegration` | `codes/unsteady/src/TimeIntegration.cpp:24` | 谱半径（`UNsSpectrum`/`UINsSpectrum`/`UTurbSpectrum`）与双时间源项共用 |
| `LusgsState`（静态 `HXVector<unique_ptr<LusgsSolver>>`） | `codes/solver/include/SolverState.h:37`、`SolverState.cpp` | 由 `UINsSolver` 注册，由 `ImplicitTaskReg` 取用 |

`nscom` 一个全局里同时装着**配置、物性、守恒态、逐面临时量**，是本次审计中最强的"没有所有权边界"的证据。

---

## 4. 状态/资源所有权

按"谁分配、谁释放、谁重置"整理：

| 资源 | 所有者 | 生命周期 | 风险点 |
|---|---|---|---|
| `NsInv inv`（工作区 `prim/q/dq/flux`） | `NsInv::Init()`（`NsInvFlux.cpp:38`） | 每次 `UNsInvFlux::CalcFlux()` 重 `resize`；跨阶段保留容量 | 单例 + 逐面复用；不存在"一次通量计算的局部状态"，因此**无法在进程内并行/重入** |
| `static UNsInvFlux uNsInvFlux` | 函数级 static（`NsRhs.cpp:188`） | 进程级，跨 zone 复用 | 持有 `limiter`（`unique_ptr<NsLimiter>`）与 `hipBoundary*` 缓存；zone 切换靠 `Init()` 重绑定，缓存键用 `grid` 指针比较（`hipBoundaryGrid`） |
| `LimField* limf` / `qf1`/`qf2` | `NsLimiter`（`codes/usolver`） | 随 `UNsInvFlux` 构造 | 面索引 MRField，承载重构结果，被 `PrepareFaceValue` 读 |
| `invflux`（MRField） | `UNsInvFlux::Alloc()` | 每次 `CalcFlux()` | 面索引残差容器 |
| `unsf.q/q1/q2/dq/...` 指针 | `BindUNsFieldPointers()`（`codes/uns/src/UTimeStep.cpp:154`） | 每次绑定覆盖全局 | **多个模块各自 bind**（`UTimeStep.cpp`、`UNsInvFlux.cpp`、`UNsUpdate.cpp`、`UNsLusgs.cpp` 均调 `unsf.Init()`），谁最后调用谁决定指针内容 |
| `unsf.tempr` | `NsCalcGamaT()`（`NsRhs.cpp`） | 每 stage | 三方程无温度配置下为 `nullptr`，代码里显式判空 |
| LUSGS `dq` 交换缓冲区 | `codes/interface`（`ExchangeInterfaceDq`） | 每 sweep | 隐式扫依赖 halo，跨越 numerics/interface 两个 workstream |
| `LusgsState::uns/str` | `SolverState.cpp:37`，`SolverMap.cpp:121` `Init(nSolver)` | solver 注册时 | `static` 全局容器，`Reset()` 在 `SolverMap.cpp:135` |
| `timeIntegration.coeff/sp1/sp2/sc*` | `TimeIntegration::Init()` ← `InitFlowSimuGlobal()` | 进程级一次 | `UUnsteady` 内部又调一次 `timeIntegration.Init()`（`UUnsteady.cpp:34`）与 `CalcResCoef/CalcSrcCoeff`，存在二次写入 |
| `ctrl.rkStage` | `RkStageScope`（`TimeIntegral.cpp:135`，RAII） | 单 stage | 唯一一处"用 RAII 表达时间阶段所有权"的正面例子 |
| `gcom.*` | `PrepareFaceValue()` | 单面 | 时间模块（`NsLusgs::CalcFaceEigenValue`）也读同一个 `gcom`，离开面循环后其值无定义 |

---

## 5. 现有测试覆盖（只描述存在性，不声明通过）

### 5.1 `tests/`（GoogleTest，CTest 注册）

| 测试目标 | 覆盖的生产文件 | 覆盖层次 |
|---|---|---|
| `inv_flux_policy_test` | `codes/ns/src/NsInvFluxPolicy.cpp` | 纯函数（Roe 平均、熵修正、AUSM split、Roe 耗散、ULP 质量） |
| `bc_policy_test` | `NsBcPolicy.cpp` + `NsInvFluxPolicy.cpp` | 纯函数 |
| `lusgs_policy_test` | `NsLusgsPolicy.cpp` | 纯函数 |
| `timestep_policy_test` | `TimeStepPolicy.cpp` | 纯函数 |
| `vis_flux_policy_test` | `NsVisFluxPolicy.cpp` | 纯函数 |
| `sutherland_policy_test` | `Sutherland.h` | 纯函数 |
| `grad_policy_test` | `codes/usolver/src/UGradPolicy.cpp` | 纯函数 |
| `limiter_policy_test` / `limiter_range_test` | `ULimiterPolicy.cpp` / `ULimiterRange.h` | 纯函数 |
| `residual_deposit_test` | `UResidualPolicy.h` | 纯函数（符号/boundary-only-left） |
| `turb_*_test`（13 个） | `codes/turb/src/Turb*Policy.cpp` | 纯函数 |
| `oneflow_euler_inv_flux_capability_test` | `codes/accel/include/EulerInvFluxCapability.h` | 能力真值表 |
| `oneflow_euler_runge_kutta_capability_test` | `codes/solver/include/EulerRungeKuttaCapability.h` | 能力真值表 |
| `oneflow_euler_backend_cpu_test` / `..._domain_backend_adapter_test` / `oneflow_euler_domain_state_*` / `oneflow_euler_cpu_domain_backend_test` | `ports/kunshan/oneflow_1d_hip/*`、`codes/accel/src/CpuEulerDomainBackend.cpp` | **1D 独立端口**，不是主 solver |
| `oneflow_flux_backend_bridge_test`（standalone main，非 gtest） | `CpuFluxBackend.cpp` vs 1D 端口 | 1D 周期算例逐步对比 |
| `hip_resident_state_contract_test` | `HipSkipHostRkSlotCopy` 等纯谓词 | 纯函数 |

### 5.2 `test/`（python 回归，CI 用 `test/test.py` 驱动，不经 CTest）

- `suites/cpu-serial.txt`：5 个算例；`suites/ci-serial.txt` = `@cpu-serial.txt` + `@grid-serial.txt`。
- 从各算例 `script/*.txt` 读出的实际数值配置（**关键证据**）：

| 算例 | `ischeme` | `time_integral` | `ilim` | `idualtime` |
|---|---|---|---|---|
| `plateuns2dslau2` | 12 (SLAU2) | 2 (LU-SGS) | 0 | 0 |
| `plateuns2dslau2_34950_35000` | 12 | 2 | 0 | 0 |
| `rae2822_roe_sa` | 1 (Roe) | 2 | 2 (vencat) | 0 |
| `m6wingroe_sa` | 1 | 2 | 0 | 0 |
| `turbplateuns2droe_sa` | 1 | 2 | 0 | 0 |

- `.github/workflows/linux.yml` 用 `1e-8` 跑 `suites/ci-serial.txt`；strict `1e-15` 是人工口径（`test/README.md`、`docs/plans/oneflow-development-todo.md` §7）。
- `test/baselines/residual-baseline.json`：`absolute_tolerance: 1e-15`、`nan_policy: reject`、50 步逐通道残差。

### 5.3 覆盖结论（可核验的负面事实）

- `grep -rn "UNsInvFlux|NsCalcRHS|NsInvFlux|UNsGrad|NsUpdate|UNsUpdate|TimeIntegral|NsLusgs|UNsLusgs" tests/` 只命中 `NsInvFluxPolicy.h` 与任务名往返测试 → **没有任何测试实例化空间/时间 kernel 类**。
- **显式多阶段 RK（`time_integral=1`）在 `test/` 中没有任何残差基线**；5 个算例全部是 `time_integral=2`。设备路径（`MULTI_STAGE`）的数值门禁只有 `ci/kunshan/f2-main-solver-*.sh` 的跨后端比较，不是仓库内 CTest/`test/` 基线。
- `UseCpuBatchAdapter()` 需要 `ONEFLOW_ENABLE_UNS_CPU_BATCH=1`（默认关闭）→ **批处理 CPU 路径在 CI 中从不执行**；它被 `tests/euler/flux_backend_bridge_test.cpp` 间接覆盖，但那条测试只把 `CpuFluxBackend` 与 1D 端口比较，**没有与生产逐面路径 `NsInvFlux::*` 比较**。
- `OneFLOW_GeometryCore`/`oneflow_euler_*` 的 `-ffp-contract=off` 是 GCC/Clang 写法（`tests/euler/CMakeLists.txt:21/51/83`）；MSVC 是否接受、以及收缩是否被钉住，本环境无法验证。

### 5.4 顺带查到的一处判据实现事实（需人工确认，见 §11）

`codes/solver/src/SolverState.cpp:117`：

```cpp
bool SolverState::Converge()
{
    if ( Iteration::innerSteps == 0 ) return false;
    if ( ctrl.idualtime == 0 ) return true;      // 稳态：内层恰好 1 次
    bool flag = true;                            // 双时间分支
    for ( ... ) SingleSolverSingleGridTask( kCalcUnsteadyCriterionTaskName );
    return flag;                                 // flag 从未被修改 → 恒真
}
```

而 `UnsteadyConvergence::conv`（`codes/unsteady/src/UnsteadyConvergence.cpp:60/71/80`）在整棵树上**没有任何读取者**（`grep "\.conv\b"` 只命中赋值）。结合 `MG::Run()` 的 `while (!SolverState::Converge())` 结构，结论是：**双时间的内层迭代始终只执行 1 次，`iConv`/`conv`/`maxIterSteps` 的可调语义当前不生效**。这条不改变本审计的边界图，但决定 §10 的 P5 是否存在。

---

## 6. 缺失的契约测试（按性价比排序）

每条给出"边界 + 可观测断言"，都可独立实现、不需要新建框架。

**C1. 逐面通量与批处理通量等价**（最高优先级）
- 边界：`NsInvFlux::LaxFriedrichs()`（逐面，用 `inv.gama`/每侧 gamma）vs `CpuFluxBackend::CalcEulerLaxFriedrichsFlux`（用 `state.gamma = nscom.gama_ref`）。
- 断言：同一组 face state 下 5 个分量逐面 `|f_face - f_batch| <= 1e-14 * |f|`；常数 gamma 时成立；非均匀 gamma 时**显式记录为已知差异**而不是假装相等。
- 价值：这是唯一能让 `ONEFLOW_ENABLE_UNS_CPU_BATCH=1` 可被接受的证据，也是 B1（Roe/AUSM 设备化）的 oracle 前置。
- 注意：两者的"面通量"量纲不同（§2.3），测试必须显式除以 `farea` 或统一约定。

**C2. 全部 `ISCHEME_*` 常量都产生有限通量**
- 边界：`NsInvFlux::SetPointer` 的 12 个分支 + `else` 默认分支。
- 断言：对每个常量，给定一组超声速/亚声速面状态，`inv.flux[0..5)` 全部 `std::isfinite`；`ISCHEME_HYBRIDROE(11)` 当前落到**空函数** `HybridRoe(){}`，`RoeOld(){}` 同理 → 测试会立刻暴露"选中即静默零通量"。
- 断言二：未知 `ischeme`（如 99）走 `else` 静默退化为 Roe，应改为显式契约（记录或 fail-fast），测试锁住选择哪一条。

**C3. 时间积分系数表**（几乎零成本，收益高）
- 边界：`TimeIntegration::InitCoef()` / `CalcResCoef()` / `CalcSpectrumCoeff()` / `CalcSrcCoeffBasic()`——纯算术，唯一障碍是它直接读全局 `ctrl`。
- 断言：`linearTwoStepMethods` 0..11 的 `(thet, xi, phi)` 与 `resc1..3`、`sp1/sp2`、`bsc1..3` 表驱动对比（BDF1/BDF2/CN 至少 3 个）；`linearTwoStepMethods = 12` 时**当前行为是打印错误并继续用 0 系数（退化为显式 Euler）**，测试要么锁住该行为，要么把它变成错误——这是必须由人决定的数值语义。
- 同时应覆盖 `rk_coef[i] = 1/(i+1)` 这一在 4 处重复实现的公式。

**C4. 每个 RK 阶段的残差基底契约**
- 边界：`DeviceOwnsExplicitRkBookkeeping()`（`codes/solver/src/ResidualBasePolicy.cpp`）+ `LoadResiduals`/`LoadQ` + `UNsInvFlux` 的设备残差基底。
- 现有：`HipSkipHostRkSlotCopy` 纯谓词已有测试（`tests/hip_resident_state_contract_test.cpp:265`）。
- 缺失：用纯 host 路径证明"3 阶段后 `res` 不等于单阶段累加"（`docs/plans/oneflow-development-todo.md` §5 第 2 项记录过该回归：残差拓到 `5.52e-4`）。可用一个 1 单元/2 面的合成状态，断言"每阶段不累积上一阶段残差"。

**C5. LUSGS 单扫的正确性（合成 2×2 三角块）**
- 边界：`UNsLusgs::LowerSweep`/`UpperSweep`（前扫正序、后扫倒序、`Update()` 写 `dq`）。
- 断言：2 个单元 1 个内面的解析 A 阵，一次 lower+upper 后 `dq` 等于直接求解的解析解（对标 `docs/source/numerical/lusgs_flux_increment.rst`、`tdma.rst`）。当前只有 `NsLusgsPolicy` 的标量片段有测试。
- 附带断言：`UNsLusgs::LowerSweep()` 会重建 `CellTopo::CalcC2f`，测试应固定"拓扑不被数值扫描改变"。

**C6. 梯度边界交换的一次滞后契约**
- 边界：`Grad::SwapBcGrad()`（`codes/grad/src/Grad.cpp:122`）与 `uns_grad.CalcGradHip(...)` 的 `mirrorBoundaryGradient` 决策。
- 断言：串行 2 zone 合成网格上，`SwapBcGrad` 后 ghost 梯度等于"邻居上一轮值"，且 `ONEFLOW_HIP_GRADIENT_MIRROR_ON_DEMAND=0` 与默认值产生相同 ghost 值（只影响 host 镜像时机，不影响物理）。这条直接对应 `docs/plans/oneflow-development-todo.md` §4.1 A2 被门禁否掉的"删第 4 轮界面交换"。

**C7. 状态更新的非物理回退契约**
- 边界：`NsUpdate::CalcFlowField()` + `WeekSolutionFix()`/`SolutionFix()`/`PositiveUpdate`。
- 断言：给定 `rho<=0` 或 `p<=0` 的 `q0+dq`，结果满足 `rho>0 && p>0`，且 `ONEFLOW_FAIL_ON_NONPHYSICAL=1` 时抛异常。当前只在算例层面被 `nonphysical_hits=0` 观察。

**C8. 空间算子的守恒性（sum of fluxes = 0）**
- 边界：`AddF2CField` + `DepositOrientedFaceFlux`。
- 断言：周期/闭合网格上 `Σ_cell res = 0`（机器精度）。`residual_deposit_test` 只测单面符号，未测全局守恒。

**C9. `nscom.gama1/gama2` 的 per-cell 与 per-face 一致性**
- 边界：`PrepareFaceValue`（逐面写 `nscom.gama1/2`）与 `UTimeStep`（逐面写 `nscom.q1/q2`）。
- 断言：通量与时间步两条路径读到的同一面状态一致；否则 dt 与 flux 来自不同的重构结果。当前无任何测试。

---

## 7. 共享文件（不能由单一 worker 拥有）

| 文件 | 为什么不能单属 | 建议归属 |
|---|---|---|
| `codes/global/include/Ctrl.h` + `src/Ctrl.cpp` | `ctrl.ilim`/`ischeme`/`rk_coef`/`idualtime`/`pdt` 同时被 spatial、temporal、IO、输出节奏读写 | 人工（shared foundation）；worker 只能提交 RFC |
| `codes/ns/include/NsCom.h` + `src/NsCom.cpp`（818 处引用） | 配置 + 物性 + 守恒态 + 逐面临时量混合 | 人工；拆分本身是一次 `need_decision` 的架构变更 |
| `codes/global/include/Iteration.h` | CFL、内外步、输出节奏，numerics 与 IO 共写 | 人工 |
| `codes/solver/include/SolverState.h` + `src/SolverState.cpp` | `LusgsState` 静态注册 + `Converge()` 判据 | 人工 |
| `codes/interface/*`（`ExchangeInterfaceDq`、`UpdateInterfaceData`） | halo/界面交换是 numerics 的硬依赖，但属 mesh/parallel 域 | 第三 workstream（未来） |
| `codes/multigrid/src/Multigrid.cpp` | **时间迭代驱动的实际所有者** | 人工或独立 `temporal-driver` scope |
| `codes/accel/include/FluxBackend.h`、`EulerCpuAdapter.*`、`CpuFluxBackend.*` | 后端中立 seam，被 spatial 与 1D 端口共用 | 人工 |
| `codes/solver/include/Euler*Capability.h` | 设备/显式路径的能力真值表 | 人工（E2 的落点） |
| `tests/solver/CMakeLists.txt`（336 行单文件） | 每个数值 scope 都要往里加 target → 必然冲突 | 人工一次性拆成 `tests/numerics/<scope>/CMakeLists.txt` |
| `codes/CMakeLists.txt` | 递归收集，新增 .cpp 自动进 `OneFLOW` 可执行文件；但新测试目标不会自动注册 | 人工 |
| `test/suites/*.txt`、`test/baselines/*` | 数值基线 = 人工所有 | 人工 |

---

## 8. 风险

**R1（高）显式多阶段 RK 无回归基线。** 设备主路径（`MULTI_STAGE`）与 `time_integral=1` 在 `test/` 里没有任何残差基线；一旦有人改动 `LoadQ`/`LoadResiduals`/`rk_coef`，只有昆山跨后端脚本能发现。`docs/plans/oneflow-development-todo.md` §5.2 已记录过一次真实回归（`5.52e-4`）。建议把"注册一个显式三阶段 CPU 算例"提到 B1 之前。

**R2（高）`ischeme`/`ilim` 无取值范围校验。** `SetPointer` 未知值静默退化为 Roe；`ISCHEME_HYBRIDROE=11` 映射到**空函数** → 零面通量、零残差、迭代"收敛"。算例脚本里 `ischeme` 是自由整数（`// 1 roe ... 12 slau2`），没有 schema 校验。

**R3（中）两条空间算子实现的量纲/参数不一致。** 逐面路径用每单元 gamma 且不乘面积；`FluxBackend` 路径用 `nscom.gama_ref` 且乘面积。两条路径没有交叉测试（C1）。任何"设备/批处理与 CPU 逐面一致"的结论目前都缺证据。

**R4（中）双时间内层迭代实际恒为 1 次**（§5.4），`iConv`/`conv`/`maxIterSteps` 语义失效。要么补齐，要么在能力表里显式拒绝。

**R5（中）`TimeIntegration` 系数被二次写入 + 非法值静默退化。** `InitFlowSimuGlobal()` 调一次，`UUnsteady.cpp:34` 又调一次；`linearTwoStepMethods >= 12` 打印错误后继续用 0 系数（显式 Euler）。

**R6（中）隐式扫描在数值步里重建拓扑。** `UNsLusgs::LowerSweep/UpperSweep` 每次都 `Init()` → `CellTopo::CalcC2f(faceTopo)`，把 mesh 派生数据的所有权重叠进时间步。HIP 化 LU-SGS（todo B5）会先撞到这个。

**R7（中）`ns_ctrl` 与 `ctrl` 双份控制面。** `NsCtrl : public Ctrl` 又声明了 `ilim`/`time_integral`/`rk_coef`/`pdt`/`linearTwoStepMethods` 等**同名字段**；`ns_ctrl` 在生产里几乎没用（仅 `NsBcSolver.cpp:124` 与 `HeatFluxTask.cpp:65` 读 `isowallbc`），但 `NsCtrl::Init()` 仍会跑并覆写这些影子字段。`INsCtrl`/`TurbCtrl` 同样复制了 `rk_coef` 构造。

**R8（低-中）单例工作区阻碍并行与重入。** `inv`、`gcom`、逐面 `nscom.gama1/2`、`ug.fId/lc/rc` 都是"当前面"的隐式状态 → 逐面循环不可并行、不可嵌套（这也是 HIP 只能按 batch 重写一遍的根本原因）。

**R9（低）`nscom` 字段语义漂移。** 同名 `q1/q2` 在 `TimeStep::CalcFaceInvSpec` 里是**面左右值**，在 `NsInv` 里是**守恒量数组**，在 `NsCom` 里又是另一套。改名/搬移极易引入静默错误。

**R10（低）`Relaxation(int nCycles)` 无 context 的重载仍存在**（`TimeIntegral.cpp:243`），会走非 `SimuContext` 的全局任务路径；`MG` 同时维护两条分支（`MG::PreRelaxationCycle`）。存在"新路径忘了同步旧路径"的风险，与 `ResidualBasePolicy.h` 注释里描述的同类问题一致。

---

## 9. 关于 starter 契约的处置意见

`.agentlab/contracts/{core_api,data_layout,interfaces}.md` 的三个抽象（`Field`/`Mesh`/`Comm`、`ISpatialOperator`、`ITimeIntegrator`）在 OneFLOW 中**没有对应物**：

- `Field::like()` ↔ 无；RK 工作区是 `MRField` + `FieldHome::Get*Field(grid, solverType, ...)` 按名字从全局 DataBase 取。
- `Mesh::dx()`拟合结构网格 ↔ OneFLOW 是非结构面/单元数组（`ug.xfn/yfn/zfn/vfn/farea`, `lcf/rcf`）。
- `Comm::haloExchange` ↔ `codes/interface` 的 `UploadInterfaceData`/`UpdateInterfaceData`/`DownloadInterfaceData` 任务，且**带宽/一次滞后语义**是离散的一部分（todo §4.1 A2）。
- `ITimeIntegrator::advance` ↔ `TimeIntegral` 静态函数 + 任务名 + 全局 `ctrl.rk_coef`；`order()` 无对应物。
- `[CONFIRM]` 项（符号约定、halo 归属、BC 归属、dt 裁剪）在 OneFLOW 中**已经由代码固定**：通量不含面积、沉积 left`-=`/right`+=`、dt 由 `CellTimeStepFromSpectrum` 给出且由 `CombineInvAndVisTimeStep` 谐波合并、BC 由 `UNsBcSolver` 写 ghost。

因此建议**不要**把这些文件改造成 OneFLOW 契约，而是用 RFC 提出替换契约（见 `.agentlab/RFC/numerics-2026-11-01-boundary-contracts.md`）。

---

## 10. 可并行的未来 scope 提案（文件级互斥）

前提：`codes/CMakeLists.txt` 递归收集，新增 `.cpp` 自动进入生产可执行文件；**测试目标注册是唯一冲突面**。建议人工先一次性建立：

```
tests/numerics/CMakeLists.txt            # 人工建一次
tests/numerics/spatial-flux/CMakeLists.txt
tests/numerics/spatial-grad-bc/CMakeLists.txt
tests/numerics/temporal-explicit/CMakeLists.txt
tests/numerics/temporal-implicit/CMakeLists.txt
tests/numerics/unsteady/CMakeLists.txt
```

之后每个 worker 只写自己目录下的 `CMakeLists.txt`。

| Scope | 建议 glob（写入 `scope/*.txt`） | 职责 | 前置 |
|---|---|---|---|
| **P1 `numerics-spatial-flux`** | `codes/ns/src/NsInvFlux.cpp`、`codes/ns/include/NsInvFlux.h`、`codes/ns/src/NsInvFluxPolicy.cpp`、`codes/ns/include/NsInvFluxPolicy.h`、`codes/uns/src/UNsInvFlux.cpp`、`codes/uns/include/UNsInvFlux.h`、`codes/usolver/src/ULimiter*.cpp`、`codes/usolver/include/ULimiter*`、`tests/solver/inv_flux_policy_test.cpp`、`tests/solver/limiter*_test.cpp`、`tests/numerics/spatial-flux/` | 12 个面通量 scheme + 重构/限制器；补齐 C1/C2/C8 | 需人工确认 `nscom.gama1/2` 是否可改为显式传参（否则只能加测试不改接口） |
| **P2 `numerics-spatial-grad-bc`** | `codes/grad/*`、`codes/usolver/src/UGrad.cpp`、`codes/usolver/include/UGrad*.h`、`codes/usolver/src/VisGrad.cpp`、`codes/uns/src/UNsGrad.cpp`、`codes/uns/include/UNsGrad.h`、`codes/ns/src/NsVisFlux*.cpp`、`codes/ns/src/NsBcSolver.cpp`、`codes/uns/src/UNsBcSolver.cpp`、`codes/uns/src/UNsVisFlux.cpp`、`tests/solver/{grad,vis_flux,bc}_policy_test.cpp`、`tests/numerics/spatial-grad-bc/` | 梯度、粘性通量、边界 ghost；补齐 C6 | 与 P1 在 `tests/solver/CMakeLists.txt` 冲突 → 走 §10 的 `tests/numerics/` 拆分 |
| **P3 `numerics-temporal-explicit`** | `codes/solver/src/TimeIntegral.cpp`、`codes/solver/include/TimeIntegral.h`、`codes/ns/src/TimeStep.cpp`、`codes/ns/include/TimeStep*`、`codes/uns/src/UTimeStep.cpp`、`codes/ns/src/NsUpdate.cpp`、`codes/uns/src/UNsUpdate.cpp`、`codes/usolver/src/ULhs.cpp`、`tests/solver/timestep_policy_test.cpp`、`tests/numerics/temporal-explicit/` | RK 阶段序列、局部 dt、显式状态更新；补齐 C4/C7 | 需要人工先在 `test/suites/` 注册显式算例（R1） |
| **P4 `numerics-temporal-implicit`** | `codes/implicit/*`、`codes/ns/src/NsLusgs.cpp`、`codes/uns/src/UNsLusgs.cpp`、`codes/ns/src/NsSpectrum.cpp`、`codes/uns/src/UNsSpectrum.cpp`、`tests/solver/lusgs_policy_test.cpp`、`tests/numerics/temporal-implicit/` | LU-SGS 扫掠与谱半径；补齐 C5 | 依赖 `codes/interface` 的 dq 交换 => 需人工确认 halo 归属 |
| **P5 `numerics-unsteady`** | `codes/unsteady/*`、`codes/usolver/src/UUnsteady.cpp`、`codes/uns/src/UNsUnsteady.cpp`、`tests/numerics/unsteady/` | 双时间、历史层、线性两步系数；补齐 C3 | **先答 §11 的问题**，否则本 scope 应标注为 `deferred` |

**不建议**并行的（保持人工或单线）：
- `codes/global/src/Ctrl.cpp`、`codes/ns/src/NsCom.cpp`（R7/R9，属 shared foundation，改动即需 RFC）。
- `codes/multigrid/src/Multigrid.cpp`（迭代驱动所有权，见 §2.1）。
- `codes/accel/*` 的 `FluxBackend`/capability 表（E2 的落点，跨 workstream）。
- `.agentlab/contracts/*`、`.agentlab/verification/*`、`test/baselines/*`。

---

## 11. 一个需要人工确认的具体问题

> **问题：双时间（`ctrl.idualtime == 1`）在当前树里是要"补齐成可用能力"，还是"显式声明不支持"？**

证据：`SolverState::Converge()` 在内层第 1 次之后恒返回 `true`（`codes/solver/src/SolverState.cpp:117-131`，局部 `flag` 从未被修改）；`UnsteadyConvergence::conv` 全树无读取者；5 个回归算例全部 `idualtime = 0`；`docs/plans/oneflow-development-todo.md` §1 把"双时间"列在"尚未覆盖"。

**建议：先声明不支持并 fail-fast（把 `idualtime == 1` 加入 `EulerRungeKuttaCapability`/启动校验的拒绝原因），等 B5（LU-SGS 设备化）之后再作为独立工作项补齐。**

理由：现状是"看起来有双时间，实际内层只走 1 次"，属于静默的数值语义缺失，比明确拒绝更危险；而一旦声明支持，就必须同时交付收敛判据的契约测试（C3 + 新的内层判据测试）与至少一个 `test/` 算例，成本远超本次 audit 的范围。若 Answer 是"要补齐"，则 §10 的 P5 从 `deferred` 变为 `todo`，并需要先决定内层判据（`iConv` 三选一的语义）与 `maxIterSteps` 的所有权。

（其余需要人确认但**不阻塞**本报告的事项，已并入 §10 的前置列与 §9 的契约替换建议。）

---

## 12. Gate 结果 / 提交 / Repro

### Gate

两次运行（交付前 / SHA 回填后复核），命令相同：

```
$ cd <checkout root> && GATE_MODE=audit \
    /Users/luql/Downloads/agent/repo/open-research-agents/platform/scripts/gate.sh numerics

# run 1（交付前）
== scope
scope: numerics OK (2 changed files)
repro record: .agentlab/repro/numerics-20261101T111218Z.json
----
GATE PASS (mode: audit, read-only)
EXIT=0

# run 2（回填 SHA 后复核；3 个改动文件 = 报告 + RFC + run 1 的 repro 记录）
== scope
scope: numerics OK (3 changed files)
repro record: .agentlab/repro/numerics-20261101T111252Z.json
----
GATE PASS (mode: audit, read-only)
EXIT=0
```

- 审计交付提交：`e936b8191`（`numerics: audit A2 numerics boundary map and contract RFC`）。报告正文无法预写自身 SHA，因此 SHA 与最终 gate 复核记录由随后的一个仅记录证据的提交回填（即 `git log -2 -- .agentlab/reports/numerics.md` 中的后一条）。
- 分支：`agent/numerics`（未 push、未 merge、未触碰 `master`/`dev`）。
- gate 运行时 `git_commit` = `df8676fc24a35439d355eff20a6fdca71b342313`。
- 注意：gate 只跑 `scope_check` + `repro_snapshot`，audit 模式**不构建、不跑测试**，因此 `GATE PASS` 只表示"只读审计的改动范围合法"，不表示任何数值正确性。

### Repro record

- commit / dirty：基线 `df8676fc2`；审计期间工作树只新增 `.agentlab/` 下文件。
- machine and hardware：本地开发机（macOS），**未构建、未运行任何测试或算例**。
- compiler and flags：未使用。
- commands and environment：仅只读检索（`git log/diff/status`、`grep`、`read`）。
- input data/version：源码树即基线提交；未引入外部数据。
- result files under `.agentlab/repro/`：本次 gate 运行生成的那一个 `.agentlab/repro/numerics-<UTC timestamp>.json`（由 `gate.sh` 在 audit 模式下自动生成；文件名带运行时刻，故此处不写死）。本报告**不提出**任何精度/收敛/性能结论，因此该记录只承载 gate 元数据，不含算例字段。其中 `host` 字段已按仓库 `AGENTS.md` "Publication boundary"（提交产物不得含主机名）改为占位说明；其余字段为脚本原样输出。
- domain fields（`.agentlab/repro_fields.yaml`）：不适用（N/A）——无 accuracy/performance 结论，故 `mesh`/`numerics`/`solver`/`timing` 全为空。

### Files changed

- `.agentlab/reports/numerics.md`（本文件，新建）
- `.agentlab/RFC/numerics-2026-11-01-boundary-contracts.md`（新建，OneFLOW 专用 numerics 契约提案）
- `.agentlab/repro/numerics-<UTC timestamp>.json`（gate 自动生成；`host` 已脱敏）
- 无其他改动（`scope/numerics.txt` 只允许上述路径）。

命名说明：任务书示例为 `.agentlab/RFC/<日期>-numerics-<主题>.md`，但 `scope/numerics.txt` 的 glob 是 `.agentlab/RFC/numerics-*`。为同时满足 `scope_check` 与"日期+主题"要求，采用 `numerics-<日期>-<主题>.md`。

---

## 13. Known limitations

1. **未运行任何构建/测试**，因此本报告不含"某测试通过"的声明；对测试的引用一律是"存在/不存在"。
2. 未审查 `codes/ins`、`codes/uins`、`codes/scalar`、`codes/uturb`、`codes/chem*` 的 numerics 细节（它们复用同一批 policy 与任务名，但不在 A2 目录清单内）；`codes/turb` 只统计了测试存在性。
3. `codes/accel` 的 HIP 具体 kernel 实现（`ports/kunshan`、`codes/accel/*Hip*`）只作为边界提及，未审计其数值正确性。
4. §5.4 / §8 R4 关于双时间内层只迭代 1 次的结论基于源码静态阅读，未用运行时 trace 验证（本环境不可构建）。
5. §5.3 "显式 RK 无基线"的结论基于 5 个算例 `script/*.txt` 的配置阅读；若 CI 之外存在其他算例集，本结论需相应收窄。
