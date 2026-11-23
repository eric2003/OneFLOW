# RFC: OneFLOW numerics 边界契约（替换 starter 的 spatial/temporal 契约）

Proposer module: numerics
Date: 2026-11-01
Status: proposed | approved | rejected (only a human changes this)

配套审计：`.agentlab/reports/numerics.md`（同一提交）。
基线提交：`df8676fc2`（`agent/numerics` == `dev`）。

---

## Background and problem

`.agentlab/contracts/` 里的 `core_api.md`（`Field`/`Mesh`/`Comm`）、`interfaces.md`
（`ISpatialOperator`/`ITimeIntegrator`）、`data_layout.md`（全 TODO）是 starter 模板。
`contracts/README.md` 已声明它们不是 OneFLOW 契约。本次审计确认：这三个抽象在
OneFLOW 主 solver 中**一个都不存在**：

| starter 抽象 | OneFLOW 现状（源码证据） |
|---|---|
| `Field`（自带存储 + `like()`） | `MRField = Marray<Real>` = `HXVector<HXVector<Real>>`（`codes/array/include/HXArray.h:29`、`Marray.h`）。工作区不是"申请一块 like 存储"，而是按名字从全局 DataBase 取：`FieldHome::GetResidualField/GetStageSolutionField/GetUpdateDeltaField(grid, solverType, k)`（`codes/uns/src/UTimeStep.cpp:154` 起） |
| `Mesh`（`nx/ny/nz/dx/periodic`） | 非结构面/单元数组：`ug.lcf/rcf`（左右单元）、`ug.xfn/yfn/zfn/vfn/farea`（面法向/网格速度/面积）、`ug.xcc/ycc/zcc/cvol`（单元中心/体积）。无 `dx()` 概念 |
| `Comm::haloExchange` | 任务三元组 `UPLOAD_INTERFACE_DATA` → `UPDATE_INTERFACE_DATA` → `DOWNLOAD_INTERFACE_DATA`（`codes/solver/src/SolverImp.cpp:29` `CommInterfaceData()`），LU-SGS 另有 `EXCHANGE_INTERFACE_DQ`（`codes/interface/src/InterfaceTaskReg.cpp:324`） |
| `ISpatialOperator::residual(U, t, R)` | 逐面 kernel `(this->*invFluxPointer)()`，输入来自全局 `inv.prim1/2`、`gcom.*`，输出 `inv.flux[]`（`codes/uns/src/UNsInvFlux.cpp:875 PrepareFaceValue` / `:903 UpdateFaceInvFlux`） |
| `ITimeIntegrator::advance` / `order()` | 静态函数 + 任务名序列（`codes/solver/src/TimeIntegral.cpp`）；无 `order()` |
| `[CONFIRM]` 符号/归一化/halo/BC/dt | 已被代码固定（见 Proposal 第 0 条），但这些事实从未写进任何契约文件 |

后果（审计 §7 的共享文件清单）：`Ctrl`、`NsCom`、`Iteration` 同时被空间、时间、
IO 读写，任何 worker 改动它们都会跨 workstream；而"每阶段残差基底"这类真正的
数值契约只以注释形式存在于 `codes/solver/include/ResidualBasePolicy.h`。

本 RFC 提出**替代契约**，明确"什么是稳定的、可以被 worker 依赖的"，并把
"什么必须由人工拥有"写死。**不提议重写求解器，也不提议引入 starter 的抽象。**

---

## Proposal (specific to interface signatures)

### 0. 先冻结既有事实（这些已经由代码决定，不要在 worker 改动中变更）

1. **数据布局**：`MRField` 为方程主序，每个方程分量是一段长度 `nTCell` 的连续
   `Real`；`nTCell = nCells + nBFaces`，ghost 单元紧接管有单元之后（依据：
   `NsSolverImp.cpp` 中 `ug.nTCell - ug.nCells == ug.nBFaces` 的校验）。
   限制：`MRField` 不携带 stride 元数据，调用方必须自行保证 `GetNEqu() >= nEqu`。
2. **守恒量顺序**：`q = [rho, rho*u, rho*v, rho*w, rho*E]`；湍流/温度等附加方程
   从 `nscom.nBEqu` 起追加（`IDX::IR..IDX::IP` = 0..4，`codes/ns/include/NsIdx.h`）。
3. **面通量量纲**：`NsInvFlux::*` 产出的是**未乘面积**的通量；面积乘与有向沉积
   由驱动层完成（`UpdateFaceInvFlux` × `gcom.farea`，`AddF2CField` →
   `DepositOrientedFaceFlux`：left `-=`、right `+=`、boundary 面只 left）。
   `FluxBackend` 路径产出的面通量**已乘面积** —— 这是两条实现之间的既有差异，
   在 §Proposal 3 里被显式建模而不是隐藏。
4. **时间推进语义**：`ctrl.rk_coef[i] = 1/(i+1)`（`codes/global/src/Ctrl.cpp:85`），
   每阶段 `ctrl.lhscoef = rk_coef[i]`；`ULhs::CalcLHS` 在 `MULTI_STAGE` 下只做
   `res *= dt/cvol*lhscoef`，**不组装矩阵**。只有 `LUSGS` 有隐式谱
   （`UNsSpectrum::CalcImplicitSpectrum`）。
5. **局部时间步**：`dt = CFL*Vol/lambda_inv`，再与粘性 dt 谐波合并
   （`TimeStepPolicy.h`：`CombineInvAndVisTimeStep`）。
6. **迭代所有权**：外层/内层循环与收敛判据在
   `codes/multigrid/src/Multigrid.cpp`（`MG::Run`）+ `codes/solver/src/SolverState.cpp`
   （`SolverState::Converge`），**不在 numerics 目录**，也不在 `codes/time/`
   （后者只有 `TimeSpan`/`TimeTest` 墙钟工具）。

### 1. 空间算子契约：`FaceFluxKernel`

取代 `ISpatialOperator`。目标：让"一次面通量计算"从全局单例变成可测试的纯函数，
同时**不改变现有数值结果**。

```cpp
// codes/ns/include/NsFaceFluxContract.h  （提案；新建）
// 单位：prim = [rho, u, v, w, p]；gama = 比热比；flux = 物理通量（未乘面积）
struct FaceKernelInput
{
    int nEquations;        // 3 或 5（含湍流标量时 >5，由 nBEqu 追加）
    Real nx, ny, nz;       // 单位面法向，指向 left -> right
    Real meshVelocityN;    // vfn
    const Real * primLeft; // 长度 nEquations
    const Real * primRight;
    Real gamaLeft, gamaRight;
};

struct FaceKernelOutput
{
    Real flux[ 8 ];        // 未乘面积；只写前 nEquations 项
};

// 前置条件：输入有限、rho>0、p>0；nEquations in {3,5}
// 后置条件：输出有限；不修改输入；不读任何全局量；可重入
// 失败语义：返回 false 且不承诺输出内容（调用方决定 fail-fast 或回退）
typedef bool ( *FaceFluxKernelFn )( const FaceKernelInput &, FaceKernelOutput & );
```

要求：
- 每个 `ISCHEME_*`（1..12）都必须映射到一个**非空**实现；当前
  `ISCHEME_HYBRIDROE(11)` → `HybridRoe(){}`、`RoeOld(){}` 为空函数，必须在契约里
  明确为"未实现 → 拒绝"，而不是"静默零通量"。
- 未知 `ischeme` 的处置（当前静默退化为 Roe）必须成为**可检查的选择**：契约要求
  初始化期校验一次，而不是每面隐式退化。
- 现有 `NsInvFlux` 可以保留为适配层（`FaceKernelInput → inv.prim1/2 + gcom → 调
  kernel → inv.flux`），但新代码与所有新测试只依赖 `FaceFluxKernelFn`。

### 2. 驱动契约：`ZoneInviscidFluxDriver`

取代 `UNsInvFlux` 的隐式逐面契约。只描述"循环 + 沉积"，数值全部委托给第 1 条。

```cpp
struct InviscidFluxBatchRequest
{
    int nCells, nFaces, nBoundaryFaces;   // 必须满足 nGhost == nBoundaryFaces
    const int * leftCell, * rightCell;
    const Real * xNormal, * yNormal, * zNormal, * meshVelocityN, * faceArea;
    // 面索引、方程主序，长度 nFaces * nEquations
    const Real * primLeft, * primRight;
    const Real * gamaCell;                 // 单元索引；可选（常数 gamma 时为 nullptr）
    bool multiplyByArea;                   // 见第 0 条第 3 项
};

// 后置条件：residual[equation][cell] +=/-= flux；Σ_cell residual == 0（闭合网格）
void AccumulateInviscidResidual(
    const InviscidFluxBatchRequest &, FaceFluxKernelFn,
    FaceKernelOutput * fluxOut /*可为 nullptr*/,
    Real * residual /*方程主序，nCells 项/方程*/ );
```

这条契约把"逐面 CPU 循环""CPU batch""HIP batch"三种实现统一成同一个可观测语义
（每个面一个 flux、面间互相独立），并给出唯一的守恒断言（`Σ residual == 0`）。
现有 `FluxBackend`/`EulerCpuAdapter` 保留为 `FaceFluxKernelFn` 的一个实现。

### 3. 时间算子契约：`TemporalScheme`

取代 `ITimeIntegrator`。描述"阶段序列 + 每阶段做什么"，不描述物理。

```cpp
enum class TemporalKind { ExplicitMultiStage = 1, Lusgs = 2, Simple = 3 };

struct StageRequest
{
    TemporalKind kind;
    int stageCount;             // ExplicitMultiStage: rk_coef.size()
    int stage;                  // 0-based；Lusgs/Simple 为 0
    Real coefficient;           // 该阶段的 ctrl.lhscoef
    Real physicalTimeStep;      // ctrl.pdt（双时间时为物理步长）
};

// 前置条件：stageCount > 0；coefficient 来自已校验的系数表
// 后置条件：只修改该阶段允许的字段（见下）；不重建拓扑；不发起 MPI 集合操作
// 阶段允许的副作用（取代"靠注释约定"）：
//   stage 内：LOAD_RESIDUALS / UPDATE_RESIDUALS / CALC_LHS / UPDATE_FLOWFIELD / CALC_BOUNDARY
//   阶段末：UPDATE_INTERFACE_DATA（仅多 zone 时，见 docs/plans/oneflow-development-todo.md §4.1 A2）
bool RunStage( const StageRequest &, StageHost & host );
```

要求：
- 阶段序列与系数由**一个**所有者提供。当前 `rk_coef` 的构造在
  `Ctrl.cpp`、`NsCtrl.cpp`、`INsCtrl.cpp`、`TurbCtrl.cpp` **复制了 4 份**；契约要求
  单一实现 + 一张表（对应 `docs/plans/oneflow-development-todo.md` §4.5 E2）。
- `TimeIntegration::InitCoef()` 的 12 组 `(thet, xi, phi)` 与派生的
  `resc1..3`/`sp1/sp2`/`bsc1..3` 属于契约的一部分：必须表驱动 + 单测，且
  `linearTwoStepMethods` 越界时**不得静默退化为显式 Euler**（当前行为是打印后继续）。
- 双时间内层判据必须先有契约才能算"支持"：当前 `SolverState::Converge()` 恒返回
  `true`、`UnsteadyConvergence::conv` 无读取者。在补齐之前，`idualtime == 1` 应在
  初始化期显式拒绝（fail-fast），而不是静默按 1 次内层迭代运行。

### 4. 能力表契约（单一真值表）

现有 `EulerInvFluxCapability`（`codes/accel/include/EulerInvFluxCapability.h`）与
`EulerRungeKuttaCapability`（`codes/solver/include/EulerRungeKuttaCapability.h`）
已经具备"请求 → 决策 + 稳定拒绝原因"的正确形状，且都有 gtest。契约要求把它们
合并为**一张**（solver × scheme × limiter × viscous × time-integral × grid-level ×
halo × dualtime）真值表，并由生产路径在初始化期调用一次，任何显式请求不支持的项
一律 fail-fast（当前 `HipProductionFeatureExplicitlyEnabled` 已实现该模式，应推广到
`ischeme`/`ilim`/`idualtime`）。

### 5. 明确**不**引入的东西

- 不引入 `Field::like()`、`Mesh::dx()`、`Comm::haloExchange`。
- 不引入虚基类 + 工厂来替换 `TimeIntegral`；任务名调度是现有 runtime 契约
  （`TaskRegister`/`MessageMap`），替换它属于 runtime 工作流，不在本 RFC 范围。
- 不在本 RFC 中把 `Ctrl`/`NsCom` 拆开（那是一次独立的 `need_decision` 架构变更）；
  本 RFC 只声明"新代码不得再向 `NsCom` 增加逐面临时量"。

---

## Affected modules and existing code

| 模块/文件 | 影响 |
|---|---|
| `codes/ns/include/NsInvFlux.h`、`src/NsInvFlux.cpp` | 新增 `FaceFluxKernelFn` 适配；删除/显式拒绝空 scheme（`HybridRoe`、`RoeOld`） |
| `codes/ns/{include,src}/NsInvFluxPolicy.*`、`TimeStepPolicy.*`、`NsBcPolicy.*`、`NsLusgsPolicy.*` | 不变（已是纯函数，是契约的参考实现） |
| `codes/uns/src/UNsInvFlux.cpp` | 逐面/CPU batch/HIP batch 三条路径收敛到 `AccumulateInviscidResidual` 语义 |
| `codes/accel/{include,src}/CpuFluxBackend.*`、`FluxBackend.h`、`EulerCpuAdapter.*` | 作为 `FaceFluxKernelFn` 的实现之一；`multiplyByArea` 差异显式化 |
| `codes/solver/src/TimeIntegral.cpp`、`TimeIntegral.h` | 引入 `StageRequest`；`rk_coef` 单一来源 |
| `codes/global/src/Ctrl.cpp`、`codes/ns/src/NsCtrl.cpp`、`codes/ins/src/INsCtrl.cpp`、`codes/turb/src/TurbCtrl.cpp` | 去掉 3 处 `rk_coef` 复制（人工拥有，需单独 RFC/decision） |
| `codes/unsteady/src/TimeIntegration.cpp` | 系数表驱动 + 越界 fail-fast（行为变更，需人工批准） |
| `codes/solver/src/SolverState.cpp`（`Converge()`） | 双时间语义（§11 问题） |
| `codes/multigrid/src/Multigrid.cpp` | 仍持有迭代循环；契约要求它只通过 `RunStage` 进入 numerics |
| `tests/solver/*`、`tests/numerics/*` | 新增契约测试；`tests/solver/CMakeLists.txt` 建议拆分为 `tests/numerics/<scope>/CMakeLists.txt` |

---

## Alternatives

1. **照搬 starter 契约**（实现 `ISpatialOperator`/`ITimeIntegrator`/`Field`）。
   否决：OneFLOW 的数据是 `MRField` + 全局 DataBase 字段名，halo 交换是"上传/更新/
   下载"三阶段任务且带一次滞后语义；照搬会把 3 份真实语义压缩成一个不存在的抽象，
   并且需要一次全量重写，风险远大于收益。
2. **什么都不做，只写审计报告**。
   部分可行：审计报告已给出边界图。但 `contracts/` 会继续是"TODO + starter"，
   下一个 worker 仍会以 `ISpatialOperator` 为起点误判边界。RFC 的成本只有一份文档。
3. **只补测试、不立契约**。
   可行但会退化：C1（逐面 vs 批处理等价）与 C2（scheme 覆盖）都需要一个"面通量
   是什么"的定义，否则测试只能照抄实现、无法判断对错。
4. **先拆 `Ctrl`/`NsCom`（globals 化）再立契约**。
   顺序相反。`docs/plans/oneflow-development-todo.md` §4.5 E1（去 ambient globals）
   已列为正确方向，但如果先拆 globals，契约测试就失去稳定的被测量；先立契约
   （纯函数 + 表驱动）可以让拆分变成"实现替换"，测试保持不变。建议 E1 在本 RFC
   的契约被批准之后再排期。

---

## Impact on verification cases

1. **不改变任何现有基线。** 本 RFC 要求所有改造在数值上与 `test/baselines/residual-baseline.json`
   一致（`absolute_tolerance: 1e-15`，strict）。任何改变 `test/cfd/*/autotest/*.dat`
   的改动都属另一个 `need_decision`。
2. **建议新增（人工执行，不在 worker scope）**：
   - 一个**显式三阶段 RK**（`time_integral=1`、`rk_stage=3`）的 CPU 算例，进入
     `test/suites/` 并生成残差基线。理由：当前 5 个算例全部 `time_integral=2`，
     显式路径在 `test/` 无任何基线（审计 §5.3 / 风险 R1），而设备主路径正是显式 RK。
   - 把 `flux_backend_bridge_test`（standalone main，当前不在 CTest）注册进 CTest
     或明确标注为手工测试，避免"看起来有覆盖但从未运行"。
3. **`.agentlab/verification/run_all.sh` 目前恒返回 2（UNVERIFIED）**，因此
   `GATE_MODE=integration` / `release` 必然拒绝。建议在批准本 RFC 后，把
   `test/test.py` + `suites/ci-serial.txt` + `baselines/residual-baseline.json`
   接成 `run_all.sh numerics`（参数化 tolerance，默认 strict `1e-15`）。
   在此之前，numerics 的 gate 只能停在 `audit`。
4. **契约测试的落点**（对应用审计 §6 的 C1–C9）：
   `tests/numerics/spatial-flux/`（C1/C2/C8）、`spatial-grad-bc/`（C6）、
   `temporal-explicit/`（C4/C7/C3）、`temporal-implicit/`（C5）。
   这些测试应只链接纯函数与新的契约头，不链接 `UNsInvFlux` 的全局状态，
   以保证与 code 布局无关。
