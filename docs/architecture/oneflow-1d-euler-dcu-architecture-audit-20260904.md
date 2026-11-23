# OneFLOW 1D Euler DCU 架构审计与接入方案

更新时间：2026-09-04  
范围：OneFLOW 主 solver、现有 accelerator seam、`ports/kunshan/oneflow_1d_hip` standalone port  
审计性质：源码与 Git 结构审计；本次未提交新的昆山作业，也未新增目标节点运行证据

## 1. 结论先行

当前仓库已经具备三块可复用基础：

1. 主工程有进程级 `AccelRuntime`，可以选择 CPU/HIP，并保存 rank/device 映射信息。
2. `codes/accel` 有 backend-neutral 的内存拷贝接口和一个 `FluxBackend` 视图 seam。
3. standalone 1D Euler 已经验证了 CPU/HIP 的 `CreateState → Upload → Advance → Download` 生命周期，以及 HIP 持久 device state、`NoTrace`/`FullTrace` 和 MPI 实验路径。

但这三块目前没有组成主 solver 的 GPU vertical slice。准确边界是：

> 当前是“主工程 runtime + scalar 验证 seam + standalone 1D Euler backend”，不是完整 OneFLOW NS solver 的 GPU backend。

因此下一步不应把 standalone 文件直接接进 `UNsInvFlux`，也不应只把一个 HIP kernel 塞进现有逐面循环。正确的最小接入单位应是以 `solverIndex + zone + gridLevel` 为生命周期边界的、backend-neutral 的 Euler/NS execution adapter；它在既有 field/task/time-integral 调度的明确边界被调用，并负责和 host field、interface、restart、residual/output 协作。

## 2. 审计基线

本次审计时：

- 工作树干净；`HEAD`/`origin/master` 为 `729d435f`。
- `upstream/master` 为 `92991805`；本地 fork master 与 upstream 已分叉，`git rev-list --left-right --count upstream/master...HEAD` 为 `2 6`。
- standalone 代码集中在 `ports/kunshan/oneflow_1d_hip/`，主工程 accelerator 代码集中在 `codes/accel/`。
- `CUDA` 和 `KOKKOS` 仍只能按接口/构建占位处理；本报告不把 HIP/DCU 证据外推给它们。

## 3. OneFLOW 主线真实调用链

```text
codes/main/src/OneFlow.cpp
  → Simulation::Run
  → SimuImp::PreProcess
      → SetUpParallelEnvironment
      → ReadControlInfo
      → InitializeAccelRuntime(rank, size)
  → SimuImp::RunSimu
      → simu_state.Init
      → ConstructSystemMap
      → FieldSimu
          → InitFlowSimuGlobal
          → MultiBlock::LoadGridAndBuildLink
          → SolverMap::CreateSolvers
          → InitializeSolver / INIT_FLOWFIELD
          → MultigridSolve
              → TimeIntegral::Relaxation
                  → RungeKutta or Lusgs
                      → LOAD_Q
                      → CALC_TIME_STEP
                      → UPDATE_RESIDUALS
                          → inviscid + viscous + source terms
                      → CALC_LHS
                      → UPDATE_FLOWFIELD
                      → CALC_BOUNDARY
                      → interface/postprocess/residual/output tasks
  → SimuImp::PostProcess
      → FinalizeAccelRuntime
      → HXFinalize
```

关键源码证据：

- 根 CMake 只把 `ThirdParty`、`codes`、`tests`、`gui` 纳入主工程：`CMakeLists.txt:66-71`。
- 主工程通过递归源码收集构建 `OneFLOW`，在 HIP 开启时追加 HIP 源码：`codes/CMakeLists.txt:296-351`。
- runtime 的初始化和释放已经位于主模拟生命周期：`codes/main/src/SimuImp.cpp:51-72,123-131`。
- field solver 主线从 `FieldSimu` 进入 grid、solver、multigrid：`codes/global/src/FieldSimu.cpp:37-45`。
- 多重网格的外层迭代和内层求解位于 `codes/multigrid/src/Multigrid.cpp:111-169,339-387`。
- RK 时间积分的实际 task 序列位于 `codes/solver/src/TimeIntegral.cpp:81-110`。

## 4. 现有组件的职责与缺口

| 组件 | 当前事实 | 可复用部分 | 不能直接承担的职责 |
|---|---|---|---|
| `AccelRuntime` | 选择 backend，初始化设备，提供 allocate/copy/synchronize | 进程级 runtime、rank/device context | solver state、zone/grid 所有权、物理 kernel、restart |
| `AccelBackend` | CPU/HIP 的运行时抽象；CPU 和 HIP 注册存在 | 内存、拷贝、同步、设备选择 | Euler/NS 的 `Advance` 语义 |
| `FluxBackend` | 只有 face flux 和 face-to-cell residual 两个操作 | backend-neutral view 的形状 | 多 RK stage、边界、梯度、粘性、源项、LHS、持久 state |
| `FieldSolver` | toy scalar path 使用 `FluxBackend`；HIP 分支在操作级别拷贝 host 数据 | 可作为“task 中调用 backend”的样板 | 不能代表 NS/Euler 主线；目前会把任意 accelerator 视为 HIP flux backend |
| `UNsInvFlux` | 按 face 循环，依赖全局 `ug/unsf/nscom`，完成重构、Riemann flux、写回 `invflux` | 物理语义和 host oracle | 不是 persistent device execution；当前是逐 face CPU 路径 |
| `UNsUpdate` | 按 cell 循环执行正性修复和状态写回 | 正性/失败语义 | 不适合作为 device state 的热路径 |
| standalone `EulerBackend` | 已有 `CreateState/Upload/Advance/Download`；HIP state 保留 device buffers | 1D Euler vertical slice 的契约和验证方法 | 尚未绑定 OneFLOW `MRField`、zone、grid level、task、restart、输出 |
| `Restart`/output | 以 `MRField`/`DataBook`/host 文件为边界 | 可以定义下载与失效时机 | 不能直接读取 device pointer |
| MPI/interface | 既有 interface 使用 host `DataBook` 和 host field；`CommInterfaceData` 是 task 序列 | host staging fallback、既有 rank/zone 语义 | 当前没有 GPU-aware MPI 或 device halo contract |

### 4.1 runtime 已接入，但 solver 还未接入

`AccelRuntime` 在主 solver 初始化阶段已经被调用，说明设备选择不需要另起一套 main 或 MPI 入口。但 `FieldSimu`/`Multigrid`/`TimeIntegral` 没有根据 `AccelRuntime` 创建或推进 Euler/NS device state。

`SimpleSimu::ToyModelSimu` 中的 scalar path 会根据 `IsAccelerator()` 选择 `HipFluxBackend`（`codes/main/src/SimpleSimu.cpp:33-49`；实际 scalar 操作在 `codes/scalar/src/FieldSolver.cpp:50-58,246-343`）。这只能证明 accelerator seam 在 toy model 中被调用，不能证明 `FieldSimu` 的 NS 主线已经使用 HIP。

另外，当前选择逻辑使用 `IsAccelerator()` 而不是具体 `AccelBackendKind`。当未来 CUDA/KOKKOS 注册后，这会把“任意 accelerator”错误地映射到 HIP-specific flux adapter；接入前必须改成 capability/backend-kind 选择。

### 4.2 standalone port 与主工程的边界

standalone CMake 在 `ports/kunshan/oneflow_1d_hip/CMakeLists.txt:1-60`，而根 CMake 没有 `add_subdirectory(ports/...)`。因此其 `OneDEulerBackend`、`OneDEulerPersistent.hip`、`OneDEulerMpi.*` 不会随主 `OneFLOW` target 自动编译，也不会被 `FieldSimu` 或 `UNsSolver` 调用。

standalone contract 的核心语义是正确的：

```cpp
CreateState(problem)
Upload(state, hostState)
Advance(state, steps, options)
Download(state, hostState)
```

但主工程的真实状态不是单个 `nx` 数组，而是由 `MRField`、ghost/interface field、cell/face metric、`ZoneState::zid`、`GridState::gridLevel` 和 `SolverState::solverIndex` 共同决定。因此不能把 standalone `EulerState` 当作全局 singleton。

## 5. 主 solver 的物理和数据流缺口

当前 NS RK stage 的真实顺序是：

```text
LOAD_RESIDUALS
  → UPDATE_RESIDUALS
      → NsCalcInvFlux
          → gradient / limiter / face reconstruction
          → per-face inviscid flux
          → AddF2CField(res, invflux)
      → NsCalcVisFlux
      → NsCalcSrcFlux
  → CALC_LHS
  → UPDATE_FLOWFIELD
      → per-cell primitive/state conversion
      → positivity repair when needed
  → CALC_BOUNDARY
```

对应证据：

- `codes/multigrid/src/Multigrid.cpp:216-245` 负责 residual load/update 和多重网格 residual 语义。
- `codes/ns/src/NsRhs.cpp:54-67,92-143` 负责 inviscid、viscous、source 组合。
- `codes/uns/src/UNsInvFlux.cpp:123-215` 负责初始化、重构、逐面通量和 residual 写回。
- `codes/uns/src/UNsVisFlux.cpp:90-150` 仍是逐面粘性通量路径。
- `codes/uns/src/UNsUpdate.cpp:46-146` 仍是逐 cell 更新及正性修复。
- `codes/lhs/src/LhsTaskReg.cpp:38-44` 与 `codes/usolver/src/ULhs.cpp:47-65` 负责 LHS 缩放。

这意味着 standalone 一维 Rusanov + SSP-RK3 只覆盖“固定三守恒量、固定一维几何、无粘性/源项、有限边界语义”的小问题。它可作为第一条 Euler vertical slice，但不能直接替换完整 `NsCalcRHS` 或 `UNsUpdate`。

## 6. 建议的 backend-neutral 接入方案

### 6.1 分离两个 contract

保留现有 `AccelBackend` 作为进程级设备资源 contract；新增或抽取一个 solver-domain contract，不把 `CreateState/Advance` 塞进 `AccelBackend`：

```cpp
class EulerExecutionState;

class EulerExecutionBackend
{
public:
    virtual ~EulerExecutionBackend() = default;

    virtual std::unique_ptr<EulerExecutionState> CreateState(
        const EulerProblemView & problem ) = 0;
    virtual void Upload(
        EulerExecutionState & state,
        const HostFieldView & fields ) = 0;
    virtual void Advance(
        EulerExecutionState & state,
        const EulerAdvanceOptions & options ) = 0;
    virtual void Download(
        const EulerExecutionState & state,
        HostFieldView & fields,
        DownloadReason reason ) = 0;
    virtual void Invalidate(EulerExecutionState & state) = 0;
};
```

`EulerProblemView` 不应暴露 `MRField` 内部实现；它应明确给出守恒量数量、cell/face 数、布局、边界/halo、几何量、`dt/dx/volume`、scheme 和 trace/stats 能力。

### 6.2 state 所有权与 key

首版 state manager 的 key 应至少包含：

```text
(solverIndex, localZoneId, gridLevel, backend kind)
```

state manager 的生命周期应覆盖一次完整 `FieldSimu → MultigridSolve`，而不是覆盖一次 face 函数调用。设备 state 的 owning layer 建议挂在 solver/backend integration manager，而不是 `UNsInvFlux` 的临时对象；原因是 `UNsInvFlux` 当前每次 task 都 `new/delete`，且 `CalcFlux()` 内部还会 allocate/deallocate `invflux`。

首个集成范围应严格限制为：

- `NS_SOLVER`、一个 local zone、`gridLevel == 0`；
- 明确的显式 RK/Euler path；
- 先支持 CPU backend 和 HIP/DCU backend；
- 先保留 host boundary/interface/MPI staging；
- 多重网格、粘性、化学源项、湍流和 LUSGS 不在第一版 Euler fast path 中伪装支持。

### 6.3 推荐调用边界

```text
INIT_FLOWFIELD 完成、host q/ghost/interface 已有效
  → backend state CreateState
  → Upload 一次

TimeIntegral::RungeKutta 的显式 Euler fast path
  → Advance(steps/stages, NoTrace 或 FullTrace)

需要 residual/convergence/checkpoint/visualization/force
  → device reduction 或按需 Download 标量/字段

ReadRestart 修改 host q
  → Invalidate + Upload，或明确的 RestoreState

退出 FieldSimu / solver map 释放
  → state manager destroy
```

最合适的 first hook 是 `TimeIntegral::RungeKutta` 的 solver-aware 分支，因为这里已经拥有完整 stage 顺序；不建议把 `Advance` 分散到 `UNsInvFlux`、`ULhs`、`UNsUpdate` 各自的 task 中，否则每个 task 都会重新面对 host/device 同步和 state 轮换问题。

但这个 fast path 必须由能力判断保护：只有 backend 明确声明支持当前 scheme、边界、方程数、viscous/source 配置和 grid layout 时才启用；否则回退现有 CPU task 序列。`ONEFLOW_ACCEL_BACKEND=HIP` 只能表示 runtime 选择 HIP，不能自动宣称 NS Euler fast path 可用。

### 6.4 host field、halo、restart 和输出策略

- `NoTrace`：device 保持主 state，stage 中间量和临时 buffer 不回传；只回传 residual/physical-check 标量以及输出所需最终字段。
- `FullTrace`：只用于小规模 correctness；允许回传 stage state、face state、flux、residual，不能作为性能路径。
- boundary：第一版可在每个 macro step 前后由 host 维护，或为一维 Euler 实现独立 device boundary kernel；必须固定“ghost 先更新还是 flux 前更新”的语义。
- interface/MPI：先使用既有 host field/interface exchange。device state 在 exchange 前做必要 D2H，exchange 后 H2D 或更新 device halo；统计这些搬运，不能把 host staging 隐藏在 kernel 时间里。
- restart：`ReadRestart` 后 host `q` 是新事实源，必须使旧 device state 失效并重新 upload；`DumpRestart` 前必须完成 download。当前 restart 读写围绕 `MRField`/`DataBook`，见 `codes/restart/src/Restart.cpp:74-166` 和 `codes/restart/src/RestartTaskReg.cpp:74-140`。
- output/visualization：继续由 CPU 负责文件格式和报告生成；只在输出边界下载所需 field。

## 7. CPU/HIP/CUDA/KOKKOS 能力边界

| 能力 | CPU | HIP/DCU | CUDA | KOKKOS |
|---|---|---|---|---|
| 主工程默认执行 | 已有 | runtime 可选 | legacy source path | 未启用 |
| `AccelBackend` 注册 | 已有 | 已有 | 当前无注册实现 | 当前无注册实现 |
| scalar `FluxBackend` | 已有 | 已有，当前操作级 host staging | 未验证 | 未验证 |
| standalone 1D Euler contract | 已有 | 已有 standalone HIP state | 无 | 无 |
| 主 NS/Euler solver 集成 | 现有 CPU 主线 | 未集成 | 未验证 | 未验证 |
| device-persistent main-solver state | 不需要 | 设计目标 | 未设计 | 未设计 |
| GPU-aware MPI | 不适用 | 未在此审计中验证 | 未验证 | 未验证 |

特别注意：`codes/CMakeLists.txt:120-215` 的 CUDA 分支负责 legacy CUDA source collection，但 `AccelBackendKind::CUDA` 当前没有对应 `Register`。所以“允许配置 CUDA”与“CUDA backend 可运行”是两个不同结论。

## 8. 分阶段实现与门槛

### 阶段 0：本报告

- [x] 主 solver 调用链和 task 边界
- [x] standalone 与主工程连接缺口
- [x] case/field/restart/residual/MPI 所有权边界
- [x] backend 能力矩阵
- [x] first hook 和不应改动的边界

### 阶段 1：CPU backend-neutral Euler slice

1. 从 standalone contract 提取不含 HIP 的 domain contract。
2. 用主工程的 `MRField`/grid view 构造固定一维 Euler problem。
3. 只接 CPU adapter，保持现有 NS path 默认不变。
4. 增加 direct main-solver test：同一 case 通过主工程 runner 执行，并检查 final state、residual、正性、边界和 restart。

门槛：失败必须返回非零；不能用空 CTest 或 skipped test 作为通过；CPU backend 与独立 oracle 不能只是同一段实现自比。

### 阶段 2：HIP/DCU persistent state

1. 复用 standalone HIP state 的 buffer 轮换和 `NoTrace`/`FullTrace` 语义。
2. 将 state key 绑定到 solver/zone/grid level。
3. 先做一维显式 Euler，无粘性、无源项、无 LUSGS。
4. 记录 allocation、kernel launch、sync、H2D、D2H、device reduction。

门槛：目标 Kunshan 节点上完成 runtime、kernel、one-step numerical、multi-step lifecycle、NoTrace/FullTrace 和最终输出验证；同时检查 Slurm State 与 workload ExitCode。

### 阶段 3：halo/MPI 和 solver 功能扩展

1. 抽象 `PackHalo/ExchangeHalo/UnpackHalo`，先 host staging。
2. 验证 rank/device mapping 和多卡结果。
3. 再处理 device-aware MPI/overlap。
4. 逐项加入 viscous、source、limiter、multigrid 和 restart fast path。

每增加一种物理或 solver 模式，都必须重新建立该模式自己的 CPU oracle、物理检查和输出语义，不从一维 Rusanov 结果推断完整 NS 正确性。

## 9. 代码分流建议

### 可以进入 upstream 的候选

- backend-neutral domain contract 和清晰的 CPU 实现；
- 不改变默认 CPU 行为的主 solver adapter seam；
- 非设备专有的 field view、state lifecycle、错误传播和测试契约；
- 不含集群私有信息的架构/验证文档。

### 先留在 fork DCU 分支的候选

- HIP/DCU kernel 和 DTK-specific CMake glue；
- Kunshan 目标架构、节点能力探测和主 solver HIP adapter；
- host-staged MPI 实验、GPU-aware MPI 探测和性能调优；
- 仍未经过主工程接口稳定性评审的 state manager。

### 只保留本地/CI artifact 的内容

- 原始 Slurm 日志、作业号、节点名、私有路径和环境 dump；
- 仅用于一次性能采样的临时脚本；
- 仍属于 standalone 实验、且没有主 solver 调用方的 benchmark 入口。

## 10. 下一步执行顺序

```text
评审本报告
  → 定义主工程 EulerProblem/FieldView/State contract
  → CPU adapter 接主 solver、主工程 runner 和 restart
  → CPU oracle + physical invariants + CTest/direct-run 门禁
  → Kunshan low-cost probe
  → HIP persistent state 接入
  → NoTrace/FullTrace 与完整输出验证
  → host-staged MPI/halo
  → 再考虑 viscous/source/multigrid/WENO5
```

在阶段 1 的 CPU 主工程 runner 通过前，不继续扩展 WENO5，也不把 standalone benchmark 的性能数字写成 OneFLOW 主 solver 性能。

