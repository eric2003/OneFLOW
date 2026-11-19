# OneFLOW 长期分阶段开发计划（accelerator 主线）

> 2026-11-19 核对。本文维护完整路线、稳定编号和阶段退出门。
> [Living TODO](oneflow-development-todo.md) 维护当前执行切片。
> 集群流程见 [Kunshan CI](../../ci/kunshan/README.md)。

## 0. 结论先行

**根本目标：保证正确性前提下做加速。** 数值、物理与状态/通信语义的正确性是优化的
前置条件；每个优化先明确契约、通过相应正确性门禁，再以同口径、可复现的端到端
测量确认收益。阶段优先级与性能目标均服从这一原则。方程、模型、容差或参考基线的
变化须按 `AGENTS.md` 的数值审查流程独立说明与验证，不能作为性能优化的隐含代价。

OneFLOW 的加速主线已经证明了两件事，长期计划的形状由它们决定：

1. **执行模型比 kernel 重要一到两个数量级。**
   1D standalone 端口：`advance` 内循环（设备驻留、每步 2 次同步）
   CPU 1 线程 vs 1 DCU 为 `102.8x/137.2x/139.1x/155.2x`（65k→4M cells）；
   同一套 kernel 走 trace/host-staged 路径只有 `2.0–2.6x`；
   应用级 lifecycle（含网格、初始化、IO、MPI）相对 32 rank CPU 为
   `3.85–9.70x`，4 卡/4M 为 `25.55x`。
   三者差两个数量级，差别全在"谁拥有状态、每步回几次 host"。
2. **3D 主求解器的完整 HIP `Advance` 尚缺，Phase B 已开始独立验证。** A.1 已收口，
   历史同节点 paired 为 14.62 ms/step（8 zone / 8 rank / 1 DCU / 1000 步）；
   不能直接与不同节点、步数或启动口径下的较低墙钟比。短程 rocprof 显示 host stage
   计时不能当 GPU 时间。H1 已完成归账；H2 已在两节点两卡复验。A.4 仍需先归因同步和传输成本。

因此长期计划的主轴不是继续抠 kernel，而是：

> **Phase A 把 3D 改成设备驻留执行（device-resident Advance）→
> Phase B 补齐数值能力（scheme/limiter/粘性/SA/隐式/multigrid）→
> Phase C 真实算例与 IO → Phase D 多卡多节点 →
> Phase E 架构通用化 → Phase F 产品化与发布。**

每个阶段都有"准入门"（oracle/门禁）和"退出门"（可发布证据），见 §3。
完整任务编号和验收条件以本计划为准。[Living TODO](oneflow-development-todo.md) 只维护当前切片和近期验收。

## 1. 代码架构现状（长期计划的事实基础）

### 1.1 分层与规模

| 层 | 目录 | 规模（cpp/hip + h） | 职责 | 设备状态 |
|---|---|---:|---|---|
| 入口/调度 | `codes/main`（1.9k）、`codes/solver`（3.8k）、`codes/task`（3.9k）、`codes/register`（1.3k） | ~11k | `Simulation → SimuImp → FieldSimu → Multigrid → TimeIntegral`，task 字符串注册与派发 | host |
| 数据/字段 | `codes/database`（2.8k）、`codes/global`（4.8k）、`codes/adt`、`codes/array` | ~9k | `DataBase`/`DataBook`、`MRField`、`FieldHome`/`FieldManager`、变量名映射 | host（设备侧另有影子缓冲） |
| 网格/几何 | `codes/geometry`（18.4k）、`codes/zone`（1.1k）、`codes/partition`（1.5k）、`codes/cgns`（8.7k） | ~30k | 拓扑、面-单元关系、BC 记录、METIS 分区、CGNS IO、wall distance | host |
| 物理/离散 | `codes/ns`（4.9k）、`codes/uns`（4.7k）、`codes/ins`（3.0k）、`codes/uins`（4.0k）、`codes/turb`（3.9k）、`codes/uturb`（3.6k）、`codes/lhs`、`codes/grad`（+`usolver`/`uns` 侧） | ~28k | NS/INS/湍流三族 solver 的 RHS、边界、LUSGS、梯度、通量 | 受限 3D inviscid 链已 device |
| 时间/迭代 | `codes/time`、`codes/multigrid`（0.9k）、`codes/unsteady`（0.9k）、`codes/implicit` | ~3k | 显式 RK、LUSGS、多重网格、双时间 | host |
| 加速器 | `codes/accel`（10.6k）、`codes/cuda`（0.36k） | ~11k | `AccelRuntime/AccelBackend/DeviceBuffer`、`FluxBackend`、`EulerDomain`、`ResidentStateContract`、`StageProfiler`、HIP kernel | **唯一 device 拥有者** |
| 后处理/IO | `codes/restart`（0.6k）、`codes/visual`（3.4k）、`codes/force`（1.0k）、`codes/residual`（0.8k）、`codes/io`（2.9k）、`codes/postprocess`（0.9k） | ~10k | restart、可视化、气动力/Cp、残差、文件 | host |
| 独立端口 | `ports/kunshan/oneflow_1d_hip`、`codes/scalar`（8.9k） | — | 1D Euler 驻留执行样板；scalar toy path | device |
| 网格生成 | `codes/special`（11.2k） | 11k | 块结构网格生成、曲线/曲面 | host，**不在加速范围** |
| 测试/CI | `tests/`、`test/`、`ci/kunshan/`、`tools/oneflow_gui` | — | GoogleTest 契约、五算例回归、多区 trace/XB/稳定性 runner、GUI | — |

### 1.2 真实执行链

```text
Simulation::Run
  → SimuImp::PreProcess            # MPI/设备/控制参数
  → FieldSimu
      → LoadGrid / Partition / CreateSolvers / INIT_FLOWFIELD
      → MultigridSolve             # 外层迭代 + grid level
          → TimeIntegral::Relaxation
              → RungeKutta / Lusgs / Simple   # host task 序列
                  LOAD_Q → CALC_TIME_STEP
                  [每个 stage] LOAD_RESIDUALS → UPDATE_RESIDUALS
                               → CALC_LHS → UPDATE_FLOWFIELD → CALC_BOUNDARY
                               → UPDATE_INTERFACE_DATA
                  post: UPLOAD/UPDATE/DOWNLOAD_INTERFACE, 残差/气动力/输出
  → SimuImp::PostProcess           # 收尾、restart、teardown
```

`UPDATE_RESIDUALS` 内部：`NsCalcRHS = inviscid flux（梯度/limiter/重构/Riemann）+ viscous flux + source`；
task 通过全局 `ug / unsf / nscom / ctrl / ZoneState::zid` 隐式定位到"当前 zone/grid"。

### 1.3 三个并存的"世界"

| 世界 | 触发 | 状态所有者 | 数值路径 | 用途 |
|---|---|---|---|---|
| legacy host | 默认 | 全局 `MRField` | 逐面/逐 cell host 循环 | CPU oracle、全部未迁移能力 |
| CPU batch | `ONEFLOW_ENABLE_UNS_CPU_BATCH=1` | 同上 + `EulerDomain` CPU state | 批量 adapter 调 `CpuFluxBackend` | 批量路径的等价性验证 |
| HIP device | `ONEFLOW_ACCEL_BACKEND=HIP` + capability | `Ns3DDeviceState` + `ResidentStateContract` | 每 stage 一批 kernel，host 调度 | 受限 3D inviscid 快路径 |

三者靠 `EulerInvFluxCapability`、`EulerRungeKuttaCapability`、`HipProductionPolicy`
准入，靠 3 步多区 trace（1e-10）+ 200 步跨后端残差（1e-12）+ CPU 五算例（1e-8/1e-15）
保证等价。

### 1.4 已有 seam（可复用的骨架）

| seam | 位置 | 现状 |
|---|---|---|
| `AccelRuntime` / `AccelBackend` / `DeviceBuffer` | `codes/accel` | 成熟：设备选择、拷贝、同步、RAII |
| `EulerDomainBackend::{CreateState,Upload,Advance,Download}` + `StateRegistry` | `codes/accel/include/EulerDomain.h` | **3D 只接了 CPU**；HIP 没有 `Advance` |
| `EulerRungeKuttaCapability`（含 `stageCallback`） | `codes/solver` | 已完成接入点：`TimeIntegral::RungeKutta(context)` 里 `backend.Advance(state,1,options)` |
| `EulerInvFluxCapability` / `HipProductionPolicy` | `codes/accel` | 逐 stage 准入 + 回退 |
| `ResidentStateContract` | `codes/accel` | 设备/主机数据所有权与代际（本轮新增 host-mirror currency） |
| `StageProfiler` | `codes/accel` | opt-in 分段计时（含 host 子阶段） |
| 1D standalone 驻留端口 | `ports/kunshan/oneflow_1d_hip` | 已证明驻留执行模型 |
| `FluxBackend`（face flux / residual） | `codes/accel` | 只覆盖 face 级操作，不足以承载 stage |

### 1.5 结构性问题（长期计划要逐个消除）

1. **ambient globals 决定热路径**：`ug/unsf/nscom/ctrl/Iteration/ZoneState::zid`
   让"哪个 zone 的哪份 device state 是当前的"靠隐式约定；本轮 device-authoritative
   工作反复踩到（初始化上传、`ug.Init()` 重绑定、contract 失效）。
2. **host 调度 + 逐 stage kernel**：每个 stage 至少一次 kernel fence、一次 task 派发、
   若干次 host 侧整场循环；这正是 3D 固定开销占比高的原因（在 31.6 ms/step 口径下
   测得 kernel 只占约 1.5 ms；A.1 收口后总时间降到 15.6，这个结构性问题仍在）。
3. **capability 判断分散**：`HipProductionPolicy` 的开关 + 各文件里的 `UseHip*`
   守卫 + 两个 capability 枚举；能力扩展时会组合爆炸。
4. **双实现与证据覆盖不完整**：已有 Roe 与 Barth 单点数值契约，但 Roe 缺多区证据，
   Barth 整链数值未过门禁；host oracle、单点契约与整链 trace 都需保留。
5. **fixture 与多卡证据缺口**：8 zone M6 基网格为 294,912 cell；原生加密与 METIS 分区已生成 ≥1M 网格。
   D.1 待固定正式验收用例的版本、读入和接口证据。
   `partition_type=0`(GenerateLink) 仍有缺口，H2 已在两节点两卡配置复验，其他多卡配置按各自门禁验证。
6. **口径/文档债**：F.1 3D 正式报告已于 2026-11-02 发布；新优化不能沿用旧 CPU
   或跨节点数字拼接加速比。当前源码的门禁与性能仍需绑定同一修订的可复验证据。

## 2. 阶段划分

依赖关系：

```text
Phase A（执行模型）
   ├── Phase B（数值能力）── Phase C（真实算例/IO）
   └── Phase D（多卡）──────┘
                 Phase E（架构通用化）, 与 A/B 并行、按痛点触发
                        Phase F（产品化/发布），每个阶段末都要做一遍
```

### Phase A — 3D 设备驻留执行（显式执行模型）

目标：让 3D 显式路径从"host 调度 + 逐 stage kernel"变成"一个宏步一次进入设备"，
把 host 往返和 fence 从每步几十次降到 ≤2 次。

| 子项 | 内容 | 预期 | 证据/入口 |
|---|---|---|---|
| A.1 | 固定开销优化已完成：按需镜像、设备残差基底、跳 host 历史层与账本、延迟同步。基准须 warmup | 历史同节点 paired **14.62 ms/step**；历史归账见 Living TODO 的证据索引 | `UNsInvFlux::CalcInvFace`、`FieldTaskReg::LoadQ/LoadResiduals`、`ResidualBasePolicy` |
| A.2 | 原每步一次交换方案暂缓；保留 legacy 滞后，不能直接删 post 轮（200 步 XB `5.54e-7`）。聚合与减轮次待 A.4 后重估 | 单卡通信开关差值约 2.2 ms，不是可兑现收益保证 | `InterfaceTask`、`UNsGrad::CalcGradHip`、`ExchangeHipInterfaceQ` |
| A.3 | stage kernel 融合：gradient→recon→flux→lhs→update→boundary 在设备内串联，减少 host 同步 | 已降权；旧 −4/−6 ms 估算失效，待配对测量 | `HipFluxBackend` 各 kernel + `TimeIntegral` |
| A.4 | **HIP `EulerDomainBackend::Advance`**：把 RK stage 循环搬进 backend，host 负责 INIT/输出/门控；宏步末按需 Download | H1 已排除 MG/派发大项；先定位同步/传输成本再量化收益，旧 9–10 ms 估算撤回 | `TimeIntegral::RungeKutta(context)`、`EulerDomain.h`、`HipSmoke` 既有 `Advance` 样板 |

退出门：
- 8 zone/8 rank/1 DCU 1000 步 paired **≤14 ms/step**（相对 10-25 的 31.6 再快 2.2x，
  相对 10-24 的 38.3 快 2.7x）；
- 3 步多区 trace（1e-10）、200 步跨后端（1e-12）、CPU 五算例（1e-8/1e-15）、
  50 步稳定性全绿；
- 每宏步 host 往返次数与 H2D/D2H 字节数进入 `StageProfiler` 报告。

### Phase B — 数值能力（决定这条后端是否有实用价值）

当前区分三条路径：

- **无粘 HIP batch**。Lax-Friedrichs 和 Roe 是候选通量，还要满足运行时、网格、
  方程数和接口交换等能力守卫。组合不受支持时，显式强制 batch 会报错，自动选择时可回落。
- **设备梯度与显式驻留执行**。设备梯度要求无粘模式、HIP batch，以及关闭限制器或使用
  Barth；多区还要求设备 interface。历史 Barth 对照差约 `1.1e-6`，超过当时的
  `1e-12` 门槛。关闭限制器仍是该验收基线的配置。
- **独立的粘性与 SA hook**。`UNsVisFlux::CalcVisFluxHip` 要求五方程、VIS_STD、
  非等温壁和必要字段。条件不满足时走主机粘性循环。
  `test/cfd/m6wingroe_sa` 使用 Roe、LU-SGS 和 SA；部分 hook 可命中，整算例设备
  端到端尚未验收。已核对状态（`cc1dd4b3a`，单 rank/单 DCU）：全关与仅开粘性各
  400 步的 `aero.dat`、`res.dat` 逐字节一致；全开 hook、收敛后气动力和多区仍未验收。
  写回门禁以 `SELFCHECK result=PASS|FAIL` 判定：先统计非有限值，再比相对误差，
  并同时看面归一化差与绝对差；只有相对误差的旧日志不能证明通过。
  当前切片与证据见 [Living TODO](oneflow-development-todo.md)。

| 子项 | 内容 | 依赖 | 关键文件 |
|---|---|---|---|
| B.1 | 无粘通量 scheme：LF → Roe / HLLC / AUSM 系列（device），LF 保留作 oracle | A | `UNsInvFlux`、`NsInvFlux`、`HipFluxBackend` flux kernel |
| B.2 | limiter：Barth / Venkat（cell 邻域极值），解除 `ilim=0` | B.1 | `UNsLimiter`/`NsLimiter`、recon kernel |
| B.3 | 层流粘性：`UNsVisFlux`（面粘性通量、T 梯度、壁面热流） | B.1/B.2 | `UNsVisFlux`、gradient（多方程）、BC 热流 |
| B.4 | 湍流：SA（1 方程）优先，SST 可选；`turb` solver 与主 solver 的 task 交互 | B.3 | `codes/turb`、`codes/uturb`、`system/turb` 字段 |
| B.5 | 隐式：LU-SGS 设备化（下/上扫、`dq` 交换、与 halo 的依赖），再评估双时间 | A.4 | `NsLusgs`、`ULhs`、`LhsTaskReg` |
| B.6 | 多重网格：restrict/prolong 的设备化，或显式声明"只加速 finest grid，粗网格留 host" | A.4 | `codes/multigrid` |

准入门（每个子项都必须先过）：
1. 单点数值契约（scheme/limiter 的 CPU vs device 单元测试，≤1e-14）；
2. 3 步多区 trace（1e-10）；
3. 该能力对应的物理门禁（正性、守恒、边界语义）；
4. 200 步跨后端 + 50 步稳定性。

退出门：`m6wingroe_sa`（Roe+LU-SGS+SA）在设备上跑通并通过同一组门禁。
这个算例在本计划中记为 **B.7**，是 Phase B 的退出条件，做完 B.1–B.6 不等于 B 结束。
之后才允许把 "HIP 支持的算例范围" 写进正式文档。

### Phase C — 真实算例与 IO（把加速变成可交付能力）

| 子项 | 内容 | 说明 |
|---|---|---|
| C.1 | IO/restart：设备输出或"整场一次性 D2H + 异步写盘" | 现在每输出步 `NsRefreshHostStateMirror` 整场回拉，长跑/多输出会吃掉收益 |
| C.2 | 气动力/Cp/热流：设备侧面量归约，或维持边界 cell 级小回拉 | 现在每步一次边界 cell 镜像，可接受但需量化 |
| C.3 | 残差/收敛判据设备归约 | 已在设备侧算 metrics，需扩到多 zone/多 solver 聚合 |
| C.4 | 内存预算与大网格 | 记录每 100 万 cell 的 device 常驻字节数，给出单卡可解规模 |
| C.5 | 定常/非定常的真实工况回归 | 与 `test/suites` 对齐，形成"设备算例集" |

退出门：≥1M cell 的单卡算例端到端可跑（含 restart 等价、气动力对比、可视化），
并有对外可引用的 paired 性能数字。

### Phase D — 多卡与多节点（结论的前提是 fixture）

| 子项 | 内容 |
|---|---|
| D.1 | **fixture**：从已有 M6 网格阶梯固定 ≥1M cell 的验收用例，核对读入与接口；另跟踪 `GenerateLink` 路径 |
| D.2 | 设备 halo：GPU-aware MPI（1D 已验证 `1.39/1.31/1.19/1.13x`）在 3D 的复现；DeviceDirect 在 packed layout 上重测 |
| D.3 | 映射与负载：1 rank = 1 zone 的约束是否长期成立；需要 zone 切分/多 zone per rank 时的 state 与 halo 语义 |
| D.4 | 评测口径：固定 per-rank cell 数、同一 repeats、1 DCU 槽 = 8 核；给出 8:1、32:1 的合法性说明 |

退出门：至少 4 卡、per-rank 恒定 cell 数的端到端数字，且与同节点 CPU 同核数配对。

### Phase E — 架构通用化（与 A/B 并行，按痛点触发）

| 子项 | 内容 |
|---|---|
| E.1 | 热路径去 ambient globals：把 zone/grid/solver 身份显式传进 stage（先 device 路径，再 host） |
| E.2 | capability 表驱动：`solver × scheme × limiter × viscous × time-integral × grid-level × halo` 一张表 + 单测，替换散落守卫 |
| E.3 | `EulerDomain` 泛化到 NS：3/5 方程已支持，补源项/粘性/能量扩展与 view 元数据 |
| E.4 | backend-kind 分派：清掉 `IsAccelerator()` 当 HIP 用的地方，为 CUDA/Kokkos 留真实入口 |
| E.5 | 清理 toy/scalar 路径与失效文档，保持 CPU oracle 单一 |

退出门：新增一种 scheme 或一种物理项，只需"新 kernel + 表项 + 单测"，不需要改
`TimeIntegral`、`UNsInvFlux`、`InterfaceTask` 的骨架。

### Phase F — 产品化与发布（每个阶段末循环一次）

| 子项 | 内容 |
|---|---|
| F.1 | 3D 正式性能报告（Markdown + HTML 同 stem）：paired、同节点、同 repeats、含环境口径 |
| F.2 | CI 矩阵：CPU normal/strict、DCU 单卡、DCU 多卡、nightly 长跑；PR 分支自跑门禁 |
| F.3 | 文档：`docs/plans/oneflow-development-todo.md` 当前切片、`ci/kunshan/README.md`、Sphinx 教程 |
| F.4 | 交付：GUI/示例算例、发行说明、对外能力矩阵（哪些 solver/scheme 在设备上已验证） |

## 2.1 编号和验收拆分

父项用 `A.1`，子项用 `A.1.1`。编号在 Living TODO 中原样复用。每个父项的子项和对应阶段退出门均保留独立证据。旧记录的 `B3` 对应 `B.3`，`B-exit` 对应 `B.7`。H1/H2 是历史调查编号。

### Phase A

- **A.1 按需镜像、设备残差归零和历史层裁剪**
  - **A.1.1** 核对末 stage 梯度镜像、设备残差基底和跳过主机历史层的条件。
  - **A.1.2** 复现同配置 paired 基线，记录逐步字节与同步次数。
- **A.2 保持 halo 滞后并优化交换**
  - **A.2.1** 记录各 RK stage 的 halo 交换方向、次数、状态时间层和字节数。
  - **A.2.2** 先锁定阶段前后交换语义，再测试聚合与轮次削减。
- **A.3 融合设备 stage**
  - **A.3.1** 用逐 stage 基线记录 kernel launch、fence 和中间字段读写。
  - **A.3.2** 从梯度、重构和通量开始，按 trace 和 paired 结果决定是否扩展。
- **A.4 在后端执行一个完整 RK 宏步**
  - **A.4.1** 归因宏步尾状态回读、清零与 RHS 队列等待。
  - **A.4.2** 保留物理状态检查及 halo 滞后，下沉 RK 阶段循环。
  - **A.4.3** 以 trace、200 步跨后端、CPU 五例和同配置 paired 验收。

### Phase B

- **B.1 无粘通量**
  - **B.1.1** 对 Roe 复核单面契约、单 zone 证据和多区差异；保留 LF 参考。
  - **B.1.2** Roe 补多区 trace 和 CPU 五例，再分别评估 HLLC 与 AUSM。
- **B.2 限制器**
  - **B.2.1** 定位 Barth 邻域极值与 φ 的首次差异，检查接口和边界 stencil。
  - **B.2.2** 对齐 Barth stencil、边界和浮点口径，再评估 Venkat。
- **B.3 层流粘性**
  - **B.3.1** 验证方程平面 stride 和求解器写回，记录首次差异。验收以写回字段的
    逐面自检 `SELFCHECK result=PASS|FAIL` 为准：先判非有限值，再比相对误差，
    并同时报告绝对差；近零面不能只凭相对数字追查。
  - **B.3.2** 复跑 M6 全关/仅粘性 400 步，再验收收敛后的 CL/CD/CD_SF。
  - **B.3.3** 补壁面热流、非支持选项的回落和多区边界契约。
- **B.4 SA 湍流**
  - **B.4.1** 固定涡粘度、源项和迎风通量各自的生产入口及字段快照。
  - **B.4.2** 分别验收涡粘度、源项和通量，再检查与主求解器的耦合；SST 独立评估。
- **B.5 隐式 LU-SGS**
  - **B.5.1** 对照主机下扫、上扫的顺序和依赖，记录 dq 交换点。
  - **B.5.2** 验证下扫、上扫和 dq halo，再评估双时间；依赖 A.4。
- **B.6 多重网格**
  - **B.6.1** 记录细网格与粗网格的字段所有权及 restrict/prolong 路径。
  - **B.6.2** 验证 restrict/prolong，或确定 finest-only 能力边界。
- **B.7 旗舰退出门**
  - **B.7.1** 旗舰 M6 通过 B 类数值门禁。
  - **B.7.2** 记录设备命中、CPU oracle、适用范围及复现命令。

### Phase C

- **C.1 输出与 restart**
  - **C.1.1** 记录 q 与 restart 历史层的回拉和落盘边界。
  - **C.1.2** 验证历史层恢复，量化长跑和高频输出代价。
- **C.2 气动力与壁面输出**
  - **C.2.1** 核对壁面 owner/ghost、梯度与热流的消费阶段。
  - **C.2.2** 以 CL/CD/Cp/热流验证镜像完整性，再对比设备归约成本。
- **C.3 残差与收敛判据**
  - **C.3.1** 核对单 zone 设备 metrics 的参考值及数据来源。
  - **C.3.2** 验证多 zone、多 solver 的归约和停止语义。
- **C.4 内存预算**
  - **C.4.1** 按几何、状态、halo 缓冲和临时空间拆分设备内存。
  - **C.4.2** 统计每百万 cell 常驻字节与峰值，给出可解规模。
- **C.5 算例集**
  - **C.5.1** 从既有 CPU 套件选择设备支持的工况及参考输出。
  - **C.5.2** 覆盖定常与非定常，记录可视化输出和同基准性能。

### Phase D

- **D.1 分区 fixture**
  - **D.1.1** 选定并固定 ≥1M M6 fixture 的版本与分区路径。
  - **D.1.2** 验证拓扑、守恒、读入和接口；GenerateLink 缺口单列证据。
- **D.2 设备 halo**
  - **D.2.1** 以 host-staged 路径建立接口数值和传输基线。
  - **D.2.2** 分开核对 host-staged、GPU-aware MPI 和 packed DeviceDirect 的正确性与性能。
- **D.3 rank/zone 映射**
  - **D.3.1** 核对每 rank 的 zone、cell、halo 大小和设备绑定。
  - **D.3.2** 先验收一对一的限制与负载，再决定多 zone 支持。
- **D.4 扩展性口径**
  - **D.4.1** 固定同节点对照和启动扣除方法，记录采样次数。
  - **D.4.2** 固定每 rank cell、节点、repeats 与 CPU 核数，再归档 paired 数据。

### Phase E

- **E.1 显式上下文**
  - **E.1.1** 列出热路径直接使用的 zone、grid、solver 全局状态。
  - **E.1.2** 先梳理设备侧的 zone/grid/solver 所有权，再迁移主机热路径。
- **E.2 能力表**
  - **E.2.1** 盘点现有 capability 枚举、回落条件和拒绝原因。
  - **E.2.2** 集中能力组合和稳定拒绝原因，再覆盖允许与拒绝矩阵。
- **E.3 NS domain**
  - **E.3.1** 盘点 NS 粘性、源项和能量所需的 view 及所有权。
  - **E.3.2** 扩展源项、粘性与能量字段及生命周期。
- **E.4 后端分派**
  - **E.4.1** 找出把任意加速器视为 HIP 的调用点。
  - **E.4.2** 用显式 backend-kind 选路径，分别验证 CPU/HIP 与新增后端。
- **E.5 清理重复路径**
  - **E.5.1** 查清 toy/scalar 路径的调用者和参考用途。
  - **E.5.2** 核对调用和行为，再移除闲置路径与失效说明。

### Phase F

- **F.1 正式性能报告**
  - **F.1.1** 固定报告里的硬件、配置、启动口径和完整原始证据。
  - **F.1.2** 同节点、同 repeats 的 paired 结果及 Markdown/HTML 报告。
- **F.2 CI 矩阵**
  - **F.2.1** 明确分支触发规则及精确 SHA 的证据要求。
  - **F.2.2** 先落实 CPU normal/strict 和单卡，再加入多卡及 nightly。
- **F.3 能力文档**
  - **F.3.1** 核对 capability、README、CI 和教程中的配置入口。
  - **F.3.2** 同步构建入口、教程和已验证能力范围。
- **F.4 发布交付**
  - **F.4.1** 列出可发布的 GUI、示例、输入和验证证据。
  - **F.4.2** 联跑 GUI、示例和发行清单，注明数据与证据版本。

### 已验收里程碑

- **A.1** 已收口。历史证据为 8 zone / 8 rank / 1 DCU 的 1000 步 paired 14.62 ms/step。
- **F.1** 3D 正式报告已发布。后续版本按相同口径重新归档。
- H2 两节点两卡复验属于 D.2 的历史证据，D.2 的其他传输模式和配置仍按各自门禁验证。

## 3. 里程碑与判定口径

| 里程碑 | 目标（8 zone/8 rank/1 DCU，1000 步，paired） | 能力 | 判定 |
|---|---|---|---|
| M0（已达成，2026-10-25） | 31.6 ms/step（≈6.4x @8:1） | 稳态 + inviscid + LF + ilim=0 | trace 24/0/0、XB 2.33e-18、CPU 5/5、PR 已合并 |
| M1（Phase A） | ≤14 ms/step（≈14–16x） | 同上，但宏步驻留 | A 退出门全部通过 |
| M2（Phase B） | 保留 M1 性能的 70% 以上 | Roe/HLLC、limiter、层流粘性、SA | `m6wingroe_sa` 设备端通过门禁 |
| M3（Phase C） | 1M cell 单卡可交付 | + IO/restart/气动力 | 设备算例集 + 正式报告 |
| M4（Phase D） | 4 卡 ≥3x 单卡（per-rank 恒定） | 多卡 halo | 合法 fixture + 配对数据 |
| M5（Phase E/F） | 不追性能 | 通用化 + 发布 | 新能力只需表项 + kernel |

口径纪律（沿用 `AGENTS.md`）：
- 只比同节点、同 repeats、同 fixture、同 rank/资源；`lifecycle_*_ms` 是 repeats 总和；
- 不以 100 步或 kernel-only 作生产结论；1 DCU 槽 = 8 核；
- 探索性单作业数据只进 todo，不进正式报告。

## 4. 风险与对策

| 风险 | 影响 | 对策 |
|---|---|---|
| 双实现漂移（host oracle 与 device 结果缓慢分叉） | 数值不可信 | 每个能力先做单元级数值契约；trace/XB/五算例是准入门而不是事后验证 |
| ambient globals 阻碍驻留执行 | Phase A.4 无法干净落地 | E.1 提前做 device 路径的显式上下文；host 侧渐进 |
| ≥1M fixture 的验收基线未固定 | Phase D 证据难复跑 | 固定现有网格阶梯的一套 fixture；单独记录 GenerateLink 路径缺口 |
| 隐式/multigrid 复杂度高 | 工期失控 | 先把显式驻留做到 14 ms；隐式与 multigrid 各自独立立项、独立门禁 |
| capability 组合爆炸 | 维护成本 | E.2 表驱动 + 单测；每个组合要么有门禁要么默认拒绝 |
| 性能口径漂移（跨会话跨节点比较） | 结论反复 | 只发 paired；报告里写清节点/会话限制 |
| 上游同步与 PR 规模 | 评审不可行 | 每轮三端同步 master、按"backend/stage/CI"小切片提 PR（#227 的教训） |

## 5. 与既有文档的关系

| 文档 | 已有结论 | 本文新增 |
|---|---|---|
| `docs/plans/oneflow-development-todo.md` | 当前进行中的任务与验收证据 | 阶段取舍以本文为准，勾选与证据以 todo 为准，编号与本文一致 |
| `doc/architecture/oneflow-1d-euler-dcu-architecture-audit-20260904.md` | backend-neutral 接入方案、state key、调用边界 | 把该方案升格为长期阶段主轴（Phase A.4） |
| `docs/plans/oneflow-euler-optimization-plan.md` | "stateful NoTrace 才是性能路径"（trace 2x vs 驻留 100x+） | 在 3D 复现同一结论：执行模型优先于 kernel |
| `docs/runbooks/oneflow-dcu-porting-playbook-20260903.md` | 移植方法论（契约→或修→门禁→benchmark） | 作为每个 Phase 子项的验收流程模板 |
| `docs/handoff/oneflow-current-handoff.md` | 环境、固定路径、启动顺序、坑（living，每次会话就地更新） | 其 §4 任务指针对应 Living TODO 当前切片 |

## 6. 需要拍板的开放问题

1. **目标规模**：这条后端的"完成"定义是 1M cell 单卡、还是 100M cell 多卡？
   （决定 Phase C/D 的投入比例）
2. **旗舰算例优先还是性能优先**：原计划先达到 Phase A 性能目标。当前已记录 M6 气动力偏差，
   先验收 B.3 的正确性，再按 A.4 的归因结果决定性能工作。
3. **隐式与多重网格是否在设备上做**：若只加速显式 finest grid，可省掉 Phase B.5/B.6
   的大半工作，但要接受"定常收敛用 LUSGS 就没加速"的现实。
4. **zone 映射**：长期是否维持 1 rank = 1 zone？若需要 zone 切分，Phase D 的工作量
   翻倍（state、halo、load balance 全变）。
5. **CUDA/Kokkos 是否真要做**：当前不做。当前继续走 HIP 生产路径。E.4 决定 CUDA/Kokkos 的后端入口，
   `EulerDomain` 继续保持后端类型独立。

## 7. 当前切片入口

当前执行顺序和验收证据见 [Living TODO](oneflow-development-todo.md)。
