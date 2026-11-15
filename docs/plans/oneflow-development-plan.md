# OneFLOW 长期分阶段开发计划（accelerator 主线）

> 最后核对：2026-11-11；本轮起点 `1d283c31a`，求解器证据修订 `6f7a49e3d`。
> 本文是长期路线图（living），按季度/阶段维护；
> **执行进度、勾选与证据以
> [`oneflow-development-todo.md`](oneflow-development-todo.md) 文首
> §4 的 Phase A–F 勾选树为准**（编号与本文 §2 一致），本文只说明取舍与阶段边界；
> 证据摘要记在 todo，交接文档只保留任务指针与操作约束，原始数据留在集群证据目录。
> 规则以根目录 `AGENTS.md` 为准；集群流程见 `ci/kunshan/README.md`。

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
2. **3D 主求解器仍在 Phase A，完整 HIP `Advance` 尚缺。** A1 已收口，
   历史同节点 paired 为 14.62 ms/step（8 zone / 8 rank / 1 DCU / 1000 步）；
   不能直接与不同节点、步数或启动口径下的较低墙钟比。短程 rocprof 显示 host stage
   计时不能当 GPU 时间。H1 已完成归账与 CI 收尾；下一开发项为 H2 多卡 VMFault，
   之后按同步/传输成本评估 A4，不能把 task 派发当作数毫秒瓶颈。

因此长期计划的主轴不是继续抠 kernel，而是：

> **Phase A 把 3D 改成设备驻留执行（device-resident Advance）→
> Phase B 补齐数值能力（scheme/limiter/粘性/SA/隐式/multigrid）→
> Phase C 真实算例与 IO → Phase D 多卡多节点 →
> Phase E 架构通用化 → Phase F 产品化与发布。**

每个阶段都有"准入门"（oracle/门禁）和"退出门"（可发布证据），见 §3。
完成/待办 checkbox 统一维护在 [`oneflow-development-todo.md`](oneflow-development-todo.md)
§0.1；本计划不复制勾选树。当前顺序为 **H2 → A4 → A2′/A3**，具体验收以 todo §5 为准。

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
   测得 kernel 只占约 1.5 ms；A1 收口后总时间降到 15.6，这个结构性问题仍在）。
3. **capability 判断分散**：`HipProductionPolicy` 的开关 + 各文件里的 `UseHip*`
   守卫 + 两个 capability 枚举；能力扩展时会组合爆炸。
4. **双实现与证据覆盖不完整**：已有 Roe 与 Barth 单点数值契约，但 Roe 缺多区证据，
   Barth 整链数值未过门禁；host oracle、单点契约与整链 trace 都需保留。
5. **fixture 与多卡证据缺口**：固定登记的 fixture 是 8 zone m6wing（294,912 cell）；
   16/32 zone 派生配方与 manifest 校验已跑通，但尚缺 ≥1M 的合法 fixture。
   `partition_type=0`(GenerateLink) 仍有缺口，多卡另受 H2 VMFault 阻塞。
6. **口径/文档债**：F1 3D 正式报告已于 2026-11-02 发布；新优化不能沿用旧 CPU
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

### Phase A — 3D 设备驻留执行（下一个迭代，最高优先）

目标：让 3D 显式路径从"host 调度 + 逐 stage kernel"变成"一个宏步一次进入设备"，
把 host 往返和 fence 从每步几十次降到 ≤2 次。

| 子项 | 内容 | 预期 | 证据/入口 |
|---|---|---|---|
| A1 | 固定开销优化已完成：按需镜像、设备残差基底、跳 host 历史层与账本、延迟同步。基准须 warmup | 历史同节点 paired **14.62 ms/step**；本轮归账的新口径见 todo §5.2 | `UNsInvFlux::CalcInvFace`、`FieldTaskReg::LoadQ/LoadResiduals`、`ResidualBasePolicy` |
| A2 | 原每步一次交换方案暂缓；保留 legacy 滞后，不能直接删 post 轮（200 步 XB `5.54e-7`）。A2′ 聚合/减轮次待 H1/H2/A4 后重估 | 单卡通信开关差值约 2.2 ms，不是可兑现收益保证 | `InterfaceTask`、`UNsGrad::CalcGradHip`、`ExchangeHipInterfaceQ` |
| A3 | stage kernel 融合：gradient→recon→flux→lhs→update→boundary 在设备内串联，减少 host 同步 | 已降权；旧 −4/−6 ms 估算失效，待配对测量 | `HipFluxBackend` 各 kernel + `TimeIntegral` |
| A4 | **HIP `EulerDomainBackend::Advance`**：把 RK stage 循环搬进 backend，host 负责 INIT/输出/门控；宏步末按需 Download | H1 已排除 MG/派发大项；先定位同步/传输成本再量化收益，旧 9–10 ms 估算撤回 | `TimeIntegral::RungeKutta(context)`、`EulerDomain.h`、`HipSmoke` 既有 `Advance` 样板 |

退出门：
- 8 zone/8 rank/1 DCU 1000 步 paired **≤14 ms/step**（相对 10-25 的 31.6 再快 2.2x，
  相对 10-24 的 38.3 快 2.7x）；
- 3 步多区 trace（1e-10）、200 步跨后端（1e-12）、CPU 五算例（1e-8/1e-15）、
  50 步稳定性全绿；
- 每宏步 host 往返次数与 H2D/D2H 字节数进入 `StageProfiler` 报告。

### Phase B — 数值能力（决定这条后端是否有实用价值）

现状：设备路径接受 `MULTI_STAGE + (Lax-Friedrichs 或 Roe) + ilim∈{0, Barth} + vismodel=0`。
Barth 的 φ 已在设备上算，但数值还差 `1.1e-6` 未过 1e-12 门禁，因此默认仍是 `ilim=0`；
Venkat / HLLC / AUSM 未实现，显式请求直接 fail-fast。
仓库自带的旗舰算例 `test/cfd/m6wingroe_sa` 是 `LU-SGS + Roe + SA`，**完全在 host**。

| 子项 | 内容 | 依赖 | 关键文件 |
|---|---|---|---|
| B1 | 无粘通量 scheme：LF → Roe / HLLC / AUSM 系列（device），LF 保留作 oracle | A | `UNsInvFlux`、`NsInvFlux`、`HipFluxBackend` flux kernel |
| B2 | limiter：Barth / Venkat（cell 邻域极值），解除 `ilim=0` | B1 | `UNsLimiter`/`NsLimiter`、recon kernel |
| B3 | 层流粘性：`UNsVisFlux`（面粘性通量、T 梯度、壁面热流） | B1/B2 | `UNsVisFlux`、gradient（多方程）、BC 热流 |
| B4 | 湍流：SA（1 方程）优先，SST 可选；`turb` solver 与主 solver 的 task 交互 | B3 | `codes/turb`、`codes/uturb`、`system/turb` 字段 |
| B5 | 隐式：LU-SGS 设备化（下/上扫、`dq` 交换、与 halo 的依赖），再评估双时间 | A4 | `NsLusgs`、`ULhs`、`LhsTaskReg` |
| B6 | 多重网格：restrict/prolong 的设备化，或显式声明"只加速 finest grid，粗网格留 host" | A4 | `codes/multigrid` |

准入门（每个子项都必须先过）：
1. 单点数值契约（scheme/limiter 的 CPU vs device 单元测试，≤1e-14）；
2. 3 步多区 trace（1e-10）；
3. 该能力对应的物理门禁（正性、守恒、边界语义）；
4. 200 步跨后端 + 50 步稳定性。

退出门：`m6wingroe_sa`（Roe+LU-SGS+SA）在设备上跑通并通过同一组门禁。
这个算例在 todo 里记为 **B-exit**，是 Phase B 的退出条件，做完 B1–B6 不等于 B 结束。
之后才允许把 "HIP 支持的算例范围" 写进正式文档。

### Phase C — 真实算例与 IO（把加速变成可交付能力）

| 子项 | 内容 | 说明 |
|---|---|---|
| C1 | IO/restart：设备输出或"整场一次性 D2H + 异步写盘" | 现在每输出步 `NsRefreshHostStateMirror` 整场回拉，长跑/多输出会吃掉收益 |
| C2 | 气动力/Cp/热流：设备侧面量归约，或维持边界 cell 级小回拉 | 现在每步一次边界 cell 镜像，可接受但需量化 |
| C3 | 残差/收敛判据设备归约 | 已在设备侧算 metrics，需扩到多 zone/多 solver 聚合 |
| C4 | 内存预算与大网格 | 记录每 100 万 cell 的 device 常驻字节数，给出单卡可解规模 |
| C5 | 定常/非定常的真实工况回归 | 与 `test/suites` 对齐，形成"设备算例集" |

退出门：≥1M cell 的单卡算例端到端可跑（含 restart 等价、气动力对比、可视化），
并有对外可引用的 paired 性能数字。

### Phase D — 多卡与多节点（结论的前提是 fixture）

| 子项 | 内容 |
|---|---|
| D1 | **fixture**：拿到 ≥1M cell 且带 `.skl` 的合法分区网格，或修 `GenerateLink`（`codes/geometry/src/FaceTopo.cpp` 的 `SetNewFace2Node`）建表 |
| D2 | 设备 halo：GPU-aware MPI（1D 已验证 `1.39/1.31/1.19/1.13x`）在 3D 的复现；DeviceDirect 在 packed layout 上重测 |
| D3 | 映射与负载：1 rank = 1 zone 的约束是否长期成立；需要 zone 切分/多 zone per rank 时的 state 与 halo 语义 |
| D4 | 评测口径：固定 per-rank cell 数、同一 repeats、1 DCU 槽 = 8 核；给出 8:1、32:1 的合法性说明 |

退出门：至少 4 卡、per-rank 恒定 cell 数的端到端数字，且与同节点 CPU 同核数配对。

### Phase E — 架构通用化（与 A/B 并行，按痛点触发）

| 子项 | 内容 |
|---|---|
| E1 | 热路径去 ambient globals：把 zone/grid/solver 身份显式传进 stage（先 device 路径，再 host） |
| E2 | capability 表驱动：`solver × scheme × limiter × viscous × time-integral × grid-level × halo` 一张表 + 单测，替换散落守卫 |
| E3 | `EulerDomain` 泛化到 NS：3/5 方程已支持，补源项/粘性/能量扩展与 view 元数据 |
| E4 | backend-kind 分派：清掉 `IsAccelerator()` 当 HIP 用的地方，为 CUDA/Kokkos 留真实入口 |
| E5 | 清理 toy/scalar 路径与失效文档，保持 CPU oracle 单一 |

退出门：新增一种 scheme 或一种物理项，只需"新 kernel + 表项 + 单测"，不需要改
`TimeIntegral`、`UNsInvFlux`、`InterfaceTask` 的骨架。

### Phase F — 产品化与发布（每个阶段末循环一次）

| 子项 | 内容 |
|---|---|
| F1 | 3D 正式性能报告（Markdown + HTML 同 stem）：paired、同节点、同 repeats、含环境口径 |
| F2 | CI 矩阵：CPU normal/strict、DCU 单卡、DCU 多卡、nightly 长跑；PR 分支自跑门禁 |
| F3 | 文档：`doc/plans/oneflow-development-todo.md` 摘要与 §2 现状、`ci/kunshan/README.md`、Sphinx 教程 |
| F4 | 交付：GUI/示例算例、发行说明、对外能力矩阵（哪些 solver/scheme 在设备上已验证） |

## 2.1 阶段与 todo 条目 id 的对应

todo §4 直接使用本节编号，不再有第二套字母表。历史编号的对应只在需要读旧记录时使用：

| 旧编号（历史文档） | 现在 |
|---|---|
| `F7.1` | A1.1（梯度镜像只在最后 stage） |
| `F7.2` | A1.2（设备残差基底归零） |
| `F7.3` | A2（每步一次界面交换） |
| `F7.4` | A3（stage kernel 融合） |
| `F7.5` | A4（HIP `Advance`） |
| `F8.1`–`F8.6` | B1–B6 |
| `F9.1` | B-exit（旗舰算例） |
| `F9.2`–`F9.3` | C1–C5 |
| `G4`–`G6` | D1–D3 |
| `H2`–`H6` | E2、E1、E3、E4、E5 |
| `G7`、`I1`–`I4` | F1–F4 |

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
| ambient globals 阻碍驻留执行 | Phase A4 无法干净落地 | E1 提前做 device 路径的显式上下文；host 侧渐进 |
| fixture 缺（≥1M 合法分区） | Phase D 无法收口 | 现在就去谈 owner；同时评估修 `GenerateLink` 的代价 |
| 隐式/multigrid 复杂度高 | 工期失控 | 先把显式驻留做到 14 ms；隐式与 multigrid 各自独立立项、独立门禁 |
| capability 组合爆炸 | 维护成本 | E2 表驱动 + 单测；每个组合要么有门禁要么默认拒绝 |
| 性能口径漂移（跨会话跨节点比较） | 结论反复 | 只发 paired；报告里写清节点/会话限制 |
| 上游同步与 PR 规模 | 评审不可行 | 每轮三端同步 master、按"backend/stage/CI"小切片提 PR（#227 的教训） |

## 5. 与既有文档的关系

| 文档 | 已有结论 | 本文新增 |
|---|---|---|
| `doc/plans/oneflow-development-todo.md` | 前置层 6/6；阶段父项 2/28（A1、F1），Roe 为 Partial | 阶段取舍以本文为准，勾选与证据以 todo 为准，编号与本文一致 |
| `doc/architecture/oneflow-1d-euler-dcu-architecture-audit-20260904.md` | backend-neutral 接入方案、state key、调用边界 | 把该方案升格为长期阶段主轴（Phase A4） |
| `doc/plans/oneflow-euler-optimization-plan.md` | "stateful NoTrace 才是性能路径"（trace 2x vs 驻留 100x+） | 在 3D 复现同一结论：执行模型优先于 kernel |
| `doc/runbooks/oneflow-dcu-porting-playbook-20260903.md` | 移植方法论（契约→或修→门禁→benchmark） | 作为每个 Phase 子项的验收流程模板 |
| `doc/handoff/oneflow-current-handoff.md` | 环境、固定路径、启动顺序、坑（living，每次会话就地更新） | 其 §4 任务指针与本文 Phase A 对应 |

## 6. 需要拍板的开放问题

1. **目标规模**：这条后端的"完成"定义是 1M cell 单卡、还是 100M cell 多卡？
   （决定 Phase C/D 的投入比例）
2. **旗舰算例优先还是性能优先**：2026-11-01 定为先把 Phase A 做到 ≤14 ms/step，
   再扩 Phase B。B1/B2 可独立评估契约，但不挤占当前 H2 → A4 主线。
3. **隐式与多重网格是否在设备上做**：若只加速显式 finest grid，可省掉 Phase B5/B6
   的大半工作，但要接受"定常收敛用 LUSGS 就没加速"的现实。
4. **zone 映射**：长期是否维持 1 rank = 1 zone？若需要 zone 切分，Phase D 的工作量
   翻倍（state、halo、load balance 全变）。
5. **CUDA/Kokkos 是否真要做**：当前不做。F7 继续走 HIP 生产路径，但不要把 HIP
   类型写进 `EulerDomain` seam。真的要做 CUDA 时再补 H5，不回头阻塞 F7。

## 7. 下一步顺序（具体 checkbox 与验收见 todo §5）

1. **H1 已收口**：集群归账、正确性门禁与 CI 测试链接修复的构建验证已过；详见 todo §5.2。
   MG 控制域与任务规划/派发仅约 0.042 ms/step，不支持靠重写外循环获得数毫秒收益。
2. **H2 多卡 VMFault**：最小复现并修复 `GradientCopyBoundaryKernel` 相关故障；
   修复与复验完成前暂停多卡结论。D1 ≥1M 合法 fixture 仍是独立的证据缺口。
3. **A4 循环下沉**：先量清宏步尾状态读取/清零与 RHS 队列等待，再选择减少同步/传输的
   最小切片，保持每步物理状态检查与 halo lag；通过 trace、跨后端、CPU normal/strict、
   稳定性及 paired；之后再评估 A2′/A3。

（F1 3D 正式报告已于 2026-11-02 发布，不再列为待办。）
