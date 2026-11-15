# OneFLOW 开发待办与衔接（living document）

> 最后核对：2026-11-14（扩展性/网格能力轮）；工作分支 dev @ `727c1ed75`。
> 数值门禁基准仍为 `6f7a49e3d`；`d41650729` 修 H2 根因，其后 5 个提交为接口/诊断改动。

## 交接摘要与整体路线图

> **本文只放三样东西**：当前状态（§1、§3）、Phase A–F 台账与证据（§4）、下一轮任务（§5）。
> 逐日轮次记录不进本文（见文末「归档说明」）。
> 长期方案的取舍见 [`oneflow-development-plan.md`](oneflow-development-plan.md)；
> 环境、固定路径与坑见 [`../handoff/oneflow-current-handoff.md`](../handoff/oneflow-current-handoff.md)；
> 分支与门禁规则见仓库 `AGENTS.md`；集群流程见 [`../../ci/kunshan/README.md`](../../ci/kunshan/README.md)
> 与 [`../runbooks/oneflow-kunshan-workspace-runbook.md`](../runbooks/oneflow-kunshan-workspace-runbook.md)。
>
> 编号约定：一律用 Phase A–F 的 id（`A1`、`B3`、`D1`…），不要再引入第二套字母表。

### 0. 一句话现状

**根本目标：保证正确性前提下做加速。** 所有任务按“明确正确性契约 → 通过相应门禁
→ 同口径端到端性能验证”推进；勾选完成与性能收益结论须同时满足该项适用的验收条件。
完整原则见长期计划 §0，数值口径变更按 `AGENTS.md` 审查。

受限无粘 3D 显式 RK 路径的 **q 已设备权威，A1 固定开销优化已收口**；
完整 HIP `Advance` 尚未实现。**H1、H2 均已完成**。Roe 只有单点/单 zone 证据，
Barth 未过数值门禁。H2 的根因是延迟状态标志的惰性分配（空指针进 `atomicExch`），与卡数无关。

**扩展性已定性（2026-11-14）**：由每 rank 的 cell 数决定，不由卡数决定。求解器核心可扩（1→32 卡 2.5–6.0×），
交换不可扩。**注意：早期把「加卡后墙钟最低点落在 2 卡」写成「净最优在 2 卡」是错误的**——
每卡吞吐随卡数单调下降（22.9 → 17.1 → 11.0 kcell/ms/卡，1/2/4 卡）说明加卡效率在退化，
2 卡只是退化曲线上的拐点，不是推荐配置；在更多卡上落后是缺陷，不是规律。
**下一主线是网格能力**：仓库无加密/导出能力，M6 弱扩展与 D1 都卡在这里，需自写原生六面体加密（§5.5）。

### 0.1 完成与待办总览（checkbox）

`[x]` 表示已有记录支持该项的限定范围，不代表当前 HEAD 已重跑；
`[ ]` 包含 Planned / Partial / Blocked / Deferred，具体证据与验收见 §4、§5。
H1/H2 是本轮任务，D1a 是子项，均不额外计入 28 个阶段父项。

- [x] 前置层：1D 端口、加速器基础层、状态契约、CPU 垂直切片、受限 3D HIP 路径（§4.7）。
- [x] A1：设备权威 q、按需镜像、跳过无读者的 host 账本、延迟同步等固定开销优化。
- [x] Roe 已实现并完成单面契约 3/3、单 zone 200 步验证（B1.0–B1.2；B1 父项未完成）。
- [x] Barth 单点契约与逐点差异定位（B2 部分完成；未通过整链数值门禁）。
- [x] D1a：结构化等块分区，9/9 单测及固定 rank 配对收益记录。
- [x] F1：3D 正式性能报告已发布（2026-11-02，Markdown + HTML）。
- [x] **H1（已完成）**：五对归账差额 0.34%–2.14%；下一刀见 §5.2。
- [x] **H2（已完成，范围见 §5.2）**：根因是延迟状态标志的惰性分配，与卡数无关；
      多卡数值路径已复验（HIP vs CPU 残差逐位一致）。**完整门禁尚未跑**，扩展性数字仍待重测。
- [x] **扩展性定性（2026-11-14）**：每 rank cell 数决定扩展性；核心可扩、交换不可扩；
      加卡后每卡吞吐单调下降（拐点约在 2 卡，但那是**待修的退化**，不是推荐配置）；弱扩展迭代层持平（§5.5）
- [x] **网格加密能力（2026-11-14 完成）**：`ci/kunshan/refine_ofl.py` 六面体 2×2×2 加密；
      M6 → 2,359,296 cell，体积守恒 5e-15，单 zone XB 3.07e-17，64-zone XB 全绿（§5.6）
- [ ] **A4**：HIP `Advance` 接管 RK 循环；通过数值门禁并测量端到端收益。
- [ ] A2′ / A3（降权）；原 A2 每步一次交换方案暂缓，不直接删除 post 轮。
- [ ] B1：补 Roe 多区 trace 与对应 CPU 五算例证据；HLLC/AUSM 后续单列能力。
- [ ] B2：修复 Barth `1.113e-06` 偏差并过门禁；Venkat 尚未实现。
- [ ] B3–B6 / B-exit：粘性、SA、LU-SGS、多重网格及旗舰算例设备验证。
- [ ] C1–C5：IO/restart、气动力、归约、内存预算与真实工况的完整设备交付。
- [ ] D1 / D1b / D2–D4：≥1M 合法 fixture、控制 halo 的弱扩展、设备 halo 与映射验证。
- [ ] E1–E5：显式上下文、能力表、NS domain、后端入口及旧路径清理。
- [ ] F2–F4：CI 矩阵、完整能力文档、示例与发行交付。

### 1. 状态快照

| 项 | 值 |
|---|---|
| refs | 工作分支 dev @ `727c1ed75`。数值门禁基准仍为 `6f7a49e3d`。`d41650729`（H2 根因）、`d44e0253b`、`970508e7c`、`2006b0a84`、`d9e6a14db`、`727c1ed75` 为本轮提交，均已在集群构建；XB 门禁逐位不变。上游基线 `f36f0c6b3` 已合入；接手仍须核对实际 HEAD |
| upstream PR | #159 / #160 / #203 / #227 / #238 / #253 全部 MERGED（#227 = 3D HIP 功能线；#238 = 昆山工作区 runbook；#253 = upstream 同步与收口） |
| 集群口径 | Kunshan DTK 26.04 / `gfx906` / 1 DCU；CPU 基线在 `kshcnormal`；**1 DCU 槽 = 8 核** |
| 性能基线 | 历史 8 zone / 8 rank / 1 DCU、1000 步 paired：**14.62 ms/step**（2026-11-07 同节点交错 5 对中位）。CPU 288.20 来自 2026-11-02 的 F1 报告。F1 当次配对 CPU 288.20 / HIP 15.59 = 18.5×；2026-11-14 复测同作业配对 CPU 282.5 / HIP 15.84 = **17.8×** |
| 扩展性（2026-11-14 实测） | 每 rank cell 数决定扩展性，不由卡数：M6 强扩展 1/2/4/8/16/32 卡 = 11.70/7.60/9.87/18.73/24.08/30.62 ms/步（墙钟最低点在 2 卡，但每卡吞吐随卡数单调下降，属待修退化，**不是**推荐配置）；立方体 884k 强扩展同形状；立方体弱扩展（110,592 cell/rank）2→32 卡迭代 13.78→13.68。求解器核心 1→32 卡 M6 2.5×、立方体 6.0×，差额全在交换 |
| 门禁基线 | `6f7a49e3d`：HIP/计时契约 22/22；3 步多区 trace PASS、resident `24/0/0`；200 步 XB residual `2.331035e-18` / aero `0.0`；CPU normal/strict 各 5/5（判据 1e-8 / 1e-15 未变）；50 步 HIP 稳定性及五对 1000 步归账 PASS。全量 GoogleTest 未重跑 |
| 加速器构建 | 可选、后端无关：`ONEFLOW_ACCEL_SDK_ROOT`（+ 可选 `_ARCH`/`_FLAGS`），种类由 `ONEFLOW_ENABLE_HIP/CUDA` 或前缀内容推断；一键 `scripts/build-accel.sh`。不请求则完全不生效（CPU/CI 不受影响） |
| 已知红灯 | **H2 已解（`d41650729`）**：根因不是多卡，而是 `ONEFLOW_HIP_DEFER_STATUS=1` 下 `stateUpdateStatus` 惰性分配被 defer 分支完全跳过，`Data()` 返回空指针进 `atomicExch`；表现为 `libgalaxyhip`/`hsa-runtime64` 里 signal 6，回溯指向 `GradientCopyBoundaryKernel`（该内核只是恰好在队列里）。**单卡也会触发**。B2 的 `1.113e-06` 未过门禁。`ActionMapReadFileTest` 曾在全量并发下偶发 |
| 能力边界 | `NS_SOLVER` / 5 方程 / **Lax–Friedrichs 或 Roe** / `ilim=0`（Barth 已实现但数值未过门禁，默认不放行）/ `vismodel=0` / 单网格层 / 稳态显式多阶段 RK / `nGhostCells==nBFaces`；其它配置自动回退 host，显式请求不支持项 fail-fast。**多 zone/rank 可通过配置检查，多卡数值路径已复验**；多卡扩展性要求每 rank 固定 cell 数（§4.8b），1 rank 多 zone 更慢（§4.1 A1.8） |
| 尚未覆盖 | Venkat、粘性/RANS、LUSGS、multigrid、双时间；IO/restart/气动力/可视化仍在 host；多卡数值路径已通。**缺 ≥1M 合法分区 fixture（D1）的唯一原因是该仓库无网格加密/导出能力**（`IsSupportedConversion` 中 `OneFLOW` 只作为 target），需自写原生六面体加密，见 §5.5 |

### 1.1 工作环境与固定路径（不要另建目录）

这一节是"agent 不要到处乱放"的唯一出处；路径一律写成可解析的形式，具体绝对路径
不写进仓库（发布边界）。接手时按 [`doc/handoff/oneflow-current-handoff.md`](../handoff/oneflow-current-handoff.md)
执行。

| 位置 | 路径 | 用途 | 规则 |
|---|---|---|---|
| 本地仓库 | 本 checkout 根（`git rev-parse --show-toplevel`） | 源码、文档、提交 | 只在 `dev` 提交；`master` 只读，与 upstream 同步；不建第二份 clone |
| 本地临时 | `/tmp` | 一次性脚本、打包、对比 | 不入库；结论写回文档 |
| 集群接入 | scnet-hpc 技能（按当前安装位置解析入口） | 查询、提交、上传、探测 | 按交接现状只用 OpenAPI；不手搓 ssh/scp 循环，不静默切换后端 |
| 昆山工作区 | **`$HOME/oneflow/`（唯一）** | 见下方固定布局 | 不新建第二个工作区 |
| 源码 | `$HOME/oneflow/src/OneFLOW` + `src/REVISION` | `git archive` → upload → 作业内覆盖解包 | 不建 `OneFLOW-<hash>` 目录；不保留旧源码残留 |
| 构建 | `$HOME/oneflow/build/app`（`build/port`、`build/platform` 预留） | 固定目标、增量编译 | 不建 `build/<jobid>`；日常不重新 configure，必要时按 runbook 的 SDK 接口处理 |
| fixture | `$HOME/oneflow/fixtures/m6-8zone-cpu8/` | 唯一登记的固定分区算例 | 只读；其它 zone 数按 runbook §1b 生成并校验 |
| 证据 | `$HOME/oneflow/evidence/<category>/<jobid>/` | `gate` / `paired-ab` / `cpu-regression` / `probe` / `fixture` / `build` | 只增不改；大二进制提取指标后可删 |
| 快照 | `$HOME/oneflow/snapshots/*.tgz` | 保留最近 5 个 | `cleanup.sh` 轮转 |
| 日志 | `$HOME/oneflow/logs/` | 作业 stdout/stderr | 保留 14 天，不写 `$HOME` |
| 弃用 | `$HOME/oneflow/_attic/<date>/` | 待删物 | 观察一轮后清（`cleanup.sh`，7 天） |

**禁止清单**：在 `$HOME` 根或 `softwares/projects/` 下新建 OneFLOW 目录/日志/一次性脚本；
按 job 号或轮次新建 build/run 目录；把绝对私有路径、主机名、账号、作业号写进提交文档；
每次会话新建 dated handoff；绕开技能直接改远端状态。

作业与环境常量（`mpirun` 启动、`scripts/env.sh`、HCA、`OMPI_MCA_coll`、cmake/DTK 路径、
1 DCU 槽 = 8 核、基准口径）见 runbook 与 living handoff，本节不重复。

### 2. 回顾（2026-08 → 2026-10-25）

| 时间 | 阶段 | 交付 | 关键证据 |
|---|---|---|---|
| 2026-08 | accelerator 框架与 CPU 基线 | `codes/accel` 抽象、CPU backend、五算例回归、昆山 CI | 框架交付报告；当时明确"只有框架，无数值 kernel" |
| 2026-09-13 | 两条线融合（A–E 层） | upstream 基线 + Phase 1–3 + EulerDomain/StateRegistry + CPU vertical slice | bridge 2.22e-16；根 CTest 163/163；五算例 5/5 |
| 2026-09-20~24 | 3D DCU vertical slice（F1–F3） | root HIP 编译/smoke/contract、3D m6 1-step/3-step trace、fixed-CFL 50 步、梯度切片 | trace 3e-15 级；50 步三路 PASS |
| 2026-10-09~19 | 多 zone 与 halo | METIS `Partition`、多 zone device state/boundary、GPU-aware MPI probe、policy 收口 | 1D 双节点 4 rank/4 DCU hash 一致；`1.39/1.31/1.19/1.13x` |
| 2026-10-20~23 | host 链优化 | cell-CSR gradient、packed gradient/q 交换、接口字段裁剪、壁面状态裁剪、ghost 打包 | 59.4 → 32.8 ms/step（跨会话） |
| 2026-10-24~25 | **device 权威 q（本轮）** | 契约增加 host-mirror currency、状态更新不 D2H、按需回拉、host 派生场按需重建 | 同机 paired 38.26 → 31.64 ms/step；PR #227 合并 |

性能演进（8 zone / 8 rank / 1 DCU，1000 步）：`59.4 → 38.3（同机旧行为）→ 31.6 ms/step`。
注意 2026-10-24 记录的 `32.8` 与 2026-10-25 的 `38.3` 不是同一节点，跨会话只比 paired 差。

### 3. 全链路审计（当前）

| 链路段 | 状态 | 备注 |
|---|---|---|
| 网格读入 / 分区 / wall distance | host | 只摊在启动阶段 |
| solver/state 创建、init q | host + 一次性 H2D | 8 zone 1000 步时占 3% 左右 |
| **时间推进调度（RK stage 循环）** | **host task 序列** | 每 stage 派发 + fence，是剩余固定开销的主因（§4 F7） |
| gradient / reconstruction / flux+residual / state update / 本地边界 / 本地时间步 | device | 每 stage 一批 kernel，q 不回写 |
| q 与 gradient 的 interface 交换 | device pack + host MPI | 每步 4 轮（3 stage + 1 post），§4 A2 |
| limiter / LUSGS / multigrid / 双时间 | limiter = device（Barth）/ 其余 host | Barth 的 φ 由设备 `ReconstructionLimiterMode::DeviceCell` 计算（主机侧重复的 `CalcLimiter()` 已删）；数值仍差 1.1e-6 未过门禁，默认 `ilim=0`；Venkat 未实现，LUSGS/multigrid/双时间回退 host |
| 粘性 / 湍流 / SA | host | 能力外，回退；旗舰算例正在此列 |
| 残差诊断 / 正性检查 | device 归约标量 | 单 zone 已闭环，多 solver 聚合待做 |
| restart / 可视化 / 气动力 / Cp | host | 按需 D2H（气动力只回拉边界相邻 cell） |
| 分区与 halo 拓扑 | host + 设备端打包 | DeviceDirect 在 packed layout 上未重测 |

### 4. 路线图与执行进度（Phase A–F）

编号与 [`oneflow-development-plan.md`](oneflow-development-plan.md) §2 完全一致：**A1–A4 / B1–B6 /
C1–C5 / D1–D4 / E1–E5 / F1–F4**。本文只记勾选与证据，阶段的取舍、准入门与退出门见计划。
不再另起一套字母表；早期文档里的 `F7.x`/`F8.x`/`F9.x`/`G4`/`H3` 等编号已在各项后就地标注对应关系。

规则：**只勾有证据的项**；父项只有在子项全勾时才勾；`[ ]` 表示未开始或未通过门禁。

#### 4.1 Phase A — 3D 设备驻留执行（1/4）

- [x] **A1 固定开销**（原 `F7.1`/`F7.2` 另加三项死拷贝；同节点 29.9 → 26.7 → 21.0 → 18.5 → 17.3 → 15.6 → **14.62** ms/step（A1.7c 收口）；F1 报告已发布）
  - [x] A1.1 梯度主机镜像只在最后一个 RK 阶段做。气动力与 Cp 在步末才读，中间 stage 的镜像没人看。同节点 paired 每 stage 29.703/30.017 → 最后 stage 26.660/26.779 ms/step（−10.5%）；单 rank `gradient_mirror` 3.20 → 0.93。`ONEFLOW_HIP_GRADIENT_MIRROR_ON_DEMAND=0` 恢复每 stage
  - [x] A1.2 RK 基线状态与残差基底不再从主机搬。第一次只删 `LOAD_Q`/`LOAD_RESIDUALS` 的 host 写入、没补设备侧基线，残差在上一 stage 上继续叠加，200 步跨后端拓到 `5.52e-4`；把“不填 host 槽位”和“设备残差归零”绑成同一个谓词（`DeviceOwnsExplicitRkBookkeeping()`）后对齐。paired 26.91 → 20.97 ms/step（−22.1%）。多重网格、粗网格、双时间、粘性仍填槽。`ONEFLOW_HIP_SKIP_HOST_RK_SLOT_COPY=0` 恢复
  - [x] A1.3 设备权威路径不再搬 host 历史层。`UPDATE_UNSTEADY_FLOW` 每步把 `q2←q1←q` 各做一次全场拷贝；双时间被 capability 排除，restart 写盘时又从设备场重建，这些层在推进中无人读。跳过条件是 `NsDeviceAuthoritativeFlowEnabled()` + 显式定常 NS，`ONEFLOW_HIP_SKIP_HOST_HISTORY_SHIFT=0` 恢复。同节点 paired 21.50 → 18.50 ms/step（−14%）；restart auth-on `5.03e-15`、auth-off `4.84e-15` 均 PASS，trace `24/0/0`，200 步跨后端 `2.331035e-18`/`0.0`，CPU normal/strict 5/5
  - [x] A1.5 设备单网格路径不再写 host 残差账本。`ZERO_RESIDUALS`（命名残差清零）与 `STORE_RHS`（拷进 `FIELD_RHS`）只服务多重网格缺陷与 host 显式更新；设备单网格路径两者都不读，却是每步两次全场 host 写。判据统一为 `HostResidualBookkeepingNeeded()`（多重网格配置或设备不拥有账本时保留）。8 rank 门禁 3 发 `17.27 / 17.25 / 18.14 ms/step`，trace `24/0/0`、200 步 XB `2.331035e-18`/`0.0`、CPU normal/strict 5/5
  - [x] A1.6 分项归因修正。任务级 profile 里 `POST_PROCESS` 的 `mean_ms` 曾显示每步 7.7 ms，实为**周期性可视化输出**摊到每一步的均值：总 2193 ms / 300 步，把 `nVisualSave` 设成不输出后该任务降到 0.055 ms/step。因此不能把它当作每步固定开销，也不能据此设计优化。教训：`mean_ms` 必须和 `count × 发生结构` 一起看。
  - [x] A1.7 每步固定基数的边界测量（2026-11-04）。**关通信**、`rk_stage=1`、8 rank/1 DCU 下仍有 `8.93 ms/step`；该基数与 cell 数、zone 数、rank 数都基本无关（1/2/4 zone per rank 的 1-stage 分别为 8.93 / 7.61 / 8.65）。同时测得：关闭通信只省 **2.2 ms**（23.97→14.28 首发、16.51→14.34 热态），所以“通信主导”不成立，且 A2 的收益上限就是这 2.2 ms。这条把 A2/A3/A4 的目标从“通信”改成“每步基数”。**注意：本条的 8.93 由 A1.7b 修正过**——100 步/300 步会把 1.4–1.6 s 启动摊进 ms/step；按 1000 步斜率基数 ≈5.5 ms，归因见 A1.7b/A1.7c
  - [x] A1.7b 固定基数的第一次实测（2026-11-07，昆山 1 DCU / 8 rank / 8 zone，m6wing，1000 步）：
      `rk_stage=3` 15.58 ms/step，与 F1 报告的 15.59 同一口径复现（合入 upstream `f2fc6ea2e` 后）；
      `rk_stage=1`（关界面通信）6.88 ms/step，所以每多一个 stage 的边际成本约 **4.35 ms**。
      **短基准会把它算错**：同一配置每次启动约 **1.4–1.6 s**（网格/分区/首次分配），
      100 步口径会把启动摊成 +14 ms/step、200 步约 +7 ms/step，A1.7 原来的 8.93 属于这一类读数，
      下次引用基数必须写步数并做启动扣除。**队列 drain 不是基数**：
      `ONEFLOW_HIP_DEFER_SYNC=1` 在 1000 步只把 15.58 降到 15.25（−2.1%），1 stage 时 6.88→6.39。
      结论：基数 ≈ **5.5 ms/step**（1 stage，已扣启动），下一步要用 `macro_*` 分段把它归到具体段
  - [x] A1.7c 第一版端到端优化（2026-11-07，1 DCU / 8 rank / 1000 步，同节点交错 3 发）：
      基线 **15.57**（15.50/15.57/15.62）→ 三个开关全开 **14.55**（14.54/14.68/14.55），
      **−6.5%**，样本无重叠。门禁全绿：trace `24/0/0` + `max_abs_diff=0`、200 步 XB PASS。
      改动：残差缩放并入状态更新入口（`ONEFLOW_HIP_FUSE_SCALE_UPDATE`）、非物理标志留在设备上
      逐步检查（`ONEFLOW_HIP_DEFER_STATUS`，每步少 4 次 4 字节阻塞 D2H）、drain 只在需要回读处
      （`ONEFLOW_HIP_DEFER_SYNC`）、inflow 值不变时不再每 stage 重传（原来它成了延迟队列的阻塞点，
      `bnd_static_upload` 曾显示 4.2 ms/step）。**三个开关已改为默认开**（`=0` 回退），
      另加 `ONEFLOW_HIP_MAX_LAUNCHES_PER_DRAIN`（默认 32）防止多 rank 下队列溢出 abort。
      收口复跑（2026-11-08）：gate 三发 14.75 / 15.08 / 15.78（该节点，跨度 1.0），
      trace `24/0/0`、XB PASS、CPU 五算例 5/5 + 5/5；性能口径仍以同节点 paired 的 14.62 为准
  - [x] A1.8 多 zone/rank 的可行性（2026-11-04）。`Zone::CheckZoneRankCount` 只查 `nZones % nProc`，所以 16 zone / 8 rank 能通过配置检查；但运行期死在 `HIP gradient face pack has no resident gradient`——该 rank 还没有推进过那个 zone。修法：设备侧梯度打包在缺 resident gradient 时**回退主机镜像**（此时主机场是权威初值），不再抛异常。修后 `zpr=2/4` 能跑完，但**更慢**：16.45（1 zone/rank）→ 17.62（2）→ 25.24（4）。结论：**“一 rank 多 zone”不是解法**，它不减少 halo，反而按 zone 重复每步工作
  - [x] A1.4 基准 warmup 规则。首个 `hipMalloc` 与冷 case 目录会污染短基准：同一配置冷跑 `27.0`、热跑 `18.5–21.1 ms/step`。规则落在交接 §5，不是代码改动
- [ ] **A2 每步一次界面交换**（原 `F7.3`，Deferred；候选方案见 §5.2 A2′）：跳过 post 轮使 200 步跨后端从 `2.33e-18` 偏到 `5.54e-7`（aero `4.0e-5`），不能直接删轮次。单卡关掉全部通信的历史墙钟差值约 **2.2 ms/step**；旧 1.1 ms 收益上限缺乏完整归账，不再作为验收预期。多卡数据须等 H2 修复后复验
- [ ] **A3 stage kernel 融合**（原 `F7.4`，Deferred）：gradient → reconstruction → flux → lhs → update → boundary 在设备内串联。旧的 −4/−6 ms 收益估算已失效；host stage 计时不是 GPU 执行时间。待 H1/H2/A4 后再按配对实验评估
- [ ] **A4 HIP `EulerDomainBackend::Advance`**（原 `F7.5`）：RK 三阶段循环进 backend，只在宏步边界回主机（输出、门控）。设备基态留在设备上，不要把 A1.2 跳过的 host 槽位拷回来。**收益待验证**：H1 已测得 MG/任务规划与派发仅约 0.042 ms/step；A4 的价值应来自减少同步与状态往返，旧 9–10 ms 估算不再作为预期

#### 4.2 Phase B — 数值能力（0/6，退出门见 B-exit）

- [ ] **B1** 无粘通量 scheme（Partial）：**Roe 已设备化并过单面契约**；多区证据尚缺。HLLC/AUSM 未开始，
    需要时按同一路径另开子项；LF 保留作 oracle
  - [x] **B1.0** scheme 透传：`nscom.ischeme` 经 `DeviceFluxScheme()` 显式映射进内核，
    替掉原先硬编码的 `1`。**命名陷阱保留作记录**：内核私有编号 `1`=Rusanov、`2`=Roe，
    与求解器的 `ISCHEME_ROE = 1`、`ISCHEME_LAX_FRIEDRICHSS = 5` 是两套编号，映射函数
    就是为了不让两边按名字互相推断；未实现的 scheme 抛 `invalid_argument`，
    capability 增加 `roeScheme` 位后接受 Roe
  - [x] **B1.1** 设备侧 Roe：`HipFluxBackend.hip` 的 `RoeFaceFlux` 按 `codes/ns` 的
    `NsInvFluxPolicy` 逐项对齐（Roe 平均态、Harten 熵修正、绝对特征值与耗散波系数、
    法向跳跃与耦合项），由 `scheme == 2` 分支调用；LF 分支保持逐位不变
  - [x] **B1.2** 单点数值契约：`oneflow_device_roe_contract_test` 3/3；单 zone 200 步
    `res_max_rel = 0.0`，验证范围限于单点和单 zone
  - [ ] **B1.3** 补齐 3 步多区 trace 与 CPU 五算例的 Roe 口径；完成前不勾 B1 父项，
    不将单 zone 结果写成生产验证结论
- [ ] **B2** limiter：Barth / Venkat（cell 邻域极值），解除 `ilim=0`
  - 已落地：设备侧 Barth 因子的**契约测试** `oneflow_device_limiter_contract_test`
    （5/5，覆盖范围内、超上界、低于下界、反向法向、斜面），参考实现是主机
    `BarthFunction` + `ReconstructFaceIncrement` + `FaceSideDirectionDot`
  - **关键结论（2026-11-02）**：**φ 不能在主机算**。设备权威 q 下主机梯度镜像只覆盖
    边界相邻单元、且只在有读者时才更新；用它算 φ 得到陈旧因子。实测 `ilim=1`
    200 步残差偏离 `8.53e-05`，逐 step 累积
  - **已实现（2026-11-03/04）**：设备侧 Barth φ —— 每 cell 一个线程、5 个方程在寄存器、
    邻接表走两遍（邻域极值 + 逐面取 min）；`ReconstructionLimiterMode::DeviceCell`
    显式声明"因子由后端算"；重构按主机语义 `min(φ_l, φ_r)` 施加。同时**删掉主机侧
    重复的 `CalcLimiter()`**（每 stage 270 ms）改为只填中性值
  - **当前数值状态**：`ilim=0` CPU vs HIP `0.0`（逐位一致）；`ilim=1` 200 步差
    **`1.113e-06`**，未过 1e-12 门禁。已排除：面上混合规则、ghost φ、主干预。
    **下一个怀疑点是主机 `LimField::CalcFaceValue` 的正性回退**：主机按面用
    `ckfun`（`NsCheckFunction`）整组丢弃受限值，设备是逐方程判定候选后回退单元值，
    两者不等价，量级与实测偏差相符
  - **端到端（单 zone，1 核 vs 1 DCU，200 步，关输出）**：`ilim=0` 1621.9 vs 24.25
    ms/step（66.9×）；`ilim=1` 2617.3 vs 63.67（41.1×）。设备限制器仍带来约
    39 ms/step 代价，其中 `recon_limiter_kernel` 仅 0.03 ms/调用，其余未定位
    （怀疑是每 stage 的 `FillNeutralLimiterField` 全场填充）。注意这是 **1 核 vs 1 卡**
    口径，不能与 8:1 槽的 15.6 ms/step 直接比
  - **B2.0 逐点比对已完成（2026-11-08，单 zone / 1 步 / `ilim=1`）**：新增
    `ONEFLOW_UNS_LIMITER_DUMP=<dir>`，把主机与设备的 φ 和面值各写一份二进制后再比：
    - `phi`：1,474,560 个值中 **6,665 个不一致（0.45%）**，形态是主机 0 / 设备 1（幅度 1.0），
      即"该单元是否取极值"判定不同；
    - `qf1/qf2`：各 4,481,280 个值中 **40,562 / 38,415 个不一致（≈0.9%）**，
      最大差 7.49e-01，集中在 eq=1（u 分量），且差异面与 φ 差异单元对应。
    结论：**不是整体公式错，而是少量单元/面的判定不同**，与"整组丢弃 vs 逐方程回退"或
    "邻域 stencil 不同（边界/ghost 面是否参与 min/max）"两类假设都相容。
  - 待做（B2.1）：按差异面清单核对 Barth 邻域构造（设备 `BarthLimiterKernel` 走 cell→face CSR
    且按 `boundaryMask`/`boundaryOperation` 跳过边界；主机 `Limiter::CalcLimiterScalar` 的邻域
    需逐项对照），再决定是补整组判据还是对齐 stencil；之后过 200 步跨后端与 CPU 五算例，
    再压 `FillNeutralLimiterField`，最后才是 Venkat（`ilim=2`）
- [ ] **B3** 层流粘性：面粘性通量、温度梯度、壁面热流。先做 E1（去 globals）
- [ ] **B4** 湍流：SA（1 方程）；SST 可选项，不单列
- [ ] **B5** 隐式：LU-SGS 设备化（下/上扫、`dq` 交换、与 halo 的依赖），再评估双时间。依赖 A4
- [ ] **B6** 多重网格：restrict/prolong 设备化，或明确声明只加速 finest grid。依赖 A4
- [ ] **B-exit**（原 `F9.1`）`test/cfd/m6wingroe_sa`（Roe + LU-SGS + SA）在设备上通过同一组门禁。这是 Phase B 的退出门，做完 B1–B6 不等于 B 结束

#### 4.3 Phase C — 真实算例与 IO（0/5）

- [ ] **C1** IO / restart：设备侧输出，或整场一次性 D2H + 异步写盘。设备权威 q 下的 restart 历史层已在 2026-11-01 修好（见 §5 第 3 项），但长跑与多输出仍缺量化
- [ ] **C2** 气动力 / Cp / 热流：设备侧面量归约，或维持边界单元级小回拉并量化代价
- [ ] **C3** 残差与收敛判据的设备归约，扩到多 zone、多 solver 聚合（目前只有单 zone 证据）
- [ ] **C4** 内存预算：记录每 100 万 cell 的常驻字节数，给出单卡可解规模
- [ ] **C5** 定常/非定常真实工况回归，与 `test/suites` 对齐，形成“设备算例集”

#### 4.4 Phase D — 多卡与多节点（0/4）

- [ ] **D1** 合法 fixture（原 `G4`）：≥1M cell 且带 `.skl` 的分区网格，或修 `GenerateLink`（`FaceTopo::SetNewFace2Node`）。外部依赖，越早启动越好
- [ ] **D2** 3D 设备 halo：GPU-aware MPI 与 DeviceDirect 在 packed layout 上重测（原 `G5`；1D 已验证 `1.39/1.31/1.19/1.13x`）
- [ ] **D3** rank/zone 映射与负载（原 `G6`）：1 rank = 1 zone 是否长期成立；需要 zone 切分时 state 与 halo 语义全变
- [ ] **D4** 评测口径：固定 per-rank cell 数、同一 repeats、1 DCU 槽 = 8 核，并说明 8:1 与 32:1 的合法性

#### 4.5 Phase E — 架构通用化（0/5；E1 已完成的部分见前置层）

- [ ] **E1** 热路径去 ambient globals（原 `H3`）：把 zone/grid/solver 身份显式传进 stage，先设备路径再主机
- [ ] **E2** capability 表驱动（原 `H2`）：solver × scheme × limiter × viscous × time-integral × grid-level × halo 一张表 + 单测。已有 `EulerInvFluxCapability` / `EulerRungeKuttaCapability` 两个枚举与稳定拒绝原因作起点
- [ ] **E3** `EulerDomain` 泛化到 NS（原 `H4`）：补源项、粘性、能量扩展与 view 元数据；1D 与主 solver 共用 `Advance` 语义并入这一项
- [ ] **E4** backend-kind 真实入口（原 `H5`）：清掉把 `IsAccelerator()` 当 HIP 用的地方，为 CUDA/Kokkos 留真实入口。构建期 SDK 接口已可选（`ONEFLOW_ACCEL_SDK_ROOT` + `--arch sm_80`），运行期入口仍是保留位，无 NVIDIA 机器实测
- [ ] **E5** 清理 toy/scalar 路径与失效文档（原 `H6`），保持 CPU oracle 单一

#### 4.6 Phase F — 产品化与发布（1/4）

- [x] **F1** 3D 正式性能报告（原 `G7`/`I1`）：同节点、同 repeats、paired。
      2026-11-02 发布，Markdown 与 HTML 同 stem：
      [`…-current.md`](../reports/performance/oneflow-3d-euler-performance-current.md) /
      [`…-current.html`](../reports/performance/oneflow-3d-euler-performance-current.html)。
      8 核 CPU 288.20 vs 1 DCU 15.59 ms/step（中位数，3 样本交替），**18.5×**；
      口径为 1000 步、50 步 warmup、关闭中间输出（输出摊薄实测 0.8 ms/step）
- [ ] **F2** CI 矩阵（原 `I2`）：CPU normal/strict、DCU 单卡、DCU 多卡、nightly 长跑；PR 分支自跑门禁
- [ ] **F3** 文档与能力矩阵（原 `I3`）：哪些 solver/scheme 在设备上已验证
- [ ] **F4** GUI / 示例算例 / 发行说明（原 `I4`）

#### 4.7 前置层（已完成，不计入阶段进度）

这些是加速主线开工前的地基。早期用 A–E 字母编号，和 Phase A–E 撞名，现在只留名字。

| 前置层 | 证据 |
|---|---|
| 运行基线 | `master`/`dev` 三端同步与分支模型走通；`sync-master.sh` macOS/Linux 双兼容；可选、后端无关的加速器构建接口（`cmake/OneFLOWAccelSDK.cmake` + `scripts/build-accel.sh`，昆山全新目录 configure 28 s） |
| standalone 1D 端口 | CPU/HIP lifecycle、NoTrace/FullTrace、单卡 HIP contract、1/2/4/8/64 rank 与 1/4/8 DCU 的 hash 与分解、应用级 `3.85–9.70x`（对 32 rank CPU）、4 卡/4M `25.55x` |
| accel substrate | `AccelRuntime`/`AccelBackend`/`DeviceBuffer`、`FluxBackend` 与 backend-neutral views、backend-kind 分派 |
| domain contract | `EulerDomain` views、`StateRegistry` 生命周期、`INIT_FLOWFIELD`/restart hook、`ResidentStateContract`（owner token、代际、ghost epoch、镜像新鲜度） |
| CPU 垂直切片 | CPU adapter、registry 生产化、RK seam 与回退、50 步输出级 oracle、逐面 trace、物理与守恒门禁 |
| DCU 受限路径 | root HIP 编译/烟测/契约、3D m6 精度门禁、fixed-CFL 50 步稳定、七段设备链、设备权威 q、多 zone 与 200 步跨后端 |

#### 4.8 顺序与插空规则

1. 已定：先收口 Phase A，再扩 Phase B。**当前顺序为 H2 → A4 → A2′/A3**；原 A2 已降级、不排期。H1/H2 的验收见 §5.2，低于 14 ms 的单个样本不等于 Phase A 完成。
2. B1（Roe）与 B2（limiter）是 B-exit 的前置；B3/B4 依赖它们；B5/B6 依赖 A4。
3. E1（去 globals）应在 B3 之前，最好也在 A4 之前。CUDA/Kokkos 当前不做，但不要把 HIP 类型写进 `EulerDomain` seam。
4. D1（fixture）是外部依赖，不阻塞 A/B。
5. F1 已发布；后续性能报告只有在新配对证据到位后更新，Markdown/HTML 同步。
6. D2/D3 依赖 H2 修复与 D1；前置满足前只做只读评估，不发布多卡性能数字。

#### 4.8b 并行扩展实测（强扩展 2026-11-08；弱扩展 CPU/DCU）

**CPU 强扩展实测（2026-11-08，`a02r4n11`，m6wing 294,912 cell，1–32 核）**

| cores | ms/step | 加速比 | 效率 | 每翻倍 |
|---:|---:|---:|---:|---:|
| 1 | 1580.32 | 1.00 | 100% | — |
| 2 | 892.54 | 1.77 | 88% | 1.77× |
| 4 | 474.06 | 3.33 | 83% | 1.88× |
| 8 | 279.32 | 5.66 | 71% | 1.70× |
| 16 | 142.46 | 11.09 | 69% | 1.96× |
| 32 | 73.70 | 21.44 | 67% | 1.93× |

拟合 `T ∝ n^-0.85`（六点误差 <5%），外推 512 核约 38% 效率。关通信对比
（`a02r4n14`，51 步两点法）：halo 占 8 核 30%、16 核 16%、32 核 **13%**。
32 核分项（`cpu32`）：`rhs_inviscid` 55.64（其中 gradient/reconstruction 各 20.98）、
界面约 11.5、其余约 6.8 → 每 cell 每 stage 约 0.76 µs（合理值 20–50 ns）。
**弱扩展阻塞已解除（2026-11-08）**：`codes/special/src/Transfinite.cpp` 的三维
transfinite 插值有两处复制粘贴错误——k 方向混合用了 j 向的边值（`dk*v0` 应为 `dk*w0`），
`uv` 角点项最后一项乘了 `ck` 应为 `cj`。这会把内部节点放到错误位置，单元数越大反转越多：


| Dim（每边节点） | 单元数 | 负体积单元 | 占比 |
|---:|---:|---:|---:|
| 3 | 8 | 0 | 0% |
| 9 | 512 | 95 | 18.6% |
| 17 | 4,096 | 921 | 22.5% |
| 33 | 32,768 | 7,734 | 23.6% |
| 65 | 262,144 | 62,816 | 24.0% |

`minimal_cube` 与 `two_connected_cubes` 两种布局在 Dim=33 下给出**完全相同**的 7,734，
说明问题在生成器的单元连接/体积约定（随规模放大才显形），不是布局方向写错。
修复（`172500ca1`）后重测：Dim=9/17/33/65 **负体积单元全部为 0**，可正常求解。

**弱扩展实测（2026-11-08，每核/每 rank 固定 ≈8,000 cell）**

注意：单核 cube 是单 zone，`nIFaces = 0`，不参与 halo 交换，因此 1 核点不能当弱扩展基线；
下面以 8 核点为基准。

| cores | 每边节点 | cells | CPU ms/step（3 发中位） | CPU µs/cell | 弱扩展效率 |
|---:|---:|---:|---:|---:|---:|
| 8 | 41 | 64,000 | 58.33 | 0.911 | 100% |
| 16 | 51 | 125,000 | 60.42 | 0.483 | 96.5% |
| 32 | 65 | 262,144 | 66.30 | 0.253 | **88.0%** |

DCU 历史同口径记录（同一套 cube 与分区，1/2/4 卡 × 8 rank/卡）：
**2026-11-10 的 H2 VMFault 使该组多卡结果当前不可复现；下表保留作差异分析，
不作为当前能力或性能结论。**

| ranks / DCU | cells | HIP ms/step（3 发中位） | HIP µs/cell | 弱扩展效率 | 对 CPU 加速 |
|---:|---:|---:|---:|---:|---:|
| 8 / 1 | 64,000 | 5.06 | 0.079 | 100% | 11.5× |
| 16 / 2 | 125,000 | 6.14 | 0.049 | 82.4% | 9.8× |
| 32 / 4 | 262,144 | 7.54 | 0.029 | **67.1%** | 8.8× |

历史读数曾给出 CPU 88%、DCU 67% 的弱扩展效率，但 halo、节点占用与频率未控制，
且 DCU 多卡现受 H2 阻塞。表中的 µs/cell 是全局 cell 数归一化的墙钟，
随资源增加而下降不证明每 rank 计算效率改善；不能据此认定扩展性或归因接口瓶颈。
单发噪声可达 3–4 倍（32 rank/4 卡三个样本 5.87 / 23.32 / 7.54），必须取中位并同作业内比较。
**但这组数字有一个未控制变量：每 rank 的 halo 并不相同，见下。**

**匹配性复核（2026-11-09，等块设计，CPU 8→128 rank）**

为了让每个 rank 拿到相同的 20³ 块，按 40³ / 40×40×80 / 40×40×160 / 80³ / 80×80×160
生成 cube 并切 8/16/32/64/128 zone（1/2/4 节点 × 32 核，1000 步 ×3，斜率扣启动）：

| ranks | cells/rank | nIFaces/rank（min..max，mean） | ms/step（中位斜率） | 弱扩展 |
|---:|---|---|---:|---:|
| 8 | 8000 | 1212..1494，**1324** | 57.9 | 100% |
| 16 | 7989..8010 | 1254..2214，**1647** | 58.4 | **99.1%** |
| 32 | 7977..8027 | 1281..2459，**1786** | 61.2 | 94.6% |
| 64 | 7983..8025 | 1339..3026，**2200** | 67.8 | 85.4% |
| 128 | 7980..8028 | 1285..3165，**2319** | 78.9 | 73.4% |

读数：**cell 平衡是 ±0.3%，但每 rank 的接口面数涨了 1.75×**——METIS 没有给出等块
（2×2×8 本应是 20³）。所以旧表里的 88% 混进了分区 halo 增长。

**D1a 结构化等块分区（2026-11-09 已实现）**：`partition_type = 2` 按 cell-count 分位切
a×b×c 块（`codes/partition/{include,src}/StructuredPartition.*`，9/9 单测）。同一套
cube 重测（1/2/4 节点，1000 步 ×3）：

| ranks | cells/rank | nIFaces/rank（min..max，mean） | 块排列 |
|---:|---|---:|---|
| 8 | 8000 | **1200..1200，1200** | 2×2×2 |
| 16 | 8000 | 1200..1600，1400 | 2×2×4 |
| 32 | 8000 | 1200..1600，1500 | 2×2×8 |
| 64 | 8000 | 1200..2400，1800 | 4×4×4 |
| 128 | 8000 | 1200..2400，1900 | 4×4×8 |

现在 cell 严格平衡、接口数完全由块的位置决定（角块 3 面=1200，边块 4 面=1600，
面内 5 面=2000，体块 6 面=2400），不再是分区器的随机结果。

**但弱扩展仍不能直接用 8 当基线**，两个原因：
1. **角块效应**：8 块全是角块（3 面=1200），域变大后内部块占比上升，每 rank halo 向
   2400 饱和；这是开放边界下弱扩展的固有几何效应，不是代码缺陷。
2. **睿频**：8/16 rank 只占节点一半核，频率高于 32 rank 全核；本作业 8→16 时间
   +22%，远大于 halo +17% 能解释的部分。
所以 **99–100% 只能用“整节点 + 同 halo”的配置来主张**（例如 32→64→128 或更长的
4×4×4→8×8×8），并且表里必须带 `nIFaces/rank`。本日的 64/128 三个样本跨度 21–31%
（集群同时有其他作业），时间斜率还不足以下结论；再用专门分配重测。

**固定 rank 的端到端收益（2026-11-09，32 rank / 1 节点，同作业 3 对交错，1000 步）**：
METIS 平均 `nIFaces` 1786 → 结构化 1500（−16%）；每步 **61.32 → 58.62 ms
（中位，−4.4%，3/3 对为正）**。样本：METIS 66.10 / 61.32 / 60.87，结构化
58.61 / 59.03 / 58.62。所以 D1a 不只是测量工具：固定 rank 下它直接省约 4%，
因为步长由 halo 最大的 rank 决定，而结构化块把最坏 rank 的 halo 从 2459 降到 1600。

**接口开关对照（2026-11-09，启动无关，回答“瓶颈在哪”）**：同一套结构化
8,000 cell/rank 分区，整节点 32/64/128 rank，比较 `ONEFLOW_DIAG_NO_COMM=1` 与正常，
用分项和（排除了 1.4–1.6 s 启动；含嵌套计时，只作相对量）：

| ranks | 计算段（comm off） | 接口项（on − off） | 平均 `nIFaces`/rank |
|---:|---:|---:|---:|
| 32 | 139.5 / 140.6 | **13.5** | 1500 |
| 64 | 138.8 / 140.7 | **22.9** | 1800 |
| 128 | 140.2 / 141.1 | **27.9** | 1900 |

读数：
- **关掉接口后三段完全持平**（139–141），每 rank 计算路径的扩展是平的——不是计算、
  也不是整节点之间的睿频；
- **接口项 13.5 → 27.9（+106%），而面数只涨 27%**，多出来的是跨节点延迟；
- 8/16 rank 只占半个节点、睿频更高，不能和整节点点混在一起比（这也是早先
  8→16 “+22%” 的来源）。

结论（**2026-11-10 修正**）：上面那张表只说明"关掉接口后计算段不变"，不等于接口是
主瓶颈。短程 rocprof 与独立墙钟数据提示 host 编排/等待值得优先排查，
但不能据此把 GPU、接口和 host 百分比拼成不重叠的时间账。
接口 recv 等待的增长与邻居到达偏差相容，尚不能排除其它因素。这是 H1 的排期依据；
本轮 H1 归账结果见 §5.2，接口（A2′）和融合（A3）仍待结合等待归因评估。

DCU 侧同一套等块分区还未重跑；H2 修复并复验前，历史 67% 的结论暂停。

**端到端并行效率实测（2026-11-10，两套同时跑）**

CPU 强扩展（m6wing 294,912 cell，METIS，1→32 rank，1 步 warm + 200 步 ×2 取小、扣启动）：

| ranks | ms/step | 加速比 | 并行效率 |
|---:|---:|---:|---:|
| 1 | 1602.4 | 1.00 | 100% |
| 2 | 862.7 | 1.86 | 93.0% |
| 4 | 466.1 | 3.44 | 86.0% |
| 8 | 275.1 | 5.82 | 72.8% |
| 16 | 140.8 | 11.38 | 71.1% |
| 32 | 76.2 | 21.02 | **65.7%** |

与 2026-11-08 基线（1/2/4/8/16/32 = 1580/893/474/279/142/74 ms，32 核 67%）在同一噪声内：
**这一轮的改动没有提升 CPU 端到端并行效率**（分区器那 4.4% 是固定 rank 的收益，不是斜率）。

DCU 单卡内扩展（m6wing，1/2/4/8 rank 同一张卡，同样的扣启动方法）：

| ranks | ms/step | 加速比 | 并行效率 |
|---:|---:|---:|---:|
| 1 | 14.36 | 1.00 | 100% |
| 2 | 10.21 | 1.41 | 70.5% |
| 4 | 8.82 | 1.63 | 40.8% |
| 8 | 7.43 | **1.93** | **24.1%** |

读数：**同一张 DCU 上把 rank 从 1 加到 8 快 1.93×**。这支持优先调查 host 工作与等待，
但同时改变了分区和调用结构，不能单凭该表证明 host 编排的独占比例。
H1 要补齐同配置计时，区分 host 执行、设备等待与 MPI 等待。

**kernel 真实耗时（2026-11-10，rocprof，DTK 26.04，1 rank / 294,912 cell / 2 步）**

| kernel | calls | 总 ms | 平均 µs | 占比 |
|---|---:|---:|---:|---:|
| GradientCell | 6 | 0.672 | 112 | 19.0% |
| ReconstructFaceEuler | 6 | 0.654 | 109 | 18.5% |
| GatherFaceFluxEuler | 6 | 0.460 | 77 | 13.0% |
| UnpackBoundary | 68 | 0.371 | 5.5 | 10.5% |
| PackBoundary | 70 | 0.305 | 4.4 | 8.6% |
| CalcAndScatterPrimitiveEuler | 6 | 0.271 | 45 | 7.7% |
| UpdateBoundaryPrimitive | 14 | 0.176 | 12.6 | 5.0% |
| 其余（gradient unpack、谱半径、LHS、state update…） | | 0.62 | | 17% |

**2 步合计 3.54 ms → 1.77 ms/step 的 kernel dispatch 时间和**。与独立的 1 rank
墙钟 14.36 ms/step 相除约 12%；假设 8 rank 聚合 kernel 成本不变，再除以 7.43 得约 24%。
后者未经 8 rank 同配置 profile 验证，两者都不是直接测得的 GPU 忙碌率；
剩余时间还可能包含拷贝、等待、通信和初始化，不能全记作 host 计算。

这推翻了此前"device kernel 4.46 ms（33%）"的读数：那个数是 host 侧入队+等待，不是
GPU 执行时间。**因此 kernel 融合（A3）不是主杠杆；要动的是 host 侧的每步编排**。
另外 `PackBoundary`/`UnpackBoundary` 每步约 35+34 次小 kernel、合计占 kernel 时间 19%，
它们同时也是 host 调用点，随 host 侧一起处理。

#### 4.9 整体进度

| 阶段 | 勾选/总数 | 状态 |
|---|---|---|
| Phase A 设备驻留执行 | 1/4 | 🟨 A1 完成（**14.62 ms/step**，默认配置），A2/A3/A4 在途 |
| Phase B 数值能力 | 0/6（另有 B-exit） | 🟨 B1 Partial：Roe 单点/单 zone 已验证，多区证据待补；B2 在途，差 1.1e-6 |
| Phase C 真实算例与 IO | 0/5 | ⬜ |
| Phase D 多卡 | 0/4 | 🟥 H2 VMFault 阻塞；另缺 D1。历史弱扩展记录不可作为当前结论 |
| Phase E 架构通用化 | 0/5 | ⬜ |
| Phase F 产品化 | 1/4 | 🟨 F1 报告已发布 |
| **阶段合计** | **2/28** | A1、F1 收口；B1 按多区证据缺口改回 Partial，不抹去已完成子项 |
| 前置层 | 6/6 | ✅ 见 §4.7 |

### 5. 下一轮任务（按序，含验收）

> 本节只列"接下来一轮"的工作，用勾选表示状态；能力台账与证据长期放在 §4，不在这里加新能力项。
> 每轮结束就地更新勾选，不新增日期章节。

#### 5.1 P0 — 先把基线恢复到可引用

- [x] **G0 门禁复跑**（合入 upstream `f2fc6ea2e` 后）。2026-11-08 全绿：
      trace `24/0/0` + `max_abs_diff=0`、200 步 XB `2.331035e-18`/`aero 0.0`、
      CPU 五算例 normal 5/5 + strict 5/5、1000 步 rc=0
- [x] **A1.7+ 固定基数已归因**（详见 §4.1 A1.7b/A1.7c）。结论：短基准把 1.4–1.6 s 的启动
      摊进 ms/step，旧的 8.93 属于这类读数；按 1000 步斜率，每步成本 ≈ 同步往返次数 × 单次开销，
      与 cell 数无关。据此做的第一版改动见 A1.7c（默认生效，14.62 ms/step）
- [x] **A1.7d 默认开关复核**（2026-11-09，m6wing 8 zone / 8 rank / 1 DCU，1000 步、3 对交错）：
      `ONEFLOW_HIP_DEFER_SYNC` on/off = **13.220 / 13.651 ms/step（中位，−3.2%，3/3 对为正）**，
      默认（on）正确；同节点绝对值 13.2 低于 14.62，但跨节点/跨日不可直接比。
      分项（同作业，defer off）：`gradient` 2.28、`reconstruction` 1.95、`flux_assembly` 1.83，
      都是 host 侧计时，不能当 GPU 执行时间。旧“kernel 合计 4.46”解释已由 §4.8b 撤回；
      `macro_calc_lhs` 0.57、`macro_time_step` 0.28。启动总计约 1.6 s，1000 步口径摊薄约 1.6 ms/step

**已收口的前置修正**（稳定结论，数字见 handoff §4）

- [x] **扩展性口径**：固定网格强扩展与固定每 rank 工作量的弱扩展分别报告，
      不由单个效率值推断瓶颈；多卡能力评估采用**每 rank 固定 cell 数**（≈8,000），
      同时记录 halo 与节点占用。单 zone 算例
      `nIFaces = 0`，**不能**当弱扩展基线。
- [x] **网格生成**：`Transfinite.cpp` 的 3D 插值曾有两处复制粘贴错误（k 方向取错边值、
      `uv` 项系数错），导致网格放大后约 24% 单元负体积。已修，现在 `minimal_cube` 可
      任意放大到 65³（262,144 cell）并正常求解。
- [x] **分区配方**：`gridObj = 3` 写 `script/grid.txt`，`npartproc` 写 `script/part.txt`
      （后加载者生效）；求解算例必须补 `grid/bcRegionMap.txt`；`time_integral` 必须显式
      设 1。三条缺一都会伪装成"扩展性问题"或崩溃。完整配方见 runbook §1b。
- [x] **合入 upstream `68b1e5f6a` 后门禁复跑**：trace `24/0/0`、200 步 XB
      `2.331035e-18`/`aero 0.0`、CPU 五算例 normal 5/5 + strict 5/5。

#### 5.2 P1 — Phase A 收口（目标 ≤14 ms/step）

H1 与 H2 均已完成，下一开发顺序为 **A4 → A2′/A3**。H1 的数据推翻了
“MG 外层循环/任务派发占掉数毫秒”的假设；A4 应围绕设备所有权与同步往返设计，
不再沿用仅删除派发可到 9–10 ms/step 的估算。H2 的根因是延迟状态标志的惰性分配，
与多卡/zone 数无关；历史扩展性数字仍须用合法 fixture 重测。

- [x] **H1 host 侧每步账（Verified）**：新增仅在
      `ONEFLOW_HOST_ACCOUNTING=1` + `ONEFLOW_STAGE_BREAKDOWN_FILE` 下启用的
      `ScopedHostTimer`；MG 根域内输出 inclusive/exclusive `.host.tsv`，不新增同步，
      不改数值调用顺序。校验 rank 完整性、实际迭代次数、exclusive 耗时闭合与启动扣除。
      修复旧 `residual_zero` / `host_mirror_refresh` 漏登记，类别表长度由初始化列表推导。
      验证入口：两步 `scripts/build-app.sh`；CTest `-R '^StageProfiler\.|^HIP\.'`；
      标准 `f2-main-solver-multizone-trace.sh`、`f2-main-solver-cross-backend-residual.sh`
      与 `scripts/cpu-regression.slurm`；`summarize-host-accounting.py` 逐对核验。
      **求解器修订 `6f7a49e3d` 的证据**：HIP/计时 CTest 22/22；trace `24/0/0`；
      200 步 residual `2.331035e-18` / aero `0.0`；CPU normal/strict 各 5/5；
      50 步 HIP 稳定性通过，日志无 dropped 记录。本地计时 6/6、task 58/58、
      Python 归账与 rank 校验各 4/4；macOS CLI 因既有可执行目录定位为空，未进入 CFD。
      Linux Build 在该修订暴露几何测试库漏链 `StructuredPartition.cpp`；
      `fb246bb3a` 仅补测试库依赖，未改求解器源码。本地相关测试 21/21、4/4、9/9；该修订的 Linux Build 已通过。
      数值门禁/归账对应 `6f7a49e3d`，未将测试 CMake 修复后的构建冒充新一轮 CFD 验证。

**归账口径与结果（探索性数据，不替换 F1 正式报告）**

m6wing 294,912 cell / 8 zone / 8 rank / 1 DCU，同一分配；独立 50 步预热，
每对 1 步与 1000 步以 `(T1000−T1)/999` 扣启动，5 对 accounting on/off 交错，
CFL 0.01、LF、ilim=0、vismodel=0、显式三阶段 RK、关闭中间输出。
每对都校验 8 个 rank 的自身耗时之和，并与最长 MG 耗时对照墙钟；各 rank 时间不相加。

| 对次 | 扣启动墙钟 ms/step | 最长 rank 的 MG ms/step | 差额 |
|---:|---:|---:|---:|
| 1 | 11.5459 | 11.5847 | 0.34% |
| 2 | 11.5340 | 11.6239 | 0.78% |
| 3 | 12.0296 | 11.7915 | 1.98% |
| 4 | 11.3561 | 11.5987 | 2.14% |
| 5 | 11.6173 | 11.6969 | 0.69% |

墙钟中位 **11.55 ms/step**；差额均低于 10%。accounting on/off 的配对墙钟差中位
−0.94%，范围 −9.53%～+0.84%，不据此宣称计时加速或零开销。
每对临界 rank 的 exclusive 分项中位（**各分项中位不能相加重建总时间**）：
`macro_update_residuals` 5.10、`resident_macrostep` 自身 2.03、`interface_mpi` 1.96、
`task_execute` 自身 1.58、`interface_pack` 0.50 ms/step。
MG 控制域 + `task_plan` + `task_dispatch` 在每对内先求和再取中位仅 **0.042 ms/step**
（约 0.36%）。`task_execute` 包含处理函数和等待，不是纯命令派发成本。
MPI trace 的已覆盖调用中，稳态等待主要在 `recv_char`；该工具未覆盖 `MPI_Bcast`，
不能将缺少某个 kind 解读为所有集合操作均不存在。

**下一刀的候选位置**：H2 修复后，先对
`NsResidentMacroStep` 尾部 `HipFluxBackend::CheckResidentStatus()` 的阻塞状态读取/清零，
以及 RHS 中 `UNsGrad` 的 `swap_download` 等待做设备时间与 host 等待对账。
宏步自身的 2.03 ms 包含尾部检查，不能直接认定全部可省；
`swap_download` 也会承接延迟队列工作。保持每步非物理状态检查与 post 轮 lag 语义，
再评估合并状态传输/减少同步，不先重写 MG/task 骨架或删检查。
原始日志、各 rank TSV、五个 `summary*.json`、`timing.tsv`、分析脚本与结果留在本轮
`evidence/gate/`、`evidence/cpu-regression/`、`evidence/probe/` 的对应目录。

- [x] **H2 多卡 VMFault（已完成，`d41650729`）**：**根因不是多卡，而是 `ONEFLOW_HIP_DEFER_STATUS=1`（默认）**
      下 `stateUpdateStatus` 缓冲区惰性分配；deferStatus 分支完全跳过 `CopyFromHost`，
      `Size()==0`、`Data()` 返回 **nullptr** 传给内核里的 `atomicExch` → 设备内存故障。
      **单卡也会触发**；表现为 `libgalaxyhip`/`hsa-runtime64` 里 signal 6，
      回溯指向恰好在队列里的 `GradientCopyBoundaryKernel`。
      次因：最小复现算例（2 zone m6）CFL 过高（`cflst=0.1 → cfled=10 over 100 steps`），
      step 2/13 就产生非物理状态；但因为主因把错误路径变成了 VMFault abort，掩盖了真实问题。
      **修复**：首次调用（`Size()<1`）仍然分配并上传零标志；后续调用保留设备侧标志直到
      `CheckResidentStatus()` 读取并清零。
      **验证**（固定 CFL=0.01，2 zone m6，50 步）：
      - 单卡 HIP：PASS
      - 多卡 HIP（2 节点/2 DCU）：**PASS**，无 VMFault
      - CPU 2 rank：PASS
      - HIP vs CPU 残差：**逐位一致**（`evidence/probe/127626632`）
      **门禁 Verified（`d41650729`，2026-11-11）**：gate 作业 `127628802`
      trace `24/0/0` + `max_abs_diff=0`、200 步 XB `2.331035e-18`/`aero 0.0`、
      3×1000 步 rc=0（14.53/14.58/14.58 ms/step）；
      CPU 回归作业 `127628910` normal 5/5 + strict 5/5。
      **扩展性数据仍需重测**：原历史 16/32 zone 的 1/2/4 卡数据可能落在同一缺陷路径上，
      不能直接引用；须用合法 fixture（D1）在该修订上重跑弱扩展。

- [x] **A1.9 梯度 unpack 路径的两处重复工作**（2026-11-08）：ghost-face→cell 映射按
      (grid, ghost level) 缓存（`006f9ff39`+`6174b21bd`），unpack 的 ghost id 列表按指针复用
      上传（`bbecf8241`）。门禁全绿（trace `24/0/0`、XB `2.331035e-18`）。
      500 步分项：`gradient` 5.27 → 4.87、`swap_download` 5.04 → 4.55；
      但 1000 步墙钟在噪声内（14.61–15.56），**尚未换成端到端收益**。
      **2026-11-08 复核（1/2/4 卡，`e07r2n03`，1000 步 ×2）**：8 zone/1 卡 17.37/14.61，
      16/2 卡 10.64/10.76，32/4 卡 10.27/10.21 → 1→2 卡 1.37×、2→4 卡 1.04×，
      扩展性没有变化。分项 `swap_download`：8 zone 4.45、16 zone 2.60、32 zone 1.09。
      结论修正：之前的"4.5 ms 是每 stage 的 D2H"判断不成立——设备 unpack 路径本身不拷贝，
      这个计时器是**延迟队列在下一个阻塞点的等待**，随 rank 数下降正好说明它承接的是
      排队工作量。所以"安全桶 5–6 ms"的估计过高，实际只兑现约 0.4 ms（分项），
      端到端在噪声内。剩下的差距是结构性的（A4 把循环下到设备），不是缓存能解决的
- [ ] **A3 stage 融合，第一段（Deferred）**：`gradient → reconstruction → flux` 串进设备，
      stage 内减少 host 同步。旧 −4/−6 ms 估算及 4.46 ms 的“设备 kernel”读数已撤回：
      前者基于旧同步口径，后者包含 host 入队与等待。rocprof 的 1.77 ms 是全部 dispatch
      的每步时间和，不只是这三个 kernel；融合收益待 H1/H2/A4 后重估。
      验收：trace `24/0/0`、XB `2.331035e-18`/`0.0`、同节点交错 paired 有可测收益
- [ ] **A2′ 接口交换（Deferred）**：单卡墙钟通信开关 A/B 差值曾为 **2.1–2.2 ms/step**，
      约占该次墙钟 16%，不能与其它配置的 kernel 比例相加。多卡历史接口分项受 H2
      与等待归因影响，暂不据其排期。Isend/Irecv 尝试已否，见 handoff 坑 30。
      2026-11-08 已验证两条路：
      - **GPU-aware 设备直发（2026-11-15 复核：旧"不可用"结论作废）**：
        旧记录说 2 rank 探针 PASS 但求解器 8/32 rank abort（rc=134）。
        **那次 abort 在 H2 修复 `d41650729` 之前，根因是 `DEFER_STATUS` 下状态缓冲
        空指针 VMFault，与 GPU-aware 无关**（同类 abort 在 host 路径同样触发）。
        `probe/128116476` 在当前构建重测：2 rank 设备指针 **pass=1**（UCX 默认传输）；
        显式加 GPUDirect 提示反而失败；**4/8 rank 的 FAIL 是探针程序写死
        `requires exactly two MPI ranks` 所致，不是 MPI 失败**。
        因此 **8 rank/卡 生产密度是否可用，至今无证据**——需要写一个多 rank 探针。
        开关仍保持关闭，理由改为"未验证"而非"已验证不可用"。
      - 减轮次或聚合消息仅作为待验证方案；MPI 计数未显示消息量随 rank 增长，
        尚无证据证明聚合或 allgather 能改善端到端性能
      - **轮次 4→3 已再次验证不可行**：`ONEFLOW_SKIP_POST_INTERFACE=1` 后
        200 步 XB `5.540033e-07`、aero `4.0e-05`（fail），但 1000 步从 14.62 降到
        **13.73**（−0.9 ms）。即第 4 轮买 0.9 ms，却承载着 legacy 的一步滞后，
        要削它必须找到"同一滞后、少一次往返"的等价写法，不能直接删
- [ ] **A4 HIP `EulerDomainBackend::Advance`（下一主线）**：RK 阶段循环进 backend，只在宏步边界回主机。
      H2 已解除多卡阻塞，H1 已给出下沉边界证据；保持状态检查与 halo lag，
      以减少实际同步/传输验证收益。验收须包含完整数值门禁与同节点 paired。
      **前置（必做）**：先在 `d41650729` 上跑完 H2 遗留的完整门禁（HIP 契约 + CPU 五算例 + trace + XB），
      否则 A4 的数值回归没有干净基线可对照。
- **A2 已降级，不排期**：单卡关全部通信的差值约 2.2 ms；可兑现收益待归账重估

**接口同步点优化（`d44e0253b`，2026-11-11）**：两处结构性 host 往返已收：
同 rank 接口改为设备直传（pack → D2D copy → unpack，去掉 D2H + std::copy + H2D
以及长度握手的 MPI 自环），宏步末尾的状态读取改为“请求一次、由首个已排空的接口
pack 消费”，不再为 4 字节单独排空队列。
**结果**：XB 200 步 `2.331035e-18` / aero `0.0` 逐位不变；1000 步三发
`14.39 / 14.42 / 14.42` ms/step，对比 `d41650729` 的 `14.53 / 14.58 / 14.58`，
**约 −1.1%（样本更集中，但量级接近噪声）**。
**原因（重要）**：同 rank 直传只对“一 rank 多 zone”生效；生产配置是 8 rank × 1 zone，
**没有任何同 rank 端口**，因此该项在生产口径不产生收益。状态读取搭车只去掉若干次
排空中的一次。
**尚未定位**：关通信只省 2.2 ms、GPU 只忙 1.6 ms，但扣除启动后每步仍有约
5–8 ms host 侧工作没有归因；实测 `UpdateBoundaryPrimitive` 之后的停顿约 4.0 ms/步
（1 rank/8 zone 口径），而 RK stage 内部的 kernel 间隙仅 4.8 µs（已接近下限，
**A3 融合无收益**）。下一刀应给 `NsResidentMacroStep` 结束段与 post-process 加计时，
先定位这段 host 工作，而不是继续在接口上微调。

#### 5.3 P1′ — Phase B 并行（与 A 不抢同一段代码）

- [ ] **B2 收口**：按 B2.0 的差异清单（φ 0.45%、面值 0.9%）核对 Barth 邻域 stencil，
      定位 `1.113e-06`；若确认是主机 `ckfun` 的整体丢弃语义，设备侧补同样的按面判定。
      过 200 步跨后端与 CPU 五算例后，再压 `FillNeutralLimiterField`（当前 `ilim=1` 仍多约 39 ms/step）
      验收：`ilim=1` 200 步 XB ≤`1e-12` + CPU 五算例；之后才讨论放开默认
- [ ] **B1.3 补证据**：设备 Roe 目前只有单 zone 证据（契约 3/3、200 步 `res_max_rel = 0.0`），
      补 3 步多区 trace 与 CPU 五算例的 Roe 口径，否则 B1 不能写进 PR 结论
- [ ] **E1 热路径去 ambient globals**（§4.5 E1）：先设备后主机；**B3 的前置**，早点做可以避免 globals 扩散
- [ ] **B3 层流粘性**：面粘性通量、温度梯度、壁面热流（依赖 E1）
      验收：同 §4 的数值门禁 + 粘性算例 trace

#### 5.4 P2 — 外部依赖与后续

- [x] **D1a 结构化等块分区**（2026-11-09）：`partition_type = 2`，按 cell-count 分位切
      a×b×c 块（`StructuredPartition.*`，9/9 单测）。8→128 rank 每 rank `nCells` 严格
      8000，`nIFaces` 由块位置决定（1200/1600/…/2400），不再随 METIS 随机波动。
      验收前半（cell 相同）已达成；**后半（用它重跑 CPU/DCU 弱扩展）待专门分配**：
      本日样本跨度 21–31%，无法下结论，见 §4.8b
- [x] **D1b 等 halo 弱扩展（部分，2026-11-14）**：M6 首个固定 36,864 cell/rank 弱扩展点
      已测（8 rank/1 卡 vs 64 rank/8 卡，同作业交错中位，17.0 → 22.2 ms/步，§5.6）；
      但两侧 nIFaces/rank 不等（METIS 4.4k vs 6.3k，+43%），“扣除 halo 后每 cell 计算时间
      为常数”的验收还需 partition_type=2 的等块分区版本；只在整节点配置间比的要求未满足（64 rank 跨 2 节点）
- [x] **D1 合法 fixture**（2026-11-14 完成）：`fixtures/m6-refined-r1/`（集群）登记 2,359,296 cell
      的 M6 加密网格及 8/32/64 zone 分区（METIS，manifest 验证通过）；单 zone XB 残差 3.07e-17、
      64-zone XB 全绿、加密 vs 原 CPU 气动力一致到 0.2%。制作链路见 §5.6，工具
      `ci/kunshan/refine_ofl.py`。拿到之前 D2/D3 只做只读评估、不出数字。历史 1/2/4 卡弱扩展记录受 H2 阻塞的结论已由 H2 修复解除，
      但那些数字仍是小子域读数（每 rank 36,864 cell），不能当弱扩展证据
- [ ] **B4 SA**（依赖 B3）
- [ ] **B5 LU-SGS**（依赖 A4）
- [ ] **B6 multigrid**（依赖 A4）
- [ ] **B-exit**：`test/cfd/m6wingroe_sa`（Roe + LU-SGS + SA）在设备上过同一组门禁；
      B1–B5 全绿才可勾，这是旗舰算例离开 CPU 的唯一途径

C1–C5、D2–D4、E2–E5、F2–F4 均未启动，定义见 §4，不在本轮排期。

#### 5.5 交接 2026-11-14：扩展性已定性，网格能力是唯一缺口

**基线**：`dev @ 727c1ed75`。本轮交付的提交（按序）：

| 提交 | 内容 | 证据 |
|---|---|---|
| `d41650729` | H2 根因：`DEFER_STATUS` 下状态缓冲惰性分配被跳过 → 空指针进 `atomicExch` | trace `24/0/0`、XB `2.331035e-18`、CPU 5/5 + 5/5 |
| `d44e0253b` | 同 rank 接口改设备直传 + 状态读取搭车已有同步点 | XB 逐位不变；np8 端到端 −1.1% |
| `970508e7c` | `ONEFLOW_HIP_LOCAL_DEVICE_EXCHANGE=0` 回退开关 | — |
| `2006b0a84` | 端点拷贝改 pinned；修正 post 轮的错误注释 | np8 无变化，np1 −5% |
| `d9e6a14db` | plan 与 per-rank exchange 子集按拓扑指纹缓存 | XB 逐位不变；256 zone 时迭代 4.44 → 2.00 ms/步 |
| `727c1ed75` | `CMD::ExecuteCmd` 按命令计时（诊断用，未开 breakdown 时零开销） | 定位到 `DUMP_RESIDUAL` |

**结论一：扩展性由每 rank 的 cell 数决定，不由卡数决定。** 两套几何、两个方向一致：

| 口径 | 1 卡 | 2 卡 | 4 卡 | 8 卡 | 16 卡 | 32 卡 |
|---|---:|---:|---:|---:|---:|---:|
| M6 强扩展 端到端 ms/步 | 11.70 | **7.60** | 9.87 | 18.73 | 24.08 | 30.62 |
| M6 核心 `resident_macrostep` | 2.70 | 2.23 | 1.45 | 1.41 | 1.22 | 1.08 |
| 立方体强扩展 端到端 | 17.53 | **14.81** | 16.08 | 24.59 | 24.70 | 39.29 |
| 立方体核心 | 7.07 | 4.48 | 2.78 | 1.66 | 1.19 | 1.18 |
| 立方体弱扩展 迭代（110,592 cell/rank） | 23.07 | 13.78 | 13.93 | 14.06 | 12.69 | 13.68 |

- **求解器核心能扩展**：M6 2.5×、立方体 6.0×（1→32 卡）。
- **交换不能扩展**：它在总时间里从约一半涨到主导。**墙钟最低点落在 2 卡，但这是退化不是最优**：
  每卡吞吐随卡数单调下降（22.9 → 17.1 → 11.0 kcell/ms/卡），加卡越多单卡效率越低；
  正确目标是让每卡吞吐随卡数保持（接近）恒定，而不是接受 2 卡后的退化。
- **弱扩展在迭代层近乎持平**（立方体，2→32 卡 13.78→13.68）。
- 经验阈值：每 rank ≳5 万 cell 正常；≈1 万开始变差；<5 千明显变差。

**结论二：瓶颈不在通信带宽**（`interface_mpi` 占 7%，rank 数 ×16 只涨 25%）。

**结论三（关键）：M6 的扩展性做不了，是因为缺更大的 M6 网格，不是网格质量问题。**
仓库无细化工具；`GridConversion` 里 `OneFLOW` 只作为 target（`IsSupportedConversion`），
`.ofl` 无导出路径；集群无 gmsh。所以下一步必须是**自己写原生加密**。

**已验证的 `.ofl` 二进制布局**（读 8-zone M6 实测，字节闭合差 8）：

```
int32  nZones
int32  pid[nZones]
int32  zoneType[nZones]
每 zone:
  int32 zoneBytes, int32 ?, int32 nNodes, int32 nFaces, int32 nCells
  double xN[nNodes], yN[nNodes], zN[nNodes]
  int32 volBcType
  int32 numFaceNode[nFaces]
  int32 faceNodeMem[sum(numFaceNode)]      # 0-based
  int32 leftCells[nFaces], rightCells[nFaces]
  int32 nBFaces, bcType[nBFaces], bcNameId[nBFaces]
  int32 nIFaces
  (nIFaces>0) int32 zoneId[], localInterfaceId[], i2b[]
```

实测值（8-zone M6）：`nZones=8`、每 zone `nNodes=40,817 / nFaces=114,432 / nCells=36,864 /
nBFaces=7,680 / nIFaces=4,785`，**114,432 个面全部是四边形 → 六面体网格**；总计 294,912 cell。

**下一轮唯一主线：原生六面体加密（`refine_ofl.py`）**

1. 读 8-zone `.ofl`（上述布局，含 `zoneBytes` 校验）
2. 每个 hex 从 6 个面恢复 8 个角点 → 加 12 棱中点 + 6 面心 + 1 体心 → 8 个子 hex
3. 节点去重：棱 `(n1,n2)`、面排序 `(n1..n4)`、体心按 cell
4. 子面重编号、left/right 重映射；边界子面沿用 `bcType/bcNameId`；界面子面沿用 `zoneId`
5. 按同一格式写回（用 `zoneBytes` 自检）
6. 验证顺序：解析自检 → 单 zone 200 步与原始 M6 比残差 → 分区成 64 zone（每 rank 36,864 cell）测弱扩展

预期：294,912 × 8 = **2,359,296 cell**，仍是**一个连通算例**（不是拼接），8 个 zone 的界面连通性保持。

**本轮新踩的坑（务必不要再踩）**

1. **每步残差/气动力输出会让多卡测量失真**：fixture 的 `cfd.txt` 是 `nResSave=1`、`nForceSave=1`，
   256 rank 时 `DUMP_RESIDUAL` 每步 37.2 ms，占整步 93%。基准必须走 `prep-case.sh`
   （它把三个 `*Save` 设成 100000），或自己显式设。**本轮前半的所有 M6 扩展性数字都因此作废。**
2. **`sed` 改参数必须用 `[[:space:]]+`**：fixture 里是 `int    nResSave = 1;`，写成 `^int nResSave = .*` 不匹配，会静默失败两次。
3. **`--dcus` 是每节点、上限 4；`--cpus` 上限 32/节点**；`mpirun --host` 必须带槽位（`h1:32,h2:32`）。
4. **`NDEV = local_rank % 4`**（每节点 4 卡），不是 8。
5. **GPU 占用率必须按分配到的卡号读**，不能取 `/sys/class/drm/card*/` 的最大值——同节点别的租户会把读数顶到 100%。
6. **作业脚本里用 `python3` 必须先 `source scripts/env.sh`**，否则 `command not found`。
7. **`RecordTask` 是有调用者的**（`CmxTask.cpp:70` 的 `TaskDispatchTimer`）；`tasks.tsv` 只记 dispatch 级，
   命令级要 `CMD::ExecuteCmd` 里补（`727c1ed75` 已做，未开 breakdown 时零开销）。
8. **单 zone 的 `m6wing.ofl` 是旧格式**（zone 头多 2 个 int），当前 reader/writer 用 3-int 布局；
   做加密请以 `m6wing-8.ofl` 为准。
   **修正（2026-11-14）**：仓库内 `test/cfd/m6wingroe_sa/grid/m6wing.ofl` 逐字节闭合解析就是
   当前 3-int 布局（nLength 精确、EOF 零剩余）；加密工具已直接读它并通过全部自检，
   无需迁移到 8-zone 文件。

**不要做**：靠加卡提升（§4.8b 已说明）；"一 rank 多 zone"（A1.8 实测更慢）；
继续端到端盲试（曾因读数噪声或归因错误返工三次）。

#### 5.8（撤回更正）：定向合并 2N/4N **可行**；此前"严格排除"是两个实现缺陷的产物

**2026-11-15 撤回**：§5.8 此前记录的"严格 UNSAT / O 型接缝本质不可积"结论**错误**，
原因是两个生产实现缺陷（均由用户提供的回归定位）：

1. `analyse_face` 按坐标集合匹配恢复对应，丢失节点身份——坐标恰好相同的两个
   **不同**节点会被静默互换。已改为按真实 node ID 逐点对应
   （`latL[n] == T(latR[n]) + offset`，四角全验证）；集合版 `affine_ok` 已删除。
   回归：手工构造的期望映射 `0→0, 1→2, 2→1`（此前被错判为恒等）现已正确恢复。
2. 带置换并查集把可满足网络误判为 UNSAT（缩减反例 + 已知合法分配）。
   已替换为确定性求解器：每连通分量枚举 3 个根方向，BFS 沿正/逆置换传播，
   逐边检查 + 全边复核；三个种子全部失败才判 UNSAT，并记录每个种子的阻断边。
   回归：已知合法分配 PASS；边序重排不变；孤立单元处理正确；
   200 个小图与穷举完全一致。

**修正后的全新核验**（不复用旧允许表；873,216/873,216 张内部面全部解析为双射约束）：
**2N 模式 PASS、4N 模式 PASS**；第二次全网格独立复核 2N 分配
（873,216 面逐面）**0 违反**。此前 1128/1104/2256/2280"反例环"的判定作废。

状态与 N 系列记录（分开标记，不混写）：

| 网格 | 已生成 | 结构验证 | 分区 | 数值 |
|---|---|---|---|---|
| N（294,912） | ✅（仓库原始） | ✅（reader 逐位闭合） | ✅ METIS 8/16/32/64 zone + manifest | ✅（生产算例、XB 门禁历史基线） |
| 8N（2,359,296） | ✅ | ✅（计数/体积 5e-15/重读） | ✅ METIS 8/16/32/64 + manifest | ✅（单 zone XB 3.07e-17；64-zone XB 全绿；CPU 气动力 vs N 差 0.2%） |
| 64N（18,874,368） | ✅ | ✅（结构+体积守恒） | ✅ METIS 512 zone | ✅（512 rank/64 卡 1000 步 ×3 跑通，`probe/128104413`；原 VMFault 判定为 i02 节点态，非代码缺陷） |
| 2N / 4N / 16N / 32N | 生成中（拓扑可行性已 PASS） | — | — | — |

#### 5.10 进行中 2026-11-15：M6 DCU 每步等待归因（H1/H2/H3）

**目标**：回答"慢 rank 的 interface_mpi 到底在等什么"，只做归因 + 一个最小干预，
不同时启动 A4、分区器重写、kernel 融合或 GPU-aware MPI。

**已完成（只读审计 + 两个判别作业）**

1. **调用链已核实**（`NsResidentMacroStep.cpp` → `InterfaceTask.cpp` →
   `DistributedDomainPlan.cpp` → `UNsGrad.cpp` → `HipFluxBackend.hip`）：
   每 RK stage 末尾一次 `CUPDATE_INTERFACE`，跨 rank **只交换 q**
   （`ExchangeHipInterfaceQ`；host DataBook 交换被 `ONEFLOW_HIP_SKIP_HOST_INTERFACE`
   跳过，`interface_host_*` 计时恒为 0）。梯度 halo 不每 stage 跨 rank：
   接收侧每 stage 解包一份**每步只打包一次**的载荷。`NsPostprocess` 另有每步一次、
   含 q+梯度的第 4 轮，保留 legacy 一步滞后。
2. **H3（rank→DCU 映射/租户争用）已证伪**：`evidence/probe/128108949` 用 HIP 运行时
   （`hipGetDeviceCount` + `hipDeviceGetPCIBusId`）验证 `rank.sh` 的卡映射正确；
   ROCR 清掉后下标越界会直接报 `no ROCm-capable device is detected`，不会静默跑错卡。
   同批探测确认 `rocminfo` 不理会这两个变量，**不能用它判断设备视图**。
3. **分项表逐 rank 横比无效**（补强坑 32）：每 rank cell 数相同（36,864）时
   `gradient_kernel` 差 5.46 倍。原因是 `hipDeviceSynchronize` 是设备级同步，
   drain 落在哪个 rank 就替同卡全部 rank 背排队时间。

**判别作业**（`probe/128109059`，1 节点 / 8 rank / 1 卡 / 1000 步 / 三变体 ×3 发交错）：
`a` 默认、`b` `DIAG_NO_COMM=1`（仅计时）、`c` `DEFER_SYNC=0`。判据：
`b` 的 per-rank 计算段离散度塌缩 ⇒ H1 证伪；`a` vs `c` 的差异量 ⇒ 计时污染幅度；
per-rank MPI trace 的 `recv` 秒数是否覆盖 `interface_mpi` ⇒ 等待在 MPI 还是设备拷贝。

**结论（2026-11-15，`probe/128109059` + 对照 `probe/128109153`）**：

| 变体 | 1 卡 8 rk/卡 | 2 卡 4 rk/卡 |
|---|---:|---:|
| a 默认 | **15.41** ms/步 | **9.00** ms/步 |
| b 通信全断 | 11.41（离散度 1.02×） | 6.85 |
| c `DEFER_SYNC=0` | 13.95 | — |

- **H1 证伪**：通信关闭后 8 个 rank 的 `macro_update_residuals` 为 7,309–7,420 ms
  （离散度 **1.02×**），每 rank 36,864 cell 完全一致。分区没有把负担分歪。
- **H3 证伪**：卡映射正确（见 §5.10 第 2 条）。但其**近亲**成立并被量化：
  同 8 rank、同 8 zone、同分区、只把卡从 1 张变 2 张（rank/卡 8→4），
  r3..r7 的 `recv` 等待均值 **2,455 → 1,192 ms（2.06× 少）**，rank3 单消息延迟
  **1.031 → 0.412 ms（2.50× 少）**。消息数与字节数不变，所以这是**同卡进程争用
  拉长 rendezvous 往返**，不是带宽。
- **H2 成立为主因**：`recv` 占 `interface_mpi` 的 **96.3–98.8%**；1000 步
  `send` 合计 414 ms 而 `recv` 18,716 ms（**45×**）。~420 MB 走 18.7 s ≈ 22 MB/s，
  远低于卡间带宽，说明每次都在等对端到齐而非传不动。邻居数与等待**负相关**（−0.64）：
  rank3 每步 3.02 次 recv、1.031 ms/次（最慢）；rank2/6 每步 6.03 次、0.275/0.379 ms/次。
- **`DEFER_SYNC=0` 更快 1.45 ms/步**（13.95 vs 15.41）：当前 8-rk/卡配置下
  延迟同步是净负债。同步模式下 rank3 的 `recv` 仍为 3.34 s，所以等待不是计时伪影。
  **该结论已由 §5.11 的三密度筛选推翻**（三发一致地更慢），此处保留为历史读数；
  差异来源是旧作业含分项+MPI trace 诊断，新筛选完全关诊断、每配置独立全新 case。
- **每卡吞吐随资源增加而下降**：19,141 → 16,386 cell/ms/card（−14%），
  即 2 卡对同一工作只有 1.71×（理想 2×）。
- **通信成本口径更新**：默认 − 断通信 = **4.0 ms/步（26%）**（1 卡），
  取代历史"关通信省 2.2 ms"的旧口径（后者是不同配置、不同步数）。

**未定位（唯一必要实验）**：`DEFER_SYNC=0` 省下的 1.45 ms 与 MPI 等待是否重叠。
`c` 变体未落 MPI trace；补跑 `c` 的 per-rank trace 即可判定"先动 drain 还是先动交换"。

#### 5.11 阶段一交付 2026-11-15：运行配置收敛（M6 N 与 8N）

**口径**：每配置独立全新 case 目录（`prep-case.sh`），1000 步，`date` 包住 mpirun 的
墙钟 ÷ 1000；不扣启动（启动对每个配置相同，配对时抵消）、不开分项、不开 MPI trace。
配置在每发内**交错**执行，节点漂移对所有配置同样作用。三发，报告离散度。

**单节点筛选**（`probe/128110089`，N=294,912 cell / 8 zone / 8 rank / 1 zone per rank，
1 节点 4 卡，三发中位）：

| 配置 | ms/步 | 离散度 | 相对默认 | 用卡 | 总吞吐 kcell/ms | **每卡吞吐 kcell/ms/卡** |
|---|---:|---:|---:|---:|---:|---:|
| 8 rank/卡 `DEFER_SYNC=1`（现默认） | **12.900** | 1.011 | 1.00× | 1 | 22.9 | **22.9** |
| 8 rank/卡 `DEFER_SYNC=0` | 13.399 | 1.003 | 0.96× | 1 | 22.0 | 22.0 |
| 4 rank/卡 `DEFER_SYNC=1` | **8.643** | 1.013 | 1.49× | 2 | 34.1 | **17.1** |
| 4 rank/卡 `DEFER_SYNC=0` | 9.274 | 1.005 | 1.39× | 2 | 31.8 | 15.9 |
| 2 rank/卡 `DEFER_SYNC=1` | **6.685** | 1.004 | 1.93× | 4 | 44.1 | **11.0** |
| 2 rank/卡 `DEFER_SYNC=0` | 7.098 | 1.008 | 1.82× | 4 | 41.5 | 10.4 |

**重要澄清：本表未测"每卡放几个 rank 最优"**。这里 rank 总数固定为 8，
所以"rank/卡 降低"必然同时"卡数增加"：12.90 → 8.64 是**多了一张卡**的结果，
不是每卡 rank 更少的功劳。固定卡数下的每卡 rank 数另有历史结论
（§4.8b 单卡 m6wing，1/2/4/8 rank = 14.36/10.21/8.82/**7.43** ms/步），
即**8 rank/卡 是单卡墙钟最优**。两类实验不能混成一个 rank/卡 因果结论。

**两条可交付结论**：
1. **`DEFER_SYNC=0` 在三个密度下一致更慢**（每个密度三发同向，差异 3–7%，远大于
   1.3% 的离散度）。**默认保持不变**，§5.10 中"更快 1.45 ms/步"的读数作废
   （那次含诊断开销，非同一口径）。
2. **加卡能换时间，但每卡吞吐单调下降**：8→4 rank/卡（1→2 卡）让固定问题从
   12.90 降到 8.64（−33%），**每卡吞吐 22.9 → 17.1 kcell/ms/卡（−25%）**；
   8→2 rank/卡（1→4 卡）降到 6.69（−48%），每卡 11.0（−52%）。
   即用 2 倍卡换 1.49 倍速、4 倍卡换 1.93 倍速，**总吞吐上升、单卡效率下降**。
   这是用资源换时间，不是效率提升；选它只在 wall time 优先且卡预算宽裕时成立。

**跨节点**：8N 上 64 rank / 8 rank/卡实测 **29.85 ms/步**（与历史 22.32 差 34%，
跨作业噪声内不可比）；32-zone/32-rank 配置在网格读取阶段**挂起**（只打印
`prjBaseDir`，无任何 rank 完成 zone 读取），记为新缺陷，待查。
同预算的粗细粒度配对（64 rank/8 卡 vs 16 rank/8 卡，见 `probe/128110470`）在跑。

**跨节点归因（`probe/128110682`，2 节点 / 64 rank / 8 卡 / 8N / 三个交错发）**：

| 口径 | 单节点 8r/1 卡 | 跨节点 64r/8 卡 | 差 |
|---|---:|---:|---:|
| 通信 ON | 12.900 | **20.596**（20.286/23.970/20.596） | +7.70（1.60×） |
| 通信 OFF | 11.408 | **16.181**（16.181/17.565/14.335） | +4.77（**1.42×**） |
| 通信份额 | 1.492（11.6%） | 4.415（21.4%） | |

**分解**：跨节点惩罚 7.70 ms/步中，**4.77（62%）在通信完全关闭时依然存在**，
只有 2.92（38%）可归给通信。即：**同一份代码、同样 36,864 cell/rank，
跨节点在没有任何接口交换时仍比单节点慢 1.42 倍**。这部分是节点计算/资源状态，
交换层代码无法回收。**这否定了"改通信调度能拿回 1.6×"的预期**。

跨节点单发离散度 1.18–1.23，**远大于单节点的 1.003–1.013**，所以跨节点数字
只能在同一作业内用时序配对比较，不能跨作业引用（这也解释了此前 29.85 与历史
22.32 的 34% 差异：单发噪声）。

**已知缺陷（非主线，记录备查）**：32-zone/32-rank 网格在两节点上启动时挂起
（只打印 `prjBaseDir`，无 zone 读取完成）；同一网格在 `probe/128021420` 里
以 32 rank / 4 卡**单节点**跑通过 33.07 ms/步。差异是节点数与 mpirun 布局，
不是网格文件损坏。

#### 5.9 交付 2026-11-15：2N/4N/16N/32N 生成并全阶梯弱扩展已测

**网格生成链**（全部从 N 派生，保留 M6 几何、边界定义与测量家族）：
- `ci/kunshan/agglomerate2_ofl.py`：2N/4N = 8N 按方向分配合并（分配在 N 网格上求解，
  PASS；A(c) 映射与真实 face ID 驱动，面环由无向边界走线恢复，角点按出现次数
  唯一性选取、环序取走线序）；
- 16N = refine8(2N)，32N = refine8(4N)（已验证的八倍细化器）。

**每个尺寸的验证状态**（全部重读自磁盘）：

| 网格 | cell（目标） | 生成 | 结构 | 几何 | 分区 | 数值 |
|---|---|---|---|---|---|---|
| 2N | 589,824 ✓ | ✅ | ✅ 六面拓扑/无自连接/面积闭合 3.1e-15 | ✅ 体积 vs N 1.06e-12 | ✅ 16 zone | ✅（本作业内 2 卡档实测跑通） |
| 4N | 1,179,648 ✓ | ✅ | ✅ 面积闭合 4.9e-15 | ✅ 体积 1.99e-14 | ✅ 32 zone | ✅（4 卡档） |
| 16N | 4,718,592 ✓ | ✅ | ✅（refine 自检全过） | ✅ 体积 1.12e-12 | ✅ 128 zone | ✅（16 卡档） |
| 32N | 9,437,184 ✓ | ✅ | ✅ | ✅ 体积 7.99e-14 | ✅ 256 zone | ✅（32 卡档） |

BC 身份逐级保留；pinch 面直方图 {1: 大量, 2: O(pinch)}，无超过 2 的对。
2N/4N 的直切四边形角点相对细网格并集存在 g-凸包偏差（分析和为 0.67 vs V=886.8，
**带符号总量相互抵消到 1e-12**），已记录在工具提交说明中。

**M6 DCU 全阶梯弱扩展**（`probe/128098635`：36,864 cell/rank 恒定，8 rank/卡，
同作业 3 发交错中位，1000 步；64 卡档为 `probe/128104413`，另作业不可直接横比，
但同口径同 cell/rank）：

| 卡 | rank | cell | 三发 ms/步 | 中位 | 相对 1 卡 |
|---:|---:|---:|---|---:|---:|
| 1 | 8 | 294,912 | 14.66/14.75/18.66 | **14.75** | 1.00× |
| 2 | 16 | 589,824 | 15.52/15.78/16.27 | **15.78** | 1.07× |
| 4 | 32 | 1,179,648 | 17.57/17.90/18.55 | **17.90** | 1.21× |
| 8 | 64 | 2,359,296 | 21.45/22.05/22.32 | **22.05** | 1.50× |
| 16 | 128 | 4,718,592 | 33.08/33.22/34.42 | **33.22** | 2.25× |
| 32 | 256 | 9,437,184 | 35.55/38.86/38.96 | **38.86** | 2.63× |
| 64 | 512 | 18,874,368 | 61.74/51.43/56.94 | **56.94** | 3.86× |

阶梯各档跑通本身即各档网格的端到端数值验证（求解器完成 1000 步、无 NaN/非物理中断；
prep-case 口径）。弱扩展从 2 卡的 1.07× 到 64 卡的 3.86× 单调退化，16 卡（4 节点）
跳变最大（22→33 ms，+50%），与既有归因（跨节点 interface_mpi）一致；64 卡（16 节点）
对 32 卡再 +46%。

#### 5.7 支线（暂停）：gmsh 网格阶梯

**2026-11-15 暂停**（按用户决定：本轮比赛无对应模型，先回主线）。已确认与已验证的事实：

- 集群计算节点**无外网**，`pip install gmsh` 不可用；预编译 `gmsh-4.13.1-Linux64` 需要
  **glibc 2.23**，集群 CentOS 7.6 是 2.17，无法运行（`/lib64/libm.so.6: GLIBC_2.23 not found`）。
- 源码编译路线可行且已走到 66%：需 `-DCMAKE_POLICY_VERSION_MINIMUM=3.5`（cmake 4.4.3 拒绝旧
  policy）、`-DENABLE_FLTK=OFF -DENABLE_BLAS_LAPACK=OFF -DENABLE_HXT=OFF`（HXT 编译失败且本任务
  不需要四面体）。**当前失败点**：内置 METIS 的 `GKlib/string.c` 里 `strptime` 隐式声明，
  修法是加 `-DCMAKE_C_FLAGS=-D_GNU_SOURCE -DCMAKE_CXX_FLAGS=-D_GNU_SOURCE`。
- 模块目录可写：`/public/home/luql/apprepo/modules`（与 `module avail` 同一批），安装前缀定为
  `/public/home/luql/apprepo/gmsh/4.13.1`，modulefile 写 `gmsh/4.13.1`。脚本见集群
  `scripts/of-build-gmsh3.sh`。
- 网格生成与转换脚本已就绪（未成功运行）：集群 `scripts/ogrid.geo`（参数化 4 块 O 型
  transfinite + 展向拉伸，`4*nc*nn*ns` 个六面体）与 `scripts/of-fixture-gmsh.sh`
  （逐级 `gmsh -3` → 用生产 `gridObj=3` 的 gmsh→oneflow 转换出 `.ofl`）。
- 阶梯尺寸（每卡恒定 294,912 cell）：N 64×64×72、2N ×2、4N ×4、8N 128×128×144、
  16N ×16、32N ×32。

**恢复条件**：gmsh 编译通过并装成 module。届时 2N/4N/16N/32N 由 gmsh 直接生成，
不再走自研的 `refine_ofl.py --cuts`（该模式在 M6 O 型极线的 72 对折叠单元处不可用）
或 `agglomerate_ofl.py`（反向合并，晶格映射已通，bc 继承在折叠面仍失败）。

#### 5.6 交付 2026-11-14：网格加密完成，M6 DCU 弱扩展已测

**基线**：`dev @ 727c1ed75` 部署树（本轮无 C++ 改动，未重建）。工具与作业脚本：
`ci/kunshan/refine_ofl.py`（仓库）+ 集群 `scripts/of-fixture-{refine,validate,partition}.sh`、
`of-probe-weakscale.sh`。证据（集群 `evidence/`，jobid 链）：

| 步骤 | jobid | 结果 |
|---|---|---|
| 加密 + 自检 | `fixture/128020525` | `REFINE_ALL_PASS`：V=2,405,665 / F=7,123,968 / C=2,359,296 / nB=92,160，体积守恒 rel 5.0e-15，重读全结构校验过 |
| 分区 + manifest | `fixture/128020818` | METIS（type 1）8/32/64 zone，`MULTIZONE_FIXTURE_PASS`×3；64 zone：cells/rank 36,458–37,076（±0.9%），nIFaces/rank 4,211–8,330（均 6,292） |
| 单 zone 数值 | `fixture/128021352` | 原/加密网格 XB 残差 1.16e-17 / 3.07e-17；加密 vs 原 CPU 气动力一致 0.2%；aero 差 8e-3 为单 zone HIP 既有缺陷（见坑 41），非网格问题 |
| 弱扩展 + 64-zone XB | `probe/128021420` | XB（64 rank，残差+aero+flow）全绿后测时序，18 次运行同作业交错 |

**M6 DCU 弱扩展结论（36,864 cell/rank，8 ranks/卡，同作业交错中位，1000 步 wall/1000）**：

| 口径 | A：原网格（294,912 cell） | B：加密网格（2,359,296 cell） | 弱扩展代价 |
|---|---|---|---|
| 36,864 cell/rank | 8 rank / 1 卡：**17.0 ms/步** | 64 rank / 8 卡：**22.2 ms/步** | **+31%（总规模 ×8）** |

分项归因（rank0）：`interface_mpi` 2.12 → 7.70 ms/步（3.63×，+5.6 ms），
恰覆盖整步 +5.2 ms；计算段（rhs_inviscid/gradient/macro_update_residuals）反而下降。
**弱扩展代价全部在接口交换**，与 §5.5 强扩展定性一致（跨节点延迟 + METIS halo +43%）。

**加密网格强扩展，逐倍阶梯**（`probe/128040681`：固定 2,359,296 cell，8 ranks/卡恒定，
同作业 3 发交错取中位，1000 步）：

| 卡 | rank | cell/rank | ms/步（三发） | 中位 | 加速比 | 效率 | kcell/ms/卡 | nIFaces/rank |
|---:|---:|---:|---|---:|---:|---:|---:|---:|
| 1 | 8 | 294,912 | 90.03 / 90.22 / 99.01 | **90.22** | 1.00× | 100% | 26.2 | 18.0–23.8k |
| 2 | 16 | 147,456 | 49.29 / 50.01 / 56.33 | **50.01** | 1.80× | 90% | 47.2 | 10.5–15.3k |
| 4 | 32 | 73,728 | 32.61 / 36.34 / 39.15 | **36.34** | 2.48× | 62% | 64.9 | 5.9–12.0k |
| 8 | 64 | 36,864 | 22.83 / 25.36 / 26.42 | **25.36** | 3.56× | 44% | 93.0 | 4.2–8.3k |

每发之间的跨度 10–18%（同一作业同一批节点），所以效率只能读作量级：
**2 卡 90% → 4 卡 62% → 8 卡 44%**。4 卡（32 rank）仍是单节点满配，
8 卡首次跨节点；halo 随 rank 增长反而下降，所以下降是节点内争用 + 跨节点延迟，
与 §5.5 的分项归因一致。

另一作业（`probe/128021420`）同配置点给出 8/32/64 rank = 105.45/33.07/22.21 ms/步，
与本作业相差 15–18%，落在坑 15 的跨作业噪声范围内；**只有同作业数据可比**。

73,728/rank（4 rank/1 卡 vs 32 rank/4 卡）与 294,912/rank（1 rank 独占卡 vs 8 rank/卡）
两个系列因 A/B 两侧卡内 rank 密度不等，不能作弱扩展结论，只作参考
（17.85→33.07、22.97→105.45 ms/步）。

**下一刀候选（按归因排序）**：
0. **弱扩展阶梯的口径（已纠正）**：三维下卡数翻倍对应网格每个方向都翻倍，即 cell 数
   **×2^3 = ×8**，所以 2 卡/4 卡不存在等 cell/rank 的弱扩展点。正确的阶梯固定
   36,864 cell/rank，三级为 **1 卡 N(294,912)、8 卡 8N(2,359,296)、64 卡 64N(18,874,368)**。
   64N 已由 `fixtures/m6-refined-r2/`（对 8N 再做一次同构加密，10.6 min、峰值内存 8.6 GB）
   生成：V=19,059,265 / F=56,807,424 / C=18,874,368 / nB=368,640，`REFINE_ALL_PASS`。
   512 zone 已完成（`fixture/128042784`：METIS 4 min 20 s，峰值内存 10.7 GB，2.22 GB 网格）。
   **64 卡这一级已解锁（2026-11-15）**：原 `probe/128043090` 在 512 rank / 16 节点
   （i02 系列）触发 `ReconstructFaceEulerKernel` VMFault。复测（`probe/128103456` 3 步
   DEFER_SYNC on/off 判别 3 发 + `probe/128104413` 1000 步 ×3 发）在 h03/b13 系列节点
   全部通过、零复现，包括默认环境与 `ONEFLOW_HIP_DEFER_SYNC=0`。kernel 参数转储逐字节
   解码自洽（nFaces/totalCells 与各 buffer 字节数一致），故障地址在全部参数缓冲区之外。
   结论：**不是代码缺陷，判定为 i02 系节点状态问题**（先例：kshdnormal e12r1n04 无法
   打开 `/dev/kfd`；坑 37 同节点租户干扰）。注意 512 zone 网格只能用 512 rank
   （zone 数必须等于 rank 数）；64 rank 时每 rank 持有 8 zone，走的是坑 18 的
   多 zone/rank 路径（`probe/128043428` 已证实该配置无效）。64N 档实测：
   56.9 ms/步（中位，61.7/51.4/56.9，36,864 cell/rank 恒定口径）。
   `--cuts K`（单方向切割）保留为工具能力与诊断入口，**不是本阶梯所需**：M6 的 O 型极线有
   **72 对单元共享两个面的折叠**（`fixture/128040158` 实测，均落在 |v|≈1e-8 的薄片里，
   加密后继承为 144 对），退化单元的方向标签只能事后补全，使跨面的“切向一致”在折叠处
   失效（最后失败点：`child 1 has no canonical face`）。若以后要非整数倍增长再用。
1. **`partition_type = 2`（结构化等块）重切 64 zone**：METIS 的 nIFaces/rank +43% 直接
   抬高 interface_mpi；D1a 在 32 rank 曾收 −4.4%，64 rank 跨 2 节点时 halo 占比更大；
2. **A4**：设备驻留宏步循环，减少每步 host 往返与交换轮次（既有下一主线）；
3. **大网格用更粗的 rank/卡 packing**：22.4 vs 13.3 kcell/ms/卡说明 8 ranks/卡对大
   per-rank 网格不是最优（本口径结论，非生产建议）。

**本轮新坑**：41（单 zone HIP aero 读陈旧 host 梯度）、`--cpus`/`--memory` 为每节点值
（见 handoff §5）。

### 6. 关键开关与脚本

| 用途 | 位置 |
|---|---|
| 后端选择 | `ONEFLOW_ACCEL_BACKEND=CPU|HIP` |
| 设备权威 q（默认开） | `ONEFLOW_HIP_DEVICE_AUTHORITATIVE_Q=0` 回退旧行为；`ONEFLOW_HIP_GHOST_Q_MIRROR=1` 只把 ghost q 的 D2H 加回（bisect） |
| 其它性能开关（默认全开） | `ONEFLOW_HIP_SKIP_HOST_INTERFACE`、`ONEFLOW_HIP_GRADIENT_MIRROR_ON_DEMAND`（默认只在最后 RK stage 镜像；`=0` 恢复每 stage）、`ONEFLOW_HIP_SKIP_HOST_RK_SLOT_COPY`（默认跳过设备显式路径上的 `LOAD_Q`/`LOAD_RESIDUALS` 槽位拷贝；`=0` 恢复）、`ONEFLOW_HIP_CACHE_INTERFACE_HANDSHAKE`、`ONEFLOW_HIP_SKIP_HOST_FIELD_PACK`、`ONEFLOW_HIP_SKIP_UNUSED_INTERFACE`、`ONEFLOW_HIP_DEVICE_GRADIENT_UNPACK` |
| 本轮新增（默认全开，`=0` 回退） | `ONEFLOW_HIP_DEFER_SYNC`（drain 只在需要回读处）、`ONEFLOW_HIP_DEFER_STATUS`（非物理标志留设备，宏步末检查一次）、`ONEFLOW_HIP_FUSE_SCALE_UPDATE`（残差缩放并入状态更新入口） |
| 队列深度 | `ONEFLOW_HIP_MAX_LAUNCHES_PER_DRAIN`（默认 32；多 rank 下不设上限会在 `libgalaxyhip` 里 AQL 溢出 abort） |
| 诊断 | `ONEFLOW_STAGE_BREAKDOWN_FILE`、`ONEFLOW_UNS_STAGE_TRACE_FILE`、`ONEFLOW_UNS_HIP_RESIDENT_CONTRACT=1`、`ONEFLOW_UNS_HIP_STAGE_DIAGNOSTIC=1` |
| 多 zone fixture 制作 | runbook §1b：`script/grid.txt`（并列入 `control.txt`）→ `simutask=Grid` → 补 `grid/bcRegionMap.txt` → 导出/merge/verify manifest |
| 门禁脚本（仓库内） | `ci/kunshan/f2-main-solver-{multizone-trace,cross-backend-residual,stability}.sh`、`f2-main-solver-paired-performance.slurm` |
| 标准作业脚本（集群固定工作区） | `scripts/{env,rank,build-app,prep-case,gate,paired-ab,cpu-regression,cleanup}`，规则见 [`runbooks/oneflow-kunshan-workspace-runbook.md`](../runbooks/oneflow-kunshan-workspace-runbook.md) |
| 加速器构建（仓库内，可选） | `scripts/build-accel.sh`；参数只有 `ONEFLOW_ACCEL_SDK_ROOT`（+ 可选 `_ARCH`/`_FLAGS`），种类自动推断；不设则 CPU/CI 不受影响 |
| 分支同步 | `scripts/sync-master.sh`（macOS/Linux 双兼容，支持 `--dry-run`） |
| 基准口径 | 1000 步、50 步 warmup、`nFieldSave/nResSave/nForceSave=100000`、CFL 0.01、`ilim=0`、`ischeme=5`、`time_integral=1`、`rk_stage=3`、`vismodel=0` |

### 7. 口径纪律（不遵守的结论一律作废）

1. 只比同节点、同 fixture、同 `steps/warmup/repeats`、同 rank/资源；`lifecycle_*_ms` 是 repeats 总和。
2. 不以 100 步或 kernel-only 作生产结论；1 DCU 槽 = 8 核。
3. 探索性单作业数据只进本文档或集群 artifacts，不进正式报告。
4. 任何性能改动先过：3 步多区 trace → 200 步跨后端 → CPU 五算例（normal/strict）→ 50 步稳定性。
5. 正式报告的 Markdown 与 HTML 必须同 stem 同步更新。
6. **数值参考数据的更新纪律**：先跑 strict 量清越界范围 → **只重生成越界算例**的
   `autotest/*.dat` → 用 `build_residual_db.py --absolute-tolerance 1e-15` 重建并
   `verify_residual_db.py` 校验 → 用**不变的** strict 判据复跑 5/5。判据本身不得为了
   通过而放宽；改成"精度优先"时只更新数据，不动判据。

---

## 归档说明

本文件在 2026-11-08 精简为 **living 大纲**：只保留当前状态（§1）、全链路审计（§3）、
Phase A–F 台账与证据（§4）、下一轮任务（§5）、开关与口径（§6–§7）。

**逐日的轮次记录、一次性实验数据与当时的判断已从本文件移除**，按 `doc/README.md`
的保留策略以 git 历史为档案（`git log -p doc/plans/oneflow-development-todo.md`）。
需要复用某轮结论时，从 §4 对应条目或 `handoff/oneflow-current-handoff.md` §5 的坑列表
进入；两者只保留**当前仍成立**的结论。

#### 5.8 合并可行性已判否（2026-11-15，可执行检查器 + 反例）

`ci/kunshan/check_merge_topo.py`（不写网格，只做拓扑判定）在本地生成的 8N 网格上实测：

- **provenance 已验证**：父单元内部 **3,538,944/3,538,944** 张面满足
  `child_cell = parent*8 + k` 且两侧 octant 坐标恰好差 1 个分量。
- **跨父单元不兼容**：3,493,404 张面中 **33,374 张（0.96%）** octant 差为 0/2/3。
  最小反例：`face=92175 fine=(5,9219) parents=(0,1152) octants=(5,3) diff=3`
  （5=(1,0,1)，3=(0,1,0)，三个分量全不同——两侧父单元的局部坐标轴彼此旋转）。
- 按 `cell//36864` 猜 8 个块不能隔离这些冲突（同块内仍有 25,635 张不一致），
  所以不是单纯的块界面旋转。

**结论已撤回（2026-11-15，同一轮内更正）**：跨父单元的 raw octant 差值**不是**一般成立的
必要条件——相邻两父单元的局部坐标轴可能切向反向或置换，octant 下标只在本单元内有意义。
因此上面的 33,374 张“不兼容面”**不能**作为不可行的证据，结论改为：**原判据不足，
待坐标变换与覆盖约束验证**。已改用基于真实 face ID + 共享面角点对应 + 3×3 候选组合表
的检查器 `ci/kunshan/check_merge_dirs.py`（每张 parent_face 保持独立身份，不按
unordered(left,right) 合并）。该检查器当前会把 99.3% 的面判为“不可解析”，原因是
`hex_corners` 给出的底面环起点是**循环歧义**，corner→晶格坐标的映射必须从共享面角点
对应中恢复，而不能按固定表假设——这是下一处要修的点。在此之前对 2N 可行性不下结论。参数语义同时纠正为
`merge_factor`（brick 内含细单元数）与 `bricks_per_parent`：
`dims(2,2,1)` → 每 parent 2 个 brick、每 brick 4 个细单元 → 输出 `C/4`；
`dims(2,1,1)` → 4 个 brick、每 brick 2 个细单元 → 输出 `C/2`。

**检查器适用边界**：只判定「同一 parent 内部 / 相邻 parent 之间 octant 坐标是否
单分量相邻」，不验证粗面几何、不写网格、不检查 BC 继承；它否定的是局部方向合并，
不否定基于真实 face ID 的其它合并策略。

#### 5.12 阶段一收口：运行配置建议（2026-11-15）

**推荐运行方式**（按卡预算选，数据见 §5.11）：

| 优先级 | 场景 | 配置 | 实测 |
|---|---|---|---|
| 默认 | 现有口径、**单卡** | **8 rank/卡**，`DEFER_SYNC=1`（不变） | 12.90 ms/步，每卡 22.9 kcell/ms |
| 要时间 | 卡预算宽裕、单节点 | 4 rank/卡（即 2 卡） | 8.64 ms/步（−33%），每卡 17.1 |
| 要极限 | 卡预算很宽 | 2 rank/卡（即 4 卡） | 6.69 ms/步（−48%），每卡 11.0 |

注意：这里的"4/2 rank/卡"是**在 8 rank 总量下多买了卡**，不是"每卡少放 rank 更好"。
若卡数固定为 1，历史 §4.8b 已测 8 rank/卡 最优（1/2/4/8 rank = 14.36/10.21/8.82/7.43）。
两个问题不要混谈。

- `DEFER_SYNC=0` **不推荐**：三个密度下一致更慢（3–7%，远超 1.3% 离散度）。
- 放宽 rank/卡是**用资源换时间**，不是效率提升：每卡吞吐 22.9 → 17.1 → 11.0
  kcell/ms/卡（单调下降），选它只在 wall time 优先时成立。
- 跨节点基线：64 rank / 8 卡 / 8N = 20.60 ms/步（中位）。
  **注意：§5.11 里"62% 惩罚与通信无关"的分解不可引用** —— 它比较的是
  「1 卡 + N 网格（294,912 cell）」与「8 卡 + 8N 网格（2,359,296 cell）」，
  同时差了三件事（节点数、网格规模、启动成本，8N 网格文件 270 MB 启动远超 N 的
  ~1.5 s），用 `wall/1000` 直接相减，不是受控分解。跨节点瓶颈目前**未定性**。

**阶段二判断**：单节点通信份额 11.6%（1.49 ms/步），是交换层改动的**绝对上限**；
现实估计减半即 0.7 ms/步。对标：配置层已有 −33% / −48%。跨节点份额 21.4%
同样不可引用（同上混剔除）。**先交付配置，代码改动只在固定卡预算时才有意义**。

**数值验证状态**：本轮所有数字都是**配置变化**，未改代码、未改方程/CFL/容差/
halo lag；求解器源码与部署树 `95c991cf3` 的 accel/ns/parallel/task/uns/interface/
grad/global **逐文件一致**（`git diff` 为空）。因此沿用该修订已有门禁：
trace `24/0/0`、200 步 XB `2.331035e-18`/aero `0.0`、CPU normal/strict 5/5。
**`DEFER_SYNC` 与 rank/卡的变化不改变数值路径**，但本轮未重跑门禁，
不主张"新配置已过门禁"。

#### 5.13 负结果 2026-11-15：每 rank 块大小**不是** M6 DCU 的主因

**假设**（被本节否定）：36,864 cell/rank 落在 todo 记的"≈1 万开始变差、≳5 万正常"
的边界下沿，而立方体在 110,592 cell/rank 下弱扩展近乎持平（13.78→13.68，2→32 卡）。
若成立，整条 M6 阶梯曲线测在最不利粒度上，运行建议应改为"每 rank ≥1.5e5 cell"。

**实验**（`probe/128112206`）：固定总工作量 8N=2,359,296 cell、固定 4 卡、单节点，
只改切法；三个分区全部由单 zone 8N 网格以 `partition_type=2`（结构化等块，
cell 严格等分）现生成并用 `src_elements` 与块排列行双重校验
（`probe/128112156`：2x8x2 / 2x8x1 / 1x8x1）。每配置 3 发交错，1000 步墙钟。

| 配置 | rank | rank/卡 | cell/rank | 1000 步墙钟（中位） | 离散度 | 相对 F32 | 每卡吞吐 kcell/ms/卡 |
|---|---:|---:|---:|---:|---:|---:|---:|
| F32 | 32 | 8 | 73,728 | **71,800 ms** | 1.019 | 1.000× | 8.22 |
| F16 | 16 | 4 | 147,456 | **62,973 ms** | 1.016 | **0.877×** | 9.37 |
| F8 | 8 | 2 | 294,912 | **70,361 ms** | 1.004 | 0.980× | 8.38 |

**结论：假设否定。** 从 73,728 放大到 294,912 cell/rank（4 倍）墙钟只变 −2%，
在离散度内；中间点 147,456 快 12%，是**浅最优**而非区间切换。
即每 rank 块大小在 7e4–3e5 范围内只有 ~12% 的影响，**不足以解释阶梯退化**。
历史弱扩展曲线的解读**不需要**因粒度改写。

**顺带确立的方法论**：8N 网格（270 MB）的启动成本让小步数斜率法失效——
40 步与 1000 步的墙钟比在所有配置下都近似 4.5:1，说明 40 步几乎全被启动占据；
第一发（冷）与后续发差 20%+。**该规模只能用同分配内交错、多发的 1000 步墙钟比较**，
不能用短步数扣启动。

#### 5.14 交付 2026-11-15：单节点生产配置的完整每步账本

**配置**：生产默认（N=294,912 cell / 8 zone / 8 rank / 1 卡 / `DEFER_SYNC=1`），
1000 步，3 发交错；二进制 sha256 `33c62105…`，部署树 `95c991cf3`（含 H1 全套
host 记账）。证据 `probe/128113967`：`.host.tsv`（incl/excl）、`break-%r.tsv`、
`mpi trace`（收尾段脚本 bug 未汇总，核心账不依赖它）。

**完整性门（先过再读数）**：闭合恒等式 `mg_run.incl = Σ scope.excl` 三发全
CLOSED（11695/11878/11693 ms）；`mg_inner_iteration` count=1000 确认步数生效；
host 账本与 stage 分项在 4 个共同段上差 <0.5%（互证）；同 rank 跨发离散 <15%。

**12.90 ms/step 的去向**（rank0，启动扣减后 mg_run=11.69）：

| 段 | ms/步 | 占比 | 性质 |
|---|---:|---:|---|
| `macro_update_residuals`（3 个 RK stage 的 host 包装：gradient/recon/flux 入队+等待） | 5.29 | 41% | host |
| `task_execute` 自有池（= 每步 post 轮的 gradient pack/upload/download 3.12 + ~0.2 杂项） | 3.33 | 26% | host |
| 设备 kernel 纯执行（rocprof 既有证据） | ~1.77 | 14% | device |
| `interface_mpi`（recv 主导的 MPI 等待） | 1.37 | 11% | 通信 |
| `interface_pack`（每 stage 设备打包+D2H） | 1.08 | 8% | 通信 |
| `interface_unpack` + time step/LHS/boundary/flowfield | 0.55 | 4% | host |
| 派发/计划/任务开销（7000 次/步） | 0.04 | 0.3% | host |
| mg_run 之外（启动 + mg_run 外段） | ~1.2 | 9% | — |

（`interface_*` 是 3 个 stage 界面派发 + post 批次的子 scope；`post_upload` 2.46 +
`post_download` 0.66 与 task_execute 自有池 3.33 吻合——池的主体就是 post 轮
gradient 打包上传下载。`swap_upload/swap_download` 是嵌套计时，不可与 5.29 相加。）

**结构性发现：等待与包装是同一现象的两面。** rank0（MPI 等待 1.4–1.6，stage 包装
5.2–5.6）与 rank3（MPI 等待 4.6–6.2，stage 包装 2.9–3.7）三桶之和几乎相同
（10.0–10.5 vs 11.0–11.4）：先到交换点的 rank 把等待记在 `MPI_Recv`，后到的把
同一个等待记在 stage 包装里（等设备队列/下一次入队）。**因此"interface_mpi 只占
11%"不成立为"通信便宜"的论据**——真实 rendezvous 等待总量约等于
`interface_mpi + interface_pack + 两 rank 包装差`，在关键路径上约为 4–5 ms/步。

**对下一步的含义（A4 预算的依据）**：
- 每步可归给 host 编排/同步的池：41%（5.29）+ 26%（3.33）中可压缩的部分。
  A4（设备驻留宏步已含 RK 循环，现状 `resident_macrostep` excl 仅 0.031）要消的
  是 5.29 里的入队-等待与 3.33 里的 post 轮 host 打包，不是任务派发（仅 0.04）。
- 明确**不是** A4 目标的：设备 kernel 1.77、mg_run 外的启动、task_plan/dispatch。
- 通信轮次合并/重叠的收益上界：~2.5 ms/步（1.37 等待 + 1.08 pack），与此前
  §5.11 的 1.49 上限同量级。

#### 5.15 新发现 2026-11-15：rank 与卡的 NUMA 错配（未验证收益）

**实测拓扑**（`probe/128116197`，kshdnormal 节点）：

```
NUMA node0: cpu 0-7    node1: 8-15    node2: 16-23    node3: 24-31
card0 0000:23:00.0 numa=1   card1 0000:04:00.0 numa=0   card2 0000:26:00.0 numa=1
card3 0000:43:00.0 numa=2   card4 0000:63:00.0 numa=3
```

`--bind-to core` 把 8 个 rank 全部绑到 **CPU 8-15（NUMA node 1）**，与 Slurm 分到哪张卡
无关；本次分到 `SLURM_JOB_GPUS=1` → card1 → **NUMA node 0**。
即 **8 个 rank 的 CPU 在一个 NUMA node，显卡在另一个**，每次 H2D/D2H 跨 NUMA。

**为什么值得追**：账本（§5.14）显示每步有 3 处阻塞 D2H（`post_upload_host`
1.2–2.3 ms、`interface_pack` 1.1–1.4、`s_gradient` 内的 unpack 等待），
跨 NUMA 的拷贝延迟正是这类阻塞成本的候选解释之一。

**但当前无法做对照**：Slurm 的 cgroup 只授予核 8-15，`--cpu-set 0-7` 直接报
"not enough slots"，**rank 无法挪到卡所在的 NUMA node**。`--mem-bind local` 只约束
内存分配、不改变 CPU 位置。`probe/128116411` 因此只跑完 baseline 一发即中止
（且该节点 baseline 17.05 ms/步，与他节点 12.9 差 32%，跨节点不可比）。

**要验证需要**：申请更多 CPU（多点 NUMA node）或 `--gres=dcu:4` 拿全卡，
再按卡所在 NUMA 分布 rank。这是一个独立实验设计，不是调 mpirun 参数。

**同时确认的负面结果**：`probe/128116476` 的 4/8 rank GPU-aware FAIL 是探针程序
本身要求恰好 2 rank，**不是 MPI 失败**；2 rank 在当前构建上 pass=1。

#### 5.16 归因框架（2026-11-15，源码核实）

后续所有"加卡/加节点为什么变慢"的分析，都应基于以下已核实的代码事实，而不是猜测：

**生产默认下，跨 rank 交换在 1 卡、多卡、多节点走的是同一条 host staging 路径**
（`codes/task/src/InterfaceTask.cpp:306-324`）：
`useDeviceDirect` 需要 `ONEFLOW_ENABLE_UNS_GPU_AWARE_MPI=1` **且**
`ONEFLOW_GPU_AWARE_MPI_CONFIRMED=1`，两者默认都不设；因此每轮交换恒为
设备 pack kernel → 阻塞 `CopyToHost` → host 缓冲上的阻塞 `MPI_Sendrecv` →
阻塞 `CopyFromHost` → unpack kernel。

**推论**：卡数/节点数增加**不改变代码路径**，只改变底下 MPI 的传输距离
（同卡 SHM → 跨卡 PCIe → 跨节点 IB）。所以归因时应表述为"同一路径下传输距离变长"，
除非数据明确显示别的机制。

**每步 4 轮交换**（3 个 RK stage + 1 个 post 轮），已核实于
`NsResidentMacroStep.cpp:117`、`NsSolverImp.cpp:382`（`NsPostprocess` →`CommInterfaceData`）。
每轮每 rank 在**每个相邻 zone 对**上阻塞一次，所以邻居数为 3 的 rank 每步付 6 次阻塞调用、
邻居数为 2 的付 4 次——这解释了此前测到的"消息多反而单条延迟低"
（固定每轮成本被并行 pair 数摊薄）。

**归因时的判别口径**：
- `interface_mpi` 增长 → 传输延迟
- `interface_pack` / `post_upload_host` 增长 → 设备侧拷贝异常（**本应节点内不变**，
  若随节点数增长要单独查，说明有非预期的跨节点数据流）
- `cells/rank` 下降 → 粒度效应，与传输无关

#### 5.17 dmabuf 反编译结论（2026-11-15，REA）

**推翻了此前"海光把 dmabuf 功能阉割"的判断。** 反编译 `libamdhip64.so`（DTK 26.04）
与 `libhsa-runtime64.so.1.5.0` 的事实：

| 事实 | 证据 |
|---|---|
| `hipMemGetHandleForAddressRange` 的 dmabuf 分支**完整实现** | 0x54f840，type==1 分支在 0x5503fd，参数正确传递设备指针+size+fd 输出 |
| 它调 `dlsym("hsa_amd_portable_export_dmabuf")` | 0x83ee31 处加载符号到 .bss 槽 0xdd9aa0 |
| **失败在 dlsym 返回 NULL** | 该 HSA 库 226 个导出函数、48 个 `hsa_amd_*`，三个 dmabuf 函数全缺 |
| 空指针 → 返回 801 `hipErrorNotSupported` | 0x550a42，正是实测到的 "operation not supported" |
| 内核/驱动侧 `AMDKFD_IOC_EXPORT_DMABUF` 常量在用户态 thunk 中不存在 | 扫描全部 `_IOWR('K',nr,size)`，行 0x00–0x3f 均无 |
| **导入方向可用**：`hsa_amd_interop_map_buffer` 已导出 | DRM fd → 设备指针这条路今天就能走 |
| HIP 支持 `LD_PRELOAD` 补符号 | 其 dlsym 走全局作用域，预加载的自定义 .so 可满足 |

**所以不是"功能缺失"，是"符号缺失"** —— 一个是死路，一个可补。

**两条候选路径**：
- **(a) LD_PRELOAD shim**：自己实现 `hsa_amd_portable_export_dmabuf`。阻塞点：
  内核是否接受 `AMDKFD_IOC_EXPORT_DMABUF`（用户态 thunk 没调用，不代表内核不支持）。
  正在用一个 ioctl 探针实测判定。
- **(b) 导入方向**：用 DRM/PRIME 自行分配 + `hsa_amd_interop_map_buffer` 映射进 HIP。
  无需 hook，今天可用，但要求数据源由我方分配。

（此前 §"dmabuf 缺失"相关表述应以本节为准。）
