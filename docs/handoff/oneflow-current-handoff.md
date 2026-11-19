# OneFLOW 交接（living）

> 2026-11-19 核对。工作分支 `dev`，当前接手项为 M6 粘性路径 B.3。
> 历史测量和旧章节从 Living TODO 的证据索引进入。

> **这是唯一的交接文档，就地更新，不再按日期新建。**
> 当前切片与验收勾选在 [`../plans/oneflow-development-todo.md`](../plans/oneflow-development-todo.md)；
> 长期阶段的取舍与里程碑在 [`../plans/oneflow-development-plan.md`](../plans/oneflow-development-plan.md)；
> 规则以仓库根 `AGENTS.md` 为准；集群操作流程见
> [`../runbooks/oneflow-kunshan-workspace-runbook.md`](../runbooks/oneflow-kunshan-workspace-runbook.md)。
>
> 本文只放"接手时立刻需要"的东西：**环境与固定路径、启动顺序、常量和坑、当前任务指针**。
> 当前验收以 Living TODO 为准，完整路线和任务编号以 Phase Plan 为准。

## 0. 60 秒上手

1. 读 `AGENTS.md`（分支模型、门禁顺序、口径纪律、工作区规则）。
2. 打开 **Living TODO 当前切片**。先验证 M6 粘性通量写回和气动力；
   H2 已在两节点两卡配置复验，A4 排在本轮正确性验收之后。
   **主线优先级只以 Living TODO 当前切片 为准**；本文只记“怎么跑”和“已经踩过的坑”，
   两边不一致时改本文，不要改 todo 的台账。
3. 要动集群先读 §1（固定路径）与 runbook：**SSH 已不可用，只能走 OpenAPI 技能**，
   因此一切远端操作必须写成作业脚本提交，不能在交互 shell 里敲。
4. 收尾：勾选 todo、证据落 `evidence/<gate>/<jobid>/`、更新本文 §4 与 §5；
   环境或路径变化同步 runbook 与 本文 §1。

## 1. 环境与固定路径（不要另建目录）

| 位置 | 路径 | 用途 | 规则 |
|---|---|---|---|
| 本地仓库 | 本 checkout 根（`git rev-parse --show-toplevel`） | 源码、文档、提交 | 只在 `dev` 上提交；`master` 只读、与 upstream 同步 |
| 本地临时 | `/tmp` | 一次性脚本、打包、比较 | 不入库；需要保留的结论写进文档 |
| 集群接入 | scnet-hpc 技能（按当前安装位置解析入口） | 查询、提交、上传、探测 | 按当前记录只用 OpenAPI；不要手搓 ssh/scp 循环，也不要静默切换后端 |
| 昆山工作区 | **`$HOME/oneflow/`（唯一）** | `src/ build/ fixtures/ evidence/ snapshots/ scripts/ logs/ _attic/` | 固定布局，见下表与 runbook |
| 门禁证据 | `$HOME/oneflow/evidence/<category>/<jobid>/` | **6 个固定 category**：`gate` / `paired-ab` / `cpu-regression` / `probe` / `fixture` / `build` | 只增不改；不要按实验临时起目录名（规则见 runbook §1c）；`*.stage` 等大二进制提取指标后可删 |
| 源码快照 | `$HOME/oneflow/src/OneFLOW` + `src/REVISION` | 每轮覆盖 | 不建 `OneFLOW-<hash>` 目录 |
| 构建 | `$HOME/oneflow/build/app` | 唯一增量构建 | 日常不重新 configure；必要时按 runbook 的 SDK 接口在 DCU 分配内处理 |
| fixture | `$HOME/oneflow/fixtures/m6-8zone-cpu8/` | 唯一**登记**的分区算例（8 zone / 294,912 cell） | 只读；派生副本用 `prep-case.sh`。其它 zone 数用 runbook §1b 的配方现生成 |
| 快照归档 | `$HOME/oneflow/snapshots/*.tgz` | 保留最近 5 个 | `cleanup.sh` 轮转 |
| 日志 | `$HOME/oneflow/logs/` | 作业 stdout/stderr | 保留 14 天 |
| 弃用 | `$HOME/oneflow/_attic/<date>/` | 待删物 | 观察一轮后清（`cleanup.sh`，7 天） |

**集群接入现状（2026-11-08）**：`ssh kseshell` 返回 `Permission denied (publickey)`（密钥
失效或 IP 白名单变更），**只能用 OpenAPI**。由此带来的硬约束：

- 一切远端动作都要写成**作业脚本**用 `submit --command` 提交（建目录、解包、构建、清理都是）；
- OpenAPI 没有 `exec`，所以不能交互式 `make`，也不能 rsync；源码以
  `git archive` 打 tgz → `upload` 到 `snapshots/` → 作业内解包到 `src/OneFLOW`；
- 提交参数有内存/CPU 比例上限（`--memory` 过大直接被 sbatch 拒绝），4 核请用 ≤8G；
- 偶发 `Socket timed out` / `HTTP 502`：隔十几秒原样重试即可，不要当成作业失败。

**禁止清单**（写在这里就是为了不再重复）：

- 在 `$HOME` 根目录或 `softwares/projects/` 下新建 OneFLOW 目录、日志、一次性脚本；
- 按 job 号或轮次新建 build/run 目录（`build/<jobid>`、`runs/<name>-<jobid>`）；
- 在仓库里提交绝对私有路径、主机名、账号、作业号、原始 CI/Slurm 日志；
- 为一次会话新建 dated handoff 或多份并行计划（本文就地更新）；
- 直接跑 `scnet`/`ssh` 绕开技能去改远端状态；
- 按实验内容临时给证据目录或作业起名（`a17`/`ab2`/`attr` 这类）。作业名固定
  `of-<category>[-<purpose>]`，证据目录固定 `evidence/<category>/<jobid>/`，规则见 runbook §1c。

## 2. 昆山环境常量

| 项 | 值 | 原因 |
|---|---|---|
| 启动 | **只能 `mpirun --bind-to core --map-by slot`** | `srun` 下 MPI size=1，每个 rank 会解全部 zone |
| 环境 | 常规作业 source `scripts/env.sh` | MPI/UCX/DTK 必须显式进 PATH；rocprof 另需坑 31 的已验证 module 初始化 |
| HCA | `UCX_NET_DEVICES=mlx5_0:1`、`UCX_IB_DEVICES=mlx5_0` | 其它 HCA 的 `ibv_create_ah` 返回 EINVAL |
| MPI coll | `OMPI_MCA_coll=^hcoll,han` | 非均匀分布会崩在 `coll/han` |
| CMake | `/public/home/luql/apprepo/cmake/4.4.3/bin/cmake`（`env.sh` 已把它加入 `PATH`） | 旧的 `$HOME/tools/cmake-3.30.5` 已随家目录清理删除；模块里的 cmake 是死链 |
| DTK | `dtk-26.04`（显式加 `bin/llvm/bin/hip/lib` 等） | 不显式设会落到系统旧 ROCm |
| 资源口径 | 1 DCU 槽 = 8 核；CPU 基线在 `kshcnormal` | 跨队列比较必须换算到同一口径 |
| 基准 | 1000 步、50 步 warmup、`nFieldSave/nResSave/nForceSave=100000`、CFL 0.01 | 100 步或 kernel-only 不能作生产结论 |

## 3. 常用命令

**注意：SSH 不可用（见 §1），下面的 `sbatch` 只是作业脚本名；实际提交走 OpenAPI 技能。**

```bash
# 仓库侧
git -C <repo> status --short          # 提交前必须干净
# 集群侧脚本（都在固定工作区，由作业内执行）
source $HOME/oneflow/scripts/env.sh
$HOME/oneflow/scripts/build-app.sh                    # 增量编译（两步 make，见坑 12）
$HOME/oneflow/scripts/gate.slurm                      # 3 步 trace + 200 步 XB + 500 步分项 + 3x1000 步
$HOME/oneflow/scripts/paired-ab.slurm                 # 同二进制 paired A/B
$HOME/oneflow/scripts/cpu-regression.slurm            # CPU 五算例 normal/strict
$HOME/oneflow/scripts/cleanup.sh                      # 日志/快照/attic 轮转

# 提交（唯一方式，OpenAPI）
SKILL=$HOME/.agents/skills/scnet-hpc
python3 $SKILL/scripts/scnet.py --region 11250 submit --name of-gate \
  --command "bash /public/home/<user>/oneflow/scripts/gate.slurm" \
  --work-dir <ws>/logs --queue kshdnormal \
  --nodes 1 --cpus 8 --dcus 1 --memory 27G --walltime 01:30:00 \
  --stdout <ws>/logs/of-gate-%j.out --stderr <ws>/logs/of-gate-%j.err
python3 $SKILL/scripts/scnet.py --region 11250 --details job <jobid>
python3 $SKILL/scripts/scnet.py --region 11250 logs --path <abs path> --lines 200
python3 $SKILL/scripts/scnet.py --region 11250 --timeout 180 upload <local tgz> <ws>/snapshots
```

门禁顺序固定：3 步多区 trace（含 resident contract）→ 200 步跨后端 → CPU 五算例
normal/strict → 50 步稳定性；性能改动再加 paired。任何一步失败都不算完成。

## 4. 当前任务指针（每次会话就地更新）

**2026-11-19 接手点**。从 `origin/dev` 接续，核对当前 SHA、工作树及集群
`src/REVISION`。三端 master 对齐 `cecd66219`，已合入 dev（`e7efe7d02`）。
本轮主线是 M6 粘性设备路径正确性，阶段与证据状态见 Living TODO 的 B.3 切片。

1. 用固定 Kunshan 工作区构建此次源码快照，核对实际运行的二进制与 hook 命中，
   从单 rank、单卡 M6 全关/仅粘性开始，粘性开关为
   `ONEFLOW_ENABLE_UNS_HIP_VISFLUX`。用 `ONEFLOW_VISCOUS_SELFCHECK` 比较
   求解器写回后的逐面通量，先确认 CPU 与设备输入的时间层和边界状态一致。
2. 记录首次分叉的步、面及方程；检查壁面和内部面，并显式检查非有限值。
   通量对照一致但气动力仍分叉时，对齐同一输出步的壁面 ghost、梯度及剪切。
   host ghost 陈旧目前只是待验证线索。
3. 关掉探针，复跑 400 步对照及收敛后气动力验收。旧结果是总阻力偏差
   22%、摩擦阻力偏差 143%。收齐证据前，粘性设备路径保持未验收状态。

**本轮自检切片（`169d59892`，Linux 本机，HIP 关闭）**。该提交已推送 `origin/dev`（备份；
常规 CI 不监听 dev）。`cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
配置通过；`viscous_flux_selfcheck_test` 10/10、`oneflow_euler_domain_contract_test`
15/15、`vis_flux_policy_test` 24/24 通过；`codes/uns/src/UNsVisFlux.cpp.o`
单独编译通过。`OneFLOW_RuntimeObjects` 全量目标在本机被 `codes/cgns`
的既有问题中断（系统头被包含进 `ONEFLOW` 命名空间，`intptr_t` 未声明），
与本次改动无关，需在集群工具链上复核全量构建。本轮未运行 HIP、集群门禁
或 M6 算例，不得据此改变 B.3 的验收状态。

**本地验证与待办**。契约目录 54/54，CPU 粘性策略 24/24。
`python3.12 -m unittest tools.oneflow_gui.test_mesh_project tools.oneflow_gui.test_mesh_validation tools.oneflow_gui.test_path_utils`
计 50 项（48 通过，2 跳过）。GUI 全套 discovery 有 388 项、1 failure、13 errors、
3 skips。12 个导入错误缺 PySide6；另有一次 capture mock 返回 str 给平台 API 的
decode，以及一项 Windows 路径在 macOS 上的断言不符。这三类问题留待 GUI 验收，
目前没有确认相对合并前的归属。

本地 CPU 构建沿用 Release/MPI，HIP 关闭。启用 METIS 时重新配置需设置
`METIS_HOME_INC` 与 `METIS_HOME_LIB` 分别指向本机 METIS 头文件目录和动态库，
否则先遇到缺头文件，或在 GeometryCore 测试链接时遇到 METIS 未定义符号。
集群门禁仍需按准确源码重跑，不能沿用此前的全开/全关数据判定新版本。

**CPU 构建与 CTest 复验**。按上面的 METIS 环境重新运行
`cmake -S . -B build` 与 `cmake --build build -j 10` 后全量构建通过。
`ctest --test-dir build -j 8 --output-on-failure` 为 1367/1374，失败为
测试 385、386、388–391、394。它们在 macOS 启动时找不到 OneFLOW system 目录。
`HX_GetExeDirectory()` 的非 Windows 分支读取 Linux 的 `/proc/self/exe`，
在本机得到空路径。Linux/Kunshan 上重新跑 CPU 五例 normal/strict 和 HIP 门禁。
GitHub 常规 Linux/Windows push 工作流监听 master/hexin，dev 的推送结果仅代表备份成功。

以下为历史诊断记录。

**2026-11-17 历史诊断：粘性通量写回错位与临时改动**

症状：全 hook 开启，旗舰算例第 50 步与 host 一致，第 100 步偏离 2×，
第 110 步 `rho = -nan` 崩溃。

**历史根因**。初始写回布局错位。kernel 按
`faceFlux[eq * faceFluxStride + face]`（stride = nEqu = 5）写，
而写回按 `visfluxBuffer[equation * ug.nFaces + face]` 读——两个完全不同的布局。
结果是粘性通量被读成几乎全零，**残差里的粘性项几乎消失**。
`395cbbb03` 调整读回索引，仍保留 stride=5，方程写入区间重叠。
`aada33da2` 将 stride 改为面数，端到端效果待复验。

**为什么查了十几轮**：单面 dump（896,256 面）与 host 逐位一致，因为 dump 的是
kernel 局部变量、不经过写回缓冲。我一直盯着「6.7e-3 相对偏差」看，
而真相是量级差了几个数量级。

**修复验证**（粘性 hook 单开 150 步）：修复前第 100 步 3.31e-06（host 1.67e-06）、
第 150 步崩溃；修复后第 150 步 9.96e-07，与 host 1.00e-06 同量级。

**方法教训**（重要）：
1. 残差轨迹不是验收标准，**收敛后的气动力**才是。
2. 验证 kernel 输出 ≠ 验证数据到达消费端，必须查完整链路。
3. 三次踩同一个坑：旧二进制、快照早于代码、`--remote-path` 覆盖共享
   `job_BASE.slurm`。**投作业前先确认被测对象真的是被测对象**：
   `strings $BIN | grep -c <flag>`、确认 `job_BASE.slurm` 里的脚本名。

**仍未解决**：`interfaceMask` 未填（分区算例会错，旗舰单 zone 不受影响）；
SA `Update` 未接（审计结论：不该接，host 只有 2 FLOP/单元）。

**上一轮（2026-11-16）：三条端口全部通过设备契约验证**

本轮的核心发现是：**「契约测试全绿」不等于「移植完成」**。温度场、SA、粘性三条
port 都出现过这个缺口——kernel 实现完整、测试通过，但生产代码里没有调用者。

已实测（`evidence/build/128157757-merge`，单节点 1 DCU，`gfx906` 实机）：

| 目标 | 结果 | 设备 |
|---|---|---|
| `oneflow_device_temperature_contract_test` | 4/4 | `HIP/DCU gfx906` |
| `oneflow_device_sa_contract_test` | 12/12 | `HIP/DCU gfx906` |
| `oneflow_device_viscous_flux_contract_test` | 6/6 | `HIP/DCU gfx906` |
| `turb_sa_arithmetic_test` | 11/11 | CPU（host parity） |

**判别力已验证**（job 128158376）：把 `ViscousFluxArithmetic.h:123` 的
`2.0/3.0` 改成 `2.5/3.0` 后，粘性契约 9 项失败，其余三个目标全绿——
失败精确局限在唯一使用该系数的目标上。

**已接线（生产构建验证通过，未跑算例）**：

**`WIRED_ALL_PASS`**（job 128161160）：`configure rc=0`、**生产运行时 rc=0**、
**生产二进制 rc=0**、契约测试回归 33/33。生产运行时是唯一会编译
`UNsVisFlux.cpp` / `UTurbSrcFlux.cpp` / `UTurbInvFlux.cpp` 的目标，所以这一步
证明了四处接线的接口全部正确。

1. **粘性面通量** `UNsVisFlux::CalcFlux` → `CalcVisFluxHip()`，门控
   `ONEFLOW_ENABLE_UNS_HIP_VISFLUX`。同时修掉一个实质缺陷：原实现每次调用从 host
   重传 26 个 cell 平面（300k 单元约 **61 MB / 4.7–9.4 ms**），而设备上本来就有
   `cellState`/`cellGradient`/`temperature` 和全部几何。现改为命中 `boundStates_`
   就复用，26 个平面降到 2 个（只有两个粘度是每次输入）。
2. **SA 涡粘度** `UTurbSrcFlux::CalcVist1Equ()`，门控
   `ONEFLOW_ENABLE_UTURB_HIP_VIST`。
3. **SA 源项** `UTurbSrcFlux::CalcSrcFlux1Eq()`，门控
   `ONEFLOW_ENABLE_UTURB_HIP_SOURCE`。
4. **SA 迎风面通量** `UTurbInvFlux::CalcInvFlux()`，门控
   `ONEFLOW_ENABLE_UTURB_HIP_INVFLUX`。

**仍缺**：SA 的标量扩散面通量与标量更新两处未接。

**接线过程中抓到三个错误，其中两个是静默的**：

| 错误 | 类型 | 抓到方式 |
|---|---|---|
| `turbcom.fw_star` 应为 `fwStar` | 编译期 | 生产构建 |
| `velocityGradient[9]` **索引转置** | **静默数值错误** | 自查读 kernel |
| 面密度/速度从 `limiter` 取（越界） | **静默越界** | 自查读 `PrepareFaceValue` |

后两个会编译通过、契约测试通过、求解器跑出错误残差且不报警。共同根因：
**view 的字段装配内联在求解器循环里，契约测试直接手工构造 view，装配代码零覆盖**。

补救：抽出 `codes/uturb/include/SaStateAssembly.h` 的 `FillSaVelocityGradient`，
配 `tests/solver/sa_state_assembly_test.cpp`（4 例，含 `VelocityGradientIsNotTransposed`
和一个守卫哨兵值互异性的用例）。验证作业 `of-asm-verify` 会同时跑正确映射和
转置映射，要求后者必须失败。

**已知缺陷（本轮修掉或记录）**：

- `HipDrain` 的 `sinceDrain` 是普通 `static` 而非 `thread_local`（多 rank 同进程
  共用计数器）；同文件 1941 行的 `LocalPinnedStaging` 用了 `thread_local`，说明是遗漏。
- `CalcCellTemperature` 原来**无条件**做 7 MB 的 D2H，而全仓库没有 host 读者；
  已改为按需。
- `interfaceMask` 未填：多 zone 界面面当前按非界面处理（旗舰 1 zone 不受影响）。

**上一轮遗留**（SA 内核族本身，仍有效）：`TurbSaArithmetic.h` 是 host/device
共享标量算术（`ONEFLOW_TURB_SA_HD` 在普通主机编译器下展开为空，否则 HIP 拒绝从
`__global__` 调用 host 函数）；`TurbulenceBackend` 是 SA 专用 seam（**没有**加到
`EulerDomainBackend`/`FluxBackend`）。证据：`evidence/probe/128156945/`。

**本轮（2026-11-14 下午）交付：网格加密 + M6 DCU 弱扩展测量**。全部数字与证据链在
历史 TODO §5.6，要点：M6 单 zone `.ofl` 加密至 2,359,296 cell（体积守恒 5e-15）；METIS 分区
8/32/64 zone 带 manifest；单 zone XB 残差 3.07e-17；64-zone XB 全绿；加密 vs 原
CPU 气动力一致 0.2%。弱扩展（36,864 cell/rank，8 ranks/卡，同作业交错）：1 卡 17.0 →
8 卡 22.2 ms/步（+31%），归因全部在 `interface_mpi`（2.12→7.70 ms/步）。
新 fixture：集群 `fixtures/m6-refined-r1/`（只读登记）。

**2026-11-15 更正：定向合并 2N/4N 可行**。此前"严格判否/O 型接缝不可积"是两个
生产缺陷（坐标集合匹配丢失节点身份、DSU 误判可满足网络为 UNSAT）造成的，
均已由用户提供的回归定位并修复（`test_merge_dirs.py` 回归 + 200 小图穷举交叉验证，
`f858f9a44`）。修正后全网格核验：873,216 面全部解析为双射约束，**2N/4N 模式均
PASS**，2N 分配独立复核 0 违反。下一动作：按 PASS 的方向分配实现 2N 生成
（合并逻辑沿真实 face ID 与 A(c) 映射），随后 4N，再 refine8 派生 16N/32N。

**下一刀（按归因排序，详见 历史 TODO §5.6）**：① `partition_type=2` 重切 64 zone 压 halo
（64N 档 512 zone 的等块重切待做，若开做需先扩 `StructuredPartition` 支持非立方块排列）；
② A4 设备驻留宏步循环；③ 大 per-rank 网格的 rank/卡 packing。另有一个独立小修候选：
坑 41（单 zone HIP aero）——`UNsInvFlux::CalcInvFace` 的 `mirrorBoundaryGradient`
依赖 `packedExchange`（要求 nZones>1），单 zone 时 host 梯度镜像从不更新，aero 力积分
读到陈旧梯度；多 zone 不受影响。修复需过单 zone XB + 8-zone 门禁 + paired。

**H1 当前交付**：新增显式开启的 MG 根域 inclusive/exclusive wall accounting，
保留原 stage/task TSV，另写 `.host.tsv`；修复两个旧漏登记类别。集群 HIP/计时契约、
trace、200 步 XB、CPU normal/strict、50 步稳定性和五对 1000 步归账通过；
验收数字与口径只维护在 历史 TODO §5.2。CI 状态与修订分开记录，不把构建当数值验证。

**H2 终结（`d41650729`）**：多卡 VMFault 的根因**不是多卡**。`ONEFLOW_HIP_DEFER_STATUS=1`
（默认）下 `stateUpdateStatus` 惰性分配，而 deferStatus 分支完全跳过 `CopyFromHost`，
`Data()` 返回 nullptr 进 `atomicExch`。单卡也会触发；回溯指向 `GradientCopyBoundaryKernel`
只是因为它在队列里。同时发现原来的最小复现算例 CFL 过高（`cfled=10`），
CPU 侧 step 2/13 就非物理，因此它**不具备 oracle 资格**。固定 `CFL=0.01` 后，
多卡 HIP（2 节点/2 DCU）与 CPU 残差逐位一致。

**扩展性已定性（2026-11-14）**：由每 rank 的 cell 数决定，不由卡数决定。
M6 强扩展 1/2/4/8/16/32 卡 = 11.70/7.60/9.87/18.73/24.08/30.62 ms/步。
**墙钟最低点落在 2 卡，但这不是"最优配置"——是退化**：每卡吞吐随卡数单调下降
（22.9 → 17.1 → 11.0 kcell/ms/卡，1/2/4 卡），卡越多单卡效率越低。
把"墙钟最低"读成"推荐只用到 2 卡"是逻辑错误（2026-11-15 更正），正确目标是让
每卡吞吐随卡数保持恒定。立方体 884k 强扩展同形状；立方体弱扩展（110,592 cell/rank）
2→32 卡迭代 13.78→13.68。求解器核心 1→32 卡可扩（M6 2.5×、立方体 6.0×），
差额全在接口交换。经验阈值：每 rank ≳5 万 cell 正常，≈1 万开始变差，<5 千明显变差。
详细数据见 历史 TODO §5.5。

**历史任务：原生六面体加密（后续结果见 历史 TODO §5.6、§5.9）**

当时仓库无网格加密/导出能力（`IsSupportedConversion` 里 `OneFLOW` 只作为 target；无细化工具；
集群无 gmsh），所以 M6 弱扩展和 D1 都卡在这里。`.ofl` 二进制布局已实测确认（历史 TODO §5.5 有字段表），
M6 为**六面体**（面全为四边形），8 zone / 294,912 cell / 306,577 node。

写 `refine_ofl.py`：读 8-zone `.ofl` → 每个 hex 剖 8 个 → 棱/面/体心去重、边界与 interface 保持 → 按同格式写回。
预期 2,359,296 cell，**仍是单一连通算例**。验证顺序：解析自检 → 单 zone 200 步比残差 → 分区成 64 zone 测弱扩展。

其余已停：A4（版图已变，优先级下降）、通信优化（7%）、halo 均衡（r=−0.007 无相关）、kernel 融合（依赖阻塞）。

原始证据与分析脚本留在本轮固定 category 目录。清理可再生成的 trace/case 大文件后，
保留每 rank breakdown、汇总、日志与复现脚本；继续使用唯一 `build/app`。
本轮临时 trace/case 已清，快照/日志已轮转；旧 `_attic` 清理仍在独立的
`of-probe-cleanup` CPU 作业中，接手先查该作业与 evidence/probe 下的结果。

## 5. 已知坑与设计约束（只加有长期价值的）

**粘性探针注意事项（2026-11-19，2026-11-19 更新）**。通量探针 `ONEFLOW_VISCOUS_PROBE_FILE`
会覆盖输出，CPU 与 DCU 用不同路径，同一进程开启自检时也须避免互相覆盖。
壁面探针 `ONEFLOW_WALL_FORCE_PROBE` 按气动力调用持续写入，当前没有末步过滤或
多 rank 隔离。先在单 rank 使用，按输出步对齐。自检重新执行 CPU 粘性循环，
需核对共享状态副作用及 NaN/Inf 处理；探针运行不进入性能对照。

自检输出现在是：`SELFCHECK step= zone= faces= equ=`、首差面/最差面（含方程与
相对误差）、`nonfinite_device/host` 计数，最后一行 `SELFCHECK result=PASS|FAIL`。
判定逻辑在 `CompareViscousFaceFlux`（`codes/uns/include/ViscousFluxSelfCheck.h`），
**以 result 行判定，不要再用旧日志里的 `differing_faces=0` 推断通过**。

- 旧的相对误差判据会让 NaN 静默通过：`rel > tolerance` 对 NaN 恒为假。非有限值
  现在单列计数并直接 FAIL，回归测试 `tests/solver/viscous_flux_selfcheck_test.cpp`。
- 壁面标签只对 `bcType == SOLID_SURFACE` 的面打印；`face < nBFaces` 还包含入口、
  出口和 interface，不可当成壁面。
- `nonfinite_*` 为 0 且 `differing_faces` 非 0 时，先看首差面是否落在壁面或
  interface，再对齐该面的输入时间层。

**测量与工具链（2026-11-08 新增，全部踩过）**

12. **`make OneFLOW/fast` 不会重编 `codes/` 下的源文件**：`<dir>/fast` 是 CMake 的"仅本目录"
    目标，只重链 `main/src/OneFlow.cpp`。必须两步：`make OneFLOW_RuntimeObjects` 再
    `make OneFLOW/fast`（`scripts/build-app.sh` 已修）。验证办法：`strings bin/OneFLOW |
    grep <新符号>`，没有就是没编进去。
13. **`StageProfiler` 会静默丢弃未登记的类别**：名字必须先出现在 `StageProfiler.cpp` 的
    `Totals()` 表里，否则计时被丢掉、TSV 里看起来"该段零成本"。现在 `Flush()` 会打印
    `dropped N timer records for unregistered category ...`，看到就补表。
    H1 复验发现 `residual_zero` / `host_mirror_refresh` 两个旧类别漏登记，已补齐；
    表长改为从初始化列表推导。性能归账必须确认日志无 dropped 记录，再解释耗时。
14. **启动固定开销 1.4–1.6 s**：100 步口径把它摊成 +14 ms/step、200 步 +7。所有 ms/step
    必须写步数并扣启动（同一算例跑 1 步与 N 步取斜率）。旧的"8.93 ms 固定基数"就是这类读数。
15. **跨节点/单发噪声很大**：同一配置跨节点差 32%（13.12 vs 9.93），DCU 单发可达 3–4 倍
    （5.87 / 23.32 / 7.54）。**只比同作业内交错重复的中位**，跨节点单发一律不算数。
16. **延迟 drain 会让队列溢出**：`ONEFLOW_HIP_DEFER_SYNC=1` 下不设上限时，16/32 rank 会在
    `libgalaxyhip` 里 AQL 溢出 abort。`ONEFLOW_HIP_MAX_LAUNCHES_PER_DRAIN`（默认 32）必须保留。
17. **GPU-aware 探针通过 ≠ 求解器可用（该结论的证据基础已失效，2026-11-15 复核）**：
    旧记录称 `gpu-aware-mpi-probe` 2 rank 单节点 PASS，但求解器开
    `ONEFLOW_ENABLE_UNS_GPU_AWARE_MPI=1` 后 8/32 rank 都 abort。
    **但那次 abort 发生在 H2 修复（`d41650729`）之前，根因是 `DEFER_STATUS` 下
    状态缓冲空指针（VMFault），与 GPU-aware 无关**——同一类 abort 在 host 路径也会
    触发。所以"8/32 rank abort"不能归给 GPU-aware。
    **2026-11-15 在当前构建上重测**（`probe/128116476`）：2 rank 设备指针 **pass=1**
    通过（UCX 默认传输）。但**8 rank/卡 的生产密度至今没有测过**——仓库里的探针
    程序写死 `requires exactly two MPI ranks`，4/8 rank 会直接 rc=4 退出
    （**这是探针限制，不是 MPI 失败，不要再把它读成失败证据**）。
    待验证：多 rank/卡的设备指针路径。在此之前保持开关关闭，但**理由是"未验证"，
    不是"已验证不可用"**。

**集群作业与算例（同上）**

18. **派生算例必须显式设 `time_integral = 1`**：`test/cfd/m6wingroe_sa` 自带 `= 2`（LU-SGS），
    设备驻留链路只在 1 下启用；漏了会跑成 host LU-SGS 并伪装成"多卡扩展性差"（实测 2.6 s/step）。
    同时显式设 `rk_stage=3`、`ischeme=5`、`ilim=0`、`vismodel=0`，并删掉 `solver.txt` 的
    `TurbSolver`。
19. **分区参数分两处**：`gridObj = 3` 在 `script/grid.txt`，`npartproc`/`part_uns_file` 在
    `script/part.txt`（control.txt 里后者后加载，**覆盖**前者）。只在 grid.txt 写 `npartproc`
    会得到 `The number of partitions should be greater than 1`。
20. **求解算例必须带 `grid/bcRegionMap.txt`**：缺它会在 `UNsBcSolver::CalcBcRegion` 段错误，
    在分片导出阶段表现为 SIGSEGV。
21. **单 zone 算例 `nIFaces = 0`**，不做 halo 交换 —— 不能当弱扩展基线。弱扩展必须每核/每 rank
    固定 cell 数（≈8,000），且以多 zone 点为基准。
22. **分区作业的 case 目录要传绝对路径**：`... OneFLOW 0 .` 会静默什么都不写（rc=0 但没有输出
    文件），看起来像分区失败。
23. **cell 平衡 ≠ halo 平衡；扩展性用 `partition_type = 2`**（2026-11-09）：METIS 在等块
    设计下仍给出 1212–3165 的 per-rank `nIFaces`。新增结构化分区（D1a，`partition_type = 2`）
    按 cell-count 分位切 a×b×c，cell 严格等分、接口数由块位置决定（角块 1200、边块 1600、
    体块 2400）。固定 32 rank 同作业配对：平均 `nIFaces` 1786 → 1500，每步
    **61.32 → 58.62 ms（−4.4%，3/3 对为正）**。**扩展性表必须带 `nIFaces/rank`**。另两个坑：
    - **8 rank 不能当弱扩展基线**：只占节点一半核心，睿频更高；本作业 8→16 时间 +22%，
      远大于 halo +17%；
    - 开放边界下每 rank halo 会从角块的 1200 饱和到体块的 2400，弱扩展 99–100% 只能在
      整节点且同 halo 的配置间主张（D1b）。
    - **瓶颈已定位（2026-11-09）**：同一套等块分区、整节点 32/64/128 rank，
      `ONEFLOW_DIAG_NO_COMM=1` 后计算段完全持平（139–141，相对量），接口项
      13.5 → 22.9 → 27.9 ms（+106%），而面数只 +27% → 剩下的是跨节点延迟。
      本轮顺序已由 2026-11-10 的记录修正为 H1 → H2 → A4；A2′ 降权，见历史 TODO §5。

1. **设备权威不能从初始化就开始**：bound 状态的第一次 `FullField` 上传依赖 host 的
   ghost 区；跳过条件必须是 `HostMirrorStale(grid)`，否则第一步就偏离 CPU oracle。
2. **设备写 ghost 而未回写 host 时必须把镜像标脏**（`NoteDeviceGhostWrite(false)`），
   否则 contract 诊断"通过"而 host ghost 落后（1e-16 级），后续整场上传会把旧值灌回设备。
3. **`UploadBoundaryState` 不要按 generation 差异回灌 host 场**：设备权威时 owner
   generation 每 stage 前进，会把上一 stage 的旧场覆盖到设备。
4. **resident bcQ 复用已撤回**：solver 与 backend 两侧判据一旦不一致会上传未初始化
   buffer，3 步内 non-physical；收益只有 0.3–0.7 ms/step，不值。
5. **异常排查**：`Simulation::RunImpl` 的 `catch(...)` 会先 `MPI_Finalize` 再重抛，
   跨 rank 抛异常时 stderr 可能永远打不出来；临时在 catch 里先打印 `what()`。
6. **capability 边界**：只有 `NS_SOLVER / 5 方程 / (Lax–Friedrichs 或 Roe) / ilim∈{0, Barth} /
   vismodel=0 / 单 grid level / 稳态显式 RK / 1 rank 1 zone` 走设备；Barth 的数值还差
   `1.1e-6` 未过门禁，默认不放行；Venkat / HLLC / AUSM 未实现，显式请求 fail-fast。
   其它残缺配置回退 host。旗舰算例 `m6wingroe_sa` 是 `LU-SGS + Roe + SA`，仍在 host。
7. **短基准的假象**：`bnd_resize` 的首个 `hipMalloc` 会被平均进 ≤100 步的基准
   （10–15 ms/step 量级），任何基准先 warmup 至少一步。
8. **fixture**：固定登记的是 m6wing 8 zone（294,912 cell）；16/32 zone 按 runbook §1b
   生成并做 manifest 校验。≥1M 带 `.skl` 的合法 fixture 仍缺，多卡还须先解决 H2。
9. **加速器构建用统一 SDK 接口，别再手拼**：运行时 `ONEFLOW_ACCEL_BACKEND`
   （AUTO/CPU/HIP/CUDA/KOKKOS）与构建期 `ONEFLOW_ACCEL_SDK_ROOT`（+ 可选 `_ARCH`/`_FLAGS`）
   是两条独立的轴；**SDK 种类不单独设参数**，由 `ONEFLOW_ENABLE_HIP/CUDA` 或前缀内容推断。
   `cmake/OneFLOWAccelSDK.cmake` 从前缀派生 `ROCM_PATH`、厂商包前缀、`amd_comgr`/
   `AMDDeviceLibs`/`hsa-runtime64` 与 `CMAKE_TRY_COMPILE_PLATFORM_VARIABLES`（CUDA 分支派生
   `CUDAToolkit_ROOT` 等）；不请求就完全不生效。一键入口：
   `scripts/build-accel.sh --arch gfx906 --flags "$ONEFLOW_DTK_HIP_FLAGS"`。
   原理与故障对照见 runbook §1.2 / §4。
10. **strict（1e-15）历史越界已解决，最近记录为 5/5**：`plateuns2dslau2_34950_35000`
   已解决：声速**改用一次除法**（精度更高：与参考逐位一致率 100% vs 乘倒数 86.8%；
   护栅对物理密度无影响），并按“精度优先、可更新历史精度数据、strict 判据不变”重生成
   了**该算例**的参考曲线与基线条目（`validated_commit=7bd358c0d`）。现在
   **normal 5/5、strict 5/5（绝对 1e-15）**。其余四个算例参考未动、仍在判据内。已定位两处上游重构：
   `f7ae944c0`（双边声速）与 `6ca93f57c`（面谱半径/局部时间步）；探针数值非单调
   （中间提交会回到 `0.0`），说明多处末位扰动互相抵消。**已排除 fast-math**：
   当前代码在 `-Ofast` 与 `-Ofast -fno-fast-math` 下给出完全相同的 `1.778709301357682e-15`，
   旧代码两种设置都是 `0.0` → 差异是源码里运算结合顺序的确定性 ≤1 ulp 改变，不是编译器
   重排。最终按上述精度优先流程只更新该算例参考数据，未放宽 1e-15 判据；
   2026-11-08 复跑 normal/strict 均 5/5。本 HEAD 没有新增复跑证据。
11. **gradient host 镜像默认只在最后 RK stage（A1.1）**：`ONEFLOW_HIP_GRADIENT_MIRROR_ON_DEMAND=0`
    才恢复每 stage 镜像。不要用上一会话的 31.64 ms/step 当本轮 A 侧；本轮同节点
    对照是 29.86 → 26.72。分项文件如果路径里没有 `%r`，8 个 rank 会写同一个文件，
    读到的是其中一个 rank，不是 8 rank 均值。
12. **`LOAD_Q`/`LOAD_RESIDUALS` 不是死写（A1.2）**：前者填 host RK 基态，后者把命名残差置成
    `-FIELD_RHS`，而 fused flux 核把它当**累加基底**上传。跳任务必须同时把基底改成
    设备侧归零，两边由 `DeviceOwnsExplicitRkBookkeeping()` 一个谓词控制
    （`codes/solver/include/ResidualBasePolicy.h`）。只改一边，读数会变成残差在上一
    stage 上累加：trace 仍过，200 步 XB 会拓到 `5.52e-4`。
13. **restart 历史层是 host 派生数据，dump 时必须从设备场重建**：restart 流写
    `q1`/`q2`，读盘用 level 1 复原当前场；设备链只回填当前 `q`，历史层会停在旧值
    （曾导致 auth-on dump 与 auth-off 差 41%）。现在 `Restart::Dump` 在
    `NsDeviceAuthoritativeFlowEnabled()` 为真时先回填 `q` 再 `InitUnsteady()` 重建。
    **判据要用所有权谓词，不要用刷新函数的返回值**：镜像可能已被别的读者刷新，
    那样会静默跳过修复（踩过一次，读数逐位不变）。
14. **`mean_ms` 会被周期性输出污染，必须连 `count` 一起读**：`POST_PROCESS` 曾显示
    每步 7.7 ms，实际总 2193 ms / 300 步、只在第 0/50/100… 步输出；把 `nVisualSave`
    设成不输出后降到 0.055 ms/step。误把均值当每步开销会选错优化方向（踩过一次）。
    `ONEFLOW_STAGE_BREAKDOWN_FILE` 的 `.tsv` 是嵌套计时，不能相加求总量。
15. **GPU-aware MPI 的"更贵 48–50 vs 15.7"读数不可作为定论（2026-11-15 复核）**：
    旧记录：同节点两 rank 设备指针直达 `pass=1`，但接到真实界面交换后 48–50 ms/step，
    host staging 15.7，结论"设备指针路径实际更贵，保持关闭"。
    **该对照是在 `DEFER_SYNC`/`MAX_LAUNCHES_PER_DRAIN` 成为默认之前测的**，当时
    每 kernel 后强制 drain、队列深度无上限，与现在的排队行为不同，所以"更贵"的
    结论**需要在新默认下重测**才能引用。
    当前状态：设备指针能力本身**可用**（`probe/128116476` 2 rank pass=1），
    端到端收益未测。开关 `ONEFLOW_ENABLE_UNS_GPU_AWARE_MPI` + `..._CONFIRMED`
    保持关闭，理由是"未验证收益"，不是"已验证更慢"。
16. **固定基数必须先扣启动再归账**：旧通信开关对照约省 2.2 ms；
    `rk_stage=1` 的 8.93 ms 含短基准启动摊销，已由 todo A1.7b 修正，不能继续当稳态基数。
    单配置通信开关差值不是 A2 收益保证；H1 同配置归账已完成，结果见 历史 TODO §5.2。
17. **别用 2 个样本下结论**：8-zone 同配置两次样本 17.59 / 14.28（跨度 23%），
    16-zone 与 32-zone 的 7% 差异完全在噪声内。要比较配置，至少 5 个样本或做配对交替。
18. **`zones % ranks == 0` 只是配置层放行**：`Zone::CheckZoneRankCount` 只查整除，
    但多 zone/rank 在设备路径曾直接 Fatal（`HIP gradient face pack has no resident gradient`）。
    已改为“缺 resident gradient 时回退主机镜像”，能跑完；**但实测更慢**
    （16.45 → 17.62 → 25.24 ms/step 于 1/2/4 zone per rank），因为它不减少 halo、
    只按 zone 重复每步工作。不要把它当扩展手段。
19. **`StageProfiler` 的分类表尺寸必须等于初始化项数**：声明 76 槽只填 73 项时，
    剩余槽 `name == nullptr`，`Record()` 里 `strcmp` 解引用空指针 → 一开分项就段错误；
    而且**未注册的分类名会被静默丢弃**，据此判断“某段没执行”会得到错误结论（踩过一次）。
20. **上游把 runtime 源码移进 OBJECT 库后，HIP 要在 `tests/` 作用域再次 enable**：
    `CMAKE_HIP_LINK_EXECUTABLE` 是目录作用域变量，`codes/` 里 enable 不够，否则
    generate 阶段报错。修法见 `tests/CMakeLists.txt` 顶部。
21. **设备权威 q 下不要在主机算 φ**：主机梯度镜像只覆盖边界相邻单元、且只在有读者的
    步更新。用它算限制因子得到的是陈旧值，200 步残差偏 `8.53e-05`（100 步 `3.05e-05`，
    逐 step 累积）。限制器要么整段留 host（回退），要么把 φ 计算搬上设备，不能半搬。
    设备路径现在只在 `ilim=0` 时启用，`UseHipGradient` 里写明了这个原因。
22. **限制器在面上是"两侧取最小"，不是各用各的**：主机 `LimField::CalcFaceValue`
    用 `ModifyLimiter(phil, phir) = min(phil, phir)`（`ug.ireconface != 0`），
    并把同一个因子用在两侧重构上。设备若按左右各自的 φ 计算，200 步残差会差
    `8.5e-5`（实测）。改设备限制器时先对齐这条，再谈别的。
23. **post 步的 interface 轮不是冗余**：显式 RK 每个 stage 末尾都交换一次，看起来
    第 4 轮是重复；实测跳掉后 HIP 与 CPU oracle 在 200 步上偏离 `5.54e-7`（aero `4e-5`）。
    它参与 legacy 的一轮 lag。这类“看起来多余”的交换，改之前先跑 XB。
24. **冷启动第一次 1000 步不能当 paired 端点**：同一配置在冷 case 目录上首跑
    `27.0 ms/step`，同一作业第二次 `21.1`。跨作业/跨冷热的绝对值不可比，只有
    交替 interleaved 的 A/B 才算数；报数字要写清是不是热态。
25. **`HIP.MainSolverFluxOneCall` 的 boundary contract 已修，不要再按数值红灯查**：
   旧读数 `reconstruction 8.600e+00` 是第 4 个方程的期望 `bcQ`，设备侧读到 0；
   `stale-op 3.000e+00` 是旧夹具里的 `+3` 偏移，不是失败判据。原因是烟测仍按
   `nFaces` 步长填 `bcQ`，而 `1c470414a` 之后的上传按连续的 `nBoundaryFaces` 步长读。
   夹具必须与生产路径一样：固壁面是 face 数组前缀，`bcQ[equation] ==
   bcQ[0] + equation*nBoundaryFaces`。不满足时后端会在上传前拒绝，不再静默读错方程。
26. **`_attic` 里删构建树要留足时限**：本轮探测留下 ~3.6 GB（多份 `build/` + git
   bundle），`rm -rf` 在 15 min 的作业里被 scheduler 切断。清 `_attic` 用
   `scripts/cleanup.sh` 或单独申请 ≥30 min 的 CPU 作业，不要塞进验证作业尾部。
   另：`snapshots/` 已超过保留 5 个的策略，下次 `cleanup.sh` 轮转。
27. **`cmake` 必须在 `PATH` 上**：`env.sh` 只导出 `$CMAKE` 变量时，`cmake` 会落到系统
   2.8.12.2，`make` 自动 reconfigure 还会去用缓存里已删除的旧 cmake 路径。
28. **`POST_PROCESS` 的 `UPDATE_UNSTEADY_FLOW` 调度只值约 0.1 ms**（2026-11-09）：同二进制
   pair（跳过 dead dispatch vs 排队执行），1000 步、8 rank/1 DCU：13.49 vs 13.26 ms/step，
   `POST_PROCESS` 3.83 vs 3.93 —— 在噪声内。之前那一版开关把 handler 也关掉了，测到的
   2.0 ms 是 history shift 本身。**`POST_PROCESS` 的 3.8 ms ≈ 2.0 shift（默认已跳）
   + 1.5 post 轮 pack/unpack（`post_upload` 0.9 + `post_download` 0.7）**；要动的是
   post 轮的 upload/download，不是命令调度。
29. **rsync 打补丁时不要多传同名目录层级**：`rsync a/codes/ns/src/NsRhs.cpp host:.../codes/ns/`
   会把文件放到 `codes/ns/NsRhs.cpp`（与 `codes/ns/src/NsRhs.cpp` 并存），
   `ConstructSolutionDirTree` 递归收全 `codes/**`，于是同一个符号被编两次、链接报
   `multiple definition`。传文件要显式落到 `.../codes/ns/src/`。
   另外：**`bin/OneFLOW` 被链接失败删掉后，后续所有作业都会 `No such file or directory`**，
   先修链接再测。

30. **非阻塞接口交换实测更慢，且第二版段错误**（2026-11-10，已回退，代码在 HEAD）：
    - 消息计数（`ONEFLOW_MPI_TRACE_FILE`，32/64/128 rank、1000 步、结构化等块）：
      **每步 24 条 send + 24 条 recv、0.66–0.78 MB，三者几乎不随 rank 变**；但 `recv`
      时间 5.7 → 8.9 → 15.8 ms/步，`send` 只有 0.2–1.4 ms。
    - 所以瓶颈**不是消息条数、不是带宽、也不是"send/recv 串行"**：阻塞路径本来就是
      "pack 一个邻居→立刻发"，已经和邻居的接收重叠；改成"先全部 pack 再统一 Isend"反而
      把流水线打断了（32 rank 1000 步：61.2 → 64.3 ms/步）。
    - `recv` 随 rank 增长是**邻居到达偏差（负载不均）**：等最慢的那个邻居。杠杆是
      **减轮次（每步 4 轮）和均衡每 rank 的 halo**，不是换 Isend/Irecv。
    - 第二版（持久 databook + 先 post 全部 recv）直接段错误（rc=139），不要再基于
      "持久 DataBook 复用"改这一层。

31. **`rocprof` 需要"module + 复制工具树 + 改一行"才能用（2026-11-10，DTK 26.04）**。
    **module 是必要但不充分的**：`module load compiler/dtk/26.04` 让 `rocminfo` 正常
    （它设 `HYHAL_PATH=/opt/hyhal`，env.sh 没有），但 `tblextr.py` 仍硬编码
    `/opt/dtk-26.04/bin/rocminfo`，而节点 `/opt` 只有 `dtk-22.04.1`，所以照样
    `Profiling data corrupted`。可用配方（已在集群验证，1 rank 拿到 253 条 dispatch）：
    ```bash
    module load compiler/dtk/26.04
    cp -a $DTK/rocprofiler $HOME/oneflow/tools/rocprofiler-26.04   # 8.8 MB
    sed -i "s#/opt/dtk-26.04/bin/rocminfo#$DTK/bin/rocminfo#" \
      $HOME/oneflow/tools/rocprofiler-26.04/libexec/rocprofiler/tblextr.py
    HIP_VISIBLE_DEVICES=0 $HOME/oneflow/tools/rocprofiler-26.04/bin/rocprof \
      -o out.csv --stats env ONEFLOW_ACCEL_BACKEND=HIP OneFLOW 0 <case>
    ```
    产出 `out.csv`（每次 dispatch）+ `out.stats.csv`（按 kernel 汇总）。8 进程一起 profile
    仍会 OOM，先 profile 单 rank。
    历史失败记录（不要再试）：`$DTK/rocprofiler/bin/rocprof` 直接跑（路径坏）；
    `rocprofv2 --kernel-trace`（1 rank 3.5 GB OOM、`-d` 不落盘）；`-t`/`-d` 都不保留原始
    数据；`24.04.3` 不带分析器；`25.04.4` 硬编码 `/opt/dtk-25.04.2`。
    - 工具在 `$DTK/rocprofiler/bin/{rocprof,rocprofv2}`，**不在 PATH 上**；
    - v1 的 `libexec/rocprofiler/tblextr.py:839` 硬编码 `/opt/dtk-26.04/bin/rocminfo`，
      本集群 DTK 在 `/public/software/...`；即使改路径，节点上 `rocminfo` 也报
      `C-3000 module is NOT loaded, possibly no HCU devices`；
    - v2 `--kernel-trace` 在 1 rank / 3.5 GB cgroup 下直接 OOM（8 rank 也 OOM）；
      1 rank / 28 GB 能跑完但 `-d` 目录不落 CSV。
    - 绕路（`-t` 保留原始数据 + 打补丁的 tblextr + 跳过 metadata）能跑通解析，
      但原始 `results.txt` 里只有 memcpy、没有 kernel dispatch（profiler 在失败前
      没 flush）。
    **当前结论**：补丁后的 rocprof v1 可提供 dispatch 计时；以上失败路径保留作避坑记录。
    仅有 kernel 耗时不能回答离 roofline 多远，还需 FLOP/访存量及硬件上限证据。
    `hipEvent` 可作为独立计时方案，但本轮未新增实现或测量。

32. **stage 分项表不是 GPU 时间**（2026-11-10）：`gradient/reconstruction/flux_assembly`
    这些 host 计时器量的是"入队 + 等待"，`ONEFLOW_HIP_DEFER_SYNC` 一开一关能差 8 倍
    （`gradient_kernel` 0.18 vs 1.48 ms/step）。短程 rocprof 的全部 kernel dispatch
    时间和为 **1.77 ms/step**，不是三个 kernel 的合计，也不是同配置直接测得的忙碌率。
    **kernel 相关排期应使用设备侧计时证据**
    （配方坑 31），不要用分项表；分项表只用来分 host 段。
    **补强（2026-11-15，8 rank / 1 卡分项实测）**：`HipDrain` 的计数器是
    **每进程（每 rank）** 的 static，`DEFER_SYNC=1`（默认）下每 32 次启动才真正
    `hipDeviceSynchronize` 一次。但 `hipDeviceSynchronize` 是**设备级**同步：
    它等的是整张卡空闲，包含同卡其他 rank 的在飞工作。所以 drain 落在哪个 rank、
    落在那一步的什么位置，决定了这次排空要替同卡其他 rank 背多少排队时间。
    后果：同一作业内、每 rank cell 数完全相同（36,864）时，`gradient_kernel` 在
    8 个 rank 间差 **5.46 倍**（16.7 vs 91.1 ms/1000 步），`recon_kernel` 差 2.85 倍。
    **因此 per-rank 的 `*_kernel` 行逐 rank 横比无效**；要比较必须同时跑
    `DEFER_SYNC=0` 对照。同理 `swap_download` 是延迟队列下一个阻塞点的承接，
    不是它自己的拷贝成本。

33. **上游 `ProcessRunner.cpp` 假定 glibc ≥ 2.29**（2026-11-10 合入 upstream 后暴露）：
    `posix_spawn_file_actions_addchdir_np` 是 glibc 2.29 才有的扩展，而昆山节点是
    CentOS 7.6 / glibc 2.17；原来的 `#if defined(__GLIBC__)` 判据不够，直接用会编译失败
    （`‘posix_spawn_file_actions_addchdir_np’ was not declared`）。已改成按
    `__GLIBC__/__GLIBC_MINOR__` 版本判断，低于 2.29 走 `ENOTSUP` 分支。以后合并 upstream
    先在本集群跑一次 `OneFLOW_RuntimeObjects`，不要假设上游在 CentOS 7 上编译过。

34. **每步残差/气动力输出会让多卡测量完全失真**（2026-11-14，本轮最大坑）：
    `fixtures/m6-8zone-cpu8/case` 的 `cfd.txt` 是 `nResSave=1`、`nForceSave=1`，
    `DumpFile()` 每步对 `results/res.dat` 做 open/append/close。256 rank 时
    `DUMP_RESIDUAL` 实测 **37.2 ms/步**，占整步 93%，而且随 rank 数增长（共享存储争抢），
    看起来就像"卡越多越慢"。基准必须走 `prep-case.sh`（它把三个 `*Save` 设成 100000）。
    **本轮前半所有 M6 扩展性数字都因这个作废。**
35. **`sed` 改参数必须用 `[[:space:]]+`**：fixture 里是 `int    nResSave = 1;`，
    写成 `^int nResSave = .*` 不匹配且不报错。改完用 `grep` 确认。
36. **集群资源按节点算**：`max_cpus_per_node=32`、`max_dcus_per_node=4`；
    `--dcus` 是每节点值；`mpirun --host` 必须带槽位（`h1:32,h2:32`），否则 `-np 64` 报槽位不足；
    内存 ≤2 GB/core。`NDEV = local_rank % 4`，不是 8。
37. **GPU 占用率要按分配到的卡读**：**正确路径是
    `/sys/class/drm/cardN/device/gpu_busy_percent`**（`device/` 这一层必须有；
    写成 `cardN/gpu_busy_percent` 读不到文件，2026-11-15 实测确认）。
    不能取 `/sys/class/drm/card*/device/gpu_busy_percent` 的最大值
    会把同节点其他租户的 100% 当成自己的（本轮因此误判一次）。先看哪个卡在我们的运行期间
    从 0 起变化。
38. **作业脚本里的 `python3` 必须先 `source scripts/env.sh`**，否则 `command not found`。
39. **`RecordTask` 是有调用者的**（`CmxTask.cpp:70` 的 `TaskDispatchTimer`），`tasks.tsv` 只记
    dispatch 级；父任务包住的子命令要自己在 `CMD::ExecuteCmd` 里记（`727c1ed75` 已做，
    未开 breakdown 时零开销）。排查时先看 `.tasks.tsv` 再看 `.host.tsv`。
40. **`.ofl` 有两种写法的残留**：单 zone 的 `m6wing.ofl` 是旧格式（zone 头多 2 个 int，
    朴素解析会读出 `nNodes = 文件大小−20` 这种值）；当前 reader/writer 与 `m6wing-8.ofl`
    用 3-int 布局（`zoneBytes, ?, nNodes, nFaces, nCells`，实测见 历史 TODO §5.5）。
    做网格工具以 8-zone 文件为准。
    **修正（2026-11-14）**：仓库当前 `test/cfd/m6wingroe_sa/grid/m6wing.ofl` 逐字节闭合
    解析就是 3-int 布局（nLength 精确、EOF 零剩余、nIFaces=0）；`refine_ofl.py` 直接读它
    通过全部自检。该坑里的“旧格式”描述不适用于仓库现行副本，保留作历史记录。
41. **单 zone HIP 的 aero 力积分读陈旧 host 梯度（2026-11-14 定位）**：
    `UNsInvFlux::CalcInvFace` 里 `mirrorBoundaryGradient` 以 `packedExchange` 为前提，
    而 `packedExchange` 要求 `crossZoneGhostGradient`（`nZones>1 || nLocal>1`）。
    单 zone/单 rank 时三个开关全 false，`CalcGradHip` 既不下载也不镜像梯度 →
    host `dqdx/dqdy/dqdz` 停在初值，`AeroForceTask` 的壁面压力重构
    `wp = q[lc] + dqdx·dx + …` 退化为格心值 → HIP aero 与 CPU 差 ~8%
    （cl 8.2e-3），而残差 1e-17 逐位一致（残差不走 host 梯度）。**多 zone 不受影响**
    （镜像激活，XB aero=0.0 历史仍成立）。单 zone XB 门禁因此只能拿残差当判据；
    修复方向：`mirrorBoundaryGradient` 改为“存在 host 读者即镜像”，不依赖跨 zone 前提，
    修复后需补单 zone XB + 8-zone 门禁 + paired。
42. **OpenAPI 提交的 `--cpus`/`--memory` 是每节点值（2026-11-14）**：`--nodes 2 --cpus 64`
    会被拒（`max_cpus_per_node=32`），`--memory` 同理按节点算（128G/节点直接报
    “Memory specification can not be satisfied”）。64 rank 的作业用
    `--nodes 2 --cpus 32 --memory 60G`。上传覆盖同名文件要 `--cover`；偶发超时原样重试。
43. **弱扩展测量的 A/B 必须同卡内 rank 密度（2026-11-14）**：固定 cell/rank 的 A/B 配对
    若一侧 1 rank 独占卡、另一侧 8 ranks/卡，差异主要来自卡内分时而非扩展性
    （本轮 294,912/rank 系列因此只能作参考）。设计系列时保证两侧都是 8 ranks/卡，
    或从 16 rank/2 卡起（前一轮立方体弱扩展即从 2 卡起）。
44. **`rank.sh` 的卡映射已验证正确，不要再怀疑“跑错卡”（2026-11-15 复验）**：
    `scripts/rank.sh` 设 `HIP_VISIBLE_DEVICES=$((n/RPC))` 后 `unset ROCR_VISIBLE_DEVICES`。
    按 ROCm 文档语义这看似会退化成“全机物理卡下标”，但**HIP 运行时实测**
    （`evidence/probe/128108949`，`hipGetDeviceCount` + `hipDeviceGetPCIBusId`）表明：
    ROCR 清掉后 HIP 下标仍限于分配到的设备，下标越界直接报
    `no ROCm-capable device is detected`，**不会静默跑到别人的卡上**。
    同批探测也确认 `rocminfo` **不理会** `ROCR_*`/`HIP_VISIBLE_DEVICES`
    （五种视图都列出全部物理卡），所以**不能用 rocminfo 判断设备视图**；
    要判断就用 `hipGetDeviceCount`+`hipDeviceGetPCIBusId`（探针源码
    `scripts/gpu-view.cpp`）。另注意 `SLURM_JOB_GPUS` 是物理卡号、
    `ROCR_VISIBLE_DEVICES` 是重编号后的视图，两者数值不相等是正常的，不是缺陷。

## 6. 本文维护规则

- **不新建 dated handoff**：每次会话结束就地更新 §4（任务指针）和 §5（新坑）。
- §1/§2 只在环境或固定路径变化时改，并同步 `runbook` 与 本文 §1。
- 能力边界、性能数字、Phase A–F 进度不写进本文，避免与 todo 双份维护。

45. **bash `local` 同一语句内不能引用刚声明的变量（2026-11-15，一天内踩了四次）**：
    在 `set -u` 下，
    ```bash
    local tag="$1" dir="$EV/case-$tag"     # ✗ tag: unbound variable
    ```
    会失败——bash 把同一条 `local` 里的多个赋值视为**同时求值**，`$tag` 在展开时
    还未定义。必须分行：
    ```bash
    local tag="$1"
    local dir="$EV/case-$tag"               # ✓
    ```
    互不引用的同行赋值（`local a=$1 b=$2`）是安全的。
    **后果不只是报错**：如果这段代码在计时之后才执行，作业会以 FAILED 收场
    但**所有测量数据已经落盘且有效**——不要因为作业状态是 FAILED 就丢掉数据，
    先看 `summary.tsv` 是否完整；也不要因为数据有效就忽略这个 bug。
    排查顺序：`grep -n 'local .*=.*\$' script.sh` 找出所有同行引用。

46. **作业脚本里 `mpirun ./rank.sh` 要求 rank.sh 与 cwd 同目录（2026-11-15，踩了两次）**：
    `rank.sh` 在 `$SCRIPTS` 下，但 `mpirun` 的相对路径 `./rank.sh` 是相对**当前工作
    目录**解析的。做法照抄别处脚本时容易漏掉这行：
    ```bash
    EV="$ONEFLOW_ROOT/evidence/probe/${SLURM_JOB_ID:-manual}/<name>"
    mkdir -p "$EV"
    cp "$SCRIPTS/rank.sh" "$EV/rank.sh"; chmod +x "$EV/rank.sh"   # ← 必须
    # ... 之后 ( cd "$EV" && mpirun ... ./rank.sh ... )
    ```
    漏了会得到
    `prterun was unable to launch ... could not access an executable: ./rank.sh`，
    而且**只在 mpirun 那一步失败**，前面的检查全部通过，容易误判为"环境问题"。

47. **多个 agent 共享工作树时不要用 `git add -A`（2026-11-16，自己犯的）**：
    并行子 agent 都在同一个 checkout 里写文件，`git add -A` 会把**别人未提交的
    改动**扫进自己的提交。实际后果：一次 SA 数学优化的提交里混进了 156 行的
    头文件拆分重构，提交信息完全没提，后来者用 `git log` 查不到那次重构。
    规则：
    - 提交前先 `git status --short`，逐个文件确认归属；
    - 只 `git add <具体路径>`，不用 `-A` / `.`；
    - 如果发现已经混进去了且**已推送**，不要重写历史，改用说明性提交澄清，
      并把这次的归属记在提交信息里。
