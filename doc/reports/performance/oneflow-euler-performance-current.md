# OneFLOW 一维 Euler 回归与 CPU/DCU 性能阶段报告

> 报告日期：2026-09-02（修订版）；测试平台：昆山 Z100 / native HIP / `gfx906`。
> 本报告只记录当前工作树和本轮实际作业证据；不把 CUDA、Kokkos、4 卡 MPI
> 或完整 Navier–Stokes 求解器描述为已验证。

## 1. 执行摘要

本轮先完成高精度 correctness，再进行性能比较。单卡 HIP 端口、32 核 CPU MPI
和 4-rank/4-DCU HIP MPI 均在昆山实际计算节点完成。

- CPU/HIP FullTrace：smooth-periodic、Sod-transmissive 均通过；逐项误差和 ULP 为 0。
- 单卡 NoTrace 四规模：最终状态误差为 0，CPU/HIP checksum 一致。
- Stateful HIP 四规模：单卡 DCU 完整生命周期相对 32 核 CPU MPI 为 3.85–9.70×。
- GPU-only 持续测试：4M 网格、100 steps、20 repeats，12,000 次 kernel launch，HIP kernel 占 87.7% 左右。
- 32-rank CPU MPI：4 个规模的 1-rank/32-rank global hash 全部一致。
- 4-rank/4-DCU HIP MPI：严格 CPU/HIP 对比通过，4 个规模均通过，程序报告 visible_devices=4。
- Trace 路径仍受 D2H 主导，端到端相对串行 CPU 只有 2.04–2.55×；它是诊断路径，不是生产性能路径。

按行业常用口径，应用级 GPU 加速应优先使用同等精度、同一算例和同一迭代量下的端到端 wall-clock；kernel-only 或相对单 CPU 核的结果只能作为辅助诊断。以此标准看，当前 1-DCU stateful 相对 32 核 CPU 为 3.85–9.70×，4-DCU/4-rank MPI 在 4M 规模达到 13.10×；小规模 4 卡结果受 MPI/设备启动和 halo 固定开销影响，不能单独代表吞吐。

## 2. 回归与残差标准

### 2.1 统一回归对象

一维 Euler 推进三个守恒量 `[rho, rho*u, rho*E]`，使用一阶 Rusanov 数值通量和
SSP-RK3。CPU 是 reference backend，HIP 通过同一套 backend-neutral 接口执行。

回归同时检查 face state、numerical flux、residual、RK state 和 primitive
`rho/u/p`，并记录 absolute error、scaled-relative error、ULP、有限值和正状态。

### 2.2 残差标准

```text
abs(error) <= 1e-15 + 1e-15 * max(1, abs(reference))
```

该规则适用于普通高精度路径。Sod 激波附近若 CPU 的 `rho` 或 `p` 相邻跳变超过
2%，该点从逐点 smooth-residual 比较中屏蔽并单独计数；守恒量、物理状态和最终
状态仍然检查。

### 2.3 普通与严格模式的统一方式

普通和严格回归共用现有测试架构、同一 residual 输出路径和同一 backend 接口；差异
只在序列化精度及 tolerance 配置。当前昆山验证使用严格的 `1e-15` 判据，未维护
第二套回归体系。

## 3. 资源与测试口径

| 测试 | 分区 | 资源 | 用途 |
| --- | --- | --- | --- |
| HIP correctness / trace / stateful | `kshdnormal` | 1 node, 1 task, 8 CPU, `dcu:1`, 27G | 单卡 HIP/DCU 实际验证 |
| GPU-only + dmon | `kshdnormal` | 1 node, 1 task, 8 CPU, `dcu:1`, 27G | 只运行 HIP stateful backend，并独立采样设备 |
| CPU MPI reference/regression | kshcnormal | 1 node, 32 ranks, 32 CPU, 3569 MB/CPU | 1-rank/32-rank hash 一致性与性能对照 |
| HIP MPI regression/benchmark | kshdnormal | 1 node, 4 ranks × 8 CPU, dcu:4, 3569 MB/CPU | 4-rank/4-DCU，一卡一进程映射 |

DCU 队列按昆山 profile 的固定配比申请 dcu:4；4 个 rank 共 32 个 CPU，内存按
3569 MB/CPU 请求。CPU 对照使用 HPC-X 2.4.1 GCC 7.3.1，HIP 使用 DTK 26.04
Clang 17 和 CMake 3.25.0。

## 4. 架构整合

```text
CreateState(problem)
Upload(state, hostState)
Advance(state, steps, options)
Download(state, hostState)
```

- CPU state 持有主存状态；HIP state 持有 device buffers、stream 和 event。
- `FullTrace` 用于逐 RK stage 诊断；`NoTrace` 用于多步生产式推进。
- HIP state 创建时分配一次，初始状态 H2D 一次，最终结果 D2H 一次。
- 旧 `Step(..., EulerTrace&)` 保留为 compatibility/correctness wrapper，不作为性能主路径。

本轮新增项目文件：

- `ports/kunshan/oneflow_1d_hip/EulerGpuOnlyBenchmark.cpp`
- `ci/kunshan/euler-dcu-smoke.slurm`
- `ci/kunshan/euler-dcu-matrix.slurm`
- `ci/kunshan/euler-dcu-gpu-only.slurm`
- `ci/kunshan/euler-cpu-mpi-matrix.slurm`
- ports/kunshan/oneflow_1d_hip/OneDEulerMpi.h/.cpp/.hip：CPU/HIP MPI 生命周期与节点内设备映射。
- ports/kunshan/oneflow_1d_hip/EulerMpiRegression.cpp：同一 MPI 分解下 CPU/HIP 严格结果回归。
- ci/kunshan/euler-cpu-mpi-regression.slurm：1-rank/32-rank hash 回归。
- ci/kunshan/euler-dcu-mpi4-regression.slurm：4-rank/4-DCU 严格回归与性能矩阵。

首次 DCU 构建失败是 DTK HSA 依赖回退到系统 `system ROCm fallback path`，与 DTK HIP 头文件混用
导致的重复定义；脚本已显式指定 DTK 的 HSA CMake 配置、头文件和运行库路径后通过。
Euler 数值代码未因该问题改写。

## 5. 精度门槛

| 层次 | 案例/规模 | 结果 | 误差证据 |
| --- | --- | --- | --- |
| FullTrace | smooth-periodic | PASS | face、flux、residual、state、primitive 均为 0 |
| FullTrace | Sod-transmissive | PASS | smooth residual 通过；激波 mask 202 点 |
| NoTrace | 65,536 / 262,144 / 1,048,576 / 4,194,304 | PASS | 四规模 final max abs error = 0，checksum 一致 |
| GPU-only | 4,194,304，100 steps，20 repeats | PASS | 独立 HIP stateful 运行；设备为 Device 66a1/gfx906 |
| CPU MPI regression | 65,536 / 262,144 / 1,048,576 / 4,194,304 | PASS | 1-rank 与 32-rank final hash 全部一致，rho/p 为正 |
| HIP MPI regression | 4 ranks × 4 DCU，65,536，20 steps | PASS | max_scaled_error=0、violations=0、physical=1、CPU/HIP hash 一致 |

所有性能数据均在上述精度门槛通过后才纳入本报告。

## 6. 性能对比

### 6.0 行业口径与本报告的指标分层

行业 benchmark 通常不采用单一“GPU 加速平均值”，而是先固定比较边界：相同算例、网格规模、精度、时间步/迭代数和收敛条件，再报告端到端执行的 wall-clock。对于 CPU 多核，wall-clock 比累计 CPU time 更适合作为跨设备主指标；累计 CPU time 会随 CPU 并行度增加，不能直接与 GPU 时间等价。

| 端到端加速区间 | 工程上通常意味着 | 本报告中的对应关系 |
| --- | --- | --- |
| 1–3× | 小规模、CPU 基线较强、双精度、频繁 CPU/GPU 数据搬运或仅部分 offload。 | 当前 HIP Trace 的 2.04–2.55×；属于诊断路径结果。 |
| 3–10× | 较成熟且有工程价值的单 GPU/单加速器应用级结果，通常要求主要计算驻留设备侧。 | 当前 1 张 DCU 相对 32 核 CPU MPI 的 3.85–9.70×；属于主指标。 |
| 10–30× | 通常需要更高计算强度、更大问题规模、GPU-native 数据流和较少同步/搬运。 | 是阶段 D 优化后的合理进取目标，不是当前已达成结果。 |
| >30× | 常见于 kernel-only、GPU 对单 CPU 核、或 CPU/GPU 资源口径不对等；不能自动解释为整体应用加速。 | 当前 157.7–293.5× 仅属于 Stateful compute-path 微基准。 |

公开可比案例也说明了口径的重要性：MFC 报告单张 A100 相对完整 CPU socket 约 7×，但相对单个 CPU 核约 300×；Euler 求解器研究则报告了局部/求解过程约 10×以上的加速。LAMMPS 官方文档也明确指出，实际倍数受问题规模、算法、精度和数据传输影响。因此这些公开数字用于建立经验边界，不构成 OneFLOW 的直接横向排名。

[MFC GPU benchmark](https://arxiv.org/abs/2305.09163) · [Euler solver wall-clock example](https://onlinelibrary.wiley.com/doi/10.1155/2017/4610138) · [LAMMPS GPU performance guidance](https://docs.lammps.org/Speed_gpu.html)


### 6.1 32 核 CPU MPI 与单卡 HIP Stateful

**主性能指标：端到端 wall-clock。**两边均为 100 steps、2 repeats、1 warmup。HIP 数值为完整生命周期 `create + upload + advance + download`；CPU 数值为 32-rank MPI 作业生命周期。该表是当前最适合做资源级比较和对外沟通的口径，结果为 3.85–9.70×，落在上述 3–10× 工程区间内。

| nx | 32-rank CPU lifecycle (ms) | 1-DCU HIP lifecycle (ms) | CPU/HIP |
| ---: | ---: | ---: | ---: |
| 65,536 | 71.366 | 18.547 | 3.85× |
| 262,144 | 281.295 | 52.716 | 5.34× |
| 1,048,576 | 1,284.895 | 188.891 | 6.80× |
| 4,194,304 | 7,082.042 | 730.417 | 9.70× |

### 6.2 串行 CPU 与 HIP Trace

**辅助端到端指标：诊断路径。**Trace benchmark 使用 50 steps、2 repeats、1 warmup；CPU 是串行 reference。该表用于定位诊断路径开销，不与 32 核 CPU 表混用，也不能替代 6.1 的应用级资源对比。HIP copy ratio 随规模升高至 93.2%，解释了为什么 trace 路径不能代表生产式 GPU-resident 性能。

| nx | 串行 CPU (ms) | HIP Trace (ms) | speedup | HIP copy ratio |
| ---: | ---: | ---: | ---: | ---: |
| 65,536 | 1,156.722 | 568.178 | 2.04× | 75.8% |
| 262,144 | 4,924.655 | 1,931.599 | 2.55× | 84.0% |
| 1,048,576 | 20,190.864 | 8,933.606 | 2.26× | 91.5% |
| 4,194,304 | 93,809.070 | 38,792.485 | 2.42× | 93.2% |

### 6.3 Stateful 相对串行 CPU

**辅助微基准：纯求解阶段。**这里的 Stateful HIP 是完整的设备驻留执行路径中的 `Advance` 段，即多步 HIP kernel 连续推进；它不是单独测一个 kernel，但排除了 create、H2D 和最终 D2H。该结果用于判断 GPU 计算内核潜力，不能写成整体应用加速。

| nx | HIP Stateful 相对串行 CPU | HIP syncs | final error |
| ---: | ---: | ---: | ---: |
| 65,536 | 157.7× | 2 | 0 |
| 262,144 | 209.1× | 2 | 0 |
| 1,048,576 | 226.5× | 2 | 0 |
| 4,194,304 | 293.5× | 2 | 0 |

### 6.4 GPU-only 持续测试

| 指标 | 实测值 |
| --- | --- |
| 规模 / steps / repeats / warmup | 4,194,304 / 100 / 20 / 2 |
| GPU-only total | 约 7,457 ms |
| GPU-only kernel time | 约 6,542 ms，占 87.7% |
| kernel launches / syncs | 12,000 / 20 |
| 设备 | Device 66a1，`gfx906:sramecc-:xnack-` |


### 6.5 32-rank CPU MPI 与 4-rank/4-DCU HIP MPI

两组数据都使用 100 steps、1 warmup、2 repeats 的 benchmark lifecycle。CPU 组同时执行
1-rank 和 32-rank，回归脚本要求两者 final_hash 完全一致；DCU 组使用 4 个 MPI rank
和 4 张可见 DCU，脚本要求 visible_devices=4，并在 65,536 网格上用
EulerMpiRegression 做 CPU/HIP 严格结果比较。

| nx | 32-rank CPU lifecycle (ms) | 4-DCU HIP MPI lifecycle (ms) | CPU/HIP wall-clock |
| ---: | ---: | ---: | ---: |
| 65,536 | 40.209 | 70.214 | 0.57× |
| 262,144 | 130.277 | 82.219 | 1.58× |
| 1,048,576 | 653.098 | 121.239 | 5.39× |
| 4,194,304 | 3,629.603 | 277.217 | 13.10× |

这一表是同一节点的资源级端到端 wall-clock 对比，不是单卡结果的线性外推。4 卡路径的
fixed MPI halo、rank/device 初始化和小问题规模固定成本在 65,536 上超过计算收益；
当局部问题规模扩大后，4 卡结果显示出明显吞吐优势。4 卡结果已通过功能回归，但尚未
覆盖跨节点 MPI 或其他 accelerator backend。

### 6.6 单卡与四卡对照

下表把单卡和四卡结果放到同一端到端 wall-clock 口径下。两组均使用 32-rank CPU
lifecycle 作为分母，HIP 侧均包含 create、upload、advance 和 download；四卡额外
包含 4-rank MPI halo 通信。

| nx | 32-rank CPU (ms) | 1-DCU HIP (ms) | 1-DCU 加速 | 4-DCU HIP MPI (ms) | 4-DCU 加速 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 65,536 | 71.366 | 18.547 | 3.85× | 70.214 | 0.57× |
| 262,144 | 281.295 | 52.716 | 5.34× | 82.219 | 1.58× |
| 1,048,576 | 1,284.895 | 188.891 | 6.80× | 121.239 | 5.39× |
| 4,194,304 | 7,082.042 | 730.417 | 9.70× | 277.217 | 13.10× |

结论：

- 单卡已经在完整生命周期口径下稳定优于 32-rank CPU，规模增大时从 3.85×提升到 9.70×；
- 四卡在 65,536 上为 0.57×，主要受 MPI、rank/device 初始化和 halo 固定开销影响；
- 从 1M 开始四卡超过单卡，4M 达到 13.10×，说明多卡域分解在较大问题规模上开始摊薄固定成本；
- 四卡数据证明的是单节点 4-rank/4-DCU HIP MPI 整体 wall-clock，不是单卡结果的线性外推，也不代表跨节点扩展。

## 7. DCU 证据与限制

- 配置、编译和运行均在实际昆山 DCU 计算节点完成，不使用登录节点结果冒充设备证据。
- 程序通过 HIP device properties 报告 Device 66a1/gfx906，HIP event 记录 kernel 时间。
- 4-rank/4-DCU 作业在计算节点报告 local_ranks=4、device_index=0（rank 0 日志）和 visible_devices=4；代码按节点内 local rank 映射可见设备，并在严格模式拒绝 rank 数超过可见设备数。
- GPU-only 作业单独启动 DTK dmon；采样输出显示 HCU VRAM 使用情况。该 DTK 版本的 dmon 输出不是 HCU compute utilization，因此不能把它写成完整利用率结论。
- CPU benchmark、HIP benchmark 和 GPU-only 采样窗口已分离；GPU-only 可执行程序本身不创建 CPU backend，但 host 仍负责初始化和最终结果回传。

未验证边界：

- 当前证据覆盖单节点 4 卡/4 rank；跨节点 MPI、跨节点 GPU 通信和多节点拓扑仍未验证。
- CUDA、Kokkos 和其他未在目标节点实测的 accelerator backend 均只保留接口，不宣称通过。
- 当前 benchmark 是隔离的一维 Euler 适配器，不等同于完整 OneFLOW NS 主求解链。

## 8. 与 Euler backend 优化计划对照

| 阶段 | 状态 | 本轮证据/剩余项 |
| --- | --- | --- |
| A：状态化统一接口 | 完成 | CPU/HIP 共用 Create/Upload/Advance/Download；FullTrace 兼容保留。 |
| B：持久 HIP state | 完成 | 四规模 NoTrace 通过；stateful 生命周期 D2H 低于 10%，allocation 与 steps 解耦。 |
| C：验证与诊断解耦 | 部分完成 | GPU-only 和独立 dmon 窗口完成；GPU reduction 仍未实现，dmon 当前主要给出 VRAM 采样。 |
| D：kernel/dataflow 优化 | 未开始 | 仍有 9 launches/step；trace 路径 D2H 最高约 93.2%。 |
| E：昆山验收 | 当前单节点范围完成 | 高精度 correctness、四规模 HIP 矩阵、32 核 CPU MPI hash 回归和 4-rank/4-DCU HIP MPI 均已完成；跨节点和其他 backend 未验证。

阶段 D 的性能目标不能用本轮 trace 的 2.04–2.55× 代替；当前主指标是 stateful 完整生命周期相对 32 核 CPU 的 3.85–9.70×，这与单卡 GPU 应用常见的 3–10× 工程区间一致。157.7–293.5× 仅用于定位设备侧计算潜力，仍不代表多节点或完整生产 NS 性能。

## 9. 下一步与交付状态

1. 将本轮 MPI 代码、统一回归脚本和报告作为一个阶段性变更块审查并 commit。
2. 后续补阶段 C 的 GPU reduction：将 checksum、最大误差、finite/positive-state 检查优先收敛为 GPU scalar reduction。
3. 再进入阶段 D：移除 NoTrace 不需要的 trace buffers，评估 Face/Residual/Update 融合、launch 数和访存效率；每个优化块统一重跑 FullTrace、NoTrace 四规模、32-rank CPU 和 4-DCU MPI 回归。
4. 当前证据尚未覆盖 CUDA、Kokkos、跨节点 MPI；这些边界继续保持未宣称验证。

可复现实验入口：

```text
ci/kunshan/euler-dcu-smoke.slurm
ci/kunshan/euler-dcu-matrix.slurm
ci/kunshan/euler-dcu-gpu-only.slurm
ci/kunshan/euler-cpu-mpi-matrix.slurm
ci/kunshan/euler-cpu-mpi-regression.slurm
ci/kunshan/euler-dcu-mpi4-regression.slurm
```
