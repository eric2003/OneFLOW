# OneFLOW 一维 Euler 原生 HIP 验证记录

> 记录日期：2026-08-28
> 本文是当前一维主线的结果记录，不代表完整 Navier–Stokes GPU 化，也不替代
> `doc/reports/architecture/oneflow-gpu-backend-delivery.md` 中较早的 GPU backend 框架交付报告。

## 1. 当前结论

- 一维验证程序已经进入本地 `master` 的当前 HEAD。
- CPU 主路径保留，HIP 只作为隔离的一维验证后端；当前不把它描述为完整生产 CFD
  GPU backend。
- 郑州 BW1000 的 `gfx936` 结果已归档，仅作为历史参考，不再进入后续代码测试。
- 昆山 Z100 的目标架构是 `gfx906`，自动探测、编译、运行和逐步比对闭环已在真实计算节点通过。
- 正确性是硬门槛：CPU/HIP stage-by-stage 回归和多规模 benchmark 最终状态误差均为 0。
- 已建立 kernel、H2D、D2H、allocation 分项计时；高性能 backend 的下一步重点是数据常驻显存和减少同步/D2H。
- 已完成两类并行验收：32-rank CPU MPI 四规模 hash 回归，以及单节点 4-rank/4-DCU HIP MPI 严格回归与四规模 benchmark；后者在 4M 达到 13.10x（相对 32-rank CPU lifecycle）。

## 2. 主线范围与数值判据

一维 Euler 验证推进三个守恒量：

```text
[rho, rho*u, rho*E]
```

实现为一阶 Rusanov 数值通量和 SSP-RK3 时间推进。HIP 与 CPU 在每个 RK stage
比较：

- face state；
- numerical flux；
- residual；
- 三个守恒量 `rho`、`rho*u`、`rho*E`；
- 派生原始量 `rho`、`u`、`p`。

每个比较项同时记录：

- absolute error；
- scaled-relative error，尺度为 `max(1, abs(reference))`；
- ULP distance；
- finite/physical-state 检查。

接受条件为：

```text
abs(error) <= 1e-15 + 1e-15 * max(1, abs(reference))
```

这表示数值误差阈值，不等同于“小数点后 15 位”或“15 位有效数字”。
Sod 激波附近不把逐点残差当作平滑区参考：CPU 的 `rho` 或 `p` 相邻跳变超过
2% 的比较点进入 residual mask，并单独报告；其它平滑区域仍按同一误差判据比较。

## 3. 集群边界

| 集群 | 入口/分区 | 硬件与目标架构 | 本记录中的角色 |
| --- | --- | --- | --- |
| 郑州 BW1000 | `zzeshell` | BW1000，实测 `gfx936` | 已归档参考 |
| 昆山 Z100 | `kseshell` / `kshdnormal` | Z100，目标 `gfx906` | 当前主线，已完成一维闭环 |

两者不能互换。郑州作业号不能证明昆山 Z100/HIP 通过。

## 4. 郑州 BW1000 归档参考（不进入当前代码测试）

原始运行目录（已归档前路径）：

```text
CLUSTER_PROJECT_ROOT/OneFLOW_1d_euler_20260825
CLUSTER_PROJECT_ROOT/OneFLOW_1d_euler_auto_20260825
```

本节结果已从郑州远端实验目录移动到以下归档位置，仅保留作历史参考；后续代码测试、回归矩阵和主线结论不再使用郑州结果：

```text
CLUSTER_PROJECT_ROOT/archive/OneFLOW_1d_euler_zhengzhou_bw1000_20260825
CLUSTER_PROJECT_ROOT/archive/OneFLOW_1d_euler_zhengzhou_bw1000_auto_20260825
```

相关作业：

| 作业 | 结果 | 说明 |
| --- | --- | --- |
| `JOB_ID_OMITTED` | PASS | CPU-only 一维验证 |
| `JOB_ID_OMITTED` | FAIL | 按错误的 `gfx906` 编译，运行时报 `invalid device function` |
| `JOB_ID_OMITTED` | PASS（探针） | 节点报告 `gfx936` |
| `JOB_ID_OMITTED` | PASS | 显式 `gfx936` 编译运行 |
| `JOB_ID_OMITTED` | PASS | 未显式传架构，自动探测 `gfx936` 后 Euler 回归通过 |
| `JOB_ID_OMITTED` | PASS（探针） | 另一节点同样报告 `gfx936` |
| `JOB_ID_OMITTED` | 待排队记录 | 郑州作业，不纳入昆山主线；今天不再处理 |

通过内容：

- 平滑周期案例：PASS；
- Sod 可传递边界案例：PASS；
- 三个守恒量和派生量均通过；
- absolute error、scaled-relative error、ULP 均为 0；
- Sod 激波附近 mask 排除 202 个逐点 residual 比较点。

这些结果只证明 BW1000/`gfx936` 路径，不证明 Z100/`gfx906`。

## 5. 昆山 Z100 当前证据

独立运行目录：

```text
CLUSTER_PROJECT_ROOT/OneFLOW_1d_euler_kunshan_auto_20260825
```

相关作业：

| 作业 | 结果 | 失败/限制 |
| --- | --- | --- |
| `JOB_ID_OMITTED` | FAIL | `CXX` 残留指向不存在的 `CLUSTER_TOOLCHAIN_ROOT/gcc-11.2.0-install/bin/g++` |
| `JOB_ID_OMITTED` | FAIL | 清除 `CC/CXX` 后，CMake 在自动架构探测阶段失败；`rocminfo` 未得到架构 |
| `JOB_ID_OMITTED` | FAIL（探针） | 节点 `compute-node-omitted` 输出 `C-3000 module is NOT loaded, possibly no HCU devices` |
| `JOB_ID_OMITTED` | PASS（设备探针） | 节点 `compute-node-omitted`，DTK 26.04/HIP 6.3.26113，`rocminfo=gfx906`，HCU/KFD 正常 |
| `JOB_ID_OMITTED` | PASS | 显式 `gfx906` 构建；CPU self-test 和 HIP Euler 回归通过 |
| `JOB_ID_OMITTED` | PASS | 不传架构，自动探测得到 `CMAKE_HIP_ARCHITECTURES=gfx906`；Euler 回归通过 |
| `JOB_ID_OMITTED` | PASS | 修复编译器 warning 后重新自动探测和回归；节点 `compute-node-omitted` |

`JOB_ID_OMITTED` 的最终验证内容：

- CMake 自动得到 `CMAKE_HIP_ARCHITECTURES=gfx906`；
- C++/HIP 配置和构建退出码均为 0；
- CPU Euler self-test：PASS；
- HIP Euler：平滑周期和 Sod 可传递边界均 PASS；
- 三个守恒量、face state、flux、smooth residual、最终 state 和 primitive
  `rho/u/p` 的 absolute、scaled-relative 和 ULP 均为 0；
- Sod 激波附近 residual mask 排除 202 个比较点；
- 构建日志不再出现此前的 `void function is missing a return statement` 和
  `hipFree` 返回值未使用警告。

昆山远端项目是源码快照，没有 `.git`，因此不能提供远端 commit hash。作为替代，
已对照本地 HEAD 和远端快照的关键文件 SHA-256，以下文件完全一致：

```text
cmake/UtilFunction.cmake
ports/kunshan/oneflow_1d_hip/CMakeLists.txt
ports/kunshan/oneflow_1d_hip/OneDEuler.cpp
ports/kunshan/oneflow_1d_hip/OneDEuler.hip
ports/kunshan/oneflow_1d_hip/OneDWeno5.hip
```

## 6. 本地提交证据

merge 前的本地仓库状态：

```text
branch: master
HEAD:   f51a184aad4f83a89e7d3328e18189155da0cef8
origin/master: c97a8c9579fb47747606aea4b9d8c78742dda38d
关系：本地 master 比 origin/master 超前 8 个提交
说明：该状态对应本地 merge commit `f51a184a`
```

一维与加速相关提交（均已在本地 HEAD）：

| 时间（Asia/Shanghai） | commit | 内容 |
| --- | --- | --- |
| 2026-08-25 09:11:45 | `7462558300a40405299c136f25635c5a74d8a03c` | 原生 HIP 标量 backend |
| 2026-08-25 10:32:47 | `5e0895c82b5a56bf7af1966f2f562a1cf5eceffa` | 隔离一维 HIP upwind port |
| 2026-08-25 11:36:15 | `58ef394bc4ba18b33a638d2894fab62b6b0c6f44` | 一维可压 Euler HIP 验证 |
| 2026-08-25 12:31:25 | `aa083377798ab32d010d808c58bffd1cfbfa55b7` | 一维 Lax-WENO5 HIP replay |
| 2026-08-25 12:45:07 | `ca3f312050db5cb6cba39055dc3cf6cf3e7b4766` | WENO5 持久 device state |
| 2026-08-25 13:39:54 | `a3d72b814da3cc3c036aa2c44ed5ee027556a471` | 跨平台 backend 自动探测 |
| 2026-08-25 16:27:01 | `157cb193b85defc1a8d2e5f385b51b4f2cc4da28` | 一维验证自动 HIP 架构探测 |
| 2026-08-25 17:54:37 | `f51a184aad4f83a89e7d3328e18189155da0cef8` | 合并上游 `origin/master`，无冲突 |
| 2026-08-26 16:22:10 | `71e86282492d1b1967da81d3dc305886d355e56a` | 昆山 Z100/gfx906 自动探测 Euler 回归和 warning 修复 |

Windows 原始 WENO3 示例目录未被这些提交覆盖或替换；一维验证放在
`ports/kunshan/oneflow_1d_hip/` 隔离目录中。

## 7. HIP backend 与多架构预留

当前保持以下边界：

- CPU backend 是默认、完整保留的主路径；
- HIP backend 代码保留为一维验证和后续移植的接口位置，不宣称已完成完整
  Navier–Stokes；
- CUDA、Kokkos、多 GPU 只保留扩展入口，当前不引入 Kokkos；
- 不写死 `gfx906` 或 `gfx936`。

架构选择顺序由 `cmake/UtilFunction.cmake` 提供：

1. 显式 `CMAKE_HIP_ARCHITECTURES`；
2. 环境变量 `ONEFLOW_HIP_ARCHITECTURES`（逗号或分号列表）；
3. 目标节点 `rocminfo` 报告的可见架构；
4. `rocm_agent_enumerator` 只有唯一候选时的保守 fallback；
5. 无法探测时显式失败，不猜默认架构。

因此未来可以直接传入多架构列表，例如：

```text
CMAKE_HIP_ARCHITECTURES=gfx906;gfx936
```

## 8. 2026-08-25 历史收尾状态

- [x] 郑州 BW1000 实验目录已移动到 archive，原始内容保留；
- [x] 郑州结果降级为历史参考，不进入后续代码测试和回归矩阵；
- [x] 本地 `master` 已 merge 上游 `origin/master`，检查无冲突；
- [x] 未提交新的 Slurm 作业，未取消已有作业，未 push；
- [x] 昆山 Z100/gfx906 设备探针、显式架构和自动架构一维 Euler 回归均通过；
- [x] 修复一维 Euler 编译 warning，并在昆山自动探测路径重新回归；
- [x] 本地验证 commit `71e86282` 已创建，工作区随后保持干净；
- [ ] 不修改 Windows 原始 WENO3 路线，不改动 `codes/ns`。

## 9. 历史恢复清单（性能基线见第 10 节）

- [ ] 在 `kseshell` 对应的 Z100 节点确认 C-3000/HCU 模块已加载；
- [ ] 在独立项目目录记录远端副本 `HEAD`、分支和工作区状态；
- [ ] 让自动探测得到 `gfx906`，或显式传入 `CMAKE_HIP_ARCHITECTURES=gfx906`；
- [ ] 先运行 CPU 一维 Euler 基准；
- [ ] 再运行 HIP 一维 Euler，逐 stage 比较三守恒量、primitive、flux、residual；
- [ ] 保存作业脚本、日志、退出码和架构探针输出；
- [ ] 只有昆山回归通过后，才讨论下一次 commit。

## 10. 2026-08-28 昆山性能基线

在昆山 Z100 `gfx906` 真实计算节点上，使用与 CPU 相同的一维 Euler 周期边界
workload（50 steps，2 repeats，1 warmup）完成四个规模测试。benchmark 同时
记录 HIP event kernel 时间、H2D/D2H、显存 allocation、kernel launch/allocation
次数、最终状态误差和 checksum：

| nx | CPU total (ms) | HIP total (ms) | speedup | kernel (ms / ratio) | H2D (ms) | D2H (ms) | alloc (ms) | final max abs error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 65536 | 1130.592 | 572.878 | 1.974x | 14.369 / 2.51% | 38.063 | 394.595 | 72.728 | 0 |
| 262144 | 4668.764 | 1889.750 | 2.471x | 37.438 / 1.98% | 134.687 | 1489.551 | 91.599 | 0 |
| 1048576 | 18967.298 | 9038.656 | 2.098x | 142.828 / 1.58% | 611.474 | 7671.997 | 121.409 | 0 |
| 4194304 | 88670.222 | 37329.346 | 2.375x | 583.860 / 1.56% | 2533.934 | 32077.289 | 237.917 | 0 |

四个规模 CPU/HIP checksum 均一致。设备属性报告为 `Device 66a1`、`gfx906`，
每个 benchmark 都执行 900 次 kernel launch 和 700 次 allocation。`rocm-smi`
采样在目标 HCU index 0 看到作业进程及其显存占用；HCU 采样平均约 0.10--2.33%，
峰值 1--55%，显存利用率峰值 0--5%。由于采样周期约 1 秒而 kernel 很短，利用率
采样用于确认设备活动，kernel event 才是 kernel 时间依据；低平均利用率不能单独
推出 CPU fallback。

当前性能结论不是“GPU 没有工作”，而是：GPU kernel 已实际执行，端到端已有约
1.97--2.47x 加速，但 D2H 主导总时间（大规模约 86--88%），算力利用率因此偏低。
下一步按“正确性回归 → per-zone 持久 device state → 减少 D2H/分配 → 重新验证
正确性 → 重新测性能”推进；二维、三维暂缓。

## 11. 当前优先级

当前不先推进 WENO5。先以一维 Euler 为唯一场景完成：

1. stateful CPU/HIP backend 生命周期；
2. per-zone 持久 device state；
3. `NoTrace` 性能模式和 `FullTrace` 验证模式解耦；
4. 昆山 correctness 复测；
5. 四规模端到端、kernel、传输、allocation 和利用率复测；
6. 根据新 profile 再做 Euler kernel/dataflow 优化。

完整阶段、接口建议和验收指标见
`doc/reports/architecture/oneflow-euler-optimization-plan.md`。Euler 达到验收
条件后再恢复 WENO5；二维、三维继续后移。


## 12. 历史：32 核 CPU 与 8 rank 共享一张卡的 MPI 对比

本节保留早期共享单卡实验，作为历史性能对照；当前生产式多卡验收见第 13 节。

### 12.1 测试口径

为避免把前面的单线程 stateful 结果误认为生产级多核加速，本轮新增两个独立
Slurm 任务，两个任务都运行同一个 MPI Euler benchmark：

| 任务 | MPI 进程与资源 | 分区 | 设备 |
| --- | --- | --- | --- |
| CPU MPI | 32 ranks × 1 CPU core | `kshcnormal` | 32 核 Hygon C86 7185 |
| GPU MPI | 8 ranks × 1 CPU core | `kshdnormal` | 8 个 MPI rank 共享 1 张 Z100，`gfx906` |

两个任务使用相同的全局网格、周期边界、100 steps、2 repeats 和 1 warmup。全局
网格按 rank 做一维块划分，每个 RK stage 交换左右 3 个守恒量 halo。GPU 任务中
每个 rank 将 device 0 设为可见设备，因此这是“8 个 MPI rank 共享一张卡”的指定
对比，不是“一 rank 一卡”的生产配置。

这组测试与第 10 节旧 trace/stateful benchmark 的区别是：

- 第 10 节旧路径是单进程、逐 stage trace 路径；
- 第 10 节 stateful 路径是单进程、一张卡一个进程、不含 MPI halo；
- 本节是两个独立 MPI 任务，重点观察多进程域分解、halo 通信和共享 GPU 的综合成本。

### 12.2 实测结果

表中时间是每个 repeat 的最大 rank wall time 累加；每步耗时按 200 个计时步归一化。
`CPU/GPU` 定义为 32-rank CPU lifecycle 除以 8-rank/1-DCU GPU lifecycle，大于 1
表示 GPU 任务更快。

| global nx | CPU MPI total (ms) | CPU ms/step | GPU MPI total (ms) | GPU ms/step | CPU/GPU | GPU kernel (ms) | GPU MPI (ms) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 65,536 | 71.729 | 0.358646 | 451.152 | 2.255759 | 0.159x | 97.067 | 170.866 |
| 262,144 | 357.378 | 1.786888 | 456.412 | 2.282060 | 0.783x | 99.054 | 168.879 |
| 1,048,576 | 1,293.912 | 6.469560 | 501.767 | 2.508834 | 2.579x | 130.438 | 161.298 |
| 4,194,304 | 7,083.452 | 35.417262 | 976.969 | 4.884843 | 7.250x | 288.790 | 419.965 |

正确性和设备证据：

- 四个规模 CPU/GPU `final_hash` 完全一致；
- `min_rho` 约为 `0.899978/0.899995/0.899999/0.900000`；
- `min_pressure` 约为 `0.919978/0.919995/0.919999/0.920000`；
- 两个任务构建、运行和退出码均为 0；
- GPU 输出设备为 `Device 66a1`，架构为 `gfx906:sramecc-:xnack-`。

### 12.3 结果分析

#### 小规模：GPU 计算量不足以摊平 MPI 和共享卡固定成本

`nx=65,536` 时，32-rank CPU 只需 71.7 ms，而 GPU 任务需 451.2 ms；
`nx=262,144` 时，GPU 仍略慢。此时每个 rank 的本地网格较小，GPU kernel 本身
只有约 97--99 ms，但 GPU 任务还包含：

- 8 个 MPI rank 各自创建 HIP context 和 device buffers；
- 每个 RK stage 通过 host buffer 完成边界 D2H、MPI Sendrecv 和 halo H2D；
- 每个 stage 需要 device synchronization，避免 host MPI 读取未完成的数据；
- 8 个进程在同一张卡上竞争 kernel launch 和执行资源。

因此这些规模测到的是固定通信/同步成本，而不是 GPU 的吞吐能力。

#### 中等规模：出现交叉点

`nx=1,048,576` 时 GPU 任务为 501.8 ms，CPU 任务为 1,293.9 ms，达到 2.58x。
此时每个 rank 的本地计算量足以摊平共享单卡的固定成本，GPU kernel 约 130.4 ms，
MPI/halo 路径约 161.3 ms，GPU 仍有明显收益。

#### 大规模：GPU 胜出，但 MPI/halo 已成为重要瓶颈

`nx=4,194,304` 时 GPU 任务达到 7.25x，但 GPU MPI 路径仍不是纯 kernel 性能：

- kernel event 约 288.8 ms；
- MPI/halo 相关 host 路径约 420.0 ms；
- 每个 rank 每个 RK stage 执行 `PackBoundary + Face + Residual + Update`；
- 两个 repeat 合计 19,200 次 kernel launch、4,800 次 halo exchange、9,632 次
  device synchronization。

因此当前 4M 加速比已经受到共享卡和 host-device halo 往返限制。若未来改成一张卡
一个 MPI rank，并把 halo 交换改成 device-aware MPI 或 GPU-aware communication，
理论上可以减少这部分开销；这是待验证方向，不能直接把本节数据外推为该配置收益。

### 12.4 与单进程 stateful 结果的关系

本节 8-rank/1-DCU 结果低于第 10 节单进程 stateful 结果是预期现象，二者不是同一
资源口径：

- 单进程 stateful `Advance` 在 4M 上约为 1.15 s / 200 steps 的纯推进累计，
  kernel event 占 `99.99%`；
- 8-rank MPI 共享卡在 4M 上约为 0.98 s / 200 steps，但其每个 rank 只负责
  1/8 网格，且包含 8 个 context、MPI halo 和大量同步；
- MPI 版本的总 wall time 更短，说明域分解确实提高了可处理吞吐，但单位资源成本
  和通信效率仍需单独评估；不能简单将两个数字相除当作线性扩展效率。

更准确的判断是：

```text
单进程 stateful：验证 GPU 连续推进和 kernel 吞吐
MPI 8-rank/1-DCU：验证共享卡域分解下的综合 wall time
32-rank CPU：提供多进程 CPU 参考
```

### 12.5 当前结论和下一步

本节结果把之前“单线程 CPU 对 GPU”的结论修正为更接近并行运行的资源对比：

- 不能再用 `102--155x` 作为 CPU/GPU 生产性能结论；它只对应串行 CPU reference；
- 在 32-rank CPU vs 8-rank/1-DCU 共享卡口径下，GPU 的有效领先区间是从约 1M
  网格开始，4M 为 7.25x；
- 小规模不适合直接上 8 个 MPI rank 共享一张卡；应考虑更少 GPU rank 或一 rank
  一卡的进程映射；
- 当前 MPI GPU 路径的首要优化对象不是继续堆 kernel，而是减少 host-device halo
  往返、device synchronization 和共享卡多 context 争用；
- 任何后续优化仍需同时通过 CPU/HIP correctness、最终 hash/物理状态检查和四规模
  wall-time 复测。

完整原始输出、运行参数和可复现脚本见：

```text
doc/reports/performance/oneflow-euler-performance-20260828.md
ci/kunshan/euler-mpi-cpu-32.slurm
ci/kunshan/euler-mpi-gpu8-1.slurm
```

## 13. 当前：32-rank CPU 与 4-rank/4-DCU HIP MPI

本轮统一 MPI 验收使用单节点、32-rank CPU 和 4-rank/4-DCU HIP 两个资源口径。CPU
脚本要求 1-rank 与 32-rank 的 final hash 一致；HIP 脚本在 4 个可见 DCU 上运行
4 个 rank，并先执行严格 CPU/HIP 结果比较。

严格回归（global nx=65,536，20 steps）：

- max_scaled_error=0；
- violations=0；
- physical=1；
- CPU/HIP global hash 一致；
- 作业报告 visible_devices=4，rank/device 映射按节点内 local rank 选择设备。

端到端 benchmark（100 steps、2 repeats、1 warmup）：

| global nx | 32-rank CPU lifecycle (ms) | 4-DCU HIP MPI lifecycle (ms) | CPU/HIP wall-clock |
| ---: | ---: | ---: | ---: |
| 65,536 | 40.209 | 70.214 | 0.57x |
| 262,144 | 130.277 | 82.219 | 1.58x |
| 1,048,576 | 653.098 | 121.239 | 5.39x |
| 4,194,304 | 3,629.603 | 277.217 | 13.10x |

4 卡结果是单节点资源级端到端比较，不是单卡线性外推。小规模受 MPI、设备初始化
和 halo 固定开销影响；大规模才显示出明显吞吐优势。跨节点 MPI、CUDA、Kokkos
和完整 Navier–Stokes 主线仍未验证。

可复现入口：

- ci/kunshan/euler-cpu-mpi-regression.slurm
- ci/kunshan/euler-dcu-mpi4-regression.slurm
