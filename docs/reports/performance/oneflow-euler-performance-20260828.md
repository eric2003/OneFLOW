# OneFLOW 一维 Euler CPU/HIP 性能对比测试报告

> 测试日期：2026-08-28
> 测试平台：昆山 Z100 / native HIP / `gfx906`
> 测试范围：一维 Euler、双精度、单进程、单线程 CPU 对单张 DCU
> 报告性质：当前 backend 架构优化的阶段性实测，不是完整 OneFLOW
> Navier–Stokes 生产性能结论。

## 1. 执行摘要

本轮优化没有修改一维 Euler 的数值格式，而是重构 CPU/HIP backend 的状态和执行
生命周期，把原来每步发生的显存分配、H2D、逐 RK stage 同步和完整 trace D2H，
改为：

```text
创建持久 state
    -> 初始状态 H2D 一次
    -> 多个时间步和全部 RK stage 留在 GPU
    -> 最终状态 D2H 一次
```

主要结果如下：

- CPU/HIP 的 smooth-periodic 和 Sod-transmissive 逐 stage `FullTrace`
  正确性回归全部通过，测试项误差和 ULP 均为 0；
- 四个规模的 `NoTrace` 最终状态最大绝对误差均为 0，CPU/HIP checksum 一致；
- HIP `Advance` 相对单线程 CPU `Advance` 的加速为 `102.76--155.22x`；
- 将 HIP state 创建、H2D、`Advance` 和最终 D2H 全部计入后，以 CPU
  `Advance` 为分母的保守加速为 `82.69--144.12x`；
- 新路径每步归一化的 HIP 完整生命周期耗时，相比旧 trace 路径降低
  `43.17--60.13x`；
- 新路径 kernel event 时间占 `Advance` 墙钟时间 `99.56--99.99%`；
- 完整 HIP 生命周期中的最终 D2H 占比为 `4.36--6.11%`；
- 每个时间步仍然执行 9 个 kernel，说明本轮收益主要来自数据生命周期和同步结构，
  还不是 kernel fusion 或算子级优化。

本轮结果证明：

1. GPU kernel 确实在目标设备上执行；
2. 旧路径的低加速主要来自 trace 数据回传、反复分配和同步；
3. stateful backend 已解决第一层架构瓶颈；
4. 目前不能把 `102--155x` 描述为“单卡对四核 CPU”的生产加速，因为 CPU
   benchmark 是串行代码，Slurm 申请的 4 个 CPU 仅用于作业与编译资源。

## 2. 测试对象

### 2.1 数值问题

测试推进三个守恒量：

```text
[rho, rho*u, rho*E]
```

数值方法：

- 空间离散：一阶 Rusanov 数值通量；
- 时间推进：SSP-RK3；
- 精度：`double`；
- 编译约束：`-ffp-contract=off`；
- 性能案例：smooth-periodic；
- `gamma = 1.4`；
- `dx = 1 / nx`；
- `dt = 0.05 / nx`。

测试规模：

```text
nx = 65,536
nx = 262,144
nx = 1,048,576
nx = 4,194,304
```

### 2.2 两条 HIP 路径

| 路径 | 用途 | 每步资源行为 |
| --- | --- | --- |
| 旧 `Step + FullTrace` | 正确性和优化前性能基线 | 每步分配、2 次 H2D、每个 RK stage 完整 trace D2H |
| 新 stateful `NoTrace` | 性能路径 | state 创建时分配一次，初始 H2D 一次，最终 D2H 一次 |

旧路径不是错误实现。它保留完整的 face state、flux、residual 和 RK state，
适合逐 stage 正确性诊断，但不应作为生产性能路径。

### 2.3 测试迭代参数

| 路径 | steps | repeats | warmup | 计时步总数 |
| --- | ---: | ---: | ---: | ---: |
| 旧 trace benchmark | 50 | 2 | 1 | 100 |
| 新 stateful benchmark | 100 | 2 | 1 | 200 |

由于两条路径的 `steps` 不同，跨路径比较统一使用“每时间步耗时”。各路径内部的
原始累计值仍单独保留。

## 3. 硬件、软件和资源口径

### 3.1 计算节点环境

| 项目 | 实测配置 |
| --- | --- |
| 集群 | 昆山 |
| 分区 | `kshdnormal` |
| CPU | Hygon C86 7185，节点 32 核 |
| 加速器 | 1 张 Z100，设备属性 `Device 66a1` |
| HIP 架构 | `gfx906` |
| 单卡显存 | 约 16 GiB |
| DTK | `compiler/dtk/26.04` |
| C++/HIP 编译器 | DTK Clang 17.0.0 |
| CMake | 3.22.0-rc1 |
| 构建类型 | Release |

### 3.2 Slurm 资源

构建/正确性和性能作业均申请：

```text
nodes=1
ntasks=1
cpus-per-task=4
gres=dcu:1
mem=13G
```

必须区分“申请资源”和“程序实际并行度”：

- CPU benchmark 没有 OpenMP、MPI 或线程池，实际为一个串行 CPU 线程；
- HIP benchmark 使用一张可见 DCU；
- 其余已申请 CPU 可用于编译、作业控制和 GPU host 侧运行时；
- 因此当前对比口径是“单线程 CPU 对单张 DCU”，不是“4 核 CPU 对单卡”。

## 4. 本轮代码改动

### 4.1 统一生命周期接口

新增 backend-neutral 接口：

```cpp
CreateState(problem);
Upload(state, hostState);
Advance(state, steps, options);
Download(state, hostState);
```

CPU 和 HIP 使用相同的调用语义：

- CPU state 持有主存状态；
- HIP state 持有 device buffers、stream 和 profiling event；
- 旧 `Step(..., EulerTrace&)` 保留为 compatibility/correctness wrapper。

### 4.2 `NoTrace` 和 `FullTrace` 分离

`EulerRunMode::FullTrace`：

- 只允许一个时间步；
- 返回每个 RK stage 的 face state、flux、residual 和 state；
- 用于 CPU/HIP 逐项正确性检查；
- 允许逐 stage 同步和 D2H。

`EulerRunMode::NoTrace`：

- 支持一次 `Advance` 多个时间步；
- 不返回逐 stage trace；
- 中间状态持续保留在 device；
- 只在显式 `Download` 时回传最终状态。

### 4.3 持久 HIP state

HIP state 当前持有：

```text
base
work
scratch
left
right
flux
residual
stream
kernel timing events
```

这些资源在 `CreateState` 中建立，在多个时间步中复用，析构时释放。RK stage
通过 device pointer 轮换，不再每步重新创建整套 device buffers。

### 4.4 新增性能与诊断程序

| 目标 | 作用 |
| --- | --- |
| `oneflow_1d_euler_cpu` | CPU self-test |
| `oneflow_1d_euler_hip` | CPU/HIP `FullTrace` 正确性回归 |
| `oneflow_1d_euler_benchmark` | 旧 trace 路径的 kernel/H2D/D2H/allocation 基线 |
| `oneflow_1d_euler_stateful_benchmark` | 新生命周期的 CPU/HIP 性能和最终状态检查 |

主要新增或修改文件：

```text
OneDEulerBackend.h/.cpp
OneDEulerPersistent.hip
OneDEulerProfile.h/.hip
EulerBenchmark.cpp
EulerStatefulBenchmark.cpp
EulerMain.cpp
CMakeLists.txt
```

### 4.5 未修改的内容

本轮没有修改：

- Rusanov 数值通量公式；
- SSP-RK3 系数；
- 周期/可传递边界定义；
- 正确性 tolerance；
- WENO5；
- 二维、三维和完整 NS 主求解链；
- Windows 原始 WENO3 示例。

## 5. 正确性结果

### 5.1 `FullTrace` 逐 stage 回归

两个案例各执行 10 个时间步：

| 案例 | 边界 | 结果 |
| --- | --- | --- |
| smooth-periodic | Periodic | PASS |
| Sod-transmissive | Transmissive | PASS |

逐 stage 比较内容：

- 三个守恒量；
- left/right face state；
- numerical flux；
- smooth-region residual；
- primitive `rho/u/p`；
- finite、正密度和正压力。

接受条件：

```text
abs(error) <= 1e-15 + 1e-15 * max(1, abs(reference))
```

所有纳入比较的 absolute error、scaled-relative error 和 ULP distance 均为 0。
Sod 激波附近按既有 2% jump mask 排除 202 个逐点 residual 比较；这些位置仍由
state、flux 和物理状态检查覆盖。

### 5.2 `NoTrace` 最终状态

四个性能规模均得到：

```text
final_max_abs_error = 0
cpu_checksum == hip_checksum
```

这证明多时间步留在 GPU 后，最终状态与同一 benchmark 中的 CPU reference
一致。它不能替代 `FullTrace`，两种验证共同构成当前正确性证据。

## 6. 优化前：旧 trace 路径实测

原始累计结果如下。每个规模包含 2 次、每次 50 steps 的计时结果：

| nx | CPU total (ms) | HIP total (ms) | speedup | kernel (ms / HIP ratio) | H2D (ms) | D2H (ms) | allocation (ms) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 65,536 | 1,130.592 | 572.878 | 1.974x | 14.369 / 2.51% | 38.063 | 394.595 | 72.728 |
| 262,144 | 4,668.764 | 1,889.750 | 2.471x | 37.438 / 1.98% | 134.687 | 1,489.551 | 91.599 |
| 1,048,576 | 18,967.298 | 9,038.656 | 2.098x | 142.828 / 1.58% | 611.474 | 7,671.997 | 121.409 |
| 4,194,304 | 88,670.222 | 37,329.346 | 2.375x | 583.860 / 1.56% | 2,533.934 | 32,077.289 | 237.917 |

每个规模的计时窗口均包含：

- 900 次 kernel launch；
- 700 次 device allocation；
- 每个时间步 9 次 kernel；
- 每个时间步 7 次 allocation。

旧路径的主要瓶颈不是 kernel：

- 4M 规模 D2H 占 HIP 总时间约 85.93%；
- kernel 仅占约 1.56%；
- 完整 trace 在每个 RK stage 回传 left/right、flux、residual 和 state；
- 因此低平均 HCU 利用率与大量传输/同步相符。

## 7. 优化后：stateful `NoTrace` 实测

### 7.1 原始累计结果

每个规模包含 2 次、每次 100 steps 的计时结果：

| nx | CPU Advance (ms) | HIP Create (ms) | H2D (ms) | HIP Advance (ms) | kernel event (ms) | D2H (ms) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 65,536 | 2,194.558 | 1.537 | 2.022 | 21.357 | 21.262 | 1.623 |
| 262,144 | 9,370.527 | 1.984 | 3.262 | 68.325 | 68.207 | 3.726 |
| 1,048,576 | 38,820.703 | 2.710 | 7.267 | 279.068 | 278.916 | 13.182 |
| 4,194,304 | 178,937.295 | 4.916 | 21.784 | 1,152.799 | 1,152.721 | 62.052 |

### 7.2 两种加速比口径

`steady speedup`：

```text
CPU Advance / HIP Advance
```

`conservative lifecycle speedup`：

```text
CPU Advance / (HIP Create + H2D + HIP Advance + D2H)
```

第二个口径把 HIP state 创建和首尾传输计入，但 CPU 侧仍只有 `Advance`，因为当前
benchmark 没有输出 CPU create/upload/download 分项。因此它是比 steady speedup
更保守的辅助口径，不是严格对称的端到端应用时间。

| nx | HIP lifecycle total (ms) | steady speedup | conservative lifecycle speedup | kernel / Advance | D2H / lifecycle |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 65,536 | 26.539 | 102.76x | 82.69x | 99.56% | 6.11% |
| 262,144 | 77.296 | 137.15x | 121.23x | 99.83% | 4.82% |
| 1,048,576 | 302.227 | 139.11x | 128.45x | 99.95% | 4.36% |
| 4,194,304 | 1,241.552 | 155.22x | 144.12x | 99.99% | 5.00% |

### 7.3 kernel、launch 和同步

每个时间步仍然是：

```text
3 RK stages × (Face + Residual + Update) = 9 kernel launches
```

每个规模在两个计时 repeat 中共执行：

```text
100 steps × 2 repeats × 9 = 1,800 kernel launches
```

`EulerRunStats` 在 `Advance` 范围内记录 2 次同步，即每个计时 repeat 的 kernel
event 最终同步一次。需要注意：

- 这个数字只统计 `Advance` 中由 profiling event 产生的同步；
- `Upload` 为保证初始数据就绪，每个 repeat 还有一次 stream synchronize；
- `Download` 当前使用阻塞 D2H；
- 因此不能把 `hip_syncs=2` 解读为整个程序只有两次所有类型的同步。

## 8. 归一化前后对比

由于旧路径计时 100 个时间步、新路径计时 200 个时间步，本节按每时间步归一化。

| nx | 旧 CPU (ms/step) | 新 CPU (ms/step) | CPU 波动 | 旧 HIP total (ms/step) | 新 HIP Advance (ms/step) | 新 HIP lifecycle (ms/step) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 65,536 | 11.3059 | 10.9728 | -2.95% | 5.7288 | 0.1068 | 0.1327 |
| 262,144 | 46.6876 | 46.8526 | +0.35% | 18.8975 | 0.3416 | 0.3865 |
| 1,048,576 | 189.6730 | 194.1035 | +2.34% | 90.3866 | 1.3953 | 1.5111 |
| 4,194,304 | 886.7022 | 894.6865 | +0.90% | 373.2935 | 5.7640 | 6.2078 |

CPU 两轮每步耗时变化在约 `-2.95%--+2.34%`，说明两轮 workload 基本一致。

| nx | 旧 HIP / 新 HIP Advance | 旧 HIP / 新 HIP lifecycle | 旧/新 kernel 单步比 |
| ---: | ---: | ---: | ---: |
| 65,536 | 53.65x | 43.17x | 1.35x |
| 262,144 | 55.32x | 48.90x | 1.10x |
| 1,048,576 | 64.78x | 59.81x | 1.02x |
| 4,194,304 | 64.77x | 60.13x | 1.01x |

大规模下 kernel 本身每步只改善约 1%，而 HIP 完整路径改善约 60 倍。这是判断本轮
收益来源最重要的数据：优化收益主要来自消除重复分配、逐 stage D2H 和同步，
不是来自改变数值 kernel。

## 9. 数据传输与 allocation 分析

以下字节数根据代码中的 buffer 大小和调用次数推导，不是驱动 profiler 直接读取的
总线计数。状态大小为：

```text
3 × nx × sizeof(double)
```

| nx | 单份 state | 旧 H2D/计时窗口 | 旧 D2H/计时窗口 | 新 H2D/计时窗口 | 新 D2H/计时窗口 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 65,536 | 1.5 MiB | 300 MiB | 2.197 GiB | 3 MiB | 3 MiB |
| 262,144 | 6 MiB | 1.172 GiB | 8.789 GiB | 12 MiB | 12 MiB |
| 1,048,576 | 24 MiB | 4.688 GiB | 35.156 GiB | 48 MiB | 48 MiB |
| 4,194,304 | 96 MiB | 18.750 GiB | 140.625 GiB | 192 MiB | 192 MiB |

在各自 benchmark 的计时窗口内：

- H2D 理论字节量降低约 100 倍；
- D2H 理论字节量降低约 750 倍；
- 旧路径执行 700 次 device allocation；
- 新路径两个计时 repeat 各创建一次 state，按当前 7 个 device buffer 推导为
  14 次 device allocation；
- 新路径每个时间步的 allocation 为 0。

旧/新计时步总数并不相同，因此 `100x/750x` 描述的是两次实际 benchmark
计时窗口的总量缩减；若固定相同时间步数，stateful 路径优势还会随一次
`Advance` 中的 steps 增长。

## 10. GPU 实际执行证据

当前结论不是根据“程序链接了 HIP”或“显存有占用”得出，而是由多类证据共同支持：

1. CMake 在计算节点配置出 HIP 编译器和 `gfx906` 目标；
2. HIP kernel 由 DTK Clang 编译并成功运行；
3. HIP event 记录的 kernel 时间随 `nx` 近似线性增加；
4. 每个规模记录到符合代码结构的 kernel launch 数；
5. `rocm-smi` 在目标 HCU 上看到 benchmark 进程和随规模增长的显存占用；
6. CPU/HIP 最终状态和 checksum 一致；
7. `FullTrace` 下 GPU 各 stage 输出与 CPU reference 一致。

`rocm-smi` 的约 1 秒采样周期长于单批 kernel 时间，因此低平均利用率不能用来判断
CPU fallback。新 stateful `Advance` 中，event kernel 时间已经占墙钟时间
`99.56--99.99%`，说明该区间的主要工作确实发生在 device kernel。

## 11. 编译步骤

### 11.1 Slurm 脚本资源头

建议在昆山 DCU 分区的计算节点编译和运行：

```bash
#!/bin/bash -l
#SBATCH --partition=kshdnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=13G
#SBATCH --gres=dcu:1
#SBATCH --time=00:15:00
```

`kshdnormal` 要求至少申请一张 DCU。4 CPU 对应的 13 GiB 内存请求符合当前分区
每 CPU 内存限制。

### 11.2 环境

```bash
module purge
module load compiler/dtk/26.04

cmake_bin=CLUSTER_SOFTWARE_ROOT/compiler/cmake-3.22.0-rc1/bin/cmake
dtk_root=CLUSTER_SOFTWARE_ROOT/compiler/rocm/dtk-26.04
llvm_root="$dtk_root/llvm"

export PATH="$dtk_root/bin:$llvm_root/bin:$PATH"
export LD_LIBRARY_PATH="$dtk_root/hip/lib:$llvm_root/lib:$dtk_root/lib:$dtk_root/lib64:${LD_LIBRARY_PATH:-}"
export CC="$llvm_root/bin/clang"
export CXX="$llvm_root/bin/clang++"
```

作业脚本使用 login shell，确保 `module` 可用。不要在登录节点上把 HIP 编译成功
当成计算节点可运行证据。

### 11.3 CMake 配置

假设：

```bash
project_root=/path/to/OneFLOW
source_dir="$project_root/ports/kunshan/oneflow_1d_hip"
build_dir=/path/to/new-build-directory
```

自动架构探测配置：

```bash
unset ONEFLOW_HIP_ARCHITECTURES

"$cmake_bin" -S "$source_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_HIP_COMPILER="$CXX" \
  -DCMAKE_PREFIX_PATH="$dtk_root" \
  -DONEFLOW_1D_ENABLE_HIP=ON
```

核对配置：

```bash
grep -E 'CMAKE_HIP_ARCHITECTURES|CMAKE_(CXX|HIP)_COMPILER|CMAKE_BUILD_TYPE' \
  "$build_dir/CMakeCache.txt"
```

本轮实测缓存为：

```text
CMAKE_BUILD_TYPE=Release
CMAKE_HIP_ARCHITECTURES=gfx906
CMAKE_CXX_COMPILER=<DTK 26.04>/llvm/bin/clang++
CMAKE_HIP_COMPILER=<DTK 26.04>/llvm/bin/clang++
ONEFLOW_1D_ENABLE_HIP=ON
```

如果目标节点自动探测失败，应先确认设备和模块环境；确需显式配置时使用：

```bash
-DCMAKE_HIP_ARCHITECTURES=gfx906
```

不要把 `gfx906` 写成所有集群的全局默认。

### 11.4 编译目标

```bash
"$cmake_bin" --build "$build_dir" --parallel 4 --target \
  oneflow_1d_euler_cpu \
  oneflow_1d_euler_hip \
  oneflow_1d_euler_benchmark \
  oneflow_1d_euler_stateful_benchmark
```

本轮实测 CPU、HIP 和 stateful benchmark 均编译成功，构建与工作负载退出码为 0。

## 12. 运行步骤

### 12.1 先做 correctness

```bash
"$build_dir/oneflow_1d_euler_cpu"
"$build_dir/oneflow_1d_euler_hip"
```

验收：

```text
OneFLOW 1D Euler CPU self-test: PASS
smooth-periodic result: PASS
sod-transmissive result: PASS
OneFLOW 1D compressible Euler CPU/HIP validation: PASS
```

### 12.2 旧 trace 基线

```bash
for nx in 65536 262144 1048576 4194304; do
  "$build_dir/oneflow_1d_euler_benchmark" "$nx" 50 2 1
done
```

### 12.3 新 stateful 基线

```bash
for nx in 65536 262144 1048576 4194304; do
  "$build_dir/oneflow_1d_euler_stateful_benchmark" "$nx" 100 2 1
done
```

输出至少检查：

```text
cpu_advance_ms
hip_advance_ms
hip_kernel_ms
hip_kernel_ratio
hip_kernel_launches
hip_syncs
hip_create_ms
hip_upload_ms
hip_download_ms
hip_steady_speedup
final_max_abs_error
cpu_checksum
hip_checksum
```

### 12.4 GPU 活动采样

`rocm-smi` 采样应只覆盖 HIP 阶段，避免把前面的长时间 CPU benchmark 混入平均值。
建议记录：

```bash
rocm-smi --showuse
rocm-smi --showmemuse
rocm-smi --showpids
```

利用率采样用于确认设备活动和进程归属；kernel event 才是短 kernel 时间的主要
定量依据。

## 13. 当前结果的限制

### 13.1 旧 benchmark 的 CPU 对比不是生产级多核基线

第 6--9 节的旧 trace 和 stateful benchmark 使用串行 CPU reference：

- 单线程；
- 没有 SIMD、OpenMP 或 MPI 性能优化；
- `NoTrace` 内部仍借助 `EulerTrace` 完成 CPU step。

因此第 7 节的 `102--155x` 只说明 GPU stateful 执行结构相对串行 CPU 的收益，
不适合直接用于硬件采购、生产部署或“单卡等价多少 CPU 核”的结论。第 15 节
补充了真正的 MPI 32-rank CPU 对 8-rank/1-DCU 共享卡测试。

### 13.2 MPI 对比的适用边界

第 15 节是独立 MPI benchmark，不是完整 OneFLOW NS 主求解链：

- CPU 任务为 32 MPI ranks × 1 core；
- GPU 任务为 8 MPI ranks × 1 core，共享 1 张 Z100；
- GPU 每个 MPI rank 都建立自己的 HIP state，并将 device 0 设为可见设备；
- 这是“8 个 rank 共享单卡”的压力/争用口径，不是“一 rank 一卡”的生产配置；
- 每个 RK stage 都通过 host halo buffer 做 MPI 交换，GPU 任务包含额外的
  D2H/H2D halo 和同步成本。

因此第 15 节的结果用于回答指定资源口径下的对比，不应外推为多卡线性扩展结论。

### 13.3 不是完整 OneFLOW 求解器

当前只覆盖隔离的一维 Euler 验证程序，不包含：

- 完整 NS 主循环；
- 多 zone；
- 网格和边界数据管理；
- MPI 通信；
- IO、checkpoint 和残差归约；
- 二维、三维；
- WENO5 当前阶段优化。

### 13.4 profiling 还不完整

当前已经有总 kernel event 时间，但还缺少：

- Face/Residual/Update 分类时间；
- 有效带宽；
- occupancy；
- 寄存器使用；
- cache miss；
- H2D/D2H 实际 profiler 字节数；
- 完整生命周期的统一同步计数；
- CPU 和 HIP 严格分离的长时间利用率窗口。

### 13.5 当前 kernel 仍有明确结构开销

代码检查已经确认以下待测项，但本轮没有实施优化：

- 每步仍有 9 次 kernel launch；
- `NoTrace` 仍分配并写入 left/right trace buffers；
- `Face` kernel 的 grid 当前按 `3 × (nx + 1)` 值数量配置，但 kernel 线程按
  单个 face 工作，超出 `nx + 1` 的线程立即退出；
- Face、Residual、Update 尚未评估融合；
- 访存合并、寄存器压力和 occupancy 尚无 profiler 证据。

这些是下一阶段候选项，不应在本报告中提前宣称收益。

## 14. 阶段性结论

本轮已经完成的不是最终 kernel 优化，而是高性能 backend 的第一层架构修正：

```text
旧路径：GPU 是逐步调用的 trace 生成器
新路径：GPU 持有状态并连续推进多个时间步
```

从数据看：

- 正确性门槛保持不变；
- GPU kernel 已真实执行；
- 旧路径的主要瓶颈已经被定位并消除；
- 大规模 kernel 每步时间基本没有变化，但完整路径提高约 60 倍；
- 当前 `Advance` 已由 kernel 主导，下一步才适合进入 kernel/dataflow 优化。

在决定下一步之前，建议以本报告作为当前 Euler backend 的冻结基线。下一轮优化
必须继续按以下顺序验收：

```text
FullTrace correctness
    -> NoTrace final-state correctness
    -> 四规模性能
    -> kernel/传输/同步/利用率证据
```


## 15. 两个独立 MPI 任务：32 核 CPU vs 8 核 + 1 卡

### 15.1 任务定义

按用户指定的“两个任务提交”口径，在昆山分别提交了两个独立 Slurm 任务，使用
相同的 MPI Euler benchmark、全局网格、100 steps、2 repeats、1 warmup：

| 任务 | MPI 配置 | Slurm 分区 | 资源 | 说明 |
| --- | --- | --- | --- | --- |
| CPU MPI | 32 ranks × 1 core | `kshcnormal` | 1 node, 32 CPU cores, 110G | 纯 CPU MPI，周期 halo 交换 |
| GPU MPI | 8 ranks × 1 core | `kshdnormal` | 1 node, 8 CPU cores, 1 DCU, 27G | 8 ranks 共享 1 张 Z100，`gfx906` |

两项任务均使用 Release 构建。GPU 任务使用 DTK 26.04/Clang 17，并显式编译
`gfx906`；MPI 使用昆山可用的 HPC-X 2.4.1 GCC 7.3.1 模块。两个 benchmark
均按 `global_nx / ranks` 做一维块划分，每个 RK stage 进行左右 halo 交换。

### 15.2 性能结果

表中 `lifecycle` 是两个 repeat 的最大 rank wall time 累加；`per step` 按
200 个计时步归一化。`CPU/GPU` 是 32-rank CPU lifecycle 除以 8-rank/1-DCU
GPU lifecycle。

| global nx | CPU MPI lifecycle (ms) | CPU ms/step | GPU MPI lifecycle (ms) | GPU ms/step | CPU/GPU | GPU kernel (ms) | GPU MPI (ms) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 65,536 | 71.729 | 0.358646 | 451.152 | 2.255759 | 0.159x | 97.067 | 170.866 |
| 262,144 | 357.378 | 1.786888 | 456.412 | 2.282060 | 0.783x | 99.054 | 168.879 |
| 1,048,576 | 1,293.912 | 6.469560 | 501.767 | 2.508834 | 2.579x | 130.438 | 161.298 |
| 4,194,304 | 7,083.452 | 35.417262 | 976.969 | 4.884843 | 7.250x | 288.790 | 419.965 |

MPI 任务的最终状态证据：

- 四个规模的 CPU 和 GPU `final_hash` 完全一致；
- `min_rho` 分别约为 `0.899978/0.899995/0.899999/0.900000`；
- `min_pressure` 分别约为 `0.919978/0.919995/0.919999/0.920000`；
- 两个作业的构建、运行和退出码均为 0；
- GPU 输出设备为 `Device 66a1`，架构为 `gfx906:sramecc-:xnack-`。

### 15.3 结果解释

- 小规模时，32-rank CPU 更快：`nx=65,536` 时 GPU/CPU 约为 `6.29x`，
  `nx=262,144` 时约为 `1.28x`；此时 8 个 rank 共享单卡的 MPI、host-device
  halo 拷贝和同步成本高于 GPU 计算收益。
- `nx=1,048,576` 后 GPU 开始胜出，达到 `2.58x`；4M 时达到 `7.25x`。
- GPU 4M 任务中 kernel event 约 `288.8 ms`，MPI halo 约 `420.0 ms`，说明
  “8 ranks 共享一张卡”的主要瓶颈已从纯 kernel 转为 MPI/host-device halo 路径。
- GPU benchmark 每个 rank 每个 RK stage 执行 `PackBoundary + Face + Residual +
  Update`，因此 100 steps、2 repeats、8 ranks 共记录 `19,200` kernel launches，
  `4,800` 次 halo exchange，以及 `9,632` 次 device synchronization。
- 该结果不能和第 7 节的单 rank stateful `Advance` 速度直接混用：第 7 节不含
  MPI halo，且是一张卡一个进程；第 15 节专门测量用户指定的共享卡 MPI 配置。

### 15.4 MPI 编译与运行脚本

仓库新增两个可复现脚本：

```text
ci/kunshan/euler-mpi-cpu-32.slurm
ci/kunshan/euler-mpi-gpu8-1.slurm
```

CPU 任务关键配置：

```bash
#SBATCH --partition=kshcnormal
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=110G
mpirun -np 32 --bind-to core \
  oneflow_1d_euler_mpi_cpu_benchmark <nx> 100 2 1
```

GPU 任务关键配置：

```bash
#SBATCH --partition=kshdnormal
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=27G
#SBATCH --gres=dcu:1
export HIP_VISIBLE_DEVICES=0
export ROCR_VISIBLE_DEVICES=0
mpirun -np 8 --bind-to core \
  oneflow_1d_euler_mpi_hip_benchmark <nx> 100 2 1
```

GPU 任务需要先加载 `compiler/devtoolset/7.3.1`，再加载
`compiler/dtk/26.04` 和 `mpi/hpcx/2.4.1-gcc-7.3.1`，因为 HPC-X 模块声明了
该 compiler module 前置依赖。8 核任务申请 27G 而不是 28G，是因为昆山每核
内存上限约为 3569 MB。
