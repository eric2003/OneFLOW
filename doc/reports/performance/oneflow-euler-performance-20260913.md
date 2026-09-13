# OneFLOW 一维 Euler CPU/DCU 性能对照报告（2026-09-13）

> 报告日期：2026-09-13
> 测试平台：昆山 Z100 / native HIP / `gfx906`
> 代码版本：`fix/contract-test-cmake-path`（commit `0c5f3d88`，含 PR #149 的 CMake 路径修复）
> 本报告记录当日全部实测数据；所有作业均在昆山真实计算节点完成，原始日志归档于集群侧
> （作业号见文末）。

## 1. 执行摘要

本轮在同一代码版本上完成四组对照测量，全部作业退出码为 0：

| 对照组 | 作业 | 说明 |
|---|---|---|
| 32-rank CPU MPI | `121984941` | `kshcnormal`，HPC-X 2.4.1 / GCC 7.3.1 |
| 1-rank / 1-DCU HIP MPI | `121985687` | `kshdnormal`，DTK 26.04，共享设备模式 |
| 4-rank / 4-DCU HIP MPI | `121985687` | 同上，4 卡并行 |
| 单线程 CPU vs 单卡 DCU（非 MPI） | `121984118` | stateful benchmark，CPU 单线程 reference |

关键结论：

- 正确性：CPU 5 算例 normal（`1e-8`）与 strict（`1e-15`）全部通过；HIP
  contract test 6/6；CPU MPI 与 HIP MPI 的全局守恒量 hash 在四个规模上完全一致。
- 单卡 DCU（非 MPI 路径）相对 32-rank CPU MPI 为 **3.75×–9.56×**。
- 4-rank/4-DCU HIP MPI 相对 32-rank CPU MPI 为 **0.95×–25.43×**；4M 规模
  **25.43×**，比 2026-08-28 基线的 13.10× 提升约 2×。
- 4 卡相对 1 卡 MPI 的并行扩展：大网格 4M 为 **3.58×**（4 卡），接近线性。

## 2. 测试方法

### 2.1 代码版本

- 分支：`fix/contract-test-cmake-path`，commit `0c5f3d88`（upstream/master `080aa893` + 2 个修复提交）。
- 二进制由集群侧从该 revision 的源码包（`git archive`）构建，构建目录与运行目录隔离。

### 2.2 Slurm 资源配置

| 作业 | 分区 | 节点 | 资源 | 时长 |
|---|---|---|---|---|
| CPU 回归 + port CPU test | `kshcnormal` | `a15r1n17` | 1 node, 16 CPU, 54G | 1:00:00 |
| 单卡 DCU contract + benchmark | `kshdnormal` | `i01r2n10` | 1 node, 1 task, 8 CPU, `dcu:1`, 27G | 0:50:00 |
| 32-rank CPU MPI | `kshcnormal` | `g07r2n12` | 1 node, 32 ranks × 1 CPU, 3569 MB/CPU | 0:40:00 |
| 1/4-rank DCU MPI | `kshdnormal` | `b15r1n17` | 1 node, 4 ranks × 8 CPU, 3569 MB/CPU, `dcu:4` | 0:50:00 |

资源配比遵循 `kseshell` profile 的固定 CPU/内存/加速器比例。

### 2.3 工具链

| 角色 | 组件 | 版本 / 路径 |
|---|---|---|
| CPU 编译 | GCC | 9.3.0（module `compiler/gcc/9.3.0`） |
| CPU MPI | OpenMPI | 4.1.5（module `mpi/openmpi/gcc-9.3.0/4.1.5`） |
| CPU MPI（对照作业） | HPC-X | 2.4.1 / GCC 7.3.1（module `mpi/hpcx/2.4.1-gcc-7.3.1`） |
| CMake | CMake | 3.25.0 |
| Python | CPython | 3.8.10 |
| METIS | METIS | 5.0.1（集群侧从源码编译） |
| CGNS | CGNS | 4.2.0（`/public/software/mathlib/CGNS-4.2.0/src`） |
| DCU 编译 | DTK clang | 17.0.0（module `compiler/dtk/26.04`） |
| DCU 目标架构 | `gfx906` | `Device 66a1` |

CPU 侧依赖通过 `MPI_HOME_*`、`METIS_HOME_*`、`CGNS_HOME_*` 环境变量注入顶层
CMake；OpenMPI 4.x 需要额外 `-DCMAKE_CXX_FLAGS="-DOMPI_SKIP_MPICXX"`。
完整环境要求与坑位见 [`ci/kunshan/README.md`](../../../ci/kunshan/README.md)。

### 2.4 构建产物与命令

四个独立构建，互不复用（HIP 与 CPU 使用不同工具链）：

1. 完整 solver（CPU 回归用）：
   ```bash
   cmake -S OneFLOW -B build-solver -DCMAKE_BUILD_TYPE=Release \
     -DMPI_ENABLE=ON -DMETIS_ENABLE=ON -DCGNS_ENABLE=ON
   cmake --build build-solver --parallel 16
   ```
2. Kunshan port CPU contract test：
   ```bash
   cmake -S OneFLOW/ports/kunshan/oneflow_1d_hip -B build-port-cpu \
     -DONEFLOW_1D_ENABLE_GTEST=ON -DONEFLOW_1D_ENABLE_HIP=OFF
   ```
3. Kunshan port DCU（合同测试 + stateful benchmark）：
   `-DONEFLOW_1D_ENABLE_HIP=ON -DONEFLOW_1D_ENABLE_GTEST=ON`，DTK 工具链，
   `-DCMAKE_HIP_ARCHITECTURES=gfx906`。
4. Kunshan port DCU + MPI（对照基准）：
   `-DONEFLOW_1D_ENABLE_HIP=ON -DONEFLOW_1D_ENABLE_MPI_BENCHMARK=ON`，
   `CC=mpicc`、`CXX=mpicxx`（HPC-X 2.4.1）、`CMAKE_HIP_COMPILER=clang++`；
   `CMAKE_PREFIX_PATH` 需以**环境变量**方式导出（`amd_comgr`/`AMDDeviceLibs`
   由 HIP 语言检测的 try_compile 子项目查找）。

### 2.5 Benchmark 程序与参数

两类基准程序，均对 CPU 与 DCU 使用同一 backend-neutral 接口和同一初值：

- `oneflow_1d_euler_stateful_benchmark`（非 MPI）：`nx steps repeats warmup`，
  本报告统一 `100` steps、`2` repeats、`1` warmup。输出完整生命周期
  （create + upload + advance + download）与 CPU/HIP 最终状态最大绝对误差。
- `oneflow_1d_euler_mpi_*_benchmark`（MPI）：`global_nx steps repeats warmup`，
  同上参数。输出各 rank 最大值（lifecycle/create/upload/advance/download、
  MPI、compute、kernel）、halo 次数、设备同步次数和全局守恒量 hash。

规模：`nx = 65536 / 262144 / 1048576 / 4194304`。

### 2.6 通过判据

- 作业负载退出码为 0；
- HIP 作业额外校验 `visible_devices` 与 `local_ranks`（4-rank 作业要求
  `local_ranks=4` 且每个 rank 绑定独立 `device_index`）；
- CPU 侧运行 CPU 回归（`test/test.py` + canonical residual baseline，normal
  `1e-8` 与 strict `1e-15` 两个 profile）；
- MPI 侧校验 CPU/HIP 全局 hash 一致与物理量下界（`min_rho`、`min_pressure` > 0）。

## 3. 正确性结果

### 3.1 CPU 回归（作业 `121982636`）

- normal（`1e-8`）：5/5 通过，五个案例最大绝对残差差 4.9e-11。
- strict（`1e-15`，`ONEFLOW_RESIDUAL_TEST_OUTPUT=1`）：5/5 通过，最大绝对差
  1.1e-17。
- 覆盖案例：`plateuns2dslau2`、`plateuns2dslau2_34950_35000`、
  `turbplateuns2droe_sa`、`rae2822_roe_sa`、`m6wingroe_sa`。

### 3.2 DCU contract test（作业 `121984118`）

HIP GoogleTest 6/6、CTest `HIP.*` 6/6 通过（含设备可见性、FullTrace/NoTrace
CPU 参考对比、生命周期复用、非法请求处理）。

### 3.3 CPU/HIP 全局 hash（作业 `121984941` / `121985687`）

| nx | 32-rank CPU MPI | 1-rank DCU MPI | 4-rank DCU MPI |
|---:|---|---|---|
| 65,536 | `0xf5b7bc809f8660da` | 一致 | 一致 |
| 262,144 | `0x37cd60c26338f9de` | 一致 | 一致 |
| 1,048,576 | `0xc84eeaee66ae9395` | 一致 | 一致 |
| 4,194,304 | `0x2143c776ad2531c9` | 一致 | 一致 |

`min_rho` 与 `min_pressure` 在所有运行中均为正（0.8999…/0.9199…）。

## 4. 性能结果

所有数值为单次生命周期对 `2` repeats 求和的 wall-clock（毫秒），`100` steps，
1 warmup；MPI 行取各 rank 最大值。

### 4.1 四组对照（生命周期，ms）

| nx | 32-rank CPU MPI | 1-DCU HIP（非 MPI） | 1-DCU HIP MPI | 4-DCU HIP MPI |
|---:|---:|---:|---:|---:|
| 65,536 | 69.44 | 18.52 | 72.72 | 72.89 |
| 262,144 | 265.41 | 53.19 | 115.77 | 81.25 |
| 1,048,576 | 1,334.22 | 188.83 | 271.69 | 119.73 |
| 4,194,304 | 6,985.94 | 730.79 | 982.35 | 274.68 |

### 4.2 加速比（相对 32-rank CPU MPI）

| nx | 1-DCU HIP（非 MPI） | 1-DCU HIP MPI | 4-DCU HIP MPI | 4-DCU / 1-DCU 扩展 |
|---:|---:|---:|---:|---:|
| 65,536 | 3.75× | 0.95× | 0.95× | 1.00× |
| 262,144 | 4.99× | 2.29× | 3.27× | 1.42× |
| 1,048,576 | 7.07× | 4.91× | 11.14× | 2.27× |
| 4,194,304 | 9.56× | 7.11× | **25.43×** | **3.58×** |

### 4.3 分项（4,194,304 规模）

| 指标 | 32-rank CPU MPI | 1-rank DCU MPI | 4-rank DCU MPI |
|---|---:|---:|---:|
| lifecycle max (ms) | 6,985.94 | 982.35 | 274.68 |
| advance max (ms) | 6,960.82 | 780.85 | 248.90 |
| MPI 通信 max (ms) | 1,321.88 | 0.54 | 7.55 |
| compute max (ms) | 6,363.19 | 780.31 | 245.51 |
| kernel max (ms) | —（CPU） | 725.21 | 192.27 |
| halo exchanges | —（未记录） | 600 | 2,400 |
| device syncs | —（未记录） | 1,204 | 4,816 |

### 4.4 与 2026-08-28 基线对比（相对 32-rank CPU MPI）

| nx | 1-DCU 历史 | 1-DCU 本次 | 4-DCU 历史 | 4-DCU 本次 |
|---:|---:|---:|---:|---:|
| 65,536 | 3.85× | 3.75× | 0.57× | 0.95× |
| 262,144 | 5.34× | 4.99× | 1.58× | 3.27× |
| 1,048,576 | 6.80× | 7.07× | 5.39× | 11.14× |
| 4,194,304 | 9.70× | 9.56× | 13.10× | **25.43×** |

单卡结果与历史一致（±10% 以内）；4 卡结果整体约为历史值的 1.7–2.0×，
大网格提升最明显。历史值与本次值的差异来源（节点状态、驱动或测量窗口）
未在本文档中归因。

## 5. 分析

1. **1-rank MPI 与 4-rank MPI 的固定开销**：小规模（65,536）下 1-rank MPI
   生命周期 72.72 ms 与 4-rank 的 72.89 ms 基本持平——halo 交换与设备同步的
   固定开销占主导，卡数增加不带来收益；该规模不适合做多卡对照。
2. **多卡扩展**：262,144→4,194,304 的 4 卡扩展效率从 1.42× 提升到 3.58×，
   接近线性；4M 时 kernel 时间为每卡 192 ms（CPU 侧 compute 6.36 s），说明
   计算已完全驻留设备。
3. **MPI 与非 MPI 路径**：单卡非 MPI stateful 路径（4M lifecycle 730.79 ms）
   快于 1-rank MPI 路径（982.35 ms），差额来自每步 halo 交换与额外同步
   （MPI 路径 100 步产生 600 次 halo、1,204 次设备同步，kernel 启动次数为
   2,400 次，非 MPI 路径为 1,200 次）。生产使用单卡时应优先非 MPI 路径。
4. **相对 32 核 CPU 的工程价值区间**：本次 4-DCU 在 1M/4M 规模分别达到
   11.14×/25.43×，落在“GPU-native 数据流、较少同步”对应的区间上沿。

## 6. 复现要点

见 [`ci/kunshan/README.md`](../../../ci/kunshan/README.md) 的
“Verified CPU regression environment”与“Verified DCU/HIP environment”两节，
要点包括：GCC 9.3.0+（GCC 7 无法编译当前代码）、`-DOMPI_SKIP_MPICXX`、
METIS 需自编译、HIP 构建需以环境变量导出 `CMAKE_PREFIX_PATH`、以及设备探测
失败时排除故障节点后重提。

## 7. 原始数据位置

集群侧（`/public/home/luql/.scnet-hpc/oneflow-regression-20260913/`）：

- `artifacts/121982636/`：CPU 回归、port CPU contract test、环境记录；
- `artifacts-dcu/121984118/`：DCU contract test、单卡 stateful benchmark；
- `artifacts-cpu-mpi/121984941/`：32-rank CPU MPI benchmark；
- `artifacts-dcu-mpi/121985687/`：1-rank/4-rank DCU MPI benchmark。
