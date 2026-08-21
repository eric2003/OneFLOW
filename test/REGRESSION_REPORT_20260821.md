# OneFLOW 回归测试与代码改动报告

**报告日期：2026-08-21**
**验证平台：昆山 SCNet，CPU 分区 `kshcnormal`，GCC 7.3.1 + HPC-X MPI**
**远端隔离目录：** `<kunshan-dedicated-ci-root>/OneFLOW_backend_smoke_20260821`

## 1. 结论摘要

本轮工作的目标是为后续 GPU 移植建立 CPU 基线、搭建 backend 框架，并验证 MPI 空 rank 与真实 4 进程域分解。

结论如下：

- CPU 串行回归：3 个算例全部通过，合计 152 秒。
- 空 rank：原始单 zone 网格用 4 MPI rank 时曾触发 `ResMax::CalcMax()` 越界；修复后求解器可正常完成 50 步，退出码为 0。
- 正式 MPI 算例：M6 网格经 METIS 拆成 4 个 zone，每个 zone 73,728 cells；4 zone / 4 rank 连续运行两次均通过。
- 4 rank 两次运行的 `aero.dat`、`res.dat`、`turbres.dat`、`wallaero.dat` 完全一致，说明当前 CPU MPI 路径具有确定性。
- 当前 GPU backend 仍是框架级接入，未替换 CFD 数值 kernel；默认 backend 为 CPU。

## 2. CPU 串行回归基线

仓库定义见 [`REGRESSION.md`](REGRESSION.md) 和 [`regression_cpu.txt`](regression_cpu.txt)。

| 算例 | 覆盖内容 | 结果 | 耗时 |
|---|---|---:|---:|
| `plateuns2dslau2` | 2D 层流平板、SLAU2、黏性通量 | PASS | 3 s |
| `rae2822_roe_sa` | 2D 跨声速翼型、Roe、SA 湍流 | PASS | 10 s |
| `m6wingroe_sa` | 3D 翼型、Roe、SA 湍流 | PASS | 139 s |
| **合计** |  | **PASS** | **152 s** |

昆山串行回归作业（内部作业号，未纳入仓库）：
调度器状态：`COMPLETED`，退出码：`0:0`。

日志：

```text
logs/regression_<job-id>.out
CASE_RESULT=plateuns2dslau2 RC=0 ELAPSED_SECONDS=3
CASE_RESULT=rae2822_roe_sa RC=0 ELAPSED_SECONDS=10
CASE_RESULT=m6wingroe_sa RC=0 ELAPSED_SECONDS=139
REGRESSION_RC=0 TOTAL_ELAPSED_SECONDS=152
```

基准输出位于每个算例的 `autotest/` 目录，作为 CPU 数值参考。

## 3. 空 rank 问题与修复

### 原因

现有 CFD 测试网格都是单 zone。运行单 zone + 4 rank 时，只有 rank 0 获得 zone，其余 rank 的本地 zone 列表为空：

```cpp
pid[iZone] = (zoneStart + iZone) % Parallel::nProc;
```

旧实现无条件执行：

```cpp
ResData & t = dataList[0];
```

因此空 rank 在 `ResidualTask::Run()` → `ResMax::CalcMax()` 路径中发生越界段错误。

### 修复

修改文件：[`codes/residual/src/Residual.cpp`](../codes/residual/src/Residual.cpp:134)

- 空 `dataList` 使用安全哨兵值，不再访问 `dataList[0]`。
- 对本地最大残差执行全局 `MPI_MAX` 归约。
- 用 owner rank 同步最大值对应的 `zone/cell/坐标/体积` 元数据。
- 空全局网格时返回安全的 `-1/0` 元数据。
- 对全局 cell 数为 0 的平均残差路径增加保护，避免除零。

### 空 rank 验证

昆山 MPI 冒烟作业（内部作业号，未纳入仓库），`mpirun -np 4`，原始单 zone M6。
求解器实际退出码：`0`，完成 50 步；测试脚本当时报告失败，是因为复用了包含旧输出的目录，`aero.dat` 行数发生追加，属于测试目录清理问题，不是求解器失败。

## 4. 正式 4-zone / 4-rank MPI 算例

### 网格分区

使用仓库已有 METIS 分区工具完成分区（内部作业号未纳入仓库）。

原始 M6 网格：

- cells：294,912
- zones：1

分区结果：

- zones：4
- 每个 zone：73,728 cells
- 输出：

```text
<kunshan-dedicated-ci-root>/OneFLOW_backend_smoke_20260821/runs/m6wingroe_sa_partition_np4/grid/m6wing_np4.ofl
```

网格校验值：

```text
SHA256 0feacacb3ee49443d5aeb86210dfc156eedff85ec1ccbbb1ae09061e3e90a84f
```

### 4 rank 求解

独立算例目录：

```text
<kunshan-dedicated-ci-root>/OneFLOW_backend_smoke_20260821/runs/m6wingroe_sa_np4
```

配置：

- `mpirun -np 4`
- 每 rank 1 OpenMP thread
- 50 steps
- `ireadwdst=0`，按分区网格现场生成壁面距离

| 作业 | 状态 | 退出码 | 总耗时 |
|---|---|---:|---:|
| MPI run 1 (internal job) | COMPLETED | 0:0 | 138 s |
| MPI run 2 (internal job) | COMPLETED | 0:0 | 138 s |

日志中 4 个 rank 均显示：

```text
nTZones = 4
```

输出文件：

| 文件 | 行数 | SHA256 |
|---|---:|---|
| `aero.dat` | 67 | `ed2dddf615cea3cf3439350fccd20ccb9dab340e2a76c3d89cb0e0c14dcdc548` |
| `res.dat` | 59 | `8744bded027c1a3026005cfbfe6151f8d76dfac84835c4be0a590415aa416a99` |
| `turbres.dat` | 55 | `65c9cc103624cefcc78fbea3c8bfd74d0a31c944b5ea4711a8d9064547f30566` |
| `wallaero.dat` | 7743 | `77794a02c4ae715cab945130b65c06c30722c6770644cedd6d640769872c8a9` |

第二次运行与第一次运行逐文件 `cmp` 结果均为 0。

## 5. GPU backend 框架改动

### 新增 `codes/accel/`

新增约 791 行框架代码，主要接口包括：

- `AccelBackend`：backend 抽象、注册表、CPU/HIP/CUDA/Kokkos backend 类型枚举。
- `AccelRuntime`：读取 `ONEFLOW_ACCEL_BACKEND`、MPI rank/local rank/device 映射和生命周期管理。
- `CpuBackend`：当前默认实现，提供 host allocation/copy/synchronize。
- `DeviceBuffer`：后续 device memory RAII 封装。
- `AccelViews`：与 `MRField` 解耦的 face state、flux、connectivity、residual view。
- `FluxBackend`：未来第一个生产 GPU kernel 的接口，目前未实现数值 kernel。

当前状态：

```text
默认 backend：CPU
HIP/CUDA/KOKKOS：仅保留注册和配置入口，尚未编译对应 adapter
GPU 数值 kernel：尚未接入
```

请求未构建的 backend 时会显式报错，不会静默回退到 CPU。

### CMake 改动

[`codes/CMakeLists.txt`](../codes/CMakeLists.txt:33) 新增：

- `ONEFLOW_ACCEL_BACKEND=CPU|HIP|CUDA|KOKKOS`
- `ONEFLOW_ENABLE_MULTI_DEVICE`
- backend 编译宏和配置检查
- `CGNS_ENABLE=OFF` 时排除 CGNS 源码，便于无 CGNS 环境构建
- MPI 关闭时不再无条件注入 MPI/CGNS 宏

### Runtime 生命周期

[`codes/main/src/SimuImp.cpp`](../codes/main/src/SimuImp.cpp:124) 在 MPI 控制参数广播完成后初始化 backend runtime，在结束阶段 finalize。

[`codes/cuda/src/HybridParallel.cpp`](../codes/cuda/src/HybridParallel.cpp:20) 增加非 MPI 编译保护，并复用统一 runtime 初始化入口。

### 测试脚本改动

[`test/test.py`](test.py:54) 修复：

- 文件一方提前 EOF 时不再误判为成功。
- 进程退出码非 0 时测试失败。
- 支持通过第三个参数指定 suite 文件。
- 通过 `sys.exit()` 正确传播总测试结果。

## 6. 当前限制与后续建议

- 4-zone 与原始单-zone 的浮点结果不应要求逐字节一致；域分解改变了接口通信和归约顺序。正式 MPI baseline 应使用 `m6wing_np4.ofl` 对应的独立 baseline。
- 目前仓库没有提交 33 MB 的分区网格和 4-rank 运行产物，避免膨胀代码仓库；它们保存在昆山隔离目录。
- `FluxBackend` 目前是接口层，尚未替换实际 inviscid flux/residual kernel。
- 下一步推荐优先实现 face flux kernel：先 CPU reference adapter，再 HIP/Z100 adapter，保持 `FaceStateView`/`ResidualView` 不变。
- 回归脚本后续应增加“每次运行前清空 results/restart/log”步骤，避免追加输出污染比较结果。

## 7. 当前本地工作树

当前未提交修改：

```text
M  codes/CMakeLists.txt
M  codes/cuda/src/HybridParallel.cpp
M  codes/main/src/SimuImp.cpp
M  codes/residual/src/Residual.cpp
M  test/test.py
?? codes/accel/
?? test/REGRESSION.md
?? test/REGRESSION_REPORT_20260821.md
?? test/regression_cpu.txt
```

本地没有构建目录、运行输出或 Python 缓存；测试构建和运行均在昆山隔离目录完成。
