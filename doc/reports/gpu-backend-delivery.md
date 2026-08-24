# OneFLOW GPU Backend 框架与昆山回归测试交付报告

**首次报告：2026-08-21**

**最近更新：2026-08-24**

**代码分支：`master`**

**首轮交付提交：`5f54fc15`**

**验证平台：昆山 SCNet CPU 分区，GCC 7.3.1 + HPC-X MPI**

## 1. 执行摘要

本阶段完成了两项工作：

1. 建立可重复的 CPU 回归基线，验证单进程、空 MPI rank 和真实
   4 进程域分解。
2. 在不改变默认 CPU 数值路径的前提下，搭建统一 accelerator
   backend 框架，为后续 HIP/Z100、CUDA 和 Kokkos 适配预留接口。

验证结论：

- 3 个 CPU 串行算例全部通过，总耗时 152 秒。
- 单 zone 网格使用 4 个 MPI rank 时的空 rank 越界已经修复，修复后
  50 步运行正常结束。
- M6 四 zone / 四 rank 算例连续运行两次均通过，关键输出逐文件一致。
- GPU backend 目前仅为框架接入，尚未替换 CFD 数值 kernel，不能据此
  宣称 GPU 加速已经生效。
- 昆山 workflow 已加入仓库并完成脚本静态检查，但尚未在配置真实
  GitHub secrets 后触发端到端运行。

## 2. 工作范围

### 2.1 已完成

- 定义 3 个短时 CPU 回归算例。
- 修复测试程序的退出码传播和文件 EOF 判断。
- 修复空 rank 下的残差最大值与平均值处理。
- 验证真实四 zone / 四 rank MPI 求解。
- 新增 backend 抽象、runtime、CPU backend 和设备内存封装。
- 为 HIP、CUDA、Kokkos 和多设备映射保留扩展入口。
- 新增昆山 Slurm 回归脚本和 GitHub Actions workflow。
- 增加 SSH 授权、远端路径和 Slurm 命令前置检查。

### 2.2 本阶段未完成

- 未实现 HIP/Z100 数值 kernel。
- 未实现 CUDA 或 Kokkos adapter。
- 未将 `FluxBackend` 接入现有 inviscid flux/residual 主计算循环。
- 未进行 GPU 性能测试或 CPU/GPU 性能对比。
- 未使用真实 GitHub secrets 执行新 workflow 的端到端验证。

## 3. Backend 架构

新增目录 [`codes/accel/`](../../codes/accel/)。

| 层次 | 组件 | 当前职责 |
|---|---|---|
| Backend 抽象 | `AccelBackend` | 初始化、内存、拷贝、同步及能力查询 |
| Backend 注册 | `AccelBackendRegistry` | 按 backend 类型创建实现 |
| Runtime | `AccelRuntime` | backend 选择、MPI/local rank 和设备上下文 |
| CPU 实现 | `CpuBackend` | 当前默认的 host 内存和同步实现 |
| 内存封装 | `DeviceBuffer<T>` | 为后续设备内存提供 RAII 接口 |
| 数据视图 | `AccelViews` | face state、flux、connectivity、residual 视图 |
| 数值接口 | `FluxBackend` | 未来 inviscid flux 和 residual kernel 接口 |

当前执行关系：

```text
OneFLOW solver
    |
    +-- AccelRuntime
    |      |
    |      +-- CPU backend（当前可用）
    |      +-- HIP adapter（预留）
    |      +-- CUDA adapter（预留）
    |      +-- Kokkos adapter（预留）
    |
    +-- 现有 CPU 数值 kernel（当前实际计算路径）
```

默认 backend 为 CPU。请求尚未编译的 HIP、CUDA 或 Kokkos backend 时，
程序会明确报错，不会静默回退到 CPU。

### 3.1 配置入口

[`codes/CMakeLists.txt`](../../codes/CMakeLists.txt) 新增：

- `ONEFLOW_ACCEL_BACKEND=CPU|HIP|CUDA|KOKKOS`
- `ONEFLOW_ENABLE_MULTI_DEVICE`
- backend 名称校验和默认 backend 编译宏
- `CGNS_ENABLE=OFF` 时排除 CGNS 源码
- MPI 关闭时不再无条件注入 MPI/CGNS 并行宏

运行时可通过以下环境变量选择 backend：

```bash
export ONEFLOW_ACCEL_BACKEND=CPU
```

多设备上下文已包含：

- world rank / world size
- local rank / local size
- requested device / selected device
- device count

目前这些字段用于建立接口，尚未执行真实 GPU 绑定。

### 3.2 生命周期

[`codes/main/src/SimuImp.cpp`](../../codes/main/src/SimuImp.cpp) 在并行环境和
控制参数初始化后启动 backend runtime，并在程序结束前 finalize。

[`codes/cuda/src/HybridParallel.cpp`](../../codes/cuda/src/HybridParallel.cpp)
增加了非 MPI 编译保护，并使用同一个 backend runtime 入口。

## 4. 回归测试标准

### 4.1 CPU 基准算例

测试清单由 [`test/suites/cpu-serial.txt`](../../test/suites/cpu-serial.txt) 定义。参考结果使用
各算例已有的 `autotest/` 输出。

| 算例 | 覆盖范围 | 实测耗时 |
|---|---|---:|
| `plateuns2dslau2` | 2D 层流平板、SLAU2、黏性通量 | 3 秒 |
| `rae2822_roe_sa` | 2D 跨声速翼型、Roe、SA 湍流 | 10 秒 |
| `m6wingroe_sa` | 3D 翼型、Roe、SA 湍流 | 139 秒 |
| **合计** |  | **152 秒** |

单个算例均满足约 1–3 分钟的目标，完整串行套件约 2.5 分钟。

### 4.2 通过判据

- 求解进程退出码为 0。
- 实际输出通过现有数值容差比较。
- 文件长度不一致或一方提前 EOF 判为失败。
- 每次运行前清理 `results/`、`restart/` 和 `log/`。
- Slurm 作业同时检查调度器状态和工作负载退出码。
- MPI 基准使用独立的四 zone 参考输出，不直接复用单 zone 基准。

### 4.3 执行入口

```bash
python3 test/test.py \
  "mpirun -np 1" \
  "/path/to/OneFLOW" \
  test/suites/cpu-serial.txt
```

## 5. 昆山实测结果

### 5.1 环境与边界

| 项目 | 配置 |
|---|---|
| 集群 | 昆山 SCNet |
| 分区 | `kshcnormal` CPU 分区 |
| 编译器 | GCC 7.3.1 |
| MPI | HPC-X MPI |
| Accelerator | 本阶段未提交 GPU 作业 |
| 运行目录 | 昆山专用隔离目录 |

以下结果均是 CPU/MPI 验证，不包含 Z100/HIP 功能或性能结论。

### 5.2 三算例串行回归

| 指标 | 结果 |
|---|---|
| 算例数 | 3 |
| 通过数 | 3 |
| 总耗时 | 152 秒 |
| Slurm 状态 | `COMPLETED` |
| 退出码 | `0:0` |

日志摘要：

```text
CASE_RESULT=plateuns2dslau2 RC=0 ELAPSED_SECONDS=3
CASE_RESULT=rae2822_roe_sa RC=0 ELAPSED_SECONDS=10
CASE_RESULT=m6wingroe_sa RC=0 ELAPSED_SECONDS=139
REGRESSION_RC=0 TOTAL_ELAPSED_SECONDS=152
```

### 5.3 空 rank 问题

原有 CFD 测试网格为单 zone。使用 4 个 MPI rank 时，只有一个 rank
持有本地 zone，其余 rank 的 `dataList` 为空。旧实现无条件访问
`dataList[0]`，导致 `ResMax::CalcMax()` 越界。

[`codes/residual/src/Residual.cpp`](../../codes/residual/src/Residual.cpp) 的
修复包括：

- 空 `dataList` 使用安全哨兵值。
- 本地最大残差执行全局 `MPI_MAX` 归约。
- 最大值对应的 zone、cell、坐标和体积由 owner rank 广播。
- 无全局 cell 时返回安全元数据。
- 平均残差对全局 cell 数为 0 的情况进行保护。

修复后，原始单 zone M6 使用 `mpirun -np 4` 完成 50 步，求解器退出码
为 0。

早期测试脚本曾因复用旧目录而报告失败：输出文件采用追加方式，导致
`aero.dat` 行数与基准不同。这不是求解器失败。新的昆山 CI 脚本会在
运行前重新建立输出目录，避免旧结果污染。

### 5.4 四 zone / 四 rank MPI 回归

原始 M6 网格：

- 294,912 cells
- 1 zone

使用仓库已有 METIS 工具拆分后：

- 4 zones
- 每个 zone 73,728 cells
- 分区网格 SHA256：
  `0feacacb3ee49443d5aeb86210dfc156eedff85ec1ccbbb1ae09061e3e90a84f`

运行配置：

- `mpirun -np 4`
- 每个 rank 1 个 OpenMP thread
- 50 steps
- `ireadwdst=0`

| 运行 | 状态 | 退出码 | 耗时 |
|---|---|---:|---:|
| MPI run 1 | `COMPLETED` | `0:0` | 138 秒 |
| MPI run 2 | `COMPLETED` | `0:0` | 138 秒 |

两次运行的以下文件逐文件一致：

| 文件 | 行数 | SHA256 |
|---|---:|---|
| `aero.dat` | 67 | `ed2dddf615cea3cf3439350fccd20ccb9dab340e2a76c3d89cb0e0c14dcdc548` |
| `res.dat` | 59 | `8744bded027c1a3026005cfbfe6151f8d76dfac84835c4be0a590415aa416a99` |
| `turbres.dat` | 55 | `65c9cc103624cefcc78fbea3c8bfd74d0a31c944b5ea4711a8d9064547f30566` |
| `wallaero.dat` | 7743 | `77794a02c4ae715cab945130b65c06c30722c6770644cedd6d640769872c8a9` |

这说明当前四 zone / 四 rank CPU MPI 路径具有可重复性。由于域分解会
改变接口通信和浮点归约顺序，四 zone 结果不应要求与原始单 zone 输出
逐字节一致。

## 6. 测试程序与 CI 标准化

### 6.1 `test.py`

[`test/test.py`](../../test/test.py) 已修复：

- 文件一方提前 EOF 时不再误判为通过。
- 求解进程非零退出码会导致测试失败。
- 支持第三个命令行参数指定 suite 文件。
- 可选读取 `test/baselines/residual-baseline.json`，逐行比较
  `res.dat`/`turbres.dat` 残差曲线。
- 使用 `sys.exit()` 将最终状态传递给 shell、Slurm 和 CI。

### 6.2 残差基线数据库

残差数据库位于
[`test/baselines/residual-baseline.json`](../../test/baselines/residual-baseline.json)。
它保存三个 CPU 串行算例的完整残差曲线，而不是只保存最后一个残差。

| 内容 | 说明 |
|---|---|
| schema/version | `oneflow.residual-baseline` / v2 |
| baseline id | `cpu-serial-e15-v1` |
| 算例 | 3 个 |
| 残差文件 | 5 个 `res.full.dat`/`turbres.full.dat` |
| 每个文件 | 50 步、所有残差分量、变量名、首值、末值、最大幅值、SHA256 |
| 比较规则 | `iter/sub-iter` 和行数必须一致，残差绝对容差 `1e-15`，拒绝 NaN/Inf |
| 输出开关 | `ONEFLOW_RESIDUAL_TEST_OUTPUT=1`；默认关闭 |
| 数值格式 | `double` 的 `max_digits10`（17 位有效数字） |

旧的 `res.dat`/`turbres.dat` 和低精度数据库仍然保留，作为兼容性参考；高精度
`.full.dat` 只在显式开启环境变量时生成，JSON 是统一索引和比较入口。只有接受
新的数值基线时才重新生成 JSON：

```bash
python3 test/baselines/build_residual_db.py \
  --high-precision \
  --validated-commit <commit-that-passed-the-baseline> \
  --output test/baselines/residual-baseline-e15.json
```

普通 CI 运行只读取数据库，不生成新基线文件。每次 CI 的日志、作业号和
环境信息保存在 workflow artifact 或昆山隔离目录中。

### 6.3 昆山 workflow

新增
[`kunshan-regression.yml`](../../.github/workflows/kunshan-regression.yml)，
支持手动选择：

- `cpu-serial`
- `mpi4`
- `all`

工作流顺序：

```text
trusted branch 检查
    -> GitHub secrets 完整性检查
    -> SSH key / known_hosts 检查
    -> 昆山 SSH + Slurm preflight
    -> 创建隔离运行目录
    -> 同步源码
    -> Slurm 构建和回归
    -> 检查 State + ExitCode
    -> 下载并上传 artifacts
```

如果没有昆山 SSH 授权、固定 host key 或私有集群配置，workflow 会在
创建远端运行目录和提交作业之前中止，并给出明确错误信息。

每次运行使用：

```text
KUNSHAN_REMOTE_ROOT/runs/<github-run-id>_<attempt>/
```

源码、构建、运行输出、Slurm 日志和 artifacts 均在该目录下隔离。

### 6.4 CI 验证状态

| 项目 | 状态 |
|---|---|
| Shell 脚本语法检查 | 通过 |
| `test.py` Python 编译检查 | 通过 |
| 残差数据库结构与源文件校验 | 2026-08-24 昆山通过 |
| 三算例残差数据库正向比较 | 2026-08-24 昆山通过，最大差值 0 |
| 超容差残差负向测试 | 2026-08-24 昆山正确拦截 |
| Git diff whitespace 检查 | 通过 |
| 昆山独立手工作业 | 通过 |
| GitHub Actions + 真实 secrets 端到端运行 | 尚未执行 |

2026-08-24 的昆山验证只覆盖新增数据库机制，没有重新运行求解器。数据库
中的数值来源仍是已通过 CPU 回归的提交 `5f54fc15`。

## 7. 风险与限制

- 当前 backend 不包含实际 GPU kernel。
- `FluxBackend` 尚未连接到生产数值路径。
- 四 rank fixture 约 33 MB，保存在昆山私有目录，未提交 Git。
- GitHub workflow 依赖仓库 secrets 和昆山私有配置。
- CPU 单 zone 与 MPI 四 zone 应分别维护数值基准。
- 当前测试主要覆盖结果一致性，尚未建立性能基线和性能退化阈值。

## 8. 后续实施方案

建议继续采用“CPU 主路径 + 原生 HIP backend + 可扩展接口”的方案：

1. 先实现 CPU `FluxBackend` reference adapter，使新接口与现有 CPU
   kernel 产生一致结果。
2. 整理 face state、normal、area、connectivity 和 residual 的连续
   数据视图，减少 GPU 侧零散拷贝。
3. 以 HIP/Z100 为第一生产 backend，实现 inviscid face flux。
4. 使用本报告中的 3 个串行算例和四 rank 算例进行数值回归。
5. 加入 rank-to-device 映射，再验证单节点多 GPU。
6. 扩展 MPI halo/边界同步，进行多节点多 GPU 验证。
7. CUDA 通过相同 backend 接口实现；Kokkos 是否引入，在 HIP 路径
   稳定后根据跨平台维护成本再决定。

## 9. 交付物

| 内容 | 路径 |
|---|---|
| 完整交付报告 | `doc/reports/gpu-backend-delivery.md` |
| HTML 报告 | `doc/reports/gpu-backend-delivery.html` |
| 测试规范 | `test/README.md` |
| CPU suite | `test/suites/cpu-serial.txt` |
| Backend 框架 | `codes/accel/` |
| 昆山 workflow | `.github/workflows/kunshan-regression.yml` |
| 昆山 CI 说明 | `ci/kunshan/README.md` |
| 私有配置模板 | `ci/kunshan/config.example.sh` |

## 10. 最终结论

本阶段已建立可在昆山复现的 CPU/MPI 正确性基线，并完成不侵入现有
CPU 数值路径的 GPU backend 框架。现有代码可以作为 HIP/Z100 移植的
起点，但当前交付的含义是“基础设施与回归基线就绪”，不是“GPU 数值
移植完成”。

## 附录 A：首轮测试证据

首轮验证记录如下：

- 三个串行算例耗时分别为 3、10、139 秒，总计 152 秒。
- 空 rank 修复后，单 zone M6 使用 4 个 MPI rank 完成 50 步，退出码为 0。
- 四 zone / 四 rank M6 连续运行两次，均耗时 138 秒，Slurm 状态均为
  `COMPLETED`，退出码均为 `0:0`。
- 两次 MPI 运行的 `aero.dat`、`res.dat`、`turbres.dat` 和
  `wallaero.dat` 逐文件一致。

后续单次 CI 运行的日志、作业号和环境快照应保存在 GitHub Actions
artifact 或昆山隔离运行目录中，不再新增日期型仓库报告。
