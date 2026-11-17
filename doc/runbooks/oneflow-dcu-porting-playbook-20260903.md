# OneFLOW 到海光 DCU 的移植经验与科学计算迁移方法

**日期：** 2026-09-03
**适用范围：** CFD、PDE、有限体积/有限差分/有限元以及其他以 C/C++/Fortran 数值内核为主的科学计算代码。
**目标平台：** Kunshan Hygon DCU / DTK / HIP；本文不把 Kunshan 结论外推到 CUDA、Kokkos 或其他 DCU 集群。

## 1. 本次实机验证结论

当前 OneFLOW 一维 Euler HIP backend contract test 已在 Kunshan `kshdnormal`
分区的真实 DCU 计算节点完成验证：

| 项目 | 结果 |
|---|---|
| 资源 | 1 节点、1 task、8 CPU、`dcu:1`、27G |
| 软件环境 | DTK 26.04、Clang 17、CMake 3.22.0-rc1 |
| 目标架构 | Hygon HCU `gfx906` |
| 配置/编译 | PASS |
| HIP contract test | PASS，6/6 |
| Slurm | `COMPLETED`，退出码 `0:0` |

6 个测试覆盖 backend identity、可见设备、FullTrace、NoTrace、state reuse
和非法请求处理。CPU reference 与 HIP 结果在测试向量上逐项通过严格比较，
并检查密度、压力和有限值等物理不变量。

## 2. 这次验证暴露的工程问题

### 2.1 `module purge` 可能把工具链降级

第一次作业在配置阶段失败。脚本在 `module purge` 后调用无绝对路径的
`cmake`，计算节点最终解析到系统 CMake 2.8.12；它不支持现代 `cmake -S/-B`
调用形式，于是把 build 目录误当成 source 目录。

处理方式是：

- 在 `module purge` 前把 `cmake` 和 `ctest` 解析为绝对路径；或显式加载经过
  profile 认可的 CMake module；
- 作业日志记录实际的 CMake/编译器版本；
- 配置失败时不得继续构建或运行。

### 2.2 CTest 空测试不能视为成功

早期作业中 HIP 测试二进制实际通过了 6/6，但 CTest 过滤器使用可执行文件名，
而 `gtest_discover_tests` 生成的是 suite/test 名称，因此输出 `No tests were found!`
且退出码为 0。这不是测试通过的证据。

修复包括：

- CPU 和 HIP discovery 分别使用 `TEST_PREFIX "CPU."` 与 `TEST_PREFIX "HIP."`；
- 作业构建 CPU/HIP 两个测试目标，删除陈旧的 `*_tests.cmake` 后再做 PRE_TEST discovery；
- CTest 使用 `-R '^HIP\.'`，并校验 `100% tests passed, 0 tests failed out of 6`；
- 同时保留完整 GoogleTest 二进制汇总 `[  PASSED  ] 6 tests.`。

当前 Kunshan 实机结果已经证明修复有效：CTest 正确发现 6 个 `HIP.*` 测试并全部通过。
一般规则仍是“发现了预期数量的测试 + 进程退出码为 0 + 关键汇总可解析”；空测试结果
不能上报为 PASS。

### 2.3 CPU 和 HIP 的测试数量可能不同

CPU contract test 有 5 个测试；HIP 编译宏额外启用 device-visibility test，因此是
6 个。测试脚本不能把一个固定数量硬编码成所有 backend 的验收标准，应按 backend
定义预期数量，或解析 suite summary。当前 CPU 为 5 项、HIP 为 6 项；CTest 前缀
同时避免两个 backend 的同名测试互相覆盖。

## 3. 推荐的科学计算移植路线

不要从“把所有循环改成 HIP kernel”开始。对于 CFD 或数值求解器，推荐按以下门槛
递进，每一层都能独立证明或否定一个假设：

1. **建立 CPU oracle**：固定算例、初值、边界、时间步和迭代数；保存最终状态、
   residual、守恒量和必要的中间量。
2. **抽出 backend-neutral seam**：把 problem、state、upload、advance、download、
   trace 等生命周期从算法调用方中分离；CPU 仍是默认 reference。
3. **先做最小设备探针**：在目标计算节点验证 driver/runtime、编译器、目标架构、
   设备数量、显存分配和最小 kernel。登录节点的 import 或 package visibility 不算证据。
4. **移植单步数值内核**：优先一个时间步或一个 stage，保持数据布局和边界语义，
   逐 stage 比较 flux、residual、state，而不是只比较最终 checksum。
5. **加入物理不变量**：有限值、正密度、正压力、守恒量、边界条件和 CFL/稳定性条件。
   数值误差通过绝对误差 + 缩放相对误差表达，不能只写“保留 15 位小数”。
6. **验证 NoTrace 和 FullTrace 两条路径**：诊断路径可以回传中间数组；生产路径应
   尽量保持 state resident，减少 H2D/D2H 和同步。
7. **验证生命周期复用**：一次分配、多步 advance、一次最终回传；确认第二次运行不会
   依赖残留 device state，也不会产生隐式 host fallback。
8. **再做多卡/MPI**：先验证 local-rank 到 device 的映射、visible device count、halo
   交换和单卡对照，再测 4-rank/4-DCU 或跨节点。
9. **最后做性能优化**：先用端到端 wall-clock 定义主指标，再用 kernel、copy、allocation、
   synchronization 分项定位瓶颈。每个优化块都重跑 correctness 和 MPI regression。

## 4. 面向 CFD/PDE 代码的设计原则

### 4.1 保留一个可验证的 CPU 参考实现

CPU reference 不只是 fallback，也是迁移期间的数值 oracle。生产代码可以共享
方程、边界和初值定义，但 backend 负责执行策略。不要为了迁移先删除成熟 CPU 路径。

### 4.2 优先迁移数据流，再迁移算术

双精度 CFD 通常受访存、临时数组和同步影响。先确定 conserved/primitive/flux 的
布局、每个 stage 的读写缓冲、可在 device 常驻的数组，以及哪些输出只在诊断或最终
结果阶段需要回传。每个小 kernel 都拷回 host 时，kernel speedup 不能代表应用 speedup。

### 4.3 把 trace 设计成可选能力

FullTrace 适合回归和调试，NoTrace 适合生产推进。接口应显式区分两种模式，避免
为了记录诊断量而在生产路径中保留大量 device-to-host copy。

### 4.4 让错误在边界处失败

CreateState、Upload、Advance、Download 应检查尺寸、空指针、步数、trace 选项和设备
状态。GPU kernel 错误、设备缺失和 runtime 初始化失败必须传播成非零退出码；不能只在
日志中打印 warning。

## 5. DCU/DTK 特有的验证注意事项

- 目标架构必须来自目标节点探针或显式 profile；当前 Kunshan Z100 证据是 `gfx906`。
- 不把上游 ROCm/CUDA 的 kernel、Triton、FP8、通信或架构能力直接推断到 DTK。
- module 环境是作业的一部分；记录 module、compiler、CMake、HIP runtime 和
  `rocminfo` 的关键结果。
- Slurm 的 partition、memory 和 GRES 必须遵循 profile 固定配比；不能临时发明资源元组。
- 作业脚本必须捕获 workload rc，并以该 rc 结束；同时检查 Slurm State 和 ExitCode。
- 共享存储上的临时文件使用 `$HOME/.scnet-hpc/tmp/$SLURM_JOB_ID` 等隔离目录，
  不依赖节点本地 `/tmp` 在多卡/多节点间共享。

## 6. 证据与报告模板

每个新 backend 至少记录 ABI/构建、Runtime、Kernel、Functional、Numerical、Integration
和 Performance 七层证据：版本与架构明确、设备可见且初始化成功、最小 kernel 确认在设备
上执行、生命周期和错误处理通过、CPU oracle/物理不变量通过、NoTrace/MPI/residual 通过，
并分别给出端到端与分项时间。

公开报告只保留脱敏后的资源元组、软件版本、架构、验证层级和结论；原始日志、账号、
节点名、作业号、私有路径和凭据留在 CI artifact 或隔离运行目录。

## 7. 当前 OneFLOW 的后续边界

这次验证证明的是一维 Euler HIP backend contract 和 DTK/HIP 运行链路，不等于完整
OneFLOW Navier–Stokes GPU 移植完成。后续应按同一方法扩展到更大的算例、真实 residual
runner、MPI halo 和性能回归；CUDA、Kokkos、跨节点扩展以及未在目标节点实测的库仍应
保持“未验证”状态。


## 8. 2026-09-03 CTest 修复实测记录

修复后的单卡 DCU 作业在 Kunshan `kshdnormal` 完成：

| 项目 | 结果 |
|---|---|
| 资源 | 1 节点、1 task、8 CPU、`dcu:1`、27G |
| 工具链 | DTK 26.04、Clang 17、CMake 3.22.0-rc1 |
| 架构 | Hygon HCU `gfx906` |
| 配置/编译 | PASS |
| 直接 GoogleTest | 6/6 PASS |
| CTest | `HIP.*` 6/6 PASS，100% |
| Slurm | `COMPLETED`，退出码 `0:0` |
| CTest 总耗时 | 约 1.05 秒（作业总耗时约 22 秒） |

这次结果把“HIP 二进制通过”和“测试框架正确发现并执行”两个层面都闭合了。
