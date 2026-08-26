# OneFLOW 一维 Euler 原生 HIP 验证记录

> 记录日期：2026-08-25
> 本文是当前一维主线的结果记录，不代表完整 Navier–Stokes GPU 化，也不替代
> `doc/reports/gpu-backend-delivery.md` 中较早的 GPU backend 框架交付报告。

## 1. 当前结论

- 一维验证程序已经进入本地 `master` 的当前 HEAD。
- CPU 主路径保留，HIP 只作为隔离的一维验证后端；当前不把它描述为完整生产 CFD
  GPU backend。
- 郑州 BW1000 的 `gfx936` 结果已归档，仅作为历史参考，不再进入后续代码测试。
- 昆山 Z100 的目标架构是 `gfx906`，但目前没有完成“可用 Z100 设备节点上的
  自动探测、编译、运行、逐步比对”闭环，因此不能把郑州结果转写成昆山通过。
- 今日已完成郑州归档和本地上游 merge；不提交新的计算作业，不清理昆山目录和日志。

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
/public/home/luql/softwares/projects/OneFLOW_1d_euler_20260825
/public/home/luql/softwares/projects/OneFLOW_1d_euler_auto_20260825
```

本节结果已从郑州远端实验目录移动到以下归档位置，仅保留作历史参考；后续代码测试、回归矩阵和主线结论不再使用郑州结果：

```text
/public/home/luql/softwares/projects/archive/OneFLOW_1d_euler_zhengzhou_bw1000_20260825
/public/home/luql/softwares/projects/archive/OneFLOW_1d_euler_zhengzhou_bw1000_auto_20260825
```

相关作业：

| 作业 | 结果 | 说明 |
| --- | --- | --- |
| `777430` | PASS | CPU-only 一维验证 |
| `777499` | FAIL | 按错误的 `gfx906` 编译，运行时报 `invalid device function` |
| `777518` | PASS（探针） | 节点报告 `gfx936` |
| `777536` | PASS | 显式 `gfx936` 编译运行 |
| `777845` | PASS | 未显式传架构，自动探测 `gfx936` 后 Euler 回归通过 |
| `777894` | PASS（探针） | 另一节点同样报告 `gfx936` |
| `777864` | 待排队记录 | 郑州作业，不纳入昆山主线；今天不再处理 |

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
/public/home/luql/softwares/projects/OneFLOW_1d_euler_kunshan_auto_20260825
```

相关作业：

| 作业 | 结果 | 失败/限制 |
| --- | --- | --- |
| `119933685` | FAIL | `CXX` 残留指向不存在的 `/public/home/luql/software/gcc-11.2.0-install/bin/g++` |
| `119933929` | FAIL | 清除 `CC/CXX` 后，CMake 在自动架构探测阶段失败；`rocminfo` 未得到架构 |
| `119934042` | FAIL（探针） | 节点 `b13r2n16` 输出 `C-3000 module is NOT loaded, possibly no HCU devices` |
| `120043554` | PASS（设备探针） | 节点 `e18r1n15`，DTK 26.04/HIP 6.3.26113，`rocminfo=gfx906`，HCU/KFD 正常 |
| `120045063` | PASS | 显式 `gfx906` 构建；CPU self-test 和 HIP Euler 回归通过 |
| `120045423` | PASS | 不传架构，自动探测得到 `CMAKE_HIP_ARCHITECTURES=gfx906`；Euler 回归通过 |
| `120045915` | PASS | 修复编译器 warning 后重新自动探测和回归；节点 `b17r4n10` |

`120045915` 的最终验证内容：

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
## 8. 今日收尾状态

- [x] 郑州 BW1000 实验目录已移动到 archive，原始内容保留；
- [x] 郑州结果降级为历史参考，不进入后续代码测试和回归矩阵；
- [x] 本地 `master` 已 merge 上游 `origin/master`，检查无冲突；
- [x] 未提交新的 Slurm 作业，未取消已有作业，未 push；
- [x] 昆山 Z100/gfx906 设备探针、显式架构和自动架构一维 Euler 回归均通过；
- [x] 修复一维 Euler 编译 warning，并在昆山自动探测路径重新回归；
- [x] 本地验证 commit `71e86282` 已创建，工作区随后保持干净；
- [ ] 不修改 Windows 原始 WENO3 路线，不改动 `codes/ns`。

## 9. 恢复昆山主线时的最小待办

- [ ] 在 `kseshell` 对应的 Z100 节点确认 C-3000/HCU 模块已加载；
- [ ] 在独立项目目录记录远端副本 `HEAD`、分支和工作区状态；
- [ ] 让自动探测得到 `gfx906`，或显式传入 `CMAKE_HIP_ARCHITECTURES=gfx906`；
- [ ] 先运行 CPU 一维 Euler 基准；
- [ ] 再运行 HIP 一维 Euler，逐 stage 比较三守恒量、primitive、flux、residual；
- [ ] 保存作业脚本、日志、退出码和架构探针输出；
- [ ] 只有昆山回归通过后，才讨论下一次 commit。
