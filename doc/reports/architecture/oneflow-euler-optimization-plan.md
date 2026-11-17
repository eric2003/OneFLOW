# OneFLOW 一维 Euler 高性能 backend 优化计划

> 计划日期：2026-08-28
> 当前平台：昆山 Z100 / native HIP / `gfx906`
> 优先级：Euler 高性能 backend 优先；WENO5、二维、三维暂缓。

## 1. 目标与原则

当前阶段只把一维 Euler 场景做扎实，用它验证 OneFLOW 面向 CPU/HIP 的完整
backend 生命周期。优化对象不是单个 kernel，而是数据所有权、显存驻留、执行
调度、验证和性能测量组成的整体路径。

必须遵守：

1. 正确性是硬门槛，任何性能改动先通过 CPU/HIP 回归；
2. 高性能 backend 是目的，不能以“调用了 GPU”或“有显存占用”作为完成标准；
3. CPU 默认路径继续保留；
4. 验证模式和生产性能模式分开；
5. 当前不展开 WENO5、二维、三维、Kokkos 或完整 NS 主线；单节点多卡 MPI 纳入当前 Euler 验收。

## 2. 已锁定基线

昆山构建和 correctness 作业 `JOB_ID_OMITTED`：

- 自动目标架构：`gfx906`；
- CPU self-test：PASS；
- smooth-periodic CPU/HIP：PASS；
- Sod-transmissive CPU/HIP：PASS；
- 三守恒量、primitive、flux、smooth residual 和 RK state 误差均为 0。

昆山性能作业 `JOB_ID_OMITTED` 使用 1 个串行 CPU 线程和 1 张 DCU。Slurm 申请了
4 个 CPU 核用于作业和构建环境，但 CPU benchmark 本身没有 OpenMP/MPI。

| nx | CPU total (ms) | HIP total (ms) | speedup | kernel ratio | D2H ratio | final error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 65536 | 1130.592 | 572.878 | 1.974x | 2.51% | 68.88% | 0 |
| 262144 | 4668.764 | 1889.750 | 2.471x | 1.98% | 78.82% | 0 |
| 1048576 | 18967.298 | 9038.656 | 2.098x | 1.58% | 84.88% | 0 |
| 4194304 | 88670.222 | 37329.346 | 2.375x | 1.56% | 85.93% | 0 |

这组是旧的 trace 型 `Step` 基线，不能代表 stateful `NoTrace` 性能。

当前 `Step(host_state, trace)` 每个时间步执行 7 次 allocation、2 次 H2D，并在
每个 RK stage 把 left/right state、flux、residual 和 updated state 全部 D2H。
4M 网格下每步约回传 1.4 GiB；100 个计时时间步累计约 140 GiB，而最终三守恒量
state 只有约 96 MiB。

## 3. Euler stateful 实施结果

昆山作业 `JOB_ID_OMITTED` 已完成 stateful backend 编译、旧 `FullTrace` correctness
回归和四规模 `NoTrace` benchmark。每个规模 100 steps、2 repeats、1 warmup；
CPU benchmark 是单线程，HIP 是 1 张 DCU。

| nx | CPU advance (ms) | HIP advance (ms) | kernel event (ms) | kernel ratio | steady speedup | final error | syncs |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 65536 | 2194.558 | 21.357 | 21.262 | 99.56% | 102.76x | 0 | 2 |
| 262144 | 9370.527 | 68.325 | 68.207 | 99.83% | 137.15x | 0 | 2 |
| 1048576 | 38820.703 | 279.068 | 278.916 | 99.95% | 139.11x | 0 | 2 |
| 4194304 | 178937.295 | 1152.799 | 1152.721 | 99.99% | 155.22x | 0 | 2 |

`NoTrace` 路径每次 state 创建时分配一次 buffer，初始状态只 H2D 一次，最终状态
只 D2H 一次；100 个时间步不产生逐步 allocation 或逐 stage D2H。按完整生命周期
（create + upload + advance + download）计算，D2H 占比约 4.4--6.1%，已低于第一阶段
10% 目标。四个规模 final error 均为 0，CPU/HIP checksum 均一致。

HIP event 统计与 Advance 墙钟时间基本一致，说明 stateful benchmark 的等待时间
主要是设备 kernel 本身，而不是主机同步；每个 repeat 只有一次最终 event synchronize。

## 3. 目标架构

统一接口从“单次函数调用统一”升级为“资源和执行生命周期统一”：

```text
Create backend context
    -> Create per-zone state
    -> Upload initial state once
    -> Advance multiple RK steps on the selected backend
    -> Optional validation/checkpoint/output
    -> Download final state once
```

建议的职责边界：

```cpp
class EulerBackend
{
public:
    virtual std::unique_ptr<EulerState> CreateState(
        const EulerProblem & problem ) const = 0;
    virtual void Upload(
        EulerState & state, const double * hostState ) const = 0;
    virtual void Advance(
        EulerState & state, int steps, const EulerRunOptions & options ) const = 0;
    virtual void Download(
        const EulerState & state, double * hostState ) const = 0;
};
```

- CPU state 持有主存；
- HIP state 持有 device buffers、stream、event 和设备元数据；
- `EulerRunOptions` 控制 `NoTrace` 和 `FullTrace`；
- 现有 `Step` 保留为 compatibility/validation wrapper，不作为性能主路径。

## 4. 实施阶段

### 阶段 A：状态化统一接口

- [x] 引入 `EulerState`、`EulerRunOptions` 和多步 `Advance`；
- [x] CPU/HIP 使用相同生命周期语义；
- [x] 保留现有 `Step`，由新接口或兼容层实现；
- [x] 不改变数值公式、边界处理和 acceptance tolerance；
- [x] benchmark 将 CPU 与 HIP 阶段分开计时。

完成判据：

- 本地 CPU build/self-test 通过；
- 昆山 HIP build 通过；
- 现有 smooth/Sod stage-by-stage 回归误差不变；
- CPU 默认路径不依赖 HIP runtime。

### 阶段 B：HIP 持久 device state

- [x] 每个 case/zone 只分配一次 device buffers；
- [x] 初始状态只 H2D 一次；
- [x] 所有 RK stage 和多个时间步连续在 GPU 上执行；
- [x] 使用 device pointer 轮换复用 `base/current/next`；
- [x] `NoTrace` 模式只在最终输出时 D2H；
- [x] 删除生产路径中的逐 stage `hipEventSynchronize`；
- [x] event 只在 profiling 模式记录，不改变生产同步边界。

完成判据：

- allocation 次数与 steps 无关；
- H2D/D2H 字节数与 steps 无关，除非显式要求 checkpoint/output；
- `NoTrace` 最终状态与 CPU reference 满足既有误差阈值；
- `FullTrace` 仍能提供逐 stage correctness 证据。

### 阶段 C：验证与诊断解耦

- [x] 完整 trace 仅由 correctness 模式启用；
- [x] 性能模式不保存逐 stage trace；
- [x] NoTrace 和 FullTrace 已通过统一 EulerRunOptions 分开；
- [x] GPU-only 持续运行测试已完成，避免短 kernel 被采样周期遗漏；
- [x] CPU benchmark、HIP benchmark 和 dmon 采样窗口分离；
- [ ] checksum、最大误差、有限值和正状态检查优先用 GPU reduction，只回传标量；
- [x] 4-rank/4-DCU 单节点 MPI 严格回归已接入统一测试入口。

完成判据：

- correctness 和 performance 两种模式输出清晰区分；
- 利用率采样不包含 CPU benchmark 阶段；
- kernel 分类时间、传输字节数、同步次数和 allocation 次数可追溯；
- MPI 回归同时验证 global hash、物理状态、严格误差和 rank/device 映射。

### 阶段 D：Euler kernel 和数据流优化

只在阶段 B、C 完成后开始：

- `NoTrace` 模式不生成 left/right trace buffers；
- 评估 face flux、residual 和 update 的融合边界；
- 减少当前每步 9 次 kernel launch；
- 复用 scratch buffers，避免不必要的 device-to-device copy；
- 检查内存合并访问、寄存器压力、occupancy 和有效带宽；
- 不启用会破坏当前数值判据的 fast-math。

完成判据：

- 每项 kernel 改动单独通过 correctness；
- 有 kernel 分类 timing 证明收益；
- 没有仅凭利用率变化宣称优化成功。

### 阶段 E：昆山验收

当前 Euler 验收覆盖单卡、32 核 CPU MPI 和单节点 4 卡/4 rank HIP MPI。

正确性矩阵：

- smooth-periodic；
- Sod-transmissive；
- 逐 RK stage `FullTrace`；
- 多步 `NoTrace` 最终状态；
- 三守恒量、primitive、flux、residual、finite/positive-state；
- 4-rank/4-DCU CPU/HIP 严格对比与 global hash。

性能矩阵：

```text
nx=65536
nx=262144
nx=1048576
nx=4194304
```

当前已完成：

- [x] 32-rank CPU MPI：1-rank 与 32-rank global hash 在四个规模一致；
- [x] 4-rank/4-DCU HIP MPI：visible_devices=4，CPU/HIP 严格比较通过；
- [x] 4-rank/4-DCU HIP MPI 四规模 benchmark 通过；
- [ ] 跨节点 MPI、CUDA、Kokkos 和完整 NS 主线仍未验证。

后续每个优化块统一重跑上述 correctness、32-rank CPU 和 4-DCU MPI 回归，不再为单个微优化点单独建立一套测试体系。

## 5. WENO5 解锁条件

以下条件全部满足后，才恢复 WENO5：

- Euler stateful backend 接口稳定；
- 持久 device state 完成；
- `NoTrace` 与 `FullTrace` 均通过昆山 correctness；
- 四规模性能复测完成；
- D2H、allocation 和同步不再是主要瓶颈；
- 接口能够复用到 WENO5，而不复制一套显存生命周期管理。

二维和三维继续排在 WENO5 之后。

## 6. 当前执行顺序

```text
状态化 backend 接口
    -> HIP 持久 device state
    -> NoTrace/FullTrace 解耦
    -> 昆山单卡 correctness 与性能
    -> 32-rank CPU MPI hash 回归
    -> 4-rank/4-DCU HIP MPI 严格回归与四规模 benchmark
    -> 阶段性代码/报告 commit
    -> GPU reduction 与 kernel/dataflow 优化
    -> 再次统一 correctness 和性能验收
    -> WENO5
```
