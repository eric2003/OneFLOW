# OneFLOW 3D Euler 主求解器 DCU 性能报告（当前）

> 状态：首份 3D 主求解器性能报告，替换 `doc/plans/oneflow-development-todo.md` 中的
> 临时数字。测量日期 2026-11-02。
> 范围：受限无粘 3D 路径（五方程 / Lax–Friedrichs / `ilim=0` / `vismodel=0` /
> 稳态显式三阶段 RK / 1 rank 1 zone）。本报告不代表粘性、湍流、Roe、LU-SGS、
> 多重网格、双时间或多卡路径已完成 DCU 验证。

## 1. 结论摘要

同一节点、同一算例、交替执行的配对测量（8 zone / 8 rank）：

| 后端 | 资源 | 中位数 ms/step | 3 次样本 |
|---|---|---:|---|
| CPU | 8 ranks × 1 core | **288.20** | 285.03 / 288.78 / 288.20 |
| HIP | 1 DCU（8 ranks 共享） | **15.59** | 15.59 / 16.49 / 14.83 |

**加速比 18.5×（中位数口径），18.4×（均值口径）。**

该比值是 8:1 槽口径：1 DCU 槽按 8 核计，两侧使用同一张 8 zone 分区网格、
同样的 1000 步和相同的内迭代设置。

## 2. 测试对象与口径

### 2.1 算例

- 网格：m6wing，294,912 cells，METIS 切成 8 个 zone，每个 zone 36,864 cells；
- 每个 rank 恰好负责一个 zone（求解器强制 `zone 数 = rank 数`）；
- 设置：`ischeme=5`（Lax–Friedrichs）、`ilim=0`、`vismodel=0`、`time_integral=1`
  （多阶段显式）、`rk_stage=3`、固定 CFL `0.01`、`startStrategy=0`；
- 步数：1000，warmup 50（warmup 不计入），每个配置交替取样 3 次。

### 2.2 为什么关闭中间输出

基准把 `nFieldSave`/`nResSave`/`nForceSave` 都设为 100000，并把 `nVisualSave`
也设为 100000。实测可视化输出在默认 `nVisualSave=50` 下摊薄约
**0.8 ms/step（约 5%）**，属于周期性写盘而不是推进成本，因此从推进口径中剔除。
含输出的数字：CPU 口径不变，HIP 约 16.8–17.0 ms/step。

### 2.3 环境

| 项目 | 值 |
|---|---|
| 平台 | 昆山 DCU 队列，单节点 |
| 加速器 | Hygon DCU，`gfx906`，1 张 |
| DTK | 26.04 |
| 编译器 | GCC 16.2.0 |
| MPI | OpenMPI 5.0.11，`mpirun --bind-to core --map-by slot` |
| 加速器后端 | `ONEFLOW_ACCEL_BACKEND=HIP` |
| 相对上次基线的开关 | `ONEFLOW_HIP_SKIP_UNUSED_INTERFACE=1`（接口未用字段不打包） |

两侧在同一作业内交替执行，使用同一份二进制与同一份算例派生副本。

## 3. 结果

### 3.1 8:1 槽配对

| 样本 | CPU 8 ranks | HIP 8 ranks / 1 DCU | 比值 |
|---|---:|---:|---:|
| 1 | 285.03 | 15.589 | 18.28× |
| 2 | 288.78 | 16.494 | 17.51× |
| 3 | 288.20 | 14.826 | 19.44× |
| **中位数** | **288.20** | **15.59** | **18.49×** |
| 均值 | 287.34 | 15.64 | 18.38× |

### 3.2 并行粒度（同一张网格切分，另一会话，非同一节点）

| 并行粒度 | cell/rank | HIP ms/step | CPU ms/step | 比值 | HIP 资源 |
|---|---:|---:|---:|---:|---|
| 1 zone / 1 rank | 294,912 | 58.58 | 1646.49 | 28.1× | 1 DCU |
| 8 zone / 8 rank | 36,864 | 20.07 | 284.27 | 14.2× | 1 DCU |
| 32 zone / 32 rank | 9,216 | 17.80 | 83.88 | 4.7× | 4 DCU |

这张表用于说明扩展行为，不与 §3.1 混用：它采集于更早的会话和节点，且当时的
HIP 侧尚未包含本报告的优化。

读法：

- 单卡上把 1 个 rank 提高到 8 个 rank（同一张卡）带来 **2.9×**，说明单 rank
  喂不满一张 DCU；
- 继续切到 32 zone 只再快 11%，却需要 4 张卡，属于过度切分；
- 比值随粒度下降是分母（CPU 扩展更好）造成的，device 侧绝对时间是下降的。

### 3.3 与一维端口的对照（同一节点，同一会话）

| 口径 | 比值 |
|---|---:|
| 1D，单线程 CPU vs 1 DCU（advance） | 282× / 299×（262k / 1M cells） |
| 3D，8 核 vs 1 DCU（本报告） | **18.5×** |
| 3D，32 核 vs 1 DCU（§3.2） | 4.7× |
| 1D，应用级 lifecycle vs 32 rank CPU（历史） | 3.85–9.70× |

一维的高比值来自两点：算例本身只有通量计算（无梯度、无重构、无界面交换、
无边界模型），以及整个宏步常驻设备。三维每 cell 的算法量高一到两个数量级，
所以上限本来就低。在应用级、32 核分母这两个条件下，3D 的 4.7× 与一维的
3.85–9.70× 处于同一档。

## 4. 这个数字是怎么得来的

从首份 3D 基线（2026-10-25，`31.6 ms/step`，同配置 paired）到本报告，共四项
改动，全部经同一组数值门禁验证：

| 项 | 内容 | 同节点 paired 收益 |
|---|---|---:|
| 1 | 梯度主机镜像只在最后一个 RK 阶段做 | −10.5% |
| 2 | RK 基线与残差基底不再从主机搬，残差在设备侧归零 | −22.1% |
| 3 | 不再搬 host 历史层（`q2←q1←q` 每步全场拷贝） | −14% |
| 4 | 设备单网格路径不再写 host 残差账本（`ZERO_RESIDUALS`/`STORE_RHS`） | 约 −0.4 ms |
| 5 | RK 循环进求解器驱动（驻留宏步），策略环境变量与就绪查询缓存 | 约 −3.5 ms |
| 6 | 重构直接读取边界阶段留在设备上的壁面状态 | 约 −0.85 ms |

合计 `31.6 → 15.6 ms/step`（1.99×）。

同轮完成、但有明确负面结论的三项（已写入交接，避免重复尝试）：

| 尝试 | 结果 |
|---|---|
| 删除每步第 4 轮界面交换 | 200 步跨后端残差由 `2.33e-18` 拓到 `5.54e-7`，属算法变更，已撤回 |
| 删除 stage 之间的全卡同步 | 无收益（22.2 vs 22.6 ms/step），已撤回 |
| GPU-aware MPI 设备直发 | 探针通过，但真实界面交换慢 3×（48–50 vs 15.7 ms/step），不可用 |

## 5. 数值验证（与本报告同一天、同一 FE 提交）

| 门禁 | 结果 |
|---|---|
| 3 步多区 trace + resident contract | `MULTIZONE_TRACE_PASS`，`records=24 failures=0 max_abs_diff=0` |
| 200 步跨后端残差 | `2.331035e-18` |
| 200 步跨后端气动力 | `0.0` |
| CPU 五算例 normal (1e-8) | 5/5 |
| CPU 五算例 strict (1e-15) | 5/5 |
| restart 等价（设备权威 q 开/关） | `4.16e-15` / `5.28e-15` |

## 6. 覆盖范围与未覆盖项

已验证：五方程、无粘、Lax–Friedrichs、`ilim=0`、单网格层、稳态显式三阶段 RK、
1 rank 1 zone、8 zone 同节点、1 DCU。

未覆盖（不得从本报告外推）：

- Roe / HLLC / AUSM；limiter（`ilim≠0`）；层流粘性与 SA 湍流；
- LU-SGS 与双时间；多重网格；
- 跨节点多卡、GPU-aware 设备直发（实测更慢）；
- ≥1M cell 的合法分区网格（仓库目前只有 294,912 cells 一套，且缺 `.skl` 大网格）；
- 旗舰算例 `test/cfd/m6wingroe_sa`（Roe + LU-SGS + SA），它仍完全在 CPU 上运行。

## 7. 复现

```bash
# 同一二进制、同一节点，交替采集 CPU / HIP
bash scripts/paired-ab.slurm        # 需要设定 ONEFLOW_ACCEL_BACKEND
bash scripts/gate.slurm             # trace + 200 步跨后端 + 稳定性
bash scripts/cpu-regression.slurm   # CPU 五算例 normal/strict
```

关键开关与含义见 `doc/plans/oneflow-development-todo.md` §6 与
`doc/handoff/oneflow-current-handoff.md` §5。

## 8. 口径纪律

- 只比同节点、同算例、同 `steps/warmup/repeats`、同 rank/资源；
- 不以 100 步或 kernel-only 作生产结论；
- 1 DCU 槽按 8 核计；
- 探索性单作业数据只进待办文档，不进本报告；
- 本报告为 Markdown 源；如需 HTML 版本，按
  `doc/reports/README.md` 的命名规则同 stem 生成。
