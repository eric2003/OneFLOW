# OneFLOW 一维 Euler MPI rank/DCU 矩阵（2026-10-19）

> 范围：`ports/kunshan/oneflow_1d_hip` standalone 一维 Euler MPI 端口。  
> 这不是 3D Navier–Stokes 主求解器的 zone/MPI 性能报告。

## 1. 实验设计

该端口没有 zone。`EulerMpiBenchmark.cpp` 的
`Partition(globalNx, rank, ranks, localNx, globalOffset)` 按当前 MPI world
rank 数重新划分同一个全局 cell array；每个 rank 都执行三阶段 RK、halo
交换和局部更新。因此不同 rank 数之间必须固定全局规模和物理输入，不能复用
另一 rank 数的 local timing。

本轮固定：

| 项目 | 设置 |
| --- | --- |
| 平台 | 昆山 `kshdnormal`，单节点 |
| CPU 资源 | 32 CPU cores，`cpus-per-task=1` |
| DCU 资源 | 4 DCUs，`gfx906` |
| 网格 | `global_nx=65,536 / 262,144 / 1,048,576 / 4,194,304` |
| 生命周期 | `steps=100`、`warmup=1`、`repeats=2` |
| CPU 矩阵 | 8 / 16 / 32 MPI ranks，1 core/rank |
| HIP 矩阵 | 1 / 2 / 4 MPI ranks，1 rank/DCU |
| 比值 | `CPU lifecycle_max_ms / HIP lifecycle_max_ms` |

每一组还强制检查：

- `local_ranks == requested ranks`；
- `decomposition_active_ranks == requested ranks`；
- `global_cells_sum == global_nx` 且 `local_nx_min > 0`；
- `halo_exchanges_sum == expected_halo_exchanges`；
- HIP `visible_devices == requested DCU count`；
- CPU/HIP `final_hash` 一致，且 `decomposition=PASS`。

因此这里的“CPU 8/16/32 核”和“DCU 1/2/4 卡”都是实际参与域分解和计算
的 MPI ranks/device，不是只申请了 `cpus-per-task` 或 GRES。

## 2. 单节点结果

下表每个单元格写作：`CPU ranks / HIP ranks (DCUs): CPU/HIP`。
比值大于 1 表示对应 DCU 作业更快。

| global_nx | CPU 8 / HIP 1 (1 DCU) | CPU 8 / HIP 2 (2 DCUs) | CPU 8 / HIP 4 (4 DCUs) |
| ---: | ---: | ---: | ---: |
| 65,536 | 4.764895x | 4.628479x | 4.583053x |
| 262,144 | 18.632080x | 22.724531x | 25.483136x |
| 1,048,576 | 30.024863x | 49.808238x | 64.643977x |
| 4,194,304 | 30.652686x | 56.529423x | 106.768257x |

| global_nx | CPU 16 / HIP 1 (1 DCU) | CPU 16 / HIP 2 (2 DCUs) | CPU 16 / HIP 4 (4 DCUs) |
| ---: | ---: | ---: | ---: |
| 65,536 | 2.364475x | 2.361155x | 2.041363x |
| 262,144 | 6.695690x | 8.262989x | 8.998183x |
| 1,048,576 | 13.999722x | 23.423211x | 30.985396x |
| 4,194,304 | 14.677293x | 27.396774x | 49.789708x |

| global_nx | CPU 32 / HIP 1 (1 DCU) | CPU 32 / HIP 2 (2 DCUs) | CPU 32 / HIP 4 (4 DCUs) |
| ---: | ---: | ---: | ---: |
| 65,536 | 1.280835x | 1.293836x | 1.217872x |
| 262,144 | 2.708388x | 3.401391x | 3.697872x |
| 1,048,576 | 5.915925x | 9.802711x | 13.472202x |
| 4,194,304 | 7.591373x | 14.146178x | 25.490468x |

代表性原始 lifecycle（ms）：

| global_nx | CPU 32 ranks | HIP 1 rank | HIP 2 ranks | HIP 4 ranks |
| ---: | ---: | ---: | ---: | ---: |
| 65,536 | 71.875828 | 56.116362 | 55.552503 | 59.017556 |
| 262,144 | 256.103894 | 94.559520 | 75.293881 | 69.257096 |
| 1,048,576 | 1,501.621673 | 253.827033 | 153.184324 | 111.460744 |
| 4,194,304 | 7,284.764827 | 959.610975 | 514.963461 | 285.783883 |

## 3. 结论与边界

1. 32 CPU cores 对 1 DCU 的严格同全局规模对照为
   `1.280835x / 2.708388x / 5.915925x / 7.591373x`；
   对 4 DCUs 为 `1.217872x / 3.697872x / 13.472202x / 25.490468x`。
2. 小规模 `65,536` 已接近固定初始化/启动/halo 成本，1–4 DCU 的变化不呈现
   单调加速；不能用小规模结果外推大规模吞吐。
3. `4,194,304` 上，CPU 32 ranks 对 1/2/4 DCUs 分别为 `7.59x/14.15x/25.49x`；
   这是 standalone 一维 Euler 生命周期比值，包含初始化、通信、推进和收尾。
4. 这组数据只证明 standalone 1D rank-local partition 和 HIP MPI device mapping
   已实际运行。3D 主 solver 的 `m6wingroe_sa` 仍是 single-zone，当前 HIP capability
   要求 `ZoneState::nLocal == 1`；3D 多 rank 实验必须先准备已分区 multi-zone fixture，
   按 rank 数重建 owner/interface plan，并验证每个 rank 的非空 zone 和 solver work。

本轮为单次矩阵作业，适合观察 rank/DCU 规模规律，不替代多轮 paired timing
或完整 3D 主 solver 性能报告。
