# OneFLOW 本地 MPI smoke 证据

更新时间：2026-09-04

环境：临时加载 `openmpi/4.1.4`，使用已构建的 CPU+MPI `OneFLOW`。

```text
mpirun --allow-run-as-root --oversubscribe -np 2 \\
  /tmp/oneflow-stage2-cpu-mpi-full-build/bin/OneFLOW
```

已观察到：

- 2 个 rank 成功初始化；
- rank 1/2 正确报告 `out of 2`；
- OneFLOW compute backend 正确选择 `CPU`；
- `OneFLOW multi-device context: enabled` 正常输出；
- 默认 HybridParallel/Jacobi 入口成功启动。

该 smoke 在默认 4096×4096 Jacobi 长计算阶段停止，未将其计为 1D Euler solver、MPI halo 或 DCU 验证。
