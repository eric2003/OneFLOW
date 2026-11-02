# OneFLOW 1D Euler DCU Stage 2 构建证据

更新时间：2026-09-04

## 已确认

- CPU-only 配置使用 GCC 13.3.0、CMake 3.28.3。
- 显式提供本地 `/usr/include` CGNS 和 METIS 路径后，主 `OneFLOW` 目标开始完整编译。
- 新增 `codes/accel/src/EulerDomainStateRegistry.cpp` 已被 `ConstructSolutionDirTree` 自动纳入主目标并成功编译。
- `EulerDomain` contract 独立测试：3/3 通过。
- `EulerDomainStateRegistry` 独立测试：3/3 通过。
- 既有 CPU CTest：7/7 通过。

## 当前环境阻塞

主目标在约 68% 处于既有 `codes/special/src/MpiTest.cpp` 停止，错误为：

```text
fatal error: mpi.h: No such file or directory
```

当前环境未发现 `mpicxx`、`mpirun` 或 `mpi.h`。因此本地不能宣称完整 OneFLOW 主目标链接成功；该阻塞与本轮 Euler domain contract/registry 无关。

## 未宣称的范围

- 尚未把新增测试注册进 CTest；既有 `tests/CMakeLists.txt` 的更新受 sandbox 文件更新 helper 限制。
- 尚未接入 `INIT_FLOWFIELD`、`TimeIntegral::RungeKutta`、主 solver residual/boundary/restart。
- 尚未新增 DCU 作业或目标节点证据。
