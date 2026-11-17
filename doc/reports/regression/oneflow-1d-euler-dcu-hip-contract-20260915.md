# OneFLOW 1D Euler DCU HIP contract evidence

日期：2026-09-15
目标：昆山 Z100，`kshdnormal`，`dcu:1`
范围：standalone 1D Euler HIP adapter；不代表 3D 主 solver 已完成 DCU 接入。

## 结果

- Slurm 资源预检查通过：8 CPU、27G、`dcu:1`。
- DTK `26.04`、clang/clang++ 17、HIP 架构 `gfx906`。
- CMake 配置通过：`config=0`。
- HIP/C++ 编译与链接通过：`build=0`。
- GoogleTest：9/9 通过。
- CTest：9/9 通过；统一 helper 生成 `hardware;hip;dcu` 元数据，Kunshan 的旧版 CTest
  兼容路径使用 `hardware` label 加 `HIP` 前缀筛选，9 个测试全部通过。
- contract 包含 Rusanov、生命周期、设备可见性、非法请求，以及 WENO5 与 CPU oracle 对照。
- 运行时实际识别 Hygon DCU `Device 66a1`，架构为 `gfx906`。

## 工具链闭环

Kunshan login 环境的系统 CMake 过旧，runner 在 module 初始化后使用 CMake 3.25；同时显式加载 GCC 9.3，并将 DTK 的 `amd_comgr`、`AMDDeviceLibs` 配置目录传递给 CMake HIP try-compile。C++ 与 HIP 编译均使用 GCC toolchain 和 `libstdc++`，避免 DTK clang 自动拾取不兼容的系统 GCC 头文件。

## 本轮修复

1. HIP runtime 头文件移到 `oneflow_1d` namespace 外，消除标准库命名空间污染。
2. HIP contract target 补齐 `OneDWeno5.cpp` 和 `OneDWeno5.hip`，闭合 CPU oracle 与 HIP kernel 的链接。
3. contract 测试改为首次使用 backend 时惰性初始化，消除静态初始化顺序问题。
4. WENO5 `HipState` 补齐 `left`、`right`、`residual` 设备 buffer。
5. WENO5 RK3 使用独立 `scratch` 保存原始 base，避免 current/next 与 RK base alias。
6. runner 的测试数量 gate 改为检查完整通过摘要，不再硬编码过时的 6 tests。
7. standalone 与根工程共用 HIP contract CMake 注册 helper；根工程通过
   `ONEFLOW_ENABLE_HIP_TESTS=ON` 显式开启，普通 CPU CTest 默认不依赖 DCU。

## 阶段边界

本证据完成阶段 F 的 standalone 1D HIP contract 子项。下一步仍是：在 E6 CPU batch/oracle 验收完成后，把同一 adapter 接入 3D 主 solver，并在昆山做主 solver correctness；MPI、多卡和四规模 WENO5 性能仍未完成。
