# OneFLOW

**面向网格生成、流场计算与可复现工程流程的开源 CFD 项目。**

OneFLOW 是一个持续演进的计算流体力学（CFD）代码库。它关注的不只是求解器能否完成一次计算，也包括计算所依赖的网格如何生成和检查、算例如何组织、不同计算之间如何传递成果，以及结果如何通过测试与回归得到验证。

项目的核心是 C++20 运行时；网格算例、命令行工具和可选的 Python 桌面工作台围绕运行时提供支持。我们希望逐步建立清晰的工程边界：**算例定义可追踪，运行产物与源文件分离，网格数据可独立检查，工作流依赖显式表达，软件行为由测试验证。**

OneFLOW 仍在积极开发中。本文尽量以当前分支中的实现为依据，并明确区分可用功能、正在整合的能力和预留的扩展方向。这里的 OneFLOW 指本仓库中的开源 CFD 与科学计算项目，与同名的合同生命周期管理商业产品无关。

[English](./README.md) · [在线文档](https://oneflow-cfd.readthedocs.io/en/latest/) · [源代码](https://github.com/eric2003/OneFLOW) · [问题反馈](https://github.com/eric2003/OneFLOW/issues)

## 项目定位与长期愿景

OneFLOW 是一款面向通用非结构网格流动模拟的开源计算流体力学（CFD）软件，也是一套持续演进的科学计算引擎与工程研究平台。项目目标覆盖一维、二维和三维偏微分方程计算，服务于数值算法学习、验证和科研，并支持工程流动问题的建模、计算、检查与复现；具体可用的方程、网格类型和物理模型仍以当前版本中的求解器与算例实现为准。项目的长期方向，是将数值计算、网格处理、自动化测试、任务编排和结果检查连接成更可编程、更可验证的工作流程。

“科学计算操作系统”可以用来描述这一长期愿景：让不同计算任务、工具和硬件后端在清晰的接口与可追踪的流程中协同工作。它不是对当前已经具备完整操作系统式能力的声明。类似地，配置选项中出现 MPI、CUDA、HIP 或其他后端名称，也不等于所有后端均已达到生产成熟度，更不能单凭名称推断特定核数下的扩展性或某种硬件上的性能；使用前应核对实际实现、构建配置和对应验证证据。

## OneFLOW 正在构建什么？

- **CFD 运行时**：以 C++20 为基础，包含求解、任务、数据、几何、网格和输入输出等代码模块。
- **网格与几何工作**：提供网格生成和格式转换算例，处理网格连接关系、边界信息，并逐步增加可独立测试的网格摘要与检查工具。
- **可复现的算例验证**：通过 GoogleTest/CTest 单元测试、网格检查以及现有 CFD 数值回归，从不同层面验证代码。程序正常退出并不自动代表所有数值结果正确。
- **工程工作台**：可选的 Python/PySide6 桌面应用，用于打开项目、配置运行、编辑网格布局、组织多算例工作流，以及查看日志和运行历史。
- **计算后端演进**：通过 CMake 选项逐步整理计算后端边界。CPU 是当前常规开发与回归的基线；不同加速后端的实现和验证程度并不相同。

## 设计理念

### 1. 运行时是核心，工具围绕运行时工作

OneFLOW 的 C++ 程序负责真正的算例执行。Python 工具和桌面工作台用于发现算例、准备运行目录、组织任务和查看结果，不应成为命令行求解器的隐式依赖，也不应复制一套数值实现。

### 2. 算例定义与运行产物分离

源代码仓库保存可追踪的输入、布局、参数和测试基准；网格、结果、日志及临时文件属于一次运行产生的产物。运行工具应优先使用独立工作目录，避免在源码算例中覆盖输入文件或留下难以辨认的生成文件。

### 3. 网格数据、摘要和序列化各有边界

网格检查不应要求启动完整的图形界面，也不应把 JSON 格式绑定到运行时网格对象上。当前演进方向是先从真实网格中提取独立的摘要快照，再由独立层负责序列化和可视化。摘要可描述维数、节点/面/单元数量和包围盒等元数据，但它不是实际网格数据的替代品。

### 4. 让算例依赖显式可见

有些算例可以独立运行，有些算例则需要使用前一步生成的网格或其他产物。工作流模型用执行模式和依赖边描述这些关系，便于验证依赖、安排执行顺序和传递产物路径。部分工作流入口仍在整合中，具体可用行为以工作台文档为准。

### 5. 用测试证据描述成熟度

单元测试、网格拓扑检查、CFD 数值回归和加速器硬件测试回答的是不同问题。它们不能互相替代，也不应把 CPU 测试通过描述为 GPU 后端已验证。项目倾向于通过小步改进、明确接口和可重复测试推进架构，而不是为了重构而重构。

## 架构概览

OneFLOW 不只是一个求解器可执行文件再加上一组互不相关的脚本，而是一个逐步建立清晰边界的工程代码库。

- **数值运行时**：C++ 程序负责真正的算例执行和数值计算。命令行运行不应依赖先启动桌面工作台。
- **网格数据流水线**：区分运行时中的真实网格对象、可独立测试的稳定摘要快照，以及 JSON 等输出格式。摘要记录维数、实体数量、包围盒等元数据，并不等同于网格本身。几何、拓扑、连接关系和边界信息应通过职责明确的接口逐步开放，而不是全部塞进一个通用导出器。
- **算例与产物边界**：受版本控制的是算例定义和参考输入；生成网格、流场结果、日志和报告属于一次运行的产物。显式的输入/输出根目录使前置算例和后续算例可以交换成果，而不依赖偶然的当前工作目录。
- **工作流编排**：算例可以单独运行、按顺序运行、彼此独立地运行，也可以组成有依赖关系的工作流。依赖应明确表达由哪个算例生成、哪个算例消费的产物，并在执行前验证，而不是隐藏在人工操作步骤里。
- **工具与呈现层**：命令行工具、Python 服务和可选的 PySide6 工作台负责发现与编排任务、检查输出和呈现结果。它们应尽可能复用共同的算例/工作流契约，而不是复制求解器行为，也不应让运行时依赖 GUI 代码。
- **可测试的渐进演进**：单元测试、网格不变量检查和数值回归覆盖不同契约。重构应沿真实调用路径推进，让对象所有权和接口边界更清楚，并在需要时保持兼容；测试用于验证行为，不意味着为了重构而重写代码。

这些内容描述的是架构边界和持续推进的方向，并不表示所有工作流或网格接口都已达到最终形态。例如，算例目录和 profile 目前还不是所有旧式运行时参数的唯一权威来源；迁移边界见[网格说明](./grid/README.md)。关于运行时入口、参数配置层次、任务分发和重构约定，参见[运行时架构与开发指南](./doc/source/runtime_architecture.rst)。

## 代码与文档导航

| 路径 | 内容 |
| --- | --- |
| `codes/` | C++ 运行时、求解器、几何、网格与输入输出实现 |
| `ThirdParty/` | 构建使用的第三方依赖，包括 GoogleTest 和 JSON 相关代码 |
| `tests/` | 由 CMake/CTest 管理的 C++ 单元测试 |
| `grid/` | 网格生成和网格转换算例、配置档案与算例目录 |
| `test/` | CFD/网格回归运行器、测试清单和数值基准 |
| `tools/` | 算例检查、网格工具和可选桌面工作台 |
| `doc/` | 正式文档根目录：架构、计划、操作手册、验证证据、报告和 Sphinx 文档源文件 |
| `docs/` | 现有辅助笔记；新增项目文档统一放在 `doc/` 下 |
| `cmake/` | CMake 辅助函数 |

建议从[在线文档](https://oneflow-cfd.readthedocs.io/en/latest/)开始，再阅读仓库内的专题说明：

- [网格算例说明](./grid/README.md)
- [测试与数值回归说明](./test/README.md)
- [桌面工作台说明](./tools/oneflow_gui/README.md)
- [文档资产盘点与整合计划](./doc/plans/oneflow-documentation-consolidation-plan.md)

## 从源码构建

### 环境要求

- CMake **3.30 或更新版本**。
- 支持 C++20 的编译器。
- Git。
- 默认构建配置使用 MPI、METIS 和 CGNS 开发库；CGNS 通常还需要 HDF5。
- Python 3 可用于算例工具、回归测试和可选桌面工作台。

不同系统的依赖安装路径不完全相同。目前项目没有统一到单一的跨平台包管理方案，因此配置后应检查 CMake 输出中的头文件和库路径。以下命令是可调整的起点，而不是对所有系统都适用的固定路径。

### Windows：Visual Studio 2022

安装 Visual Studio 2022 的 **Desktop development with C++** 工作负载，以及 CMake 和 Git。还需要安装兼容的 MPI、METIS 和 CGNS 库。

当前 CMake 默认路径包括 MS-MPI 的标准 SDK 目录、`C:/dev/METIS/METIS-VS2022-STATIC` 和 `C:/dev/cgns/4.4.0`。如果库安装在其他位置，请通过 CMake 参数覆盖这些路径。

在仓库根目录使用 PowerShell：

~~~powershell
cmake -S . -B build `
  -G "Visual Studio 17 2022" -A x64 `
  -D CMAKE_INSTALL_PREFIX="C:/dev/OneFLOW" `
  -D METIS_ROOT="C:/dev/METIS/METIS-VS2022-STATIC" `
  -D CGNS_ROOT="C:/dev/cgns/4.4.0"

cmake --build build --config Debug --parallel
ctest --test-dir build -C Debug --output-on-failure
~~~

构建 Release 版本时，将构建和测试命令中的配置都改为 `Release`。Visual Studio 是多配置生成器，应通过 `--config Debug` 或 `--config Release` 选择配置，而不是依赖 `CMAKE_BUILD_TYPE`。

可执行文件通常位于 `build/bin/Debug/` 或 `build/bin/Release/`。如果依赖库安装在其他目录，请相应调整 `METIS_ROOT`、`CGNS_ROOT` 以及 MPI 的 include/library 配置。

### Windows：Ninja 单配置构建

如果已经准备好 Ninja 和 MSVC 开发环境，也可以使用单配置构建：

~~~powershell
cmake -S . -B build-ninja -G Ninja `
  -D CMAKE_BUILD_TYPE=Debug `
  -D CMAKE_INSTALL_PREFIX="C:/dev/OneFLOW" `
  -D METIS_ROOT="C:/dev/METIS/METIS-VS2022-STATIC" `
  -D CGNS_ROOT="C:/dev/cgns/4.4.0"

cmake --build build-ninja --parallel
ctest --test-dir build-ninja --output-on-failure
~~~

切换生成器或编译器时请使用新的构建目录；不要在已有构建目录中直接切换 CMake generator。

### Linux：GCC 与系统依赖

Debian/Ubuntu 系统可先安装常用依赖：

~~~bash
sudo apt update
sudo apt install build-essential cmake git \
  libhdf5-dev libcgns-dev libmetis-dev libmpich-dev
~~~

当前 Linux CMake 配置可通过环境变量指定 MPI、CGNS 和 METIS 的头文件及库文件路径。以下示例适用于常见的 Ubuntu/MPICH 安装布局；请根据实际系统修改：

~~~bash
export MPI_HOME_INC=/usr/include/x86_64-linux-gnu/mpich
export MPI_HOME_LIB=/usr/lib/x86_64-linux-gnu/libmpich.so
export CGNS_HOME_INC=/usr/include
export CGNS_HOME_LIB=/usr/lib/x86_64-linux-gnu/libcgns.so
export METIS_HOME_INC=/usr/include
export METIS_HOME_LIB=/usr/lib/x86_64-linux-gnu/libmetis.so

cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX="$PWD/Run"
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
~~~

如果系统将 MPI、HDF5 或其他库安装到不同目录，应调整对应路径。如果配置或链接失败，先检查 CMake 实际检测到的 include 和 library 路径，再判断是否需要修改代码。

## 运行一个算例

OneFLOW 命令行入口采用以下形式：

~~~text
OneFLOW.exe r <case-directory>
~~~

例如，在 Windows 的 Debug 构建完成后，可以从仓库根目录执行：

~~~powershell
.\build\bin\Debug\OneFLOW.exe r .\grid\minimal_rectangle
~~~

如果算例会生成网格、结果或日志，请在可写的工作目录中运行，而不要直接在受版本控制的源算例目录里执行。桌面工作台对支持的任务会创建独立运行副本。算例输入和预期输出请参阅具体算例目录及[网格算例说明](./grid/README.md)。

## 测试与验证

### C++ 单元测试

CMake 将 C++ 单元测试注册到 CTest。Visual Studio 多配置构建可运行：

~~~powershell
ctest --test-dir build -C Debug --output-on-failure
~~~

单配置构建不需要 `-C Debug`。定位失败时，除了完整 CTest，还应运行相关测试可执行文件或单个 GoogleTest 用例，以区分测试隔离问题、组件契约问题和整体回归问题。

### 网格检查与 CFD 数值回归

`test/` 目录维护回归运行器、测试清单和已接受的数值基准；`grid/` 目录维护网格算例目录和配置档案。详细命令、输出检查方式和容差策略见[测试说明](./test/README.md)和[网格说明](./grid/README.md)。

这些验证层次各有用途：

- **单元测试**：验证类、函数和模块的局部契约。
- **网格检查**：验证生成网格的维数、坐标、包围盒和拓扑等属性。
- **CFD 数值回归**：将求解结果与已接受的参考数据比较。
- **加速器验证**：需要匹配的编译器、运行时和硬件，不能由 CPU 测试结果代替。

## 可选桌面工作台

桌面工作台基于 Python 和 PySide6，负责项目与工作流的组织、运行配置和结果查看；它不是构建或运行 C++ 命令行求解器的必需组件。

~~~powershell
python -m pip install -r tools/requirements-gui.txt
python -m tools.oneflow_gui
~~~

工作台当前涵盖项目打开与创建、运行配置、网格布局编辑、测试入口、多算例/依赖工作流、运行日志和历史记录。功能仍在整合，具体行为和已知边界请查看[工作台文档](./tools/oneflow_gui/README.md)。

## 计算后端现状

CMake 提供 `ONEFLOW_ACCEL_BACKEND` 选项，取值包括 `AUTO`、`CPU`、`HIP`、`CUDA` 和 `KOKKOS`。这些名称代表配置入口或扩展方向，并不意味着各后端的实现成熟度相同。

- **CPU**：当前常规开发与回归的基线。
- **HIP**：在 HIP 软件包、编译器和目标架构均可用时，可显式启用；实际加速器验证需要匹配的硬件环境。
- **CUDA**：沿用可选的旧 CUDA 源码路径，需要显式设置 `ONEFLOW_ENABLE_CUDA=ON`，默认并未开启。
- **Kokkos**：当前是预留扩展点，尚不是已启用的移植实现。

CPU 配置示例：

~~~bash
cmake -S . -B build-cpu -D ONEFLOW_ACCEL_BACKEND=CPU
~~~

如果要构建 HIP 版本，请结合 CMake 的检测结果配置编译器和目标架构；交叉编译时尤其需要明确指定架构。只有在工具链和依赖确实可用时，才应启用对应后端。

## 参与开发

欢迎围绕真实问题提交改进，尤其是边界清晰、可复现并有测试支撑的变更。建议：

1. 说明修改影响的运行时模块或工具边界。
2. 为新行为或关键契约增加、更新测试。
3. 避免将生成文件写入受版本控制的源算例。
4. 运行相关单元测试或回归，并如实记录平台和硬件限制。
5. 当命令、参数或可观察行为变化时同步更新文档。

OneFLOW 正在持续整理网格、几何、测试和工作台的架构。我们更看重能改善真实调用路径、让行为更容易理解和验证的渐进式改进，而不是没有明确需求的大规模重写。

## 开发者与联系方式

**主要开发者：** Dr. He Xin（GitHub：[@eric2003](https://github.com/eric2003)）。

- **电子邮箱：** [fantasy_2003_@hotmail.com](mailto:fantasy_2003_@hotmail.com)
- **源代码仓库：** [github.com/eric2003/OneFLOW](https://github.com/eric2003/OneFLOW)
- **问题反馈与功能建议：** [GitHub Issues](https://github.com/eric2003/OneFLOW/issues)
- **架构与实现讨论：** 建议提交 issue，并附上相关算例、预期行为、实际行为，以及可复现的日志或测试结果，方便其他贡献者了解背景并持续跟进。

该邮箱已在 OneFLOW 早期版本的 README 中公开列出；也可以通过 GitHub Issues 反馈问题和讨论设计。

## 项目与许可

OneFLOW 由 He Xin 和贡献者共同维护。使用、分发或贡献代码前，请以仓库中的许可文件和贡献说明为准。

---

*OneFLOW 是一个持续演进的开源 CFD 项目。本文以仓库当前实现和已知限制为依据，不代表所有规划能力都已完成。*
