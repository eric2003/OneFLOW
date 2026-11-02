OneFLOW Workbench 用户手册
==========================

本手册介绍 tools/oneflow_gui 中的可选桌面工作台。Workbench 负责项目浏览、运行编排、测试入口、网格布局编辑和运行记录；OneFLOW 求解器仍作为独立进程运行，GUI 不是求解器的替代品。

本页按 hexin 当前实现编写。按钮细节、模板列表和完整服务测试命令见模块 README：
https://github.com/eric2003/OneFLOW/blob/hexin/tools/oneflow_gui/README.md
设计边界见 doc/source/workbench_architecture.rst。工作台能力随代码演进，未实现的集成点会在本页明确标出。

启动
----

在仓库根目录使用 Python 3.10 或更高版本：

.. code-block:: powershell

   python -m pip install -r tools/requirements-gui.txt
   python -m tools.oneflow_gui

首次使用时打开 **Configuration...**，设置 OneFLOW 仓库根目录、可执行文件和 Debug/Release 配置。程序会按约定位置搜索可执行文件；也可以直接输入路径或使用 **Browse...**。完整搜索顺序见 GUI README。

**Run mode** 设置保存在用户设置文件中。选择 **Multi-case run** 或 **Dependency workflow** 后打开 **Saved Workflows...**，会自动创建对应模式的新工作流草稿；**Run single case** 不会创建工作流草稿，单算例从 **Project > Open Project...** 进入。

参数统一视图
------------

打开项目后，ProjectWindow 根据 control.txt 中列出的文件顺序读取参数文件，显示 Shared、CFD、Grid 和 Other 四个分类。control.txt 只作为加载清单，不会作为普通参数显示。参数行保留原始值、声明类型、来源文件、manifest 文件顺序和文件内声明顺序。数组声明和标量解析器暂不支持的声明形式会显示在 Other / Unsupported declarations 中，并附带来源行和限制说明；该扫描不代表已识别全部历史语法。

同名参数按后声明覆盖先声明处理；有效声明下方保留所有被覆盖声明及其来源。部分已核对的参数还会显示 Runtime consumer 子项，包括消费函数、读取时机和 C++ 源文件。消费者目录目前是经过源码核对的有限清单；显示 Not audited 或未列出消费者表示尚未追踪，不表示参数无用。参数视图只读，双击参数可在配置的外部编辑器中打开原始文件。 对于已经建立公式知识关联的参数，还会显示 Formula & Knowledge 条目 ID；例如 time_integral 和 rk_stage 关联 time_integral_rk_stage_weight。该关联用于导航，不代表时间离散阶数或稳定性已经得到验证。

项目和单算例
------------

**Project > Open Project...** 可以打开任意已有目录，不要求它在内置算例目录中登记。**New Project...** 从模板创建工作副本。项目页可选择可执行文件和执行位置：

* **Run in original project directory** 会直接在项目中运行；OneFLOW 可能创建或覆盖文件。OneFLOW 源码仓库内的项目不能用此方式运行。
* **Copy to OneFLOW default workspace** 将项目复制到当前用户的数据目录后运行。
* **Copy to specified directory** 将项目复制到指定目录后运行。目标目录必须位于源码仓库之外。

普通单算例运行使用 OneFLOW.exe r <case-directory>。若项目含 workflow.json，项目运行器会按其依赖图执行；每个 case 使用隔离运行目录，依赖和输出文件由工作流执行服务检查。

多算例和依赖工作流
------------------

从 **Project > Saved Workflows...** 创建、编辑或运行已保存的工作流：

* **Multi-directory** 把目录按界面顺序传给同一个 OneFLOW 进程，例如 OneFLOW.exe r <dir-1> <dir-2>。进程内逐个运行 case；重复目录也会再次执行。
* **Dependency graph** 根据 case 依赖顺序运行，先将 case 准备到本次运行目录，再把声明的输出作为后续 case 的输入。前置 case 失败或缺少所需产物时，依赖 case 不会继续；无依赖的 case 可按图执行策略独立完成。

**Open/Create Multi-Case Workspace...** 编辑版本 2 的 .ofw 定义。旧版不含 artifact 绑定的工作区会在保存时迁移；点击 **Save and Run...** 会保存定义并打开对应的 Saved Workflows 草稿；确认可执行文件和路径后，点击 **Run workflow** 启动。单算例、顺序和独立模式使用一个 OneFLOW 进程并保留目录顺序。依赖模式会暂存各 case，并按依赖顺序执行；前置步骤失败会阻止依赖 case。依赖模式的每个 case 可用 JSON `outputs` 声明产物，用 `inputs` 将产物绑定到直接依赖 case 的参数文件和设置。产物和参数文件路径均相对其所属 case，例如：

.. code-block:: json

   [
     {"case_id": "grid", "path": "grid", "outputs": [{"id": "mesh", "path": "results/mesh.ofl", "parameter_file": "script/grid.txt", "setting": "targetGridFileName"}]},
     {"case_id": "solver", "path": "solver", "depends_on": ["grid"], "inputs": [{"artifact": "mesh", "setting": "gridFileName", "parameter_file": "script/cfd.txt"}]}
   ]

运行定义和历史记录由现有工作流仓库保存。

测试
----

首页 **Tests...** 提供三类入口：

* **Unit Tests** 在源码构建目录运行 CTest 和 Python 单元测试。
* **CFD/Grid Regression** 使用仓库已有的回归清单和脚本，并在用户工作目录准备输入。
* **Development Tests** 是可选项，默认不运行；选择后默认勾选目录清单中的所有 case，也可以打开 **Select cases...** 调整。

Development Tests 使用固定的用户数据测试工作目录，并在每批开始时重新准备该目录中的测试 case；源仓库内的内置 case 保持不动。测试进程的退出状态会记录为运行结果；这不等同于所有网格几何或数值结果都已由 GUI 独立验证。

Contract Workbench
------------------
Case 领域还提供两个文件入口：**Open case directory...** 打开算例目录；
**Open suite and references...** 打开 Suite 清单、残差基线和保留的解析参考数据。
这些操作只定位已有文件，不会启动求解器，也不会修改基线。


首页 **Contract Workbench...** 是 OneFLOW 数值知识和工程契约的统一入口。当前提供四个领域入口：

* **Formula & Knowledge**：打开现有 Verification Curves 浏览器，查看公式、原理、参考文献、文档摘录、代码草图和理论曲线。
* **Geometry**：打开 Geometry Contract Explorer，查看拓扑场景、输入/期望表和 fixture 结果。Fixture 运行不等同于执行生产 C++。
* **Case**：浏览类型化 CFD Case 契约，查看算例所属 Suite、求解器、物理模型、数值格式、输入网格、输出比较清单、残差基线和分析参考数据。
* **Runtime**：浏览进程生命周期、资源所有权、状态转换、实现代码和测试证据。

Formula、Geometry 继续使用已有领域工具；Case 和 Runtime 使用统一 Contract Catalog 的只读视图。各领域共享 Overview、Principle、Architecture、Implementation / Code、Documentation / References、Verification 等导航，但领域数据保持类型化，不要求所有领域都具有公式或曲线。

Contract Workbench 展示证据声明和验证范围，不负责执行求解器或测试，也不会因为关联了一个测试名称就显示 PASS。Case 中保留的 Blasius 参考数据如果尚未被回归程序自动比较，会明确标记为“已声明、未自动比较”。Workflow 和 Performance 入口将在对应契约和证据模型有实际内容后再加入。

网格布局编辑器
--------------

首页 **Mesh Layout Editor...** 可编辑内置单矩形、垂直中点矩形和双块平板模板。通过 **Save specification** 保存布局定义，**Update preview** 查看 SVG 预览，**Compile layout** 生成布局文件。首次打开模板时会复制到用户数据目录中的网格布局工作区；再次打开复用这份用户副本。**Load JSON...** 可以打开其他目录中的用户布局。源码树中的文件不能通过该编辑器写入。

运行记录与用户数据
------------------

**Project > Run History...** 查看运行状态、日志和记录；删除运行记录或清理历史不会删除运行产物。**Background Tasks...** 查看、停止或打开后台任务结果。退出时如果仍有活动任务，Workbench 会询问如何处理。

默认用户数据目录是 Windows 的 %USERPROFILE%\.oneflow，Linux 的 $HOME/.oneflow。设置 ONEFLOW_HOME 可将其移动到其他位置。设置、工作区定义、工作流定义、运行历史、日志、后台任务状态、网格编辑器项目和隔离运行目录都保存在用户数据树中；它与源码仓库分开。GUI 启动的 Python unittest 命令仍从源码根目录运行，以保持模块导入和仓库相对路径语义；其执行日志写入 test-runs/logs，执行证据写入 verification/execution-records.jsonl。临时测试夹具使用系统临时目录或显式隔离工作区，不应在源码树留下运行产物。

架构和当前边界
--------------

Workbench 是 OneFLOW 的可选 Python/Qt 客户端，通过参数数组启动求解器进程。case、工作流、运行记录和任务结果等服务尽量保持 Qt 无关并有独立测试。项目发现、case catalog、工作流定义、运行副本和运行历史各自有不同职责。

后续可完善项是：让工作流定义与 .ofw 定义共享一个持久化模型，减少两种格式之间的重复。详细架构约束见 doc/source/workbench_architecture.rst；GUI 模块、服务和测试索引见 tools/oneflow_gui/README.md。

参数运行时消费者追踪：第一批证据
--------------------------------

新增 parameter_consumer_catalog.py，作为 Qt 无关的运行时消费者证据目录。它不是全仓库自动扫描器，而是经过当前 hexin C++ 源码核对的一组明确记录。每条记录包含参数名、源文件、消费函数、读取表达式、读取时机、适用条件以及行为说明。

第一批核对的参数包括：

* maxSteps：Iteration::Init 必需读取；在 iexitflag == 0 时，SimuIterState::Running 使用它作为外迭代终止上限。
* time_integral：分别在 Ctrl::Init 和 NsCtrl::Init 中读取，属于不同控制初始化路径。
* mach_ref：由 NsCom::Init 读取，并与 gama_ref 一起用于计算 const_cp。
* vismodel：由 VisModel::Init 读取，后续影响黏性模型检查和壁面距离处理分支。
* ptTol：不是简单的启动时立即读取，而是在首次调用 Tolerance::GetTol 时通过 Tolerance::Init 延迟读取。
* nmg：由 MG::Init 读取；当值不大于 1 时，程序会调整前、后平滑步数。

读取时机区分为 required、optional 和 lazy。参数不在当前目录中，只表示尚未完成追踪，不能据此判定参数无用。代码中存在读取语句，也不代表该参数对所有任务都生效，更不代表其物理意义、推荐范围和数值效果已经验证。

新增独立单元测试，覆盖大小写无关查询、多处消费者、延迟读取和未知参数边界。后续逐步扩充时，必须先核对 hexin 中的真实代码路径，再补充适用模型、公式依据和验证证据；不能仅凭参数名称推断公式。
