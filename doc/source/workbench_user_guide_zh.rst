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
参数统一视图（Parameter Workbench）
----------------------------------

打开项目后，ProjectWindow 会根据 ``script/control.txt`` 指定的文件顺序
构建统一参数树。参数树分为 **Shared**、**CFD**、**Grid** 和 **Other** 四个
顶层分类，不再要求用户逐个文件寻找参数。

每个已识别的标量声明显示参数名、原始值、声明类型和来源文件。来源列同时
显示它在 control manifest 中的文件顺序，以及在文件中的声明顺序。
``control.txt`` 只负责列出参数文件及其加载顺序，不作为普通参数显示。

同名参数遵循后声明覆盖先声明的规则。树中显示最终有效值，并在其下方
保留标记为 ``[overridden]`` 的历史声明，包括原值和来源文件。因此，用户
既能看到运行时应采用的值，也能追查该值覆盖了哪些早期声明。

暂时无法归类的标量参数仍然显示在 Other 或相应任务域的 Other 分组中。
数组声明以及标量解析器暂不支持的声明语法，会列在 Other / Unsupported
declarations 下，并显示来源行和不支持原因；不会被伪装成已完整识别的参数。

当前统一参数树是只读视图。选中参数后点击 **Open source file...**，或双击
参数行，可使用配置的 VS Code、Visual Studio 或其他外部编辑器打开原始文件。
参数文件仍是唯一权威来源；SQLite 参数目录仅用于索引和保存内容快照。

注意：参数能够被解析，不代表它一定被当前求解路径使用，也不代表其单位、
适用模型和数值范围已经确认。本切片尚未完成真实运行时消费者追踪，也没有
建立参数与 Formula & Knowledge 的交叉引用。这两项将在下一切片进行，
届时会以代码读取路径、适用条件和公式依据作为证据，不根据参数名称猜测。
