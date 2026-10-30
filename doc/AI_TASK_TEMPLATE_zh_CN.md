# OneFLOW AI 辅助工程任务规范

本文规定 OneFLOW 中较大规模的 AI 辅助开发与人工工程任务的标准表达方式，使需求具备可执行性、可审查性、可验证性和可复现性。本文应与 `doc/source/architecture_review.rst` 及相关子系统路线图配合使用。

## 1. 项目上下文与基准版本

每轮任务开始前，核实并记录：

| 字段 | 要求 |
|---|---|
| 仓库 | `eric2003/OneFLOW` |
| 目标分支 | 默认 `hexin`，除非明确指定其他分支 |
| 基准版本 | 获取远程 `hexin` 最新 HEAD，记录完整 commit SHA |
| 本地验证基线 | 记录用户报告的测试结果及环境，并标记为修改前证据 |
| 开发环境 | 操作系统、编译器、工具链、构建配置及相关可选后端 |
| 工作范围 | Runtime、Mesh、数值方法、I/O、Workbench、测试、文档或跨模块主题 |

不得依据历史对话推定当前 HEAD。未经明确授权，不得以 `master` 作为开发目标。修改前的测试结果不能表述为修改后版本的验证结果。

## 2. 目标定义

每项任务应定义一个可观察、可验收的主要工程目标。目标应说明需要建立或验证的行为契约、架构边界、数值性质、可复现能力或文档能力。

不得仅以“继续重构”“优化代码”“增加测试”等活动性描述作为目标；应说明修改将降低何种风险、澄清何种职责、改善何种复用能力或补足何种验证证据。

## 3. 全局收益评估与优先级

每轮确定一个主要任务，不对无关模块平均分配精力。候选任务按以下维度评估：

1. **系统影响范围**：涉及的核心工作流、模块数量及其重要程度。
2. **风险降低程度**：能够预防的问题的严重性和发生可能性。
3. **架构杠杆效应**：对职责边界、依赖方向、复用能力和可测试性的改善。
4. **科学计算价值**：对守恒性、稳定性、精度、收敛性或物理模型有效性的影响。
5. **验证充分性**：能否通过有代表性的测试和可复现实验建立可信证据。
6. **实施成本与兼容风险**：开发工作量、接口迁移成本和潜在行为变化。

优先处理能够消除跨模块风险、建立可复用契约或显著改善核心工作流的任务。孤立的命名调整、纯格式整理和低收益辅助函数拆分，不应挤占高收益架构工作的资源。

整体优先级以 `doc/source/architecture_review.rst` 和对应子系统路线图为准。源码证据或测试结果发生实质变化时，应重新评估优先级。

## 4. 实施前调查

修改前，应基于目标版本核实：

- 实际入口、调用方、下游使用方及相关构建目标；
- 数据表示、对象所有权、可变性、生命周期和清理责任；
- 进程级、算例级、任务级及对象级状态边界；
- 正常执行、重复执行及异常失败路径；
- 现有单元测试、回归算例和 CI 平台限制；
- 相关设计文档、历史决策及接口兼容要求。

调查结论必须区分为**已证实事实、尚未验证的假设和拟议设计**。路线图或历史 AI 回复不能代替对当前源码的核实。

## 5. 实现原则

- 除非明确提出并通过测试验证行为变更，否则保持既有行为兼容性。
- 仅在建立真实契约、消除不合理依赖、明确所有权或形成有效验证边界时引入抽象。
- 优先采用职责内聚、输入输出明确、副作用可识别、失败语义清晰的增量修改。
- 根据领域需要区分语法解析、类型验证、领域语义解释、运行状态应用和实际执行。
- 明确资源所有权及释放顺序；算例级清理不得终结进程级资源。
- 对无法一次性迁移的旧接口，采用边界明确的兼容适配方式。
- 在职责调整时改善不准确的命名和注释，但不得以重命名代替架构改进。
- C++ 源文件仅使用 ASCII 字符。
- 不得仅为测试或演示方便，在 Python 或 GUI 中另行实现生产算法。

## 6. 验证要求

测试必须与所声称的结论相匹配：

| 证据层级 | 验证范围 |
|---|---|
| 单元 / 契约测试 | 局部不变量、输入输出规则和失败条件 |
| 组件测试 | 相邻生产组件之间的协作行为 |
| 集成测试 | 真实生产链路及相关状态变化和副作用 |
| 回归测试 | 代表性应用行为的兼容性 |
| 数值验证 | 一致性、守恒性、稳定性和收敛性等数学性质 |
| 平台 CI | 指定操作系统、编译器和配置下的构建与测试 |

必须说明每项测试验证了什么，以及不能证明什么。无关测试通过不能作为修改链路正确性的证据；测试数量本身不是质量指标。

数值算法应按其性质采用解析解、制造解、极限情形、守恒检查、网格或时间步加密、标准基准算例等验证方式。适用时记录误差范数、加密比、容差和理论假设。

只有在目标 commit 上实际观察到构建、测试、基准或 CI 结果后，才能报告其通过。用户报告的本地结果与工具观察到的 CI 结果必须明确区分。

## 7. 文档与科学依据

文档属于持续维护的工程成果，不应被视为开发结束后的附加工作。

### 源码注释与 API 文档

对稳定的公开接口和重要不变量使用简洁、结构化的注释。根据需要记录设计意图、参数含义及单位、所有权和生命周期、前置条件、后置条件、副作用、失败行为及兼容性约束。

仅在存在对应维护文档时使用 `@doc-algorithm`、`@doc-stage`、`@doc-theory` 等交叉引用标识。

评估使用 Doxygen 生成 API Reference 的可行性。Doxygen 与 Sphinx 的集成（例如 Breathe）应先在一个职责稳定的模块上进行试点。自动生成的 API 文档应补充 `doc/source/` 中经过编辑的正式文档，而不能取代后者。

### Sphinx 专题文档

架构设计依据、模块职责、依赖和生命周期图、数学模型、公式推导、假设条件、离散形式、边界处理、算法流程、数据流、几何与网格示意、实现映射、局限性、验证方法及已知失效模式，应保留在维护中的专题文档中。

公式必须定义符号、单位、符号约定、假设条件和适用范围。流程及依赖关系优先使用 Mermaid 或 Graphviz；几何、网格和数值结果可使用 SVG 或可复现生成的图像。图示应具有标题或图注，并说明其表达的关键关系。

### 参考文献与外部软件资料

数学或工程结论应优先引用原始论文、经典教材、标准和权威技术手册。引用用于说明理论来源和知识溯源，不等同于证明 OneFLOW 实现正确。

其他开源及商业 CFD/CAE 软件的理论手册、用户手册和开发手册，可用于术语对照、工作流比较及公开算法研究。必须区分公开理论、厂商公开声明、独立分析和 OneFLOW 已验证的实现；不得暗示掌握商业软件未公开的内部算法。

## 8. 文档、代码与测试追踪

重要算法和架构契约应建立可导航的关联链：

`理论 / 需求 -> 数学或行为契约 -> 设计文档 -> 源码接口 / 实现 -> 测试 -> 代表性验证证据`。

仅在确实改善追踪能力时引入稳定标识符。应先完成少量端到端示例，再决定是否建设更完整的元数据或自动检查系统。阶段性审查应检查失效文档链接、缺少说明的公开契约，以及名称或断言已与实现不一致的测试。

## 9. 范围控制与设计决策

每次修改均应明确包含范围、排除范围、兼容性要求、后续迁移事项，以及具有长期架构影响的决策所考虑的替代方案。

仅当决策存在实质性替代方案、长期影响或跨模块后果时，才创建 Architecture Decision Record（ADR）；不为每个局部修改创建 ADR。

## 10. 完成条件与交付记录

适用的完成条件包括：

- 实现满足目标并保持要求的兼容性；
- 测试覆盖相应层级的行为契约；
- 相关设计文档、公式和图示已更新；
- 源码注释、自动生成 API 文档与专题文档不存在语义冲突；
- 构建、测试、数值验证和 CI 结果关联到明确版本；
- 未解决风险和后续工作已明确记录；
- 已重新审视本次修改对整体架构优先级的影响。

交付记录至少包括：基准及最终 commit SHA、问题定义和系统收益、修改文件及职责变化、测试契约、实际执行命令及结果、精确版本的 CI 链接、文档和参考资料更新、遗留风险及下一项最高收益任务。

必要的集成链路或平台验证未完成时，不得标记为 Complete。统一采用 **Complete、Partial、Blocked、Deferred** 四种状态，并给出证据依据。

## 11. 标准任务请求模板

```text
REPOSITORY: eric2003/OneFLOW
BRANCH: hexin
BASE: latest remote HEAD; record full SHA
LOCAL_BASELINE: <user-reported result and environment>

PRIMARY_OBJECTIVE:
<observable engineering outcome>

SYSTEM_CONTEXT:
<affected workflow, subsystem, upstream/downstream dependencies>

PRIORITY_RATIONALE:
<system impact, risk reduction, architectural/scientific leverage>

REQUIRED_INVESTIGATION:
<entry points, ownership, lifecycle, failure paths, build/test targets>

CONSTRAINTS:
- Preserve compatibility unless explicitly changed and tested.
- Keep C++ source ASCII-only.
- Do not introduce abstractions without a concrete contract and consumer.
- Do not claim unobserved build/test/CI results.
- Follow the current hexin branch and repository conventions.

IMPLEMENTATION_SCOPE:
<planned code changes and explicit exclusions>

VERIFICATION:
<unit/component/integration/regression/numerical/platform evidence>

DOCUMENTATION:
<doc/source pages, equations, diagrams, references, API comments>

ACCEPTANCE_CRITERIA:
<objective, testable completion conditions>

DELIVERY:
<commit SHAs, changed files, test evidence, CI links, limitations, next priority>
```

当前项目优先级及子系统状态以 `doc/source/architecture_review.rst`、对应路线图和 `doc/AI_CONTEXT.md` 为准。本文规范任务表达与执行方式，不替代对目标版本源码的检查。
