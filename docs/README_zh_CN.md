# OneFLOW 文档索引

[English](README.md) | [简体中文](README_zh_CN.md) | [项目主页](../README.md) | [在线文档](https://oneflow-cfd.readthedocs.io/en/latest/)

`docs/` 是 OneFLOW 的文档根目录，包含工程记录和正式 Sphinx 手册源文件。两者用途不同：工程记录服务于实现、评审、交接、可复现操作与验证；`source/` 面向读者介绍稳定概念。并非每一份工程记录都需要发布到文档网站。

公开项目文档和面向贡献者的规范以英文版为基准。本页是中文配套入口：应保持导航、技术含义和治理规则一致，同时可以保留有价值的中文补充说明，不要求逐句翻译所有工程记录。中英文内容发生冲突时，以 [英文版](README.md) 为准，并据此修订本页。

## 1. 目录结构与保留策略

| 目录 | 内容 | 维护策略 |
|---|---|---|
| [`architecture/`](architecture/) | 稳定的架构说明、设计边界与审查结论 | 每个主题只保留一份现行文档；被取代的审查文档删除，Git 历史即档案 |
| [`plans/`](plans/) | 持续维护的路线图和长期计划 | 只保留持续维护的文档，不添加日期后缀；完成的一次性清单应删除；上游架构快照按下文映射规则处理 |
| [`handoff/`](handoff/) | 持续更新的项目交接 | 就地维护唯一的 `oneflow-current-handoff.md`，不创建带日期的交接副本 |
| [`runbooks/`](runbooks/) | 可重复执行的构建、集群、移植和运维流程 | 每个流程维护一份操作手册，就地更新 |
| [`evidence/`](evidence/) | 可复核的环境与构建事实 | 证据应绑定具体结论；结论被取代后同步清理 |
| [`reports/`](reports/) | 正式性能、回归和交付报告 | 不删除已发布报告；过时报告标记为历史版本，Markdown 与 HTML 同名文件保持同步 |
| [`source/`](source/) | Read the Docs 使用的 Sphinx 手册源文件 | 按面向读者的手册结构维护 |
| `_drafts/` | 未定稿草稿 | 评审后转入正式目录或删除，不作为长期档案 |

文档命名规则见 [DOCUMENT_NAMING.md](DOCUMENT_NAMING.md)。写作前先确定目录和文档类型；新增、删除或重命名文档后，必须同步更新本索引和受影响的引用。

## 2. 当前入口

| 主题 | 文档 |
|---|---|
| Phase A–F 路线图与稳定编号 | [开发计划](plans/oneflow-development-plan.md) |
| 当前工作切片、验收步骤与证据 | [开发 TODO](plans/oneflow-development-todo.md) |
| 持续交接：环境、固定路径、启动顺序与问题记录 | [当前交接](handoff/oneflow-current-handoff.md) |
| 昆山工作区固定路径、作业规范与清理规则 | [昆山工作区操作手册](runbooks/oneflow-kunshan-workspace-runbook.md) |
| DCU 移植方法：契约 → oracle → 门禁 → benchmark | [DCU 移植手册](runbooks/oneflow-dcu-porting-playbook-20260903.md) |
| 主 solver 加速器架构与接入边界 | [DCU 架构审计](architecture/oneflow-1d-euler-dcu-architecture-audit-20260904.md) |
| 正式报告列表 | [报告索引](reports/README.md) |

## 3. Fork 与上游文档归属

仓库中既有 fork 自行维护的文档，也有与上游共享的文档。修改前请遵循以下归属规则：

- **Fork 自行维护：** `plans/oneflow-development-*.md`、`handoff/`、本索引，以及 `DOCUMENT_NAMING.md` 中 fork 新增的内容。
- **与上游共享：** `architecture/`、`evidence/`、`reports/`、`runbooks/` 中的操作文档以及 `source/`。可以改进内容，但应避免大范围移动或重写，以免产生结构性合并冲突。确需移动时，应作为独立且明确的变更，并同步修复 `ports/`、`ci/` 等位置的引用。
- 不要在本 fork 中重新创建上游的 `doc/reports/architecture/` 目录；fork 已将相关文档映射到现有目录结构。

合并上游 `doc/reports/architecture/` 中的新文件时，按以下规则映射，并在同一个工作切片中完成移动和引用更新：

| 上游路径 | Fork 目标目录 |
|---|---|
| `reports/architecture/oneflow-dcu-porting-playbook-*.md` | `runbooks/`（保留文件名） |
| `reports/architecture/oneflow-gpu-backend-delivery.*` | `reports/delivery/`（保留文件名） |
| `reports/architecture/oneflow-euler-optimization-plan.md` | `plans/`（保留文件名） |
| `reports/architecture/<other>.md`（设计说明或快照） | `plans/`（保留文件名） |

## 4. 阅读顺序与技术入口

进行整体架构审查时，先阅读[系统架构审视](source/architecture_review.rst)，再阅读[运行时架构指南](source/runtime_architecture.rst)。之后根据任务继续阅读[网格架构](source/mesh_architecture.rst)、[网格拓扑](source/grid_topology.rst)或[Workbench 架构](source/workbench_architecture.rst)。

规范性的 [OneFLOW 持续演化与科学验证开发规范](AI_TASK_TEMPLATE.md) 定义了工作主线优先级、实现约束、科学验证与确认、参考文献、可复现绘图和交付证据；其[中文版](AI_TASK_TEMPLATE_zh_CN.md)作为配套版本维护。

需要快速恢复持续开发上下文时，请阅读 [AI 项目状态与交接页](AI_CONTEXT.md)。该页面汇总分支约束、已检查的实现基线、已确认行为、尚待验证的缺口和下一步建议。它是导航入口，不能替代针对目标提交的源码检查和验证。

现有文档的初步盘点与整合顺序见[文档整合计划](plans/oneflow-documentation-consolidation-plan.md)。计划只记录预期归属和步骤，并不意味着其中列出的迁移已经全部完成。

新增工程文档统一放在 `docs/` 下，遵循[文档命名规范](DOCUMENT_NAMING.md)，适用于 `architecture/`、`plans/`、`handoff/`、`runbooks/`、`evidence/`、`reports/` 和 `_drafts/`。不要批量重命名现有文件，以免破坏历史链接。

建议按以下顺序从系统整体阅读到实现细节：

1. 项目目标与文档边界：仓库[项目主页](../README.md)和本索引。
2. 系统架构：[运行时架构与开发指南](source/runtime_architecture.rst)。
3. 子系统架构：根据任务选择网格、拓扑、Workbench 或数值模块。
4. 算法与接口：数学定义、关键公式、数据流、代码入口和已知限制。
5. 验证证据：契约测试、集成测试、网格检查和 CFD 回归；区分已验证行为与仍待验证的假设。

架构审视页是检查清单，并不意味着其中讨论的所有能力都已实现。顶层文档应解释系统组成、分层理由、主要控制流与数据流以及风险；专题文档再解释算法、接口契约和实现位置，避免在系统概览中重复函数级细节。

## 5. 构建 HTML 手册

在仓库根目录执行：

```bash
cd docs
make html
```

生成结果位于 `docs/build/html/`。CI 和 Read the Docs 使用 `docs/source/` 作为 Sphinx 源目录。修改正式手册页面时，应检查 toctree、交叉引用和构建结果。只有计划公开发布的工程记录才应加入网站导航。

同时提交 Markdown 和 HTML 的正式报告必须使用相同文件名主体（例如 `oneflow-euler-performance-current.md` 与 `oneflow-euler-performance-current.html`），并在同一工作切片中同步更新。

## 6. 每个工作切片结束时的文档检查

1. 更新 `plans/oneflow-development-todo.md` 中的完成项、状态快照和本轮证据。
2. 就地更新 `handoff/oneflow-current-handoff.md` 的任务指针和新问题记录；不要新增带日期的交接文档，并清理已完成的一次性清单和草稿。
3. 新增、删除或重命名文档后，检查本索引、`docs/reports/README.md` 和 `AGENTS.md` 的链接。运行 `python3 scripts/check-doc-paths.py`，确保没有失效路径或过时的 `doc/` 引用。
4. 运行 `git diff --check`。不要在文档中发布原始 CI/Slurm 日志、主机名、账号标识、作业编号或私有绝对路径。

## 7. 文档层次与相关契约

文档树包含两个层次：

- **工程文档**服务于开发、评审、交接、可复现和验证，按前述生命周期目录规则维护。
- **`source/`** 是面向读者的 Sphinx 手册，应概括稳定概念并链接详细工程记录，而不是整段复制。

相关入口：

| 主题 | 文档 |
|---|---|
| Case 生命周期与增量验证契约 | [Case 生命周期契约](architecture/oneflow-case-lifecycle-contract.md) |
| 系统架构审视清单 | [架构审视](source/architecture_review.rst) |
| 运行时架构 | [运行时架构](source/runtime_architecture.rst) |
| 网格、网格拓扑与 Workbench 架构 | [网格架构](source/mesh_architecture.rst)、[网格拓扑](source/grid_topology.rst)、[Workbench 架构](source/workbench_architecture.rst) |
| AI 协作任务规范（中英） | [英文任务模板](AI_TASK_TEMPLATE.md)、[中文任务模板](AI_TASK_TEMPLATE_zh_CN.md) |
| AI 上下文与交接 | [AI 上下文](AI_CONTEXT.md) |
| 文档归属与整合计划 | [文档整合计划](plans/oneflow-documentation-consolidation-plan.md) |
| 变更预检与 Change Gate | [变更预检说明](source/change_preflight.rst) |

顶层索引回答系统包含什么、为何这样组织、主要流程如何工作以及风险在哪里；专题文档解释算法、契约和实现细节。
