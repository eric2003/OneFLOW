# OneFLOW 文档索引

[English](README.md) | [简体中文](README_zh_CN.md) | [Project README](../README.md) | [Online documentation](https://oneflow-cfd.readthedocs.io/en/latest/)

`doc/` 是 OneFLOW 的工程文档根目录。页面按**生命周期**分区：长期维护的活文档、
一次性快照、可复核证据、正式报告互不混放。

> 文档命名遵循 [DOCUMENT_NAMING.md](DOCUMENT_NAMING.md)。新增文档先选目录与
> artifact，再写内容；删除或合并旧文档必须同步更新本索引和所有引用。

## 1. 目录与保留策略

| 目录 | 放什么 | 保留策略 |
|---|---|---|
| [`architecture/`](architecture/) | 稳定的架构说明、边界审计 | 每个主题只保留**一份**现行版本；被取代的审计直接删除（git 历史即档案） |
| [`plans/`](plans/) | 活文档：路线图、长周期计划 | 只放**living**文档（无日期后缀）；一次性 checklist 完成后删除；从 upstream `reports/architecture/` 映射过来的快照例外（见 §3） |
| [`handoff/`](handoff/) | 交接（living） | **只有一份** `oneflow-current-handoff.md`，就地更新；不新建 dated handoff |
| [`runbooks/`](runbooks/) | 构建、集群、移植等操作流程 | 每个流程一份；步骤变更就地更新，不建副本 |
| [`evidence/`](evidence/) | 可复核的环境/构建事实 | 与某一轮结论绑定的证据；结论被取代后随之删除 |
| [`reports/`](reports/) | 正式报告（performance / regression / delivery） | **不删已发布的报告**；过时报告保留但标注 historical，Markdown 与 HTML 同 stem 同步 |
| [`source/`](source/) | Sphinx 教程源（Read the Docs） | 与工程文档无关，按教程需要维护 |
| `_drafts/` | 未定稿草稿 | 定稿后删除，不长期保留 |

## 2. 现行入口（从这里开始读）

| 想了解 | 文档 |
|---|---|
| 当前状态、Phase A–F 路线图与勾选进度、每轮门禁 | [plans/oneflow-development-todo.md](plans/oneflow-development-todo.md) |
| 长期阶段的取舍、里程碑与开放问题 | [plans/oneflow-development-plan.md](plans/oneflow-development-plan.md) |
| 交接（living：环境、固定路径、启动顺序、坑） | [handoff/oneflow-current-handoff.md](handoff/oneflow-current-handoff.md) |
| 昆山工作区固定路径、作业规范、清理规则 | [runbooks/oneflow-kunshan-workspace-runbook.md](runbooks/oneflow-kunshan-workspace-runbook.md) |
| DCU 移植方法论（契约→oracle→门禁→benchmark） | [runbooks/oneflow-dcu-porting-playbook-20260903.md](runbooks/oneflow-dcu-porting-playbook-20260903.md) |
| 主 solver 加速器架构审计与接入边界 | [architecture/oneflow-1d-euler-dcu-architecture-audit-20260904.md](architecture/oneflow-1d-euler-dcu-architecture-audit-20260904.md) |
| 报告列表 | [reports/README.md](reports/README.md) |

## 3. 仓库内与上游共有的文件

`doc/` 里既有 fork 维护的文档，也有 upstream 拥有的文档。动手前的判断规则：

- **可自由整理（fork 维护）**：`plans/oneflow-development-*.md`、`handoff/`、
  本索引、`DOCUMENT_NAMING.md` 中 fork 新增的章节。
- **与 upstream 共有**：`architecture/`、`evidence/`、`reports/`、`runbooks/` 下的
  报告与流程、`source/`。可以更新内容，但不要整体移动或重写，否则上游同步时会产生
  结构性冲突；确需移动时单独提交并同步修复 `ports/`、`ci/` 等处的引用。
- `doc/reports/architecture/` 只存在于 upstream（fork 已按目录规范拆到
  `plans/`、`runbooks/`、`reports/delivery/`）；不要在本仓库重新创建该目录。

合并 upstream 时该目录下新增的文件按下列固定映射落到 fork 目录，并在同一合并提交里
完成移动，避免目录被“无意中”重建：

| upstream 路径 | fork 路径 |
|---|---|
| `reports/architecture/oneflow-dcu-porting-playbook-*.md` | `runbooks/`（同名） |
| `reports/architecture/oneflow-gpu-backend-delivery.*` | `reports/delivery/`（同名） |
| `reports/architecture/oneflow-euler-optimization-plan.md` | `plans/`（同名） |
| `reports/architecture/<other>.md`（设计说明/快照） | `plans/`（同名） |


For a broad architecture review, read [system architecture review](source/architecture_review.rst) first, then the [runtime architecture guide](source/runtime_architecture.rst). Continue with [mesh architecture](source/mesh_architecture.rst), [grid topology](source/grid_topology.rst), or [Workbench architecture](source/workbench_architecture.rst), depending on the task. The review page is a checklist; it does not claim that every general capability discussed there has already been implemented.

For the normative [continuous evolution and scientific verification development specification](AI_TASK_TEMPLATE.md), use the English version or its [Chinese counterpart](AI_TASK_TEMPLATE_zh_CN.md). It defines single-workstream prioritization, implementation constraints, scientific verification and validation, references, reproducible plots, and delivery evidence.

For a quick continuation of ongoing development, read [AI project context and handoff](AI_CONTEXT.md). It summarizes branch constraints, recently inspected implementation baselines, confirmed behavior, open verification gaps, and suggested next steps. It is a navigation aid, not a substitute for source inspection and verification against the target commit.

For the initial inventory and planned consolidation of existing documents, see the [documentation consolidation plan](plans/oneflow-documentation-consolidation-plan.md). The plan records intended ownership and sequencing; it does not imply that every listed migration is complete.

New engineering documents belong under `doc/`; do not grow the top-level `docs/` directory into a second documentation system. Follow the [document naming rules](DOCUMENT_NAMING.md) for new files in `architecture/`, `plans/`, `handoff/`, `runbooks/`, `evidence/`, `reports/`, and `_drafts/`. Existing files should not be renamed in bulk because historical links may depend on their paths.

## Build the HTML documentation

From the repository root:

## 4. 构建 HTML 文档


```bash
cd doc
make html     # 输出到 doc/build/html/，CI 与 Read the Docs 使用 doc/source/
```

正式报告若同时提交 HTML，必须与 Markdown 使用相同 stem（例如
`oneflow-euler-performance-current.md` / `.html`），并在同一提交中更新。

## 5. 每轮收尾检查表

1. `plans/oneflow-development-todo.md` 的勾选、状态快照与本轮证据已更新。
2. `handoff/oneflow-current-handoff.md` 的 §4（任务指针）与 §5（新坑）已就地更新；
   没有新增 dated handoff，一次性 checklist/草稿已删除。
3. 新增/删除/改名文档后，本索引、`doc/reports/README.md`、`AGENTS.md` 的链接
   全部有效（`grep` 一遍旧文件名应无残留）。
4. `git diff --check` 通过；文档中不含原始 CI/Slurm 日志、主机名、账号、作业号、
   绝对私有路径。

## 6. 上游文档入口与阅读顺序

`doc/` 同时存放两类内容，用途不同：

- **工程文档**服务于开发、评审、交接、可复现与验证，即上文 §1 的生命周期分区。
- **`source/`** 是面向读者的 Sphinx 手册源；并非每条工程记录都应出现在公开文档站点上。
  手册应概括稳定概念并链接到详细记录，而不是复制它们。
- Sphinx 手册以英文为主；中文文档是部分主题的阅读辅助，概念、术语与导航保持一致。

从系统到实现细节的阅读顺序：

1. **项目入口与目标**：仓库 [README](../README.md) 与本页。
2. **系统架构**：[运行时架构与开发指南](source/runtime_architecture.rst)，
   覆盖运行期执行、case 生命周期、网格数据、Workbench 与验证边界。
3. **子系统架构**：按任务选择网格、拓扑、Workbench 或数值组件。
4. **算法与接口**：数学定义、重要方程、数据流、真实代码入口与当前实现边界。
5. **验证证据**：契约测试、集成测试、网格检查与 CFD 回归；
   区分已验证行为与仍待验证的假设。

相关入口：

| 想了解 | 文档 |
|---|---|
| 系统架构评审清单 | [source/architecture_review.rst](source/architecture_review.rst) |
| 运行时架构 | [source/runtime_architecture.rst](source/runtime_architecture.rst) |
| 网格架构 / 网格拓扑 / Workbench | [source/mesh_architecture.rst](source/mesh_architecture.rst)、[source/grid_topology.rst](source/grid_topology.rst)、[source/workbench_architecture.rst](source/workbench_architecture.rst) |
| AI 协作任务规格（中英） | [AI_TASK_TEMPLATE.md](AI_TASK_TEMPLATE.md)、[AI_TASK_TEMPLATE_zh_CN.md](AI_TASK_TEMPLATE_zh_CN.md) |
| AI 上下文与交接 | [AI_CONTEXT.md](AI_CONTEXT.md) |
| 文档归属与合并计划 | [plans/oneflow-documentation-consolidation-plan.md](plans/oneflow-documentation-consolidation-plan.md) |

顶层文档回答系统包含什么、为何这样分层、主要控制流与数据流如何工作、风险在哪里；
主题文档回答算法如何工作、接口承诺什么、实现在哪里。不要把函数级细节复制到系统总览中。

新的工程文档放在 doc/ 下，不要把顶层 docs/ 扩张成第二个文档系统；已有文件不要批量改名，历史链接可能依赖其路径。
