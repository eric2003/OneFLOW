# OneFLOW 文档

[English](README.md) | [简体中文](README_zh_CN.md) | [项目主页](../README.md) | [在线文档](https://oneflow-cfd.readthedocs.io/en/latest/)

`doc/` 是 OneFLOW 仓库的统一文档根目录，包含工程文档和正式 Sphinx 手册的源文件。两者用途不同：

- **工程文档**服务于开发、审查、项目交接、可复现操作和验证。
- **`source/`** 是经过整理、面向读者导航的 Sphinx 文档站点源目录。
- **并非每一份工程记录都需要发布到网站。** 正式手册应总结稳定的概念，并链接到详细记录，避免重复维护。

公开项目文档和面向贡献者的主要内容以英文为准。中文文档用于帮助读者更快理解部分主题；应保持概念、术语和导航一致，但不要求把每一份技术记录都逐句翻译。

## 目录结构

- `architecture/`：架构说明、设计边界和审查结论。
- `plans/`：持续维护的 TODO、路线图、改进计划和文档治理计划。
- `handoff/`：项目交接、阶段总结和后续工作记录。
- `runbooks/`：可重复执行的构建、测试、调试、集群和移植流程。
- `evidence/`：可复核的构建、MPI 和运行证据。
- `reports/`：按主题分类的正式报告。
- `_drafts/`：临时草稿，应定期审查、转为正式文档或清理；不作为正式发布入口。
- `source/`：Sphinx 文档站点源文件。
- `Makefile`、`make.bat`：Sphinx 构建入口。
- `DevelopManual-chinese.doc`、`UserManual-chinese.doc`：历史 Word/WPS 手册，尚待内容审查。

## 文档整合

仓库根目录和顶层 `docs/` 中仍有少量技术文档。在比较内容并检查链接之前，不要批量移动、删除或重命名这些文件。应先为每个主题确定唯一权威版本，再将稳定内容适当地纳入 Sphinx 手册。

## 从系统整体读到实现细节

文档按以下顺序组织：

1. **项目入口与目标。** 从仓库 [README](../README.md) 和本页开始，了解 OneFLOW 的定位和文档边界。
2. **系统架构。** 阅读[运行时架构与开发指南](source/runtime_architecture.rst)，理解运行时执行、算例生命周期、网格数据、工作台和验证之间的边界。
3. **子系统架构。** 根据任务进入网格、拓扑、工作台或数值模块的专题文档。
4. **算法与接口。** 解释数学定义、关键公式、数据流、实际代码入口和当前实现限制。
5. **验证证据。** 链接契约测试、集成测试、网格检查和 CFD 回归，并区分已验证行为与仍需验证的假设。

顶层文档应回答：系统由什么组成、为什么这样分层、主要控制流和数据流如何工作、目前有哪些主要风险。专题文档再解释算法如何计算、接口承诺是什么、实现位于哪里。不要在系统概览中重复函数级细节，也不要让技术笔记只列代码位置而不解释设计意图。

如果需要从整体审视架构，先读[系统架构审视](source/architecture_review.rst)，再读[运行时架构指南](source/runtime_architecture.rst)。之后可根据任务阅读[网格架构](source/mesh_architecture.rst)、[网格拓扑](source/grid_topology.rst)或[工作台架构](source/workbench_architecture.rst)。架构审视页提供的是检查框架，并不表示其中讨论的所有通用能力都已在 OneFLOW 中实现。

如果需要规范化表达较大规模的开发任务，请使用[AI 辅助工程任务规范](AI_TASK_TEMPLATE_zh_CN.md)；英文版本见 [AI task specification](AI_TASK_TEMPLATE.md)。该规范统一任务目标、全局收益评估、调查要求、实现约束、验证证据和交付条件。\n\n如果需要快速恢复持续开发的上下文，请阅读 [AI 项目状态与交接页](AI_CONTEXT.md)。该页面汇总分支约束、近期检查过的实现基线、已确认的行为、尚待验证的缺口和下一步建议。它是导航入口，不能替代源码检查以及针对目标提交的验证。

关于现有文档的初步盘点和整合顺序，请阅读[文档整合计划](plans/oneflow-documentation-consolidation-plan.md)。该计划记录文档归属和预期步骤，不代表其中列出的迁移都已完成。

新增工程文档统一放在 `doc/` 下，不再扩展顶层 `docs/` 作为第二套文档体系。新增文件应遵循[文档命名规范](DOCUMENT_NAMING.md)，适用于 `architecture/`、`plans/`、`handoff/`、`runbooks/`、`evidence/`、`reports/` 和 `_drafts/`。现有文件不批量改名，以免破坏历史链接。

## 构建 HTML 文档

在仓库根目录执行：

```bash
cd doc
make html
```

生成结果位于 `doc/build/html/`。CI 和 Read the Docs 使用 `doc/source/` 作为 Sphinx 源目录。修改正式发布的手册时，应检查相关 `toctree`、交叉引用和构建结果。只有计划向网站发布的工程记录才需要加入网站导航。
