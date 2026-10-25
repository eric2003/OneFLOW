# OneFLOW 文档

`doc/` 是 OneFLOW 的统一文档根目录。工程资料与正式发布的 Sphinx 手册都由这里管理，但二者的组织目的不同：

- 工程资料按用途和生命周期分类，服务于开发、审计、交接与验证。
- `source/` 是经过整理、面向读者导航的 Sphinx 文档站点源目录。
- 并非每一份工程记录都需要发布到网站；正式手册可以提炼稳定结论并链接到详细工程记录，避免重复维护。

## 目录结构

- `architecture/`：架构说明、设计边界和审计结论。
- `plans/`：持续维护的 TODO、路线图、优化计划和文档治理计划。
- `handoff/`：项目交接、阶段进度和工作交接记录。
- `runbooks/`：构建、测试、调试、集群及移植等可重复操作流程。
- `evidence/`：构建、MPI、节点运行等可复核事实记录。
- `reports/`：正式报告，按主题子目录组织。
- `_drafts/`：临时草稿；应定期整理、提升为正式文档或清理，不作为发布内容。
- `source/`：Sphinx 文档站点源文件。
- `Makefile`、`make.bat`：Sphinx 构建入口。
- `DevelopManual-chinese.doc`、`UserManual-chinese.doc`：历史 Word/WPS 手册，状态待内容审查。

## 文档资产整理计划

当前仓库仍有少量技术文档位于根目录和顶层 `docs/`。在完成内容比较和链接检查之前，不要批量移动、删除或重命名这些文件。按主题确定唯一权威版本，再将稳定内容纳入 Sphinx 手册。

开发者需要从整体理解代码时，先阅读 [运行时架构与开发指南](source/runtime_architecture.rst)，再进入 [网格架构](source/mesh_architecture.rst)、[网格拓扑](source/grid_topology.rst) 或 [工作台架构](source/workbench_architecture.rst)。该指南记录当前代码入口、配置边界、运行流程和重构约定，并会随真实代码变更持续修订。\n\n请先阅读 [文档资产盘点与整合计划](plans/oneflow-documentation-consolidation-plan.md)。该计划记录现有文档的初步归属、建议的目标结构和分阶段整理顺序。它是持续更新的计划，不表示其中提到的迁移已经完成。

新增工程文档统一放在 `doc/` 下，不再扩展顶层 `docs/` 作为第二套体系。新建工程文档遵循[文档命名规范](DOCUMENT_NAMING.md)，该规范适用于 `architecture/`、`plans/`、`handoff/`、`runbooks/`、`evidence/`、`reports/` 与 `_drafts/`。现有文件不批量改名，以免破坏历史链接。

## 构建 HTML 文档

在仓库根目录执行：

```bash
cd doc
make html
```

生成结果位于 `doc/build/html/`。CI 和 Read the Docs 使用 `doc/source/` 作为 Sphinx 源目录。修改发布手册时，应检查相关 `toctree`、交叉引用和构建结果；仅更新工程记录时，不必强制把该记录加入网站导航。
