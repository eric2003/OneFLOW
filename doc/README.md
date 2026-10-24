# OneFLOW 文档



`doc/` 是 OneFLOW 的统一文档根目录。按文档生命周期和用途分区，避免把计划、证据、交接和正式报告混在一起。



## 目录结构



- `architecture/`：稳定的架构说明、设计边界和审计结论。

- `plans/`：living TODO、路线图、优化计划和 porting checklist。

- `handoff/`：项目交接、会话交接和阶段进度补充。

- `runbooks/`：构建、集群、DCU porting 等操作流程。

- `evidence/`：构建、MPI、节点验证等可复核事实记录。

- `reports/`：正式交付报告，按 `performance/`、`regression/`、`delivery/` 分类。

- `_drafts/`：未完成草稿和本地交接笔记，不纳入正式发布。

- `source/`：Sphinx 文档源文件。

- `Makefile`、`make.bat`：Sphinx 构建入口。

- `DevelopManual-chinese.doc`、`UserManual-chinese.doc`：原有 Word/WPS 手册。



## 构建 HTML 文档



在仓库根目录执行：



```bash

cd doc

make html

```



生成结果位于 `doc/build/html/`。CI 和 Read the Docs 均使用 `doc/source/`

作为 Sphinx 源目录。



新增或更新工程文档时，统一放在 `doc/` 下；不要重新创建顶层 `docs/` 目录。

新建工程文档还必须遵循 [文档命名规范](DOCUMENT_NAMING.md)；该规范统一约束
`architecture/`、`plans/`、`handoff/`、`runbooks/`、`evidence/`、
`reports/` 与 `_drafts/`。现有文件不批量改名，以免破坏链接。
