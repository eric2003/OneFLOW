# OneFLOW 文档

`doc/` 是 OneFLOW 的统一文档根目录，包含历史手册、Sphinx 网站文档和项目报告。

## 目录结构

- `DevelopManual-chinese.doc`、`UserManual-chinese.doc`：原有 Word/WPS 中文手册。
- `source/`：Sphinx 文档源文件。
- `reports/`：维护中的项目交付报告和测试报告。
- `Makefile`、`make.bat`：Sphinx 构建入口。

## 构建 HTML 文档

在仓库根目录执行：

```bash
cd doc
make html
```

生成结果位于 `doc/build/html/`。CI 和 Read the Docs 均使用 `doc/source/`
作为 Sphinx 源目录。

新增或更新工程文档时，统一放在 `doc/` 下；不要重新创建顶层 `docs/` 目录。
