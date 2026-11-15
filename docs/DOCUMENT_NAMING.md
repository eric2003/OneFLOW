# OneFLOW 文档命名规范

本文是 `doc/` 下工程文档的唯一命名规则。它以
[`doc/reports/README.md`](reports/README.md) 的既有 reports 约定为基础，扩展到
architecture、plans、handoff、runbooks、evidence 和 drafts。

## 1. 目标与适用范围

命名应让读者不打开文件就能判断三件事：主题、文档用途、时间属性。规则适用于
`doc/` 下新建或新命名的工程文档；不要求批量重命名已有文件，以免破坏历史链接、
外部引用和 Sphinx 索引。

Sphinx 教程源文件 `doc/source/` 按内容层级组织，不受本规范的文件名语法约束。
根目录的 `README.md`、`DOCUMENT_NAMING.md`、Makefile、构建脚本和历史
Word/WPS 手册也属于保留文件名。

## 2. 标准语法

```text
oneflow-<scope>-<artifact>[-<qualifier>][-YYYYMMDD][-rN].<ext>
```

- 全部使用小写 ASCII 和连字符（kebab-case）；不使用空格、下划线、中文文件名或个人名。
- `scope` 是稳定主题，而不是某次任务描述，例如 `1d-euler-dcu`、
  `main-solver-dcu`、`euler`、`gpu-backend`、`development`。
- `artifact` 描述文档职能，使用受控词表。
- `qualifier` 只在需要区分同类文档时使用，例如 `current`、`baseline`、
  `addendum`、`checklist`、`playbook`。
- `YYYYMMDD` 是发布、快照或证据形成日期，不是最后一次编辑日期。
- 同一主题同日有多个不可合并的快照时，用 `-r2`、`-r3`；不使用
  `final`、`latest`、`new`、`v2` 等模糊后缀。修订历史由 Git 管理。

## 3. 目录与 artifact 对照

| 目录 | 首选 artifact | 示例 |
|---|---|---|
| `architecture/` | `architecture`、`design`、`audit` | `oneflow-1d-euler-dcu-architecture-audit-20260904.md` |
| `plans/` | `plan`、`todo`、`checklist` | `oneflow-development-todo.md` |
| `handoff/` | `handoff`、`recap`、`addendum` | `oneflow-main-solver-dcu-recap-20260916.md` |
| `runbooks/` | `runbook`、`playbook` | `oneflow-dcu-porting-playbook-20260903.md` |
| `evidence/` | `evidence`、`validation` | `oneflow-1d-euler-dcu-build-evidence-20260904.md` |
| `reports/performance/` | `performance`、`benchmark` | `oneflow-euler-performance-current.md` |
| `reports/regression/` | `regression`、`validation` | `oneflow-1d-euler-dcu-hip-contract-20260915.md` |
| `reports/delivery/` | `delivery` | `oneflow-gpu-backend-delivery.md` |
| `_drafts/` | `draft` | `oneflow-main-solver-dcu-draft-20260916.md` |

目录表达存放边界，artifact 表达读者预期；不要用目录名重复替代 artifact，例如
`reports/oneflow-euler-report-report.md`。

## 4. 时间、活文档与发布格式

### 活文档

长期维护、始终指向当前状态的文档不带日期，必要时使用 `current`：

- `oneflow-development-todo.md`
- `oneflow-development-plan.md`
- `oneflow-current-handoff.md`（交接也改为 living：只有一份，就地更新）
- `oneflow-euler-performance-current.md`

它们必须在开头写明“最后更新”日期，并在同一文件中维护当前结论；不要为每次小改动
新建 dated copy。

### 快照、交接与证据

一次性事实、审计和可复核证据必须带日期（交接不在此列，见上）：

- `oneflow-1d-euler-dcu-architecture-audit-20260904.md`
- `oneflow-1d-euler-dcu-hip-contract-20260915.md`
- `oneflow-1d-euler-dcu-mpi-smoke-evidence-20260904.md`

### Markdown 与 HTML

Markdown 是工程文档和报告的源文件。若需要提交 HTML，HTML 必须与 Markdown 使用
相同 stem：

```text
oneflow-euler-performance-current.md
oneflow-euler-performance-current.html
```

不为只有 Markdown 的 handoff、plan、evidence 强行生成 HTML。

## 4.5 保留、合并与删除

文档目录必须长期保持精简，避免"每轮新增一份"。

| 类型 | 规则 |
|---|---|
| 活文档（todo、长周期计划、本索引） | 只保留一份，就地更新；不建 dated 副本 |
| 会话交接（`handoff/`） | **只有一份** living 交接 `oneflow-current-handoff.md`，就地更新 §4/§5；不再新建 dated handoff |
| 一次性计划/清单（`plans/` 里的 checklist） | 目标达成即删除；其证据并入 todo 或报告 |
| 架构审计（`architecture/`） | 每个主题只保留现行版本；被取代时删除旧审计 |
| 证据（`evidence/`） | 只保留支撑**现行结论**的证据；结论被取代时删除 |
| 正式报告（`reports/`） | **不删除已发布的报告**；过时报告标注 historical 并保留 Markdown/HTML |
| 草稿（`_drafts/`） | 定稿后删除 |

操作要求：

1. 删除或改名必须用 `git rm` / `git mv`，并在同一提交内修复所有引用
   （`doc/README.md`、`doc/reports/README.md`、`AGENTS.md`、`ports/`、`ci/`、代码注释）。
2. 提交前用 `grep -rIn --exclude-dir=.git <旧文件名>` 确认零残留。
3. Git 历史就是档案；不需要用 `_archive/` 目录堆积副本。
4. 每轮收尾按 `doc/README.md` §5 的检查表执行。

## 5. 内容边界与禁止项

文件名和已提交文档均不得包含 branch、作者、账号、主机名、节点名、作业号、私有路径、
token 或凭据。原始 CI/Slurm 日志、环境 dump、临时运行目录留在 CI artifact 或集群
工作区；文档只保留脱敏后的资源口径、工具链版本、通过判据、关键数值和结论。

不要用文件名承载状态，例如 `pass`、`failed`、`wip`、`final`。状态应写在
正文、TODO 或 Git commit 中。

## 6. 新建与维护流程

1. 先选择目录，再选择 `scope` 和 `artifact`。
2. 判断是活文档还是快照：前者无日期，后者带 `YYYYMMDD`。
3. 新增或修改正式报告时，检查是否存在同 stem 的 HTML，并按需要同步更新。
4. 在目录 README 或 `doc/README.md` 中补入口；不要靠文件系统遍历作为导航。
5. 提交前运行 `git diff --check`；确认未提交 raw log、私有元数据或无意义的重名快照。

## 7. 迁移策略

本规范从生效后约束新增文件。已有文件只有在内容发生实质重写、链接可同时更新且收益明确时
才改名；否则保留历史名称并在 README 索引中归类。任何改名必须使用 `git mv`，并同步
修复仓库内链接。
