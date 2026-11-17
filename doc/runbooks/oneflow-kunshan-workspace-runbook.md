# 昆山工作区与作业规范（runbook）

> 适用范围：昆山 SCNet 集群上所有 OneFLOW 相关工作（CPU 分区 `kshcnormal`、
> DCU 分区 `kshdnormal`）。目标：**一个固定工作区、固定路径、固定作业脚本**，
> 不再出现"每次一个目录、按 job 号建 build、日志散落在 home 根目录"。

## 1. 固定路径

工作区根固定为集群上的 `$HOME/oneflow/`（仓库里写成 `<workspace>/`，绝对路径只出现
在集群侧工作区 README，不进仓库）。远端操作走 `$HOME/.agents/skills/scnet-hpc`
（OpenAPI 优先、SSH 次之），不要手搓 ssh/scp 循环。除既有的共享工具链前缀外，不得在
别处新建 OneFLOW 工作区：

- 旧的按日期命名的实验目录（`softwares/projects/OneFLOW_*`，2026-08/09）已在
  2026-10-25 删除——源码在 Git 历史，证据在 `doc/`，不要再按那种方式建目录。
- 家目录根下不要留 OneFLOW 日志或一次性脚本：日志进 `logs/`，一次性脚本进
  `_attic/<date>/`，可复用的提升为 `scripts/` 里的标准脚本。

```text
<workspace>/
├── src/OneFLOW/                当前源码快照（每轮覆盖，不建 OneFLOW-<hash>）
│   └── REVISION                repo / branch / commit
├── deps/metis-install/         自编译 METIS 5.0.1
├── build/app/                  root 生产构建（HIP+MPI+METIS），增量
├── build/port/                 预留：standalone 1D 端口 contract 构建
├── build/platform/             预留：CPU-only root 构建
├── fixtures/m6-8zone-cpu8/     唯一合法可分区 fixture（只读）
├── evidence/<gate>/<jobid>/    门禁证据：gate | paired-ab | cpu-regression
├── snapshots/*.tgz             源码快照，保留最近 5 个
├── scripts/                    env.sh + 标准作业脚本 + cleanup.sh
├── logs/                       作业 stdout/stderr，保留 14 天
└── _attic/<YYYYMMDD>/          待删物，观察一轮后清空
```

规则：

1. **源码**：把仓库快照覆盖到 `src/OneFLOW/`（`rsync -a --delete` 或解压覆盖），
   同步更新 `src/REVISION`。不要在 `src/` 下保留多个修订目录。
2. **构建**：`build/app` 是唯一生产构建目录，按目标固定、增量编译。
   - **不要重新 configure**：计算节点 `/opt/rocm` 的旧 `hip-targets.cmake` 会让
     新配置失败；登录节点上重新 configure 同样会失败。需要恢复时，把上一个
     `build/app` 从 `snapshots/` 或 `_attic/` 取回，或到有完整 DTK 的节点重建。
   - 用 `scripts/build-app.sh`（即 `make -C build/app -j<N> OneFLOW/fast`）。
3. **fixture**：`fixtures/m6-8zone-cpu8/` 只读。派生算例用
   `scripts/prep-case.sh <dest> <steps>` 写到 `evidence/.../<jobid>/` 或
   `$TMPDIR`，作业结束随手清掉。
4. **证据**：`evidence/<gate>/<jobid>/` 只增不改；`*.stage` 这类大二进制在
   记录完指标后可以删除，但 `result.txt`、rank breakdown、作业输出必须保留。
5. **日志**：Slurm 的 `-o/-e` 指向 `logs/`（见脚本头），不要再落到 `$HOME`；
   `scripts/cleanup.sh` 保留 14 天。
6. **废弃**：旧 build/run/源码快照不要就地留着，移入 `_attic/<date>/`，下一轮
   清空。

## 2. 作业规范

| 约定 | 值 |
|---|---|
| 作业名 | `of-gate` / `of-paired` / `of-cpureg`（与证据目录同名） |
| 启动方式 | **只能 `mpirun`**；`srun` 下 MPI size=1，每个 rank 会解全部 zone |
| 布局 | `--bind-to core --map-by slot`；8 ranks/卡用 `scripts/rank.sh` 映射 `HIP_VISIBLE_DEVICES` |
| 环境 | 只 source `scripts/env.sh`（MPI 5.0.11 / UCX 1.22.0 / gcc16 / DTK 26.04 / cmake 3.30.5）；批处理里的 `module load` 是空操作 |
| 网络 | `UCX_NET_DEVICES=mlx5_0:1`、`UCX_IB_DEVICES=mlx5_0`（其他 HCA 会 `ibv_create_ah` EINVAL） |
| MPI coll | `OMPI_MCA_coll=^hcoll,han`（非均匀分布会崩在 `coll/han`） |
| 资源口径 | 1 DCU 槽 = 8 核；CPU 基线跑在 `kshcnormal` |
| 输出 | 证据目录内 `result.txt` + 日志；文档里只写脱敏后的口径与结论 |

标准脚本：

```bash
source $WORKSPACE/scripts/env.sh
$WORKSPACE/scripts/build-app.sh                       # 增量构建
sbatch $WORKSPACE/scripts/gate.slurm                  # trace + XB + 500 步分项 + 3x1000 步
sbatch $WORKSPACE/scripts/paired-ab.slurm             # 同二进制 paired A/B
sbatch $WORKSPACE/scripts/cpu-regression.slurm        # CPU 五算例 normal/strict
$WORKSPACE/scripts/cleanup.sh                         # 日志/快照轮转 + 占用
```

## 3. 每轮开发流程（与仓库文档一致）

1. 本地 `dev` 提交 → 推送 `origin/dev`。
2. 覆盖 `src/OneFLOW/` + 更新 `src/REVISION`，`scripts/build-app.sh` 增量编译。
3. 跑门禁（顺序固定）：3 步多区 trace + resident contract → 200 步跨后端 →
   CPU 五算例 normal/strict → 50 步稳定性；性能改动再加 paired。
4. 证据落到 `evidence/<gate>/<jobid>/`，数字写回
   `doc/plans/oneflow-development-todo.md` 的对应勾选项。
5. `scripts/cleanup.sh` 收尾；`_attic` 满一轮后清空。

## 4. 常见故障对照

| 现象 | 原因 | 处理 |
|---|---|---|
| `no ROCm-capable device is detected` | 节点未分配到 DCU / 卡被占用 | 检查 `--gres=dcu:N` 与节点，重投 |
| `hip-targets.cmake` 校验失败 | 在该节点重新 configure | 不要 reconfigure；用既有 `build/app` 增量编译 |
| `ldd` 找不到 `libmpi.so.40` | 未 source `scripts/env.sh` | 所有作业脚本第一行 source |
| `srun` 下 MPI size=1 | `srun` 不做 MPI 启动 | 改用 `mpirun` |
| `bnd_resize` 报 10 ms/step | 首个 `hipMalloc` 被平均进短基准 | 基准先 warmup 至少一步 |
| 多 zone 卡死 / 数值发散的假象 | 每 rank 的 cell 数太小 | 固定 per-rank cell 数再比较 |
