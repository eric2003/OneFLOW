# 昆山工作区与作业规范（runbook）

> 适用范围：昆山 SCNet 集群上所有 OneFLOW 相关工作（CPU 分区 `kshcnormal`、
> DCU 分区 `kshdnormal`）。目标：**一个固定工作区、固定路径、固定作业脚本**，
> 不再出现"每次一个目录、按 job 号建 build、日志散落在 home 根目录"。
> 最后核对：2026-11-11；本轮修正文档冲突，未重新探测远端环境。

## 1. 固定路径

工作区根固定为集群上的 `$HOME/oneflow/`（仓库里写成 `<workspace>/`，绝对路径只出现
在集群侧工作区 README，不进仓库）。远端操作走当前安装的 scnet-hpc 技能；
按 2026-11-08 的接入记录，SSH 不可用，**只走 OpenAPI**，远端动作均提交为作业。
不要手搓 ssh/scp 循环。除既有的共享工具链前缀外，不得在
别处新建 OneFLOW 工作区：

- 旧的按日期命名的实验目录（`softwares/projects/OneFLOW_*`，2026-08/09）已在
  2026-10-25 删除——源码在 Git 历史，证据在 `doc/`，不要再按那种方式建目录。
- 家目录根下不要留 OneFLOW 日志或一次性脚本：日志进 `logs/`，一次性脚本进
  `_attic/<date>/`，可复用的提升为 `scripts/` 里的标准脚本。

```text
<workspace>/
├── src/                       当前源码与版本标记
│   ├── OneFLOW/                每轮覆盖，不建 OneFLOW-<hash>
│   └── REVISION                repo / branch / commit
├── deps/metis-install/         自编译 METIS 5.0.1
├── build/app/                  root 生产构建（HIP+MPI+METIS），增量
├── build/port/                 预留：standalone 1D 端口 contract 构建
├── build/platform/             预留：CPU-only root 构建
├── fixtures/m6-8zone-cpu8/     固定登记的 fixture（只读；派生配方见 §1b）
├── evidence/<category>/<jobid>/  证据：6 个固定 category，见 §1 第 4 条
├── snapshots/*.tgz             源码快照，保留最近 5 个
├── scripts/                    env.sh + 标准作业脚本 + cleanup.sh
├── logs/                       作业 stdout/stderr，保留 14 天
└── _attic/<YYYYMMDD>/          待删物，观察一轮后清空
```

规则：

1. **源码**：`git archive` → OpenAPI upload 到 `snapshots/` → 作业内覆盖解包到
   `src/OneFLOW/`，同步更新 `src/REVISION`；清除旧快照已删除的源码，不能只叠加解包。
   OpenAPI 不提供交互 exec，不使用 rsync；不在 `src/` 下保留多个修订目录。
2. **构建**：`build/app` 是唯一生产（HIP）构建目录，按目标固定、增量编译。
   - 日常增量：`scripts/build-app.sh`（即 `make -C build/app -j<N> OneFLOW_RuntimeObjects`
     再 `make -C build/app -j<N> OneFLOW/fast`）。**两个目标都要跑**：`<dir>/fast` 是
     CMake 的"仅本目录"目标，不会重建 `OneFLOW_RuntimeObjects` 里的对象；2026-11-07
     实测只跑 `OneFLOW/fast` 时，改过 `codes/accel/src/HipFluxBackend.hip` 也报构建成功，
     但对象和二进制都还是旧的（`strings bin/OneFLOW` 里找不到新符号才发现）。
   - 上游改动 `CMakeLists.txt` 后 `make` 会自动 reconfigure。**不要再手拼 cmake 变量**，
     统一走仓库里的加速器 SDK 接口（`cmake/OneFLOWAccelSDK.cmake` +
     `scripts/build-accel.sh`）；它把下面三条一次配齐：
     1. 作业要在 **DCU 节点**（`kshdnormal`，`dcu:1`；`/dev/kfd` 与
        `HIP_VISIBLE_DEVICES` 可见）——CPU 队列上没有设备，检测不到架构；
     2. `ROCM_PATH`：DTK 的 `hip-targets.cmake` 用
        `set(_IMPORT_PREFIX "${ROCM_PATH}/hip" "/opt/dtk/hip")` 定位
        `libgalaxyhip`，不设就落到 `/opt/rocm` → 报 “install faulty /
        not all the files it references”；
     3. `amd_comgr` / `AMDDeviceLibs` / `hsa-runtime64` 要传给 cmake，**且**用
        `CMAKE_TRY_COMPILE_PLATFORM_VARIABLES` 让 `enable_language(HIP)` 的
        try_compile 子工程继承，否则报 “Could not find amd_comgr”。

     ```bash
     source $ONEFLOW_ROOT/scripts/env.sh          # 提供 ROCM_PATH 与 ONEFLOW_DTK_HIP_FLAGS

     # 仅在确需全新配置或修复配置时，在 DCU 分配内执行；日常用 build-app.sh
     bash $ONEFLOW_SRC/scripts/build-accel.sh --sdk auto --sdk-root "$ROCM_PATH" \
       --arch gfx906 --build-dir $ONEFLOW_ROOT/build/app \
       --flags "$ONEFLOW_DTK_HIP_FLAGS" --jobs 8

     # 增量更新生产构建
     bash $ONEFLOW_SRC/scripts/build-accel.sh --sdk hip --sdk-root "$ROCM_PATH" \
       --arch gfx906 --build-dir $ONEFLOW_ROOT/build/app \
       --flags "$ONEFLOW_DTK_HIP_FLAGS" --jobs 8
     ```

     SDK 种类不需要指定：由 `ONEFLOW_ENABLE_HIP/CUDA` 或前缀内容（`llvm/bin/clang++` /
     `bin/nvcc`）推断；`--sdk-root` 省略时取 `$ROCM_PATH`/`$CUDAToolkit_ROOT`/`$CUDA_HOME`。
     `--configure-only` 只配置不编译。CUDA 侧同一接口：`--sdk-root <CUDA Toolkit> --arch sm_80`。

   - CPU-only 备用构建（不涉及 HIP 时可用，可在 CPU 节点配置）：`build/platform`。
     CPU 门禁可以指向它：

     ```bash
     export ONEFLOW_CPU_BIN=$ONEFLOW_ROOT/build/platform/bin/OneFLOW
     bash $ONEFLOW_ROOT/scripts/cpu-regression.slurm
     ```

     `cpu-regression.slurm` 默认仍用 `$ONEFLOW_BIN`（`build/app`），只在显式设置
     `ONEFLOW_CPU_BIN` 时切换，并把 `system/` 拷到被测二进制同目录。
     实测两种构建对 CPU 五算例给出**完全相同的残差**，可互为对照。
3. **fixture**：`fixtures/m6-8zone-cpu8/` 只读。派生算例用
   `scripts/prep-case.sh <dest> <steps>` 写到 `evidence/.../<jobid>/` 或
   `$TMPDIR`，作业结束随手清掉。
4. **证据**：`evidence/<category>/<jobid>/` 只增不改；`*.stage` 这类大二进制在
   记录完指标后可以删除，但 `result.txt`、rank breakdown、作业输出必须保留。
   `<category>` 只能用下面 6 个固定值，**不要按实验内容临时起名**（命名规则见 §1c）：

   | category | 放什么 |
   |---|---|
   | `gate` | 标准门禁：3 步 trace / 200 步 XB / 500 步分项 / 3×1000 步 |
   | `paired-ab` | 同节点同作业内的 A/B 交替测量 |
   | `cpu-regression` | CPU 五算例 normal/strict |
   | `probe` | 探索性测量：分项剖析、启动/扩展性扫描、单点 A/B、诊断 dump |
   | `fixture` | 分区、网格生成、manifest 导出与校验 |
   | `build` | 构建、工具链、链接问题 |

   实验身份由**作业名 + 日志名**承载（`logs/of-<category>-<purpose>-<jobid>.out`），
   不靠目录名。旧轮次留下的临时类别目录按 `evidence/probe/<jobid>/` 归位。
5. **日志**：Slurm 的 `-o/-e` 指向 `logs/`（见脚本头），不要再落到 `$HOME`；
   `scripts/cleanup.sh` 保留 14 天。
6. **废弃**：旧 build/run/源码快照不要就地留着，移入 `_attic/<date>/`，下一轮
   清空。

## 1b. 多 zone fixture 的制作（16/32 zone，用于扩展性测量）

`fixtures/m6-8zone-cpu8` 是唯一登记的 fixture；更大 zone 数必须自己划分，并走完
manifest 校验，否则数据不可用（2026-10-23 的 ad-hoc 16-zone 就是因为缺 manifest 被排除）。
下面这条配方 2026-11-08 逐条验证过，缺一步就会失败：

```bash
# 1) 划分：参数必须放 script/grid.txt，并且 control.txt 里要列出 "grid.txt";
#    否则 gridObj 读不到，objective 回落到默认 ConvertOnly，报
#    "Unsupported grid conversion: oneflow -> oneflow"。
cat > <case>/script/grid.txt <<'TXT'
int gridObj = 3;
string sourceGridFileName = "grid/m6wing.ofl";
string sourceGridType = "oneflow";
string targetGridType = "oneflow";
string topoType = "u";
int npartproc = 16;
int partition_type = 1;
string part_uns_file = "grid/m6wing-16.ofl";
TXT
sed -i '1i "grid.txt";' <case>/script/control.txt
printf 'string simutask = "Grid";\n' > <case>/script/simu.txt
mpirun -np 1 "$ONEFLOW_BIN" 0 <case>        # 约 20 s

# 2) 求解算例指向新网格，并且必须补 grid/bcRegionMap.txt
#    （缺它会在 UNsBcSolver::CalcBcRegion 段错误，导出阶段直接 SIGSEGV）
sed -i -E 's#^string gridFileName.*#string gridFileName = "grid/m6wing-16.ofl";#' <case>/script/cfd.txt
printf '4\n1 BCFarfield1\n0 BCOutflow0\n2 BCSymmetryPlane2\n3 BCWall0\n' > <case>/grid/bcRegionMap.txt
printf 'string simutask = "Solve";\n' > <case>/script/simu.txt

# 3) 导出 rank 分片 → 合并 → 验证（verify 通过才算合法 fixture）
mpirun -np 16 env ONEFLOW_ACCEL_BACKEND=CPU \
    ONEFLOW_MZ_FIXTURE_EXPORT_DIR=<shards> "$ONEFLOW_BIN" 0 <case>
python3 "$ONEFLOW_SRC/ci/kunshan/merge-multizone-fixture.py" \
    --shards <shards> --output fixture-metadata.json
python3 "$ONEFLOW_SRC/ci/kunshan/verify-multizone-fixture.py" \
    --metadata fixture-metadata.json --expected-ranks 16
```

上述参数必须以实际加载顺序为准：若 `control.txt` 在 `grid.txt` 后加载 `part.txt`，
`npartproc` / `partition_type` / `part_uns_file` 应在 `part.txt` 中同步设置，
避免后加载的旧值覆盖本次分区参数（见交接坑 19）。

2026-11-08 结果：16 zone / 16 rank，`zones=16 interfaces=52946`，
`MULTIZONE_FIXTURE_PASS`。32 zone 用同一配方（`npartproc=32`）。

**选哪种 `partition_type`（2026-11-09 起）**

| 值 | 分区器 | 用途 |
|---:|---|---|
| 1 | METIS 图划分（默认） | 生产算例；最小化割边，但**不保证每 part 形状/接口数相同** |
| 2 | 结构化等块（cell-count 分位切 a×b×c） | **扩展性测量专用**：cell 严格等分，块形状由数组维度决定，接口数可解释 |

`partition_type = 2` 在固定 32 rank 下实测比 METIS 快约 **4.4%**（每步 61.32 → 58.62 ms，
平均接口面 1786 → 1500），因为步长由 halo 最大的 rank 决定。分区日志会打印
`Structured partition: <a> x <b> x <c> blocks ... (effective lengths ...)`；
块数量按 **cell 数**的有效长度选择（不是物理盒尺寸），所以同一物理立方体在 z 向加密
4 倍时会自动切成 z 向 4 倍数量的块。扩展性表必须同时记录每 rank 的 `nIFaces`，
否则无法区分分区差异和代码扩展性。

**必须同时把算例改成显式 RK**：从 `test/cfd/m6wingroe_sa` 派生的算例带
`int time_integral = 2`（LU-SGS），设备驻留链路只在 `time_integral = 1` 下启用。
2026-11-08 就因为漏了这一行，16/32 zone 的 HIP 跑成了 host LU-SGS，报出
2.6 s/step 的假数字，被当成“多卡扩展性差”。派生算例至少显式设置：

```bash
sed -i -E 's/^int[[:space:]]+time_integral[[:space:]]*=.*/int time_integral = 1;/;
            s/^int[[:space:]]+rk_stage[[:space:]]*=.*/int rk_stage = 3;/;
            s/^int[[:space:]]+ischeme[[:space:]]*=.*/int ischeme = 5;/;
            s/^int[[:space:]]+vismodel[[:space:]]*=.*/int vismodel = 0;/' <case>/script/cfd.txt
sed -i -E '/^[[:space:]]*TurbSolver[[:space:]]*$/d' <case>/script/solver.txt
```

## 1c. 命名规范（作业名、证据目录、脚本）

沿用 `doc/DOCUMENT_NAMING.md` 的 kebab-case 原则，只是把范围延伸到集群侧：

| 对象 | 语法 | 例 |
|---|---|---|
| 作业名 | `of-<category>[-<purpose>]` | `of-gate`、`of-paired-ab`、`of-cpu-regression`、`of-probe-scaling`、`of-fixture-partition`、`of-build` |
| 证据目录 | `evidence/<category>/<jobid>/` | `evidence/probe/127303578/` |
| 作业日志 | `logs/of-<category>[-<purpose>]-<jobid>.out` | `logs/of-probe-scaling-127303578.out` |
| 一次性脚本 | `_attic/<YYYYMMDD>/`，可复用才提升为 `scripts/` | `scripts/build-app.sh` |

`<category>` 是上表的 6 个固定值；`<purpose>` 只在需要区分同一类别下的多次实验时使用，
用一两个单词（`scaling`、`startup`、`limiter`），**不要**用 `a17`、`ab2`、`attr`、`bnd`
这类缩写或序号——那些信息应该由作业名和结论文档承载。

## 2. 作业规范

| 约定 | 值 |
|---|---|
| 作业名 | `of-<category>[-<purpose>]`，固定 category 见 §1c |
| 启动方式 | **只能 `mpirun`**；`srun` 下 MPI size=1，每个 rank 会解全部 zone |
| 布局 | `--bind-to core --map-by slot`；8 ranks/卡用 `scripts/rank.sh` 映射 `HIP_VISIBLE_DEVICES` |
| 环境 | 常规作业 source `scripts/env.sh`（MPI 5.0.11 / UCX 1.22.0 / gcc16 / DTK 26.04 / cmake 4.4.3；CMake 必须在 PATH）；rocprof 另需交接坑 31 的 module 初始化 |
| 网络 | `UCX_NET_DEVICES=mlx5_0:1`、`UCX_IB_DEVICES=mlx5_0`（其他 HCA 会 `ibv_create_ah` EINVAL） |
| MPI coll | `OMPI_MCA_coll=^hcoll,han`（非均匀分布会崩在 `coll/han`） |
| 资源口径 | 1 DCU 槽 = 8 核；CPU 基线跑在 `kshcnormal` |
| 输出 | 证据目录内 `result.txt` + 日志；文档里只写脱敏后的口径与结论 |
| 提交器 | OpenAPI 提交用 `--work-dir <workspace>/logs`，`--stdout/--stderr` 也指到 `logs/` | 提交器会在 work-dir 生成 `job_*` 记账文件，指到工作区根会污染根目录 |

标准脚本：

```bash
source $WORKSPACE/scripts/env.sh
$WORKSPACE/scripts/build-app.sh                       # 增量构建
bash $WORKSPACE/scripts/gate.slurm                    # trace + XB + 500 步分项 + 3x1000 步
bash $WORKSPACE/scripts/paired-ab.slurm               # 同二进制 paired A/B
bash $WORKSPACE/scripts/cpu-regression.slurm          # CPU 五算例 normal/strict
$WORKSPACE/scripts/cleanup.sh                         # 日志/快照轮转 + 占用
```

以上是作业内调用，不是本地提交命令；提交走 OpenAPI `submit --command`，资源按
`ci/kunshan/README.md` 的标准 suite 及当前作业目的选择。

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
| `The imported target "hip::galaxyhip" references "/opt/rocm/hip/lib/libgalaxyhip..." but this file does not exist` | 没设 `ROCM_PATH`，DTK 的 `hip-targets.cmake` 前缀回落 `/opt/rocm` | `export ROCM_PATH=<DTK>` 后重新 configure（见 §1 配方） |
| `Could not find a package configuration file provided by "amd_comgr"` | `enable_language(HIP)` 的 try_compile 子工程没继承目录 | 传 `-Damd_comgr_DIR/-DAMDDeviceLibs_DIR/-Dhsa-runtime64_DIR` **且** `-DCMAKE_TRY_COMPILE_PLATFORM_VARIABLES="CMAKE_PREFIX_PATH;amd_comgr_DIR;AMDDeviceLibs_DIR"` |
| `/bin/sh: .../cmake-3.30.5-linux-x86_64/bin/cmake: No such file` | `make` 用缓存里的旧 cmake 路径，且 `cmake` 未在 `PATH` | `source scripts/env.sh`（已把 `$CMAKE` 的目录前置到 `PATH`） |
| `ldd` 找不到 `libmpi.so.40` | 未 source `scripts/env.sh` | 所有作业脚本第一行 source |
| `srun` 下 MPI size=1 | `srun` 不做 MPI 启动 | 改用 `mpirun` |
| `bnd_resize` 报 10 ms/step | 首个 `hipMalloc` 被平均进短基准 | 基准先 warmup 至少一步 |
| 多 zone 卡死 / 数值发散的假象 | 每 rank 的 cell 数太小 | 固定 per-rank cell 数再比较 |
