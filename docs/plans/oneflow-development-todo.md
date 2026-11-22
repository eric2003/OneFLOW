# OneFLOW Living TODO

本页跟踪 [Phase Plan](oneflow-development-plan.md) 的当前切片。父项和子项沿用计划编号；更细的执行步骤记在对应子项下。当前验收完成后移出本页，将复验结论和证据索引留在计划或对应报告。

## 当前状态（2026-11-22 收口）

分支 `dev`，HEAD `5ddb72f72`，工作树干净，已推送 `origin/dev`。旗舰算例 `test/cfd/m6wingroe_sa`（vismodel=3 SA / Roe ischeme=1 / LU-SGS time_integral=2 / ilim=0 / 单 zone / ivischeme=0）在昆山单 rank 已跑通**完整 GPU 通量路径 + 完整 SA 链路**。固定路径、构建命令见[当前交接](../handoff/oneflow-current-handoff.md)，架构细节见[3D device residency 设计](../architecture/oneflow-3d-device-residency-design.md)。

### 已达成（准确性）

- **气动力 aero.dat 逐字节一致**（全部 20 步，max_rel=0）；主控残差 res.dat max_rel `6.78e-21`（gate 1e-12 pass）；湍流残差 turbres 末位 `1.8e-16`；当前流场 flow.dat level1 `2.7e-15` 绝对；湍流场 turbflow 全层 `1.4e-17`；可视化 visual.dat 差异全为机器零（1e-37）。证据：`xb-dlfinal/128674985`、`xb-final-check/128676813`。
- **HIP 重跑确定性**：同 binary 同配置两次，aero/res/turbres 逐字节相同（`xb-determinism/128677591`）。
- **SA 全链路设备化后数值等价**：aero 仍 max_rel=0；主控残差与湍流残差与移植前**逐位相同**（`1.808152e-12` / `8.168e-14`）。证据：`128693499-saport`、`128694886-saport`。
- **flow.dat level2（历史槽）差 7.58e-3**：= 一步物理变化，HIP `InitUnsteady` 重建历史层 vs CPU 保留上一步；restart reader 用 level1 重建 current，**不影响续算**（`xb-history-shift/128676991`）。

### 已达成（性能，20 步单 rank 实测）

所有数字都取**同一次作业、同一节点**的 CPU/HIP 两臂对照；跨节点比较不成立（实测同配置不同节点墙钟可差 8%，见下）。

- 端到端 wall **77.09s → 35.19s（2.19x）**（`128693499-saport`，节点 b13r2n04）。
- 通量层 rhs **7.9x**（无粘 23x、粘性 3.9x）；迭代段 1.40x。
- **boundary BC device 化：~20000ms → ~90ms（~220x）**，aero 仍 0（`xb-boundary/128680850`）。
- **SA 链路 device 化：`update_residuals` 33.8s → 18.5s（1.83x）**、`calc_boundary` 2.1s → 1.1s。

### 已知不达预期（实测，勿再按旧结论推理）

- **LU-SGS 设备化比 host 慢**。同运行同节点对照（`128693499`，b13r2n04）：CPU 臂上下扫合计 **43.05s**、HIP 臂 **45.91s**，慢 6.6%；去掉无用上传后的复测（`128694886`）HIP 44.4s，仍慢于同节点 CPU 的 43.9s。根因是网格 48×48×128 的依赖层级约 221 层、每层仅约 1300 单元，每步 2×221 次内核启动的开销超过并行收益。当前按用户决定**保留启用**，是下一个要处理的大项。
- **门禁口径陷阱**：`f2-main-solver-cross-backend-residual.sh:246` 的相对差是 `|d|/max(1,|a|,|b|)`；残差量级 1e-6，故 "1e-12 门禁" 实际是 **1e-6 相对**。`residual 1.808e-12 fail` 对应真实相对差 **5.3e-7**（step20），不是"接近机器精度只差末位"。
- **`tasks.tsv` 时间是真实值的 2 倍**（`CmxTask.cpp:483` 与 `Command.cpp:295` 各记一次），`stage.tsv` 才是单次。对比两组数字前必须确认口径，否则会得出错误的加速比。
- **节点间性能波动可达 8%**：同一份代码 20 步墙钟 b13r2n04 为 35.19s、b13r3n03 为 38.06s，两臂同步变慢（CPU 77.09→83.42s），是节点差异不是代码回归。**任何加速比必须取自同一次作业的两臂**。

## 进行中 / 待收口

- [ ] **B.4.1 boundary 残差 1.808e-12**（Partial，device 已跑，aero 已 0）
  - [x] 去掉 `HipBoundarySupported` 的 `vismodel!=0` 硬拒绝；kernel 补 viscous 壁面、symmetry ghost、bcQ 覆盖所有 BC。提交 `ab7a40fa8`。
  - [x] `bnd_kernel 2.78 / bnd_ghost_copy 10.04 / boundary_device 89.19ms`，`boundary_calc_face_bc` 归零；aero max_rel=0。
  - [x] **已排除 SA 移植的影响**：SA 设备化前后 res/turbres 逐位相同，证明这 1.808e-12 来自 NS 非壁面 BC，与湍流无关。
  - [ ] **待收：对称面 ghost 反射语义（首选假设，附验证方法）**
    - 证据链：`NsBcSolver::SymmetryBc`（`codes/ns/src/NsBcSolver.cpp:376`）先把 `primt1 = prims1`，**紧接着又用 `ReflectVelocityRelativeToMovingFace` 把 primt1 的三个速度分量覆盖成反射值**；`UNsBcSolver::UpdateBc`（`codes/uns/src/UNsBcSolver.cpp:184`）写回 `q[rc] = primt1`。所以 host 的 symmetry ghost **带反射法向速度**。
    - device kernel（`codes/accel/src/HipFluxBackend.hip`，type==3 分支）按 "ghost = owner" 实现，注释引用的 "NsBcSolver::SymmetryBc sets primt1 = prims1" 是**只读了前半句**。提交 `ab7a40fa8` 的说明沿用了这个误读。
    - 影响路径：ilim=0 时 `NsLimField::BcQlQrFix`（`codes/uns/src/UNsLimiter.cpp:49`）对非壁面取 `0.5*(q[lc]+q[rc])`，ghost 直接进无粘通量→残差。
    - **验证方法**：把 device type==3 分支改成与 wall 相同的反射式（`ghost = u - 2*n*(u·n - vfn)`，rho/p 取 owner），重跑 20 步 gate；若 `residual` 从 `1.808e-12` 掉到 `1e-20` 量级即证实。
    - 若证实，同时要核对 symmetry 的 **bcQ**。已核对：`SetBc` 里 symmetry 分支不改 `updateFlag`（默认 `true`），所以 `UpdateBc` 确实写 bcQ；host 写 `ArithmeticAverage(prims1, primt1)`，device 写 `0.5*(owner+反射)`，**这两者等价**（已核对，非假设）。所以 bcQ 不是差异来源，ghost q 才是。
  - [ ] **次要项：farfield 幂函数快路不一致**。host `FarFieldPower`（`codes/ns/src/NsBcSolver.cpp:44`）的 `x*x*sqrt(x)` 快路只在 `vismodel == INVISCID` 时走；device `BoundaryPower`（`HipFluxBackend.hip:1259`）无条件走。粘性算例上 host 用 `pow(x,2.5)`、device 用 `sqrt` 形式，差 ≤1 ulp（实测 3.7e-16）。修法：给 device 加同一个 `viscousWall` 守卫。
- [ ] **B.4.2 SA 湍流 device 化**（Done，code-complete）
  - [x] 梯度（CSR 重写，去掉 O(nCells×nFaces)）、显式更新、隐式谱、边界、LU-SGS 全部接线；`TurbSaGradientView`/`TurbSaUpdateView` 由孤儿代码变为生产路径。
  - [x] 三个 "已有入口但从未调用" 的根因已修：`CalcGradient`/`Update` 曾把 host 指针直接传给 kernel。
  - [x] 集群验证：编译链接通过，HIP rc=0，数值与移植前逐位等价。
- [ ] **B.4.3 HIP 开关语义统一**（Done）
  - [x] 9 个 `ONEFLOW_ENABLE_UTURB_HIP_*` 并入 `HipProductionPolicy.h` 白名单，backend=HIP 即生效；显式 =0 仍可单独关闭。与 UNS 语义一致。

## 后续（性能提升路径）

- [ ] **B.4.4 LU-SGS（当前最大项，实测 host 更快）**：三条可选路径——(a) 默认关回 host，代码留用；(b) 红黑/多色排序把层数从 221 压到 O(1)，但破坏 Gauss-Seidel 逐位等价，需重新做门禁；(c) 换 GMRES/多重网格等可并行求解器。当前状态：保留启用，44.4s vs host 40.1s。
- [ ] **B.4.5 同步开销**：`gradient_download 206ms / v_download 116ms / v_upload 88ms / gradient_mirror 37ms`（20步累计）是为 oracle 对齐付出的回传代价，SA 常驻化后可减少。
- [ ] **B.4.6 更长步数/多 rank**：35.19s 是 20 步含固定开销的数字，需 200+ 步、多 rank 复测端到端。

## 证据索引

- 准确性 gate：`evidence/gate/xb-dlfinal/128674985`（res 6.78e-21, aero 0）
- 三项输出：`xb-final-check/128676813`、`xb-history-shift/128676991`
- 确定性：`xb-determinism/128677591`
- boundary device：`xb-boundary/128680850`（220x，res 1.8e-12 待收）
- **SA 全链路**：`128693499-saport`（首版，b13r2n04）、`128694886-saport`（去无用上传后，b13r3n03）
- 早期粘性：§5.40–§5.44；H1/H2：§5.2；网格阶梯：§5.6/§5.9

## 接手第一条命令

集群作业脚本已经在位，不用重写：

```bash
# 1) 打包当前工作树（未提交改动也要进包）
cd <repo> && SHA=$(git stash create) && \
  git archive --format=tar.gz --prefix=OneFLOW/ "$SHA" -o /tmp/of-sa.tgz

# 2) 上传（超时常见，重试；成功后核对字节数）
python3 ~/.agents/skills/scnet-hpc/scripts/scnet.py --region 11250 \
  upload /tmp/of-sa.tgz /public/home/luql/oneflow/snapshots --cover

# 3) 提交：脚本内解包→构建 RuntimeObjects+fast→跑 20 步 CPU/HIP gate
python3 ~/.agents/skills/scnet-hpc/scripts/scnet.py --region 11250 submit \
  --name of-saport --command "bash /public/home/luql/oneflow/scripts/of-saport.sh" \
  --work-dir /public/home/luql/oneflow/logs --queue kshdnormal \
  --nodes 1 --cpus 8 --dcus 1 --memory 27G --walltime 02:00:00

# 4) 等结果，evidence 落在 evidence/gate/<jobid>-saport/
python3 ~/.agents/skills/scnet-hpc/scripts/scnet.py --region 11250 wait <jobid>
```

注意：**memory 上限 27G**（32G 会被 DefMemPerCPU 拒）；`of-saport.sh` 末尾已经把 `cross-backend.tsv`、两臂 stage 表、HIP 日志尾都打进 `result.txt`，不必单独取。

[历史 Living TODO 快照](https://github.com/lql341/OneFLOW/blob/1eee55499a69e6b9bc7dadf28188f484af0fc983/docs/plans/oneflow-development-todo.md)保留原章节、测量口径和证据位置。
