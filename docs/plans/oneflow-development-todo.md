# OneFLOW Living TODO

本页跟踪 [Phase Plan](oneflow-development-plan.md) 的当前切片。父项和子项沿用计划编号；更细的执行步骤记在对应子项下。当前验收完成后移出本页，将复验结论和证据索引留在计划或对应报告。

## 当前状态（2026-11-22）

分支 `dev`，HEAD `7063145ee`，工作树干净。旗舰算例 `test/cfd/m6wingroe_sa`（vismodel=3 SA / Roe ischeme=1 / LU-SGS time_integral=2 / ilim=0 / 单 zone / ivischeme=0）在昆山单 rank 已跑通**完整 GPU 通量路径**并达到机器精度等价。固定路径、构建命令见[当前交接](../handoff/oneflow-current-handoff.md)，架构细节见[3D device residency 设计](../architecture/oneflow-3d-device-residency-design.md)。

### 已达成（准确性）

- **气动力 aero.dat 逐字节一致**（全部 20 步，max_rel=0）；主控残差 res.dat max_rel `6.78e-21`（gate 1e-12 pass）；湍流残差 turbres 末位 `1.8e-16`；当前流场 flow.dat level1 `2.7e-15` 绝对；湍流场 turbflow 全层 `1.4e-17`；可视化 visual.dat 差异全为机器零（1e-37）。证据：`xb-dlfinal/128674985`、`xb-final-check/128676813`。
- **HIP 重跑确定性**：同 binary 同配置两次，aero/res/turbres 逐字节相同（`xb-determinism/128677591`）。
- **flow.dat level2（历史槽）差 7.58e-3**：= 一步物理变化，HIP `InitUnsteady` 重建历史层 vs CPU 保留上一步；restart reader 用 level1 重建 current，**不影响续算**（`xb-history-shift/128676991`）。

### 已达成（性能，20 步单 rank 实测）

- 端到端 wall **75.16s → 61.26s（1.23x）**；迭代段顶层 **1.40x**。
- 通量层 rhs **7.9x**（无粘 23x、粘性 3.9x）。
- **boundary BC device 化：~20000ms → ~90ms（~220x）**，aero 仍 0（`xb-boundary/128680850`）。

## 进行中 / 待收口

- [ ] **B.4.1 boundary 残差 1.8e-12**（Partial，device 已跑，aero 已对齐）
  - [x] 去掉 `HipBoundarySupported` 的 `vismodel!=0` 硬拒绝；kernel 补 viscous 壁面（ghost=`-u+2u_wall`、bcQ=`(rho,0,0,0,p)`）；symmetry ghost=owner（匹配 host `primt1=prims1`）、bcQ=(owner+反射)/2；bcQ 覆盖所有 BC。提交 `ab7a40fa8`。
  - [x] `bnd_kernel 2.78 / bnd_ghost_copy 10.04 / boundary_device 89.19ms`（20步），`boundary_calc_face_bc` 归零；aero max_rel=0。
  - [ ] **residual max_rel 1.808e-12（eq6，超 1e-12 门禁 80%）待收**。已排除 gamma（nTModel=1 常 gamma）、cref（两边 `sqrt(gama_ref*pref/rref)`）、BoundaryPower（exponent=1.4 不走 2.5 近似）。下一步：逐 face dump device vs host ghost q/bcQ 定位非壁面 BC（symmetry/farfield/outflow）末位差异。
- [ ] **B.4.2 SA 湍流 device 化**（Not started，第二大头）
  - agent 查明：`update_residuals` 32840ms 中 rhs 仅 18063ms，其余 **~14777ms 两臂相同 = SA 链路（TurbCalcRHS）**，四个 `ONEFLOW_ENABLE_UTURB_HIP_{VIST,SOURCE,INVFLUX,VISFLUX}` **默认关闭**（不在自动开启白名单）。
  - LU-SGS sweep 实测仅 ~2s（≤4%），**不值得移植**；`numberOfSweeps` 依赖零初始化（`NsLusgs.cpp:37-53` 未赋值）。
  - 待办：先零代码开四个 env 试跑（`SupportsFaceFlux`/`iturb_visflux=2` 条件已满足），不行再加代码接线。预期 `update_residuals` 大幅下降。
- [ ] **B.4.3 HIP 开关语义统一**（审计完成，待修）
  - 现状三类：**白名单自动开**（`HipProductionPolicy.h:105-111`：UNS_HIP_BATCH/GRADIENT/INTERFACE/RECONSTRUCTION/STATE_UPDATE/BOUNDARY/TIME_STEP）；**需显式=1**（ALLOW_VISCOUS_HIP_BATCH/GRADIENT、UTURB_HIP_\* 全部 4 个）；**默认生效=0 关闭**（HIP_SKIP_\*）。
  - 不一致：①boundary 在白名单却曾被 vismodel 守卫拒绝（已修）；②UTURB 四个 hook 不在白名单，必须手动设（SA 链路留 host 的根因）；③UNS 自动开 vs UTURB 手动开语义不统一。待决策：是否把 UTURB 纳入同类策略。

## 后续（性能提升路径）

- [ ] **B.4.4 LU-SGS 宏步**：评估结论=收益 ≤4%，deferred。除非先补 sweep 计时行（`StageProfiler` 无 `LUSGS_*` 映射，成本被误归到 `initialization`）证明确实值得。
- [ ] **B.4.5 同步开销**：`gradient_download 199ms / v_download 116ms / v_upload 89ms / gradient_mirror 38ms`（20步累计）是为 oracle 对齐付出的回传代价，SA device 化后可减少。
- [ ] **B.4.6 更长步数/多 rank**：1.23x 是 20 步含固定开销的数字，需 200+ 步、多 rank 复测端到端。

## 证据索引

- 准确性 gate：`evidence/gate/xb-dlfinal/128674985`（res 6.78e-21, aero 0）
- 三项输出：`xb-final-check/128676813`、`xb-history-shift/128676991`
- 确定性：`xb-determinism/128677591`
- boundary device：`xb-boundary/128680850`（220x，res 1.8e-12 待收）
- 早期粘性：§5.40–§5.44；H1/H2：§5.2；网格阶梯：§5.6/§5.9

[历史 Living TODO 快照](https://github.com/lql341/OneFLOW/blob/1eee55499a69e6b9bc7dadf28188f484af0fc983/docs/plans/oneflow-development-todo.md)保留原章节、测量口径和证据位置。
