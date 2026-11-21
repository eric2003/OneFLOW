# OneFLOW Living TODO

本页跟踪 [Phase Plan](oneflow-development-plan.md) 的当前切片。父项和子项沿用计划编号；更细的执行步骤记在对应子项下。当前验收完成后移出本页，将复验结论和证据索引留在计划或对应报告。

## 当前切片

分支 `dev` HEAD `1395e3056`（含 `aada33da2`/`169d59892`/`d41650729` 及 SA/viscous hook 接线修复）。**2026-11-21 已在昆山用当前 HEAD 重跑旗舰全开 hook**：`src/REVISION`、binary sha256、五个 device hook 字符串均核验通过。实际投放集群前仍须核对被测对象。固定路径、探针限制和构建命令见[当前交接](../handoff/oneflow-current-handoff.md)。

**关键纠正（2026-11-21）**：此前所有"CD 22%/CD_SF 143%"及 SA hook rc=134 崩溃的昆山证据都来自 11-16 修复前的旧 binary（`UTurbSrcFlux.cpp` 报错行标注 "Compiled on Nov 16"），不能作为当前源码结论。证据在 `evidence/gate/128189103-vdump`（粘性 kernel 整场逐位一致）、`128169609-bisect`（旧 SA hook VMFault）、`128191406-conv`（旧对照：device 臂 rc=1 只跑 100 步、host 跑满 400 步，属伪比较）。

- [ ] **B.3.1 验证粘性通量布局和写回**（Partial，单卡写回及自检已运行）
  - [x] `aada33da2` 将设备输出 stride 改为 `nFaces`，按连续方程平面写回求解器字段。
  - [x] 共享 view 校验要求 stride ≥ nFaces；契约测试覆盖紧凑布局和重叠布局（15/15）。
  - [x] CPU 粘性策略通过 24/24；本地 CPU 全量构建通过，CTest 1,367/1,374，7 项在 macOS 的 system 目录定位处失败。
  - [x] 昆山 1 rank / 1 DCU，`cc1dd4b3a` 构建及 HIP 粘性契约 6/6。M6 自检确认 hook 命中。第 1 步共 896,256 面，仅面 3356 的能量分量被报告，绝对差约 `4.38e-40`，面归一化差 0.125327；第 2、3 步在 `1e-10` 阈值内无报告。相对判据在近零面上会放大差异，需同时看绝对差。
  - [x] `169d59892` 回应上一条：抽出 `CompareViscousFaceFlux`（`codes/uns/include/ViscousFluxSelfCheck.h`），显式统计非有限值（相对判据会让 NaN 静默通过）、按 `bcType == SOLID_SURFACE` 标注壁面，输出步/区/首差面/首差方程/最差面、`absdiff` 及末行 `SELFCHECK result=PASS|FAIL`。新测试 10/10，既有契约 15/15、粘性策略 24/24。
  - [x] 2026-11-21 当前 HEAD 在昆山实证：加固 selfcheck 命中并运行（`evidence/gate/128627452-allhooks/probe-S20/hip.log`）。前两个 RK stage 的 `result=FAIL` 仅因近零面相对判据放大（首差面 3356 `absdiff=1.2e-41`、面归一化 rel=0.125），`nonfinite_device/host=0`；后 3 个 stage `result=PASS`、逐位一致。**粘性 device kernel 整场逐位等价已确认**，与整场 dump `128189103-vdump` 一致。
- [ ] **B.3.2 复验 M6 气动力**（Partial，当前 HEAD 全开 hook 有效对照已跑，根因未收口）
  - [x] 历史 200 步隔离对照将偏差定位到粘性 hook。旧版本 400 步全开/全关的 `CD` 差 22%，`CD_SF` 差 143%——**已确认是修复前旧 binary 的伪比较**（device 臂 rc=1 仅 100 步 vs host 400 步）。
  - [x] 同一 HIP 二进制、同一算例、关闭自检与探针，全关和仅开粘性均跑满 400 步。两臂的 `aero.dat`、`res.dat` 分别逐字节一致。末步 CL=`1.4097e-1`、CD=`2.1895e-2`、CD_SF=`3.7910e-3`、res1=`2.81887e-7`。
  - [x] 2026-11-21 当前 HEAD 全开 hook 对照（`evidence/gate/128627452-allhooks`）：CPU oracle 与 DCU 五 hook 全开**两臂均 rc=0、跑满 400 步、aero 各 37 行**（有效同判据比较，非旧伪比较）。偏差量级较旧证据大幅缩小：`CD_SF` step20 仅差 `5.18e-05`（≈0.5%，stride 修复生效），但 `CL` step20 差 `9.42e-03`（≈25%）、`res1` step20 差 `9.31e-07`（≈18%）。残差轨迹形状不同（CPU 先降后升、DCU 单调降），从首个输出点即系统性偏离。
  - [x] 逐面 wall force probe（`128627452-allhooks/probe-S20`）：壁面 `u_lc` 差 ≈0.86%、梯度 `dudy` 差 ≈256（0.86%）、`viscosity` 差 0.17%，集中在壁面/对称面（face 3457/3459/3461）；q[rc] 与 q[lc] 反号对称，比较有效 → §1f 的 ghost 陈旧**不是本轮主因**（对称性证明 q[rc] 已刷新）。
  - [x] 2026-11-21 stage trace（`128629176`，8 步，`nResSave=1`，trace 文件 download 到本地流式比对，集群 verifier 因 1.5GB×2 文件 OOM）：**step1 全部 stage 逐位一致；step2 seq=1 起 `qf1/qf2/invflux/inviscid_residual/state` 同时分叉**；`inviscid_residual` first_idx=22 cpu=1.28e-09 vs dcu=-5.5e-13（符号反、量级差 3 个数量级）；`state` worst_rel 从 step2 的 0.001 累积到 step8 的 0.0074。逐残差 step1 rel=0、step2 突变 9.99e-07 并保持。
  - [x] 根因定位到**第 1 步末的状态更新/时间层推进**，非 kernel 算术。代码定位：`codes/uns/src/UNsUpdate.cpp:219-262` 的 `UpdateFlowField`——`UseDeviceAuthoritativeHipState()`（`:210-215`，`startStrategy!=2` 时启用）走 `HipFluxBackend::UpdatePrimitiveState` GPU 更新（`:248-249`，primitive q 留设备、host 按需镜像）；host 路径走逐 cell `GetUpdateField+PrepareData+UpdateFlowFieldValue`（`:265+`）。step1 两者一致，step2 读回时 device authoritative 路径产生 1e-8 差异并被 RK/LU-SGS 放大。粘性 kernel（selfcheck/vdump）与 SA hook 均逐位一致，已排除。**文档坑 41/§1f 描述的"单 zone 不镜像/ghost 陈旧"在 startStrategy=0 下被设备权威+按需镜像路径取代，不能照搬**。
  - [ ] 补 CPU 五例 normal/strict、HIP 契约与准确提交 CI 证据。收敛后气动力（非固定 400 步）单独验收。
- [ ] **B.3.3 补边界契约**（紧接 B.3.2）
  - [ ] 验证非五方程、非 VIS_STD 和等温壁按预期回落。记录壁面热流是否仍走主机，并补多区 interface 面覆盖后再评估 B.3 收口。

## 证据索引

[历史 Living TODO 快照](https://github.com/lql341/OneFLOW/blob/1eee55499a69e6b9bc7dadf28188f484af0fc983/docs/plans/oneflow-development-todo.md)保留原章节、测量口径和证据位置。粘性看 §5.40–§5.44，H1/H2 看 §5.2，网格阶梯看 §5.6/§5.9。
