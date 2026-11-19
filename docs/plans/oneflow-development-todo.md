# OneFLOW Living TODO

本页跟踪 [Phase Plan](oneflow-development-plan.md) 的当前切片。父项和子项沿用计划编号；更细的执行步骤记在对应子项下。当前验收完成后移出本页，将复验结论和证据索引留在计划或对应报告。

## 当前切片

分支 `dev`。昆山复验对应求解器提交 `cc1dd4b3a`（含 `aada33da2` 的平面写回修正），其后 `169d59892` 加固写回自检（本地验证，**尚未上传集群**）。实际投放集群前核对 `git rev-parse HEAD`、`src/REVISION` 和运行二进制。固定路径、探针限制和构建命令见[当前交接](../handoff/oneflow-current-handoff.md)。

- [ ] **B.3.1 验证粘性通量布局和写回**（Partial，单卡写回及自检已运行）
  - [x] `aada33da2` 将设备输出 stride 改为 `nFaces`，按连续方程平面写回求解器字段。
  - [x] 共享 view 校验要求 stride ≥ nFaces；契约测试覆盖紧凑布局和重叠布局（15/15）。
  - [x] CPU 粘性策略通过 24/24；本地 CPU 全量构建通过，CTest 1,367/1,374，7 项在 macOS 的 system 目录定位处失败。
  - [x] 昆山 1 rank / 1 DCU，`cc1dd4b3a` 构建及 HIP 粘性契约 6/6。M6 自检确认 hook 命中。第 1 步共 896,256 面，仅面 3356 的能量分量被报告，绝对差约 `4.38e-40`，面归一化差 0.125327；第 2、3 步在 `1e-10` 阈值内无报告。相对判据在近零面上会放大差异，需同时看绝对差。
  - [x] `169d59892` 回应上一条：抽出 `CompareViscousFaceFlux`（`codes/uns/include/ViscousFluxSelfCheck.h`），显式统计非有限值（相对判据会让 NaN 静默通过）、按 `bcType == SOLID_SURFACE` 标注壁面，输出步/区/首差面/首差方程/最差面、`absdiff` 及末行 `SELFCHECK result=PASS|FAIL`。新测试 10/10，既有契约 15/15、粘性策略 24/24；`UNsVisFlux.cpp` 单独编译通过。**该提交尚未在集群运行**。
- [ ] **B.3.2 复验 M6 气动力**（Partial，400 步隔离对照已验证）
  - [x] 历史 200 步隔离对照将偏差定位到粘性 hook。旧版本 400 步全开/全关的 `CD` 差 22%，`CD_SF` 差 143%。
  - [x] 同一 HIP 二进制、同一算例、关闭自检与探针，全关和仅开粘性均跑满 400 步。两臂的 `aero.dat`、`res.dat` 分别逐字节一致。末步 CL=`1.4097e-1`、CD=`2.1895e-2`、CD_SF=`3.7910e-3`、res1=`2.81887e-7`。
  - [ ] 下一步用同一套配置跑全开钩子对照，明确收敛判据后复核气动力；补 CPU 五例 normal/strict、HIP 契约与准确提交 CI 证据。400 步结果用于复现历史偏差，收敛判断单独记录。
- [ ] **B.3.3 补边界契约**（紧接 B.3.2）
  - [ ] 验证非五方程、非 VIS_STD 和等温壁按预期回落。记录壁面热流是否仍走主机，并补多区 interface 面覆盖后再评估 B.3 收口。

## 证据索引

[历史 Living TODO 快照](https://github.com/lql341/OneFLOW/blob/1eee55499a69e6b9bc7dadf28188f484af0fc983/docs/plans/oneflow-development-todo.md)保留原章节、测量口径和证据位置。粘性看 §5.40–§5.44，H1/H2 看 §5.2，网格阶梯看 §5.6/§5.9。
