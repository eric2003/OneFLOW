# OneFLOW Living TODO

本页跟踪 [Phase Plan](oneflow-development-plan.md) 的当前切片。父项和子项沿用计划编号；更细的执行步骤记在对应子项下。当前验收完成后移出本页，将复验结论和证据索引留在计划或对应报告。

## 当前切片

分支 `dev`。最近已核对的求解器变更为 `aada33da2`，CPU 粘性契约测试修复为 `a3d7dd611`。实际投放集群前核对 `git rev-parse HEAD`、`src/REVISION` 和运行二进制。固定路径、探针限制和构建命令见[当前交接](../handoff/oneflow-current-handoff.md)。

- [ ] **B.3.1 验证粘性通量布局和写回**（Implemented, unverified）
  - [x] `aada33da2` 将设备输出 stride 改为 `nFaces`，按连续方程平面写回求解器字段。
  - [x] 共享 view 校验要求 stride ≥ nFaces；契约测试覆盖紧凑布局和重叠布局（15/15）。
  - [x] CPU 粘性策略通过 24/24；本地 CPU 全量构建通过，CTest 1,367/1,374，7 项在 macOS 的 system 目录定位处失败。
  - [x] 写回自检可信化：抽出 `CompareViscousFaceFlux`（`codes/uns/include/ViscousFluxSelfCheck.h`），
        显式统计非有限值（相对误差判据会让 NaN 静默通过）、按 `bcType == SOLID_SURFACE` 标注壁面，
        并输出步/区/首差面/首差方程/最差面。新测试 `tests/solver/viscous_flux_selfcheck_test.cpp`
        10/10；既有 `oneflow_euler_domain_contract_test` 15/15、`vis_flux_policy_test` 24/24。
        受影响的 `UNsVisFlux.cpp` 已单独编译通过；端到端自检仍待 HIP 作业。
  - [ ] 构建准确提交的 HIP 二进制，在单 rank、单卡的 M6 开启 `ONEFLOW_ENABLE_UNS_HIP_VISFLUX` 和 `ONEFLOW_VISCOUS_SELFCHECK`，确认 hook 命中及两侧输入时间层相同。以 `SELFCHECK result=PASS|FAIL` 为准，记录首个 FAIL 的步、面、方程与对应输入。
- [ ] **B.3.2 复验 M6 气动力**（Blocked，等 B.3.1）
  - [x] 历史 200 步隔离对照将偏差定位到粘性 hook。旧版本 400 步全开/全关的 `CD` 差 22%，`CD_SF` 差 143%。
  - [ ] 关闭探针，按同配置跑全关与仅开粘性 400 步，记录 CL、CD、CD_SF、残差和停止条件。若气动力仍分叉，在同一输出步比较壁面 ghost、梯度及剪切。
  - [ ] 完成收敛后的气动力对照，并补全开、CPU 五例 normal/strict、HIP 契约与准确提交 CI 证据。400 步结果用于复现历史偏差，收敛判断单独记录。
- [ ] **B.3.3 补边界契约**（紧接 B.3.2）
  - [ ] 验证非五方程、非 VIS_STD 和等温壁按预期回落。记录壁面热流是否仍走主机，并补多区 interface 面覆盖后再评估 B.3 收口。

## 证据索引

[历史 Living TODO 快照](https://github.com/lql341/OneFLOW/blob/1eee55499a69e6b9bc7dadf28188f484af0fc983/docs/plans/oneflow-development-todo.md)保留原章节、测量口径和证据位置。粘性看 §5.40–§5.44，H1/H2 看 §5.2，网格阶梯看 §5.6/§5.9。
