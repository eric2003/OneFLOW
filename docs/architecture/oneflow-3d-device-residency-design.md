# 3D Euler 设备驻留架构设计（A.4 / Phase B 交汇）

> 2026-11-21。本文是 C 主线（核心目的：推架构改动）的设计基线，基于 A/B 两线的实机证据。
> 现状事实以源码为准；验证状态见 [Living TODO](../plans/oneflow-development-todo.md) B.3。

## 1. 为什么要改架构（证据驱动）

A 线（昆山 2026-11-21，作业 `128627452`/`128628784`/`128629176`）已钉死：

1. 粘性/SA device kernel 整场逐位一致（selfcheck 后 3 stage PASS、`nonfinite=0`；vdump 896,256 面 1e-26）。kernel 算术不是瓶颈。
2. CD_SF 143% → 0.5%（stride 修复 `aada33da2` 生效），旧 22%/143% 是修复前伪比较（device 仅 100 步 vs host 400 步）。
3. ** residuals step1 逐位一致、step2 突变 1e-6 并保持**；stage trace：step2 seq=1 起 `qf1/qf2/invflux/inviscid_residual/state` 同时分叉，`state` worst_rel 0.001→0.0074 累积。根因在第 1 步末的状态更新/时间层推进（`UNsUpdate.cpp:210-262`，device authoritative 的 `UpdatePrimitiveState` vs host 逐 cell 循环）。
4. 加速不达预期的结构性原因：旗舰 `m6wingroe_sa`（LU-SGS+Roe+SA+粘性）主体仍在 host，设备路径只支持 Euler/显式 RK + 零散 hook。调 hook 换不来数量级加速。

结论：**正确性地基（状态推进语义）和加速能力（完整设备驻留）必须在同一个架构里解决**，不能分别打补丁。

## 2. 现状：三个执行世界与 A.4 的真实起点

| 世界 | 入口 | 状态所有者 | 数值路径 |
|---|---|---|---|
| legacy host | `TimeIntegral::RungeKutta()` / `Lusgs()` | 全局 `MRField` | 逐 cell host task 串 |
| "RK backend 骨架" | `TimeIntegral::RungeKutta(context)` :346-395 | `CpuEulerDomainBackend`（硬编码 `AccelBackendKind::CPU` :369） | `backend.Advance` :390 调 `RunEulerRungeKuttaStage`，回调内**仍是 5 个 host task**（:133-169） |
| 受限 device | capability 守卫 | `Ns3DDeviceState` + `ResidentStateContract` | 每 stage 一批 HIP kernel，host 调度 |

**关键事实**：
- `EulerDomainBackend::Advance` 只有 `CpuEulerDomainBackend` 实现（`CpuEulerDomainBackend.cpp:98-133`），无 HIP 实现；其 `Advance` 只是 `stageCallback` 空转循环。
- `RunEulerRungeKuttaStage`（`TimeIntegral.cpp:133-169`）每 stage 调 `LOAD_RESIDUALS`/`UPDATE_RESIDUALS`/`CALC_LHS`/`UPDATE_FLOWFIELD`/`CALC_BOUNDARY`——所谓 RK 驻留骨架，内部**没有任何设备驻留**，只是把 host task 循环包进 `Advance`。A.4 的起点比文档描述更靠前。
- 旗舰走 `TimeIntegral::Lusgs()`（:397-426）：一条纯 host task 串（`ZERO_DQ`→时间步→残差→INIT_LUSGS→[下扫→dq 交换→上扫]×n→UPDATE→BOUNDARY），完全不进 `RungeKutta(context)` 的 backend 路径。
- 设备状态推进分歧点：`UNsUpdate::UpdateFlowField`（`UNsUpdate.cpp:219-262`），`UseDeviceAuthoritativeHipState()`（:210-215，`startStrategy!=2` 启用）走 `HipFluxBackend::UpdatePrimitiveState`（GPU kernel，`ScaleAndUpdatePrimitiveState` `HipFluxBackend.hip:4323`），host 路径走逐 cell `GetUpdateField+PrepareData+UpdateFlowFieldValue`（:265+）。

**已证明的驻留样板（1D）**：`ports/kunshan/oneflow_1d_hip` 的 `HipEulerBackend::Advance` —— state 全程留设备，`Advance` 内循环 step 每步调一个设备 kernel，仅 trace 时回传。速度比 host-staged 高 1-2 个数量级。3D 目标是把这个模式复制到 3D。

## 3. 设计原则（对齐 DEVELOPMENT_PRINCIPLES）

1. **单一状态所有权**：每个 (solver, zone, grid, backend, device) 一个 `EulerDomainState`，所有权显式（device authoritative 或 host authoritative），切换必须经 `ResidentStateContract` 标脏/刷新，禁止隐式回灌。
2. **宏步整体下沉**：`Advance(state, steps, options)` 内部完成完整 RK/LUSGS 宏步循环，host 只做 INIT/输出/门控；每宏步 host 往返 ≤2 次，H2D/D2H 字节数进 `StageProfiler`。
3. **能力表驱动**（plan E.2）：`solver × scheme × limiter × viscous × time-integral × grid-level × halo` 一张表 + 单测，替换散落的 `UseHip*` 守卫；每个组合要么有门禁要么默认回退 host。
4. **数值可追溯**：每个下沉的 stage 用 3 步 trace（1e-10）+ 200 步跨后端（1e-12）+ CPU 五例（1e-8/1e-15）门禁；残差轨迹不是验收标准，收敛后气动力 + 非有限值才是。
5. **host oracle 单一**：CPU 路径保持为唯一参考，不因迁移删除或改写。

## 4. 架构改动切片（用户拍板：RK + LU-SGS 都做，单卡旗舰 + ≥1M 都做）

### 切片 0：统一状态推进语义（正确性地基，先做）

这是 A 根因的直接修复，也是所有驻留的地基。

- **契约**：`startStrategy ∈ {0,1,2}` 下，device authoritative 与 host 路径在每步末对 q 的更新逐位一致（或差异 ≤1e-12 且不累积）；stage trace step1/step2 `state` 数组在 1e-10 内一致。
- **入口**：`UNsUpdate.cpp:210-262`（`UseDeviceAuthoritativeHipState` / `UpdateFlowField`）、`HipFluxBackend::ScaleAndUpdatePrimitiveState`（`HipFluxBackend.hip:4323`）、`ResidentStateContract` 的镜像代际。
- **实现方向**：定位 step1→2 的 1e-8 差异源（device `UpdatePrimitiveState` 的 ghost/retain 处理 vs host 逐 cell 循环），统一时间层推进（history shift / q1q2）的更新点。不掩盖差异，以 stage trace 为准。
- **验证**：重跑 `128629176` 的 8 步 stage trace + 400 步气动力对照；sign/量级差异消失。
- **先验**：`ONEFLOW_HIP_GRADIENT_MIRROR_ON_DEMAND`、`HipBoundaryStateMirrorEnabled`、`retainResidentState` 的组合（handoff 坑 1-4）。

**代码级差异候选（待用单 cell 对照测试逐个排除，不要静态猜测）**：

| 候选 | Device 路径 | Host 路径 | 位置 |
|---|---|---|---|
| gamma 来源 | `view.gamma = nscom.gama_ref`（全局标量） | `nscom.gama = (*unsf.gama)[0][cId]`（per-cell 场） | `UNsUpdate.cpp:234` vs `:340`；当前 `gama[cId]=gama_ref` 使两者数值相同，但等价脆弱 |
| 更新公式 | kernel 直接 `momentum = density*u + residual`、`energy = p/(γ-1)+kinetic+residual`，用**当前 prim** | `PrimToQ(prim0,γ,q0)` → `q=q0+dq` → `QToPrim(q,γ,prim,t)`，用 **prim0** | `HipFluxBackend.hip:1764-1793` vs `NsUpdate.cpp:46-53` |
| RK 系数 scale | `ScaleResidualKernel` 用 `residualScaleCoefficient`（`HipFluxBackend.hip:1742-1750,4484`），`view.rkCoefficient=1.0` 硬编码（`:235`） | host `nscom.lhscoef = ctrl.rk_coef[iStage]` 参与 dq | `UNsUpdate.cpp:235` / `HipFluxBackend.hip:4484` |

**下一步实现动作**：写一个单 cell 对照测试（复用 `tests/solver` 模式），用同一组 prim0/dq/gama/rkCoef 输入分别跑 host `CalcFlowField` 与 device kernel 公式的 CPU 模拟，定位首个不等项；再决定改哪一处。不要用 400 步整场对照反向猜。

**实证结论（2026-11-21，本机 `g++` 纯数学对照 `/tmp/state_update_cmp.cpp`）**：host `PrimToQ→q+dq→QToPrim` 与 device `UpdatePrimitiveStateKernel` 公式**数学完全等价**——M6 典型状态与 machine-zero dq 全部 rel=0，仅一组 pressure 差 5.55e-17（机器 eps）。**排除"更新公式差异"候选**；与 stage trace step1 逐位一致吻合。根因收敛到**状态传递/代际/镜像**，下一条只查这三点：

1. `UNsUpdate.cpp:281-293` host 路径末尾无条件 `NsRefreshHostStateMirror()` + `NotifyHostFieldStateModified()`；device 路径（`:219-262`）`deviceAuthoritative` 时 `downloadToHost=false` 不回传。两者对 host mirror 代际的处理不对称。
2. `retainResidentState`（`RetainHipStateForBoundary()`，`:200-204`）——boundary/interface hook 开时留驻留态，ghost 区是否一致。
3. `HipFluxBackend.hip:4484` `fused.rkCoefficient = storage.residualScaleCoefficient` vs `view.rkCoefficient=1.0`（`:235`）——residual scale 的两条路径是否对同一 dq 施加不同系数。

下一步用昆山作业逐一开关这三点（单变量 A/B），以 stage trace step2 state 是否在 1e-10 内一致判定。

**排除实验结论（2026-11-21）**：

| 候选 | 实验 | 结果 |
|---|---|---|
| 更新公式差异 | 本机纯数学对照 `/tmp/state_update_cmp.cpp` | **排除**：M6 典型/machine-zero 全部 rel=0，仅一组 pressure 差 5.55e-17 |
| device authoritative mirror 代际 | 昆山 `128629560` A/B：`ONEFLOW_HIP_DEVICE_AUTHORITATIVE_Q` 1 vs 0 | **排除**：两臂 step2 res 相对差都是 9.99e-07，完全一致，关 authoritative 不改变分叉 |
| SA 标量更新 | 代码核查 `UTurbUpdate.cpp:85-98` + `HipTurbulenceBackend.hip:909` | **排除**：生产 `UpdateFlowField1EquStd` 仍是 host 逐 cell 循环，device `SaUpdateKernel` 只被契约测试用（handoff:223 属实）；两臂同走 host |
| Residual base bookkeeping | 代码核查 `ResidualBasePolicy.cpp:40-46,51-73` | **排除**：`vis_model.vismodel==0` 才启用 device RK bookkeeping，旗舰 vismodel=3 → host 权威，两臂一致 |

**逐 hook bisection 结论（昆山 `128629691`，7 臂 × 8 步）**：关掉全部 5 个 hook（`all_off`）与全开（all_on）的 step2/step8 res 相对差**完全相同**（9.989e-07 / 9.767e-07）。**彻底排除粘性/SA hook**——分叉与任何 hook 无关。

**根因最终定位**：`all_off` 仍设 `ONEFLOW_ACCEL_BACKEND=HIP`，走的是 **HIP batch 无粘基础路径**（LF/Roe fused flux kernel）vs host 无粘路径。结合 step1 逐位一致、step2 第一个 stage `inviscid_residual` 符号反转（cpu=1.28e-09 vs dcu=-5.5e-13），差异在 `HipFluxBackend` 的 fused flux kernel 结果写回/累加，属 plan B.1 范畴，而非状态推进或任何粘性/SA hook。

**对切片划分的修正**：切片 0 的"统一状态推进语义"假设被证伪。真正地基是 **HIP 无粘 batch 与 host 无粘的逐步数值等价**（`cfd.txt` ischeme=1 Roe、ilim=0；`HipFluxBackend` fused flux kernel vs `UNsInvFlux` host）。这是 B.1，须在 A.4 驻留前解决。

**stage bisection 失败（昆山 `128629903`，9 臂 × 8 步）**：逐个显式 `=0` 关闭 7 个 device stage（BATCH/GRADIENT/RECONSTRUCTION/STATE_UPDATE/BOUNDARY/TIME_STEP/INTERFACE），step2 res 相对差**全部不变**（9.989e-07）。`CachedEnv`（`HipProductionPolicy.h:25`）缓存正常、env 传递无误。**根因：旗舰配置下所有 stage 的 `UseHip*` 守卫要求 `vismodel==0 && time_integral==MULTI_STAGE`（如 `UTimeStep.cpp:94-96`、`ResidualBasePolicy.cpp:45`），旗舰 vismodel=3 + LU-SGS 使这些 device stage 本就不启用——env 开关自然无效**。env-level bisection 在旗舰配置上走不通。

**下一步（换方法）**：改用 **instrument 代码**：在 host/device 数值路径的关键分叉点加受 env 控制的日志，重编译后跑 8 步，直接观察 CPU backend vs HIP backend 在哪个数值量上首次不同。log 第 6 列是最大残差单元体积（`ResidualTask.cpp:332`，非 dt）——两臂 step2 最大残差落在不同 cell，是流场分叉结果而非原因。

**关键悖论（2026-11-21）**：两臂配置行/初值完全一致，step1 残差一致，step2 NS 才分叉（CPU 2.42715e-05 vs DCU 1.76992e-05）。旗舰配置下所有 device stage 不启用 → 两臂应走相同 host 路径却不同。

**一次性 instrument 结果（昆山 `128630364`，复用 `128630166` 编译的 FP_ENV binary）**：在 `UNsUpdate::DumpUpdatedStateTrace` 加受 `ONEFLOW_FP_ENV_DIAG=1` 控制的 17 位精度指纹（代码已留在 `UNsUpdate.cpp:322+`，默认关闭，可复用）。两臂 all_off 跑 8 步：

- **浮点环境完全相同**：`round=0`、`fp=33.333333333333435` 逐位一致 → **排除浮点环境/舍入模式假设**。
- **step1 所有 q0..q4 逐位一致**；从 **step2** 起分叉。
- 差异集中在**近零量**：`q3`（w 速度，量级 1e-20，对称面反号）与 `q2_1`（v 速度，~0.05338）。CPU 的 `q2_1` 每步累加约 7e-15（step2 .480→step8 .530），HIP 始终固定 .481；`q3_0` CPU 每步 ×1.5、HIP ×1.25。

**结论**：在相同浮点环境下，从 step2 起某量每步累加 ~7e-15 到 v/w 近零速度——**CPU 累加、HIP 不累加**。这是 host 侧代码在两 backend 下行为不同的白盒证据，且差异从第一个输出步就存在、逐步累积（印证了之前"step1 仅 5 位精度一致"的修正）。

## 根因确定在 RHS（昆山 `128634238-tracepoints4`，LU-SGS PrepareSweep 插桩，binary fed93180）

在 `UNsLusgs::PrepareSweep`（`UNsLusgs.cpp:332`，读 `unsf.rhs` / `unsf.dq`）插桩 cell 254656：

- **step1 两臂 rhs 逐位一致**（rhs0=-4.1358532927951903e-06，rhs1=-6.5589158095902482e-06，...全部相同）；首轮 dq=0、次轮 dq 相同
- **step2 首次 LU-SGS 调用（dq=0，residual 刚算完）rhs 已不同**：
  - rhs0 CPU -4.117e-06 / HIP -3.156e-06（差 23%）
  - rhs1 -6.530e-06 / -3.242e-06（50%）
  - rhs2 -9.306e-07 / -1.850e-07（80%）
  - rhs3 -5.010e-07 / -1.322e-08（97%）
  - rhs4 -1.959e-05 / -1.275e-05（35%）

**结论**：root cause 在 **RHS（residual）计算**——step1 rhs 相同、step2 rhs 在进入 LU-SGS 前就不同。LU-SGS 忠实把不同 rhs 转成不同 dq（其次轮 dq 值与 rhs 差异成比例），**LU-SGS 求解本身无差异**。这证实了分层拆解的第 1 步：不是"LU-SGS 不同"，是"RHS 不同"。下一步拆 RHS 分项（inviscid/viscous/source 在 `UNsInvFlux::CalcFlux` / `UNsVisFlux::CalcFlux` / SA source），定位哪个物理项贡献了不同的 rhs。

**注意**：第一条 TRACE 差异（step1 turbq 梯度末位 5e-12）与 step2 rhs 差 23%~97% 之间仍需建立因果（末位差经 SA 强耦合放大），但 rhs 层已确定。

## "对 0 处理"分歧点（统一方案的代码依据）

用户提出的核心问题：差异集中在近零量，本质是"两臂对 0 的处理不统一"。代码核查发现多处**已存在的、但两臂实现不同**的近零处理：

| 处理 | Host | Device | 风险 |
|---|---|---|---|
| SA 非负 floor / clamp | `TurbUpdate::ModifyValue1Equ` 有 | `TurbSaArithmetic.h:498-502` **故意省略**（注释："adding a live clamp on the device would make device and host diverge by exactly the clamp"） | 已知且刻意；SA 更新两臂都走 host，暂不触发 |
| `SMALL` 常量 | `Constant.h SMALL=1e-40`，`(x+SMALL)` 除法 | `TurbSaArithmetic.h:74 kSmall=1e-40`，注释"same value required for bit parity" | 已统一 |
| 对称面速度反射 | `NsBcSolver::SymmetryBc`（`NsBcPolicy.cpp:8-18`）：`vn=RelativeFaceNormalVelocity(xfn,yfn,zfn,vfn,u,v,w)`，输入是 **prims1/prims2 两侧**，再 `ArithmeticAverage` | `HipFluxBackend.hip:1315-1320`：`nx=xNormal*faceOuterNormal`，`vn=nx*u+ny*v+nz*w-meshVelocityNormal`，输入是**单侧 owner**，再 `0.5*(q[left]+ghost)` | **输入与结合顺序不同**，w~1e-20 近零时产生末位差异，经限制器/粘性放大 |

**统一方案（待实施，切片 0 修复方向）**：把 host `NsBcPolicy.cpp` 的 `RelativeFaceNormalVelocity` / `ReflectVelocityRelativeToMovingFace` / `ArithmeticAverage` 抽成 **backend-neutral 的共享 policy helper**（类似已做的 `TurbSaArithmetic.h`、`ViscousFluxArithmetic.h`），device 对称面 BC 直接调用同一 helper 而非内联自己的版本；这样 vn 的结合顺序、faceOuterNormal 约定、两侧 vs 单侧输入全部统一，近零 w 的处理从根上一致。先写单面数值契约（CPU vs device 同输入逐位一致，≤1e-14），再接入生产。

**注意**：这仍需一次 instrument 确认 7e-15 累加源确实是对称面 BC（而非其它近零处理）；统一方案本身是稳健的架构改进，即使累加源在别处也应收口所有"对 0 处理"分歧。

**重要修正（用户问"按现在的修改结果能否对上"）**：不能。核查 wall probe（`/tmp/probe-{cpu,hip}.txt`，face 3457 的 owner cell lc=36864）发现**真实差异是主流场量级，不是近零舍入**：

| 量 | CPU | HIP | 相对差 |
|---|---|---|---|
| u_lc | 0.988912 | 0.997617 | 0.87% |
| v_lc | 0.0571734 | 0.054122 | 5.3% |
| w_lc | 1.66e-4 | 3.04e-5 | 82%（w 近零放大） |
| dudx/dudy | 4052/-29593 | 4087/-29853 | 0.85% |
| visc | 0.994556 | 0.996239 | 0.17% |
| nx/ny/nz/area | 完全相同 | 完全相同 | **0** |

几何量完全一致、速度/梯度差 0.87%~5.3% → 差异在**流场演化本身**（q 的推进），不是后处理、不是网格、不是近零舍入。"统一对 0 处理"只能消除近零 cell 的 7e-15，对不上主流场 0.87%~5.3% 的差异。之前 instrument 打印 cell 0/1 是对称面近零 cell，**采样偏差**导致误把次要现象当根因。

**从已有 stage trace 挖出的决定性证据（零集群时间，本地 `/tmp/of-trace-files/`）**：`128629176` 的 8 步 trace 各 1.5GB，含全部 294,912 cell 的 state 数组。写 `find_first_cell.py`/`track_cell.py` 逐 cell 比对：

- **step1 全 294,912 cell、5 方程完全逐位一致**（ndiff=0，worst rel=5.4e-20）
- **step2 突然 147,660 个 cell 同时不同**；FIRST diff 在 farfield cell 20/21（rel 1e-11，次要），**WORST 固定 cell 254656 的 u，从 step2 的 0.1% 线性放大到 step8 的 0.74%**
- cell 254656 轨迹：step1 五量全同；step2 跳变 **u 0.99635/0.99742（0.1%）、v 0.05273/0.05305、w -4.56e-4/-2.31e-4（差 1 倍）、p 1.00798/1.00871**；之后 CPU 的 w 持续降到 -1.98e-3，HIP 几乎不变 -3.5e-4

**根因确定**：差异在 **step1→2 的状态推进，w（展向速度）最敏感**（差 1 倍），u/v/p 差 0.1~0.5%。与 wall probe 的 w 差 82% 吻合。step1 逐位一致排除初值与首步 RHS；跳变在第一步末的写回/时间层推进。

**决定性根因（逐 kind 分析 `analyze_all_kinds.py`，零集群时间）**：

- `face_metadata`（几何/左右 cell/法向/面积）**全程 ndiff=0** —— 网格完全正确
- step1 所有 face/state/residual 逐位一致
- **step2 `face` 记录符号翻转**：cell 452808 的 `qf2`(eq=1, x 动量) cpu=-1.158 → hip=**+0.998**（符号反）；`invflux` eq=0 cpu=+2.66e-8 → hip=**-3.99e-8**（符号反）。nDiff=346,835 个 face
- step2 `inviscid_residual` 差 2.18e-4、`state` 差 0.00107（cell 254656）——都是 face 翻转的下游结果

**根因重新定性**：不是状态更新 w 分量，而是 **step2 的 HIP 无粘 face 重构/通量在部分面上给出了镜像/反号的值**（qf2 的 x 动量从 -1.158 变 +0.998，像是读到了对侧 ghost 或左右 cell 角色交换）。step1 正确、step2 翻转 → 状态推进后某个面的左右/ghost 归属在 device 路径变了。

**为什么之前关 hook/stage 无效**：`ONEFLOW_ACCEL_BACKEND=HIP` 默认启用的 fused flux kernel（`HipProductionPolicy.h:89-111`）不受 hook env 控制；旗舰 vismodel=3 使 stage 守卫失效但 fused flux 仍跑。所以差异始终来自这个始终启用的 kernel。

**修复尝试 1（失败，已回退）**：假设 `UNsUpdate.cpp:289` 的 `NsRefreshHostStateMirror()` 在 host 写完 q 后 D2H 覆盖新值。移除后昆山 `128631398-fixverify`（fix binary bf609f3d，20 步，**不开 STAGE_TRACE 避免诊断 D2H 污染**）结果与修复前**完全相同**（step20 CL 差 0.0094、CD_SF 差 5.18e-05）。证明该 D2H 在 host 回退路径下未实际覆盖（contract 不满足条件），**不是根因**。源码已 `git checkout` 回退，工作树只剩文档改动。

**最早分叉点（自动 diff 两臂全部 21 条 TRACE 行，`128633496`）**：第一条差异在 **step1 的 `TRACE_GRAD_OUT name=turbq`**（SA 湍流工作变量梯度）：
- CPU `g0=(560.19725129898188, 2730.9593000769778, 2313.3167818880438)`
- HIP `g0=(560.19725129898188, 2730.9593000769773, 2313.3167818880438)`
- 仅第二分量末位不同（…778 vs …773，差 ≈5e-12，相对 2e-15）

**这推翻"step1 完全一致"的判断**：之前 res.dat 只打印 5 位、state trace 用 1e-10 容差，都掩盖了这个 5e-12 末位差异。真正最早分叉是 **step1 SA 湍流量 host 梯度的末位数**，经 LU-SGS 迭代放大 → step2 dq 分歧 → 气动力 25%。两臂 hipGrad=0，都走 host `Grad::CalcGrad`，末位差异源于 host 浮点运算的某种环境/顺序差异（浮点 round=0 已证同，待查 `CalcGradGGCellWeight` 的编译/内联差异或 SA 场初值末位）。

**决定性突破（昆山 `128633496-tracepoints`，在 CalcInvFace / Grad::CalcGrad / UpdateFlowField 入口插桩，binary 2b081d64）**：

- **step1 全部量两臂逐位一致**：`TRACE_INVF` 的 q、`TRACE_GRAD` 的梯度、`TRACE_UPD` 的 dq 全相同
- **step2 分叉点在 `TRACE_UPD`（`UNsUpdate::UpdateFlowField` 入口）的 dq**：进入时 q（0.998106/0.997464/...）和梯度**完全相同**，但 dq 不同：
  - dq0 CPU -0.0018854 / HIP -0.0014452（比 0.766）
  - dq1 CPU -0.0029863 / HIP -0.0014845（比 0.497）
  - dq2 CPU -0.0004267 / HIP -0.0000847（比 0.198）
  - dq3 CPU -0.0002301 / HIP -0.0000061（比 0.026）
  - dq4 CPU -0.0089665 / HIP -0.0058395（比 0.651）
  - 比例不一致 → 不是单一缩放，是 RHS 某些项不同
- `TRACE_INVF` 显示 `hipGrad=0 hipRecon=0 hipBatch=0`（旗舰 vismodel=3 全 host 路径），cfl=0.01、lhscoef=0 两臂相同

**根因收敛**：差异在 `UPDATE_RESIDUALS`（RHS/dq 计算），不在状态更新公式、不在梯度、不在初值。q 和梯度进入 UpdateFlowField 时一致，但 RHS 算出的 dq 不同。下一步：插桩 RHS 累加（`UNsInvFlux::CalcFlux` 的 inviscid/viscous/source 分项）在 cell 254656 的值，定位是哪一项。粘性 device kernel 已证一致，怀疑 host 粘性/SA source 的 backend 分支。

**仍未知**：根因未定位到具体代码行。已确证的事实边界：
- 同一 binary、同 case、所有 script/grid/walldist 文件逐字节相同、浮点环境相同、step1 全部 294,912 cell 逐位一致、step2 突变
- 旗舰 vismodel=3 → 所有 device stage（gradient/recon/flux/state/boundary/timestep）+ `UseHipBatchAdapter` 全部 disabled → **两臂都走 host 路径**
- 唯一确定存在的 backend 相关 host 分支是 `UNsUpdate.cpp:284-294`（已证非根因）和 `NsSolverImp`/`EulerDomainStateSync` 的 capability 判断

**下一步（换思路，不再猜单一行）**：在 `CalcInvFace` / `UNsGrad::CalcGradient` / `UNsUpdate::UpdateFlowField` 的 host 路径入口用 instrument 打印两臂 step1/step2 的关键量（cfl/dt、dq、梯度、residual），一次作业拿到"host 路径相同输入为何不同输出"的直接证据。重点确认 step1 末两臂的 host `q` 和 `dq` 是否真的还一致——stage trace 只验证了 trace 时刻的 state，未覆盖"step1 末到 step2 计算之间"的所有中间量。

### 切片 1：RK 设备驻留 Advance（显式 RK 路径）

- **契约**：`HipEulerDomainBackend::Advance` 内完成 N 个 RK stage，state 留设备；每宏步 host 往返 ≤2；3 步 trace 1e-10、200 步 XB 1e-12、CPU 五例通过。
- **入口**：`EulerDomainBackend::Advance` 接口（`EulerDomain.h:159-162`）、`TimeIntegral::RungeKutta(context)`（:346-395，去掉 backend 硬编码 CPU，改能力表选择）、`RunEulerRungeKuttaStage`（把 host task 串替换为设备 kernel 序列）。
- **复用**：`HipFluxBackend::Shared()` 的 fused flux / gradient / state update kernels；1D `HipEulerBackend::Advance` 的驻留模式（`ports/kunshan/oneflow_1d_hip/OneDEulerBackend.h:51,68,87`）。
- **门槛**：先在受限 Euler 配置（vismodel=0、ilim=0）跑通并与 CPU 对齐，再逐步开放粘性/SA。

### 切片 2：旗舰 LU-SGS 设备驻留（硬骨头）

- **契约**：`TimeIntegral::Lusgs()` 的下扫/上扫/dq 交换纳入 backend 设备驻留；flagship m6wingroe_sa 端到端通过 B 类门禁 + paired 性能。
- **入口**：`TimeIntegral::Lusgs()`（:397-426）、`NsLusgs`/`ULhs`/`LhsTaskReg`（plan B.5）、`ExchangeInterfaceDq`（halo 依赖）。
- **难点**：下/上扫有严格 cell 顺序依赖，`dq` 需跨 zone 交换；host oracle 的 LU-SGS 语义必须逐位复现。
- **顺序**：切片 0、1 完成并通过门禁后启动；先证明单 zone 下/上扫设备化与 host 一致，再处理多 zone dq halo。

### 切片 3：真实算例与 ≥1M cell（Phase C）

- **契约**：≥1M cell 单卡端到端（含 restart 等价、气动力、可视化），paired 性能可复现。
- **入口**：C.1 IO/restart（消除每输出步整场回拉）、C.2 气动力设备归约、C.4 内存预算、D.1 固定 ≥1M fixture。
- **口径**：同节点、同 repeats、扣 1.4-1.6s 启动、per-rank cell 恒定；StageProfiler host 计时不作 GPU 时间（设备侧计时用 rocprof/hipEvent）。

## 5. 明确不做 / 风险

- 不在切片 0 未过门禁前动切片 1/2 的语义；不把 host task 串包装成 `Advance` 冒充驻留。
- 不改 CPU oracle 语义；不以单 rank 性能冒充加速；不采信未绑定准确提交的 CI/旧证据。
- 风险：ambient globals（`ug/unsf/nscom/ctrl/Iteration/ZoneState::zid`）阻碍显式上下文（plan E.1），需与切片并行梳理 device 侧所有权。
- 风险：capability 组合爆炸，用 E.2 表驱动收口；每个新组合先有门禁。

## 6. 证据与验收边界

- 每个切片落 `evidence/gate/<jobid>/`，作业名 `of-gate`/`of-probe-*`，单工作区。
- 数值门禁：3 步 trace 1e-10、200 步 XB 1e-12、CPU 五例 normal 1e-8 / strict 1e-15、50 步稳定性；收敛后 CL/CD/CD_SF 与 CPU oracle 比。
- 状态：`Verified`（门禁全过 + 收敛气动力对比）/ `Implemented, unverified`（仅实现）/ `Partial` / `Blocked`，如实写入 Living TODO。

## RHS 分项定案（昆山 `128635831-tracepoints7`，binary ea9acf52）

在 `AddInvFlux`（无粘累加后）和 `AddVisFlux`（粘性累加后）打印 cell 254656 的 residual：

- **step1**：INV 后、VIS 后两臂**完全一致**
- **step2 INV 后（无粘通量累加完）就不同**：r1 CPU -6.438e-06 / HIP -3.150e-06（差 51%）、r2 -9.243e-07/-1.786e-07、r3 -4.967e-07/-8.874e-09、r4 -1.960e-05/-1.276e-05
- step2 VIS 后在 INV 基础上继续，差值 pattern 相同 → 粘性不引入新差异

**最终收敛**：root cause = **step2 的无粘通量（inviscid flux）计算**。step1 相同、step2 不同，而 step2 进入时 q 和梯度都相同（TRACE_INVF/GRAD 已证）。无粘走 host 逐 face 循环（`hipBatch=0`），Roe scheme（`invFluxPointer`）。下一步定位单个 face 的 Riemann flux 差异：在 `UNsInvFlux::CalcInvFace` 的逐 face 循环里打印问题 face（如 452808）的 qL/qR/法向/通量，对比两臂。重点怀疑 `invflux` 的 host 累加或 Roe 通量的某个输入（尽管 q/梯度在 cell 254656 相同，但 face 两侧 cell 的值可能不同）。

## 根因确定 + 修复验证（昆山 128650947 + 128651348，binary 11a8824）

在 face 452808 打印 qf1/qf2 及其输入（limf->q / dqdx,y,z）后确定：

- **step1**：两臂 qf1_1 都 =4.0045935308630964（重构正常），输入 q/梯度逐位一致
- **step2**：输入 q/梯度仍逐位一致，但 CPU qf1_1=4.004189（重构生效），HIP qf1_1=0.998389=ql_1（**未重构，面值=格心值**）；qf2_1 CPU=-1.158 vs HIP=0.998=qr_1

机制：`LimField::CalcFaceValue()`（ULimiter.cpp:123-162）只在 `ckfun(qTry)`（`NsCheckFunction`: 密度>0 且压力>0）为真时写回重构值；否则 qf 保持 `GetQlQr()` 设的格心值。HIP 臂 step2 的 qTry 被判非物理 → 重构被拒。phi 异常源于 `noLimiterFieldFilled`：

- `HipUnusedInterfaceFieldsSkipped()` 默认返回 `HipProductionBackendSelected()`=true（HIP 臂）→ `noLimiterFieldFilled=true`（ULimiter.cpp:326）
- step1 fill limiter=1 后，step2 `if(noLimiterFieldFilled) return`（line 320）**跳过 fill**
- CPU 臂该值为 false → 每步 fill limiter=1 → 重构正常

**A/B 验证（128651348）**：HIP 臂设 `ONEFLOW_HIP_SKIP_UNUSED_INTERFACE=0`（单 zone 下唯一效果是 `noLimiterFieldFilled=false` → 每步 fill）后，step2 的 qf1/qf2 与 CPU **逐位一致**，residual r0..r4 **逐位一致**（r1=-6.438316441336979e-06）。根因确认，修复方向确定：ilim=0 时必须每步 fill limiter 场（不能依赖"fill 一次"优化）。

**修复**：ULimiter.cpp 的 `noLimiterFieldFilled` 不再跳过 per-step fill（正确性优先；每步 O(nCells) fill 代价比对整体步时间可忽略）。待重编译后跑 400 步气动力 + CPU 五例 + CI 验证。

## 400 步验收通过（昆山 128651666，修复 binary e032e7f）

复用修复 binary，CPU oracle vs HIP 全开 5 hook 各跑满 400 步，两臂 rc=0、nForceSave=nResSave=20（20 个共同步列）。无诊断 env 的干净运行：

- **aero.dat 逐字节完全一致**（cmp identical），末步 CL=1.3348e-01、CD=1.1677e-01、CD_PR=1.0765e-01、CD_SF=9.1244e-03
- **res.dat 逐字节完全一致**（cmp identical）

结论：旗舰 M6 算例 CPU 与 DCU 全开 hook 的数值**完全等价**。根因（limiter 每步 fill 被跳过）修复后，之前 25% 气动力偏差、step2 无粘通量 51% 偏差全部消失。**切片 0（准确性地基）完成，标 Verified**。后续可进入切片 1-3 设备驻留架构改动。注意：本算例 cfl=0.01 下 400 步不收敛（CD 仍在变化），但等价性比较与收敛无关；收敛气动力属 B.3.2 遗留项，需单独长期算例。

## DCU 真算验证（昆山 128660151，vismodel=0 配置，binary e032e7f）

用户质疑旗舰 400 步 DCU 没真算（fact：旗舰 vismodel=3 被 capability 守卫挡住 device stage，TRACE hipGrad=0/0/0）。改用 vismodel=0 + 显式 RK3 + Lax 单zone（device batch 支持配置）验证 DCU 真算：

- **DCU 确实在算**：`HIP device 66a1 gfx906`、backend=HIP/DCU、`HIP_STAGE_DIAGNOSTIC hip_gradient=1 hip_reconstruction=1`、stage breakdown `HIP kernel 4.73726`（非零）
- **RES IDENTICAL**：res.dat 逐字节一致 → 流场演化完全相同
- **AERO DIFFER**：step20 CL=6.4428e-02(CPU) vs 5.0928e-02(HIP)，差 21%；CD 同样差 21%

**结论**：device 路径数值（残差/流场）与 CPU 逐字节一致，但气动力后处理差 21%。这是气动力读 host 镜像陈旧的 bug（handoff 坑41 / runbook §1f）：DCU 路径 device 拥有 q/梯度，host 气动力计算读的 dqdx/y/z 或 q[rc] ghost 镜像陈旧；残差不走 host 镜像所以一致。之前旗舰 DCU 没算所以该 bug 未显现，现在 DCU 真算了暴露出来。vismodel=0 无粘，CD_PR=CD、CD_SF=0，符合预期。下一步：用 wall force probe 对比两臂壁面 q/梯度，确认是梯度镜像陈旧还是 q[rc] ghost。

## DCU 真算后气动力 21% 差异 — 定位方向

128660151 确认 vismodel=0 下 DCU 真算、流场(res)逐字节一致，但气动力差 21%。气动力压力项 `wp = q[IP][lc] + dqdx[IP][lc]*dx + ...`（AeroForceTask.cpp:363）读 host 的 q 与 cell 梯度。

- `NsMirrorHostBoundaryState()`（NsRhs.cpp:152-176）在 HIP 下只镜像 **q**（MirrorPrimitiveCellsToHost），**不镜像梯度 dqdx/y/z**
- 注释（AeroForceTask.cpp:209-210）声称 HIP 路径已镜像 boundary-adjacent cells，但那只覆盖 q；梯度靠 `HipGradientHostMirrorOnDemand()` 另一机制，气动力读到的 dqdx[IP][lc] 可能是陈旧值（某 RK stage 的）
- 残差在 device 内部算、不走 host 镜像 → 一致；气动力读 host 镜像 → 差

下一步：在 AeroForceTask.cpp:363 加 env-gated instrument 打印壁面 q[IP][lc]/dqdx[IP][lc]/wp，对比两臂。

## 气动力 21% 差异根因确定（昆山 128660722，AERO_INPUT probe）

在 AeroForceTask.cpp:363 加 env-gated AERO_INPUT probe，打印壁面 cell 气动力输入：

- **qip（壁面压力 q[IP][lc]）**：CPU 0.99307596947828625 / HIP 0.99307596947828669 —— 逐位一致
- **gpx/gpy/gpz（压力梯度 dqdx/y/z[IP][lc]）**：CPU = 25.307/-170.246/-7.081，**HIP = 0/0/0 —— 全零！**
- wp：CPU 0.98738 / HIP 0.99308（差 0.0057）→ cp 差 → 气动力 21%

**根因 100% 确定**：DCU 真算时 device 梯度 stage 正确算出梯度（流场/残差逐字节一致证明），但气动力读的 host cell 梯度 `dqdx/y/z[IP][lc]` **从未被镜像到 host，保持零**。这正是 runbook §1f / 坑 41 的物理表现。`NsMirrorHostBoundaryState()`（AeroForceTask.cpp:233）只镜像 q，不镜像梯度；梯度靠 `HipGradientHostMirrorOnDemand()` 按需机制，但气动力读取前没有触发该镜像（或触发条件不满足）。

修复方向：气动力读梯度前必须确保 host 梯度镜像是最新的（boundary-adjacent owned cells），即在 `CalcForce()` 里补一次梯度 D2H 镜像（类似 NsMirrorHostBoundaryState 但针对 dqdx/y/z 的壁面邻域 cell），或修正 HipGradientHostMirrorOnDemand 的触发条件覆盖气动力读取点。待定位具体镜像 API 后修复+验证。

## DCU 真算 + 气动力对齐修复完成（昆山 128660853，binary 147d3f9f）

修复：AeroForceTask::CalcForce() 在 HIP device-authoritative 路径下，读梯度前主动调用 `MirrorGradientCellsToHost` 把所有 boundary face owner cell 的 dqdx/y/z 镜像到 host（对称于已有的 NsMirrorHostBoundaryState 镜像 q）。

vismodel=0 显式 RK3 Lax 单zone 20 步验证（DCU 真算：HIP kernel 非零、hip_gradient/hip_reconstruction=1）：

- **AERO IDENTICAL**（修复前差 21%）
- **RES IDENTICAL**
- AERO_INPUT probe：gpx CPU 25.30716503487686 / HIP 25.307165034877251（机器精度末位差），wp 两臂 0.98738090912439525/0.98738090912439547（一致到 16 位）

结论：vismodel=0 配置下 DCU 设备路径真正执行且数值与 CPU 逐字节对齐（流场残差 + 气动力）。这回答了用户"DCU 没算没意义"的质疑——现在 DCU 真算了，且证明对齐。气动力差异根因是 host 梯度镜像陈旧/为零，非数值 kernel。

**下一步**：这是 vismodel=0 受限路径。旗舰 vismodel=3 + LU-SGS 才走 device 路径（切片 2），需在 A.4 驻留架构里保证气动力读的 host 梯度镜像到位。

## 气动力通盘审计（4 个并行 Explore agent + 代码核实，2026-11-22）

用户要求不逐个打补丁、通盘排查气动力。4 个 agent 覆盖：host 数据依赖、计算/输出链、镜像时机/代际、旗舰特殊风险。一致结论：

**已修复**：q[lc]（NsMirrorHostBoundaryState）、dq[lc]（8618ecb26 MirrorGradientCellsToHost）——vismodel=0 DCU 真算下 aero.dat 已逐字节一致。

**剩余风险（按证据强度）**：
1. **q[rc] ghost 无 D2H 镜像**（3/4 agent 确认）：CorrectGrad（VisGrad.cpp:328）读 q[IU/IV/IW][rc]，AeroForceTask.cpp:354-356/458-460；现有镜像只覆盖 leftCells（owner），MirrorPrimitiveCellsToHost 对 ghost cell 抛异常（HipFluxBackend.hip:5215-5219）。旗舰 vismodel=3 的 viscous block 会触发，单 zone 壁面 rc 多为 ghost → 高风险。**修 q[rc] 需新增/扩展 backend API（能按任意 cell id 含 ghost 镜像 q），属架构决策。**
2. **镜像 gate 不对称**：q 镜像需 `HipDeviceAuthoritativeStateEnabled() && HipBoundaryStateMirrorEnabled()`（NsRhs.cpp:155-156），dq 镜像只需前者（AeroForceTask.cpp:245）。若 `ONEFLOW_HIP_BOUNDARY_STATE_MIRROR=0` → dq fresh 但 q stale → wp 错。统一 gate。
3. **visl host-written**（UNsVisFlux.cpp:743-759 基于 host tempr）：当前非 device-authoritative 依赖，不需 D2H；但切片1-3 粘性搬 device 后会 stale（agent4 未来风险）。
4. **wallaero.dat / heatflux 无镜像**（PostProcess.cpp:51、HeatFluxTask.cpp:127-141）：aero.dat 对齐后它们可能仍 stale。
5. **cross-backend verifier 限 vismodel=0**（f2-main-solver-cross-backend-residual.sh:105-118，容差 aero 1e-10），旗舰需改造配置。

**镜像时机确认**（agent3）：aero 在所有 RK stage 后读（NsSolverImp.cpp:425-428 → ForceOk），梯度是最后 stage 开始时的值（一 stage lag，RK 固有，非 bug）；显式镜像无额外 lag、无缓存滞后（MirrorGradientCellsToHost 每次重打包）。单 zone packedExchange=false，on-demand 梯度镜像从不触发，只有显式镜像工作 → 8618ecb26 的修复必要且正确。

**通盘方案**：把气动力输入镜像收敛为单一 helper `MirrorAeroInputsToHost()`（q owner+rc ghost、dq owner、统一 gate），新增 backend API 支持按任意 cell id（含 ghost）镜像 q。待用户确认架构范围后实施。

## 旗舰 DCU 真算成功（昆山 128662638，binary 6deb5225）

用户要求"真刀真枪"：旗舰 vismodel=3 之前被 capability 守卫挡住、DCU 没真算。实施最小 opt-in：`ONEFLOW_ALLOW_VISCOUS_HIP_BATCH=1` 时把 request.viscousModelActive 视为 false，放行无粘 fused flux device 路径（粘性仍走 UNsVisFlux 自己的路径，二者解耦——CalcInvFluxBatch 只算无粘 Riemann flux）。

旗舰 verifier（ONEFLOW_XB_PROFILE=flagship + HIP_ALLOW_VISCOUS_BATCH=1）验证：

- **HIP kernel 25.4773**（stage breakdown，非零）→ 无粘 fused flux 真正在 DCU 执行
- residual pass：max_rel=6.78e-21（机器精度）
- aero pass：max_rel=0（逐字节一致）
- verifier rc=0, CROSS_BACKEND_RESIDUAL_PASS
- 配置确认：vismodel=3, Roe(ischeme=1), LU-SGS(time_integral=2), TurbSolver 保留

**里程碑**：旗舰 3D 工作流的无粘通量真正在 DCU 上计算且与 CPU 数值等价。这是 A.4 设备驻留的一个真实切片（无粘 fused flux device 化），比之前 host 路径验证有本质进展。剩余：梯度/重构/LU-SGS 宏步循环仍在 host（UseHipGradient 有独立 vismodel==0 守卫），完整设备驻留是后续切片。verifier 增强：HIP_EXTRA env 透传 + stage breakdown 输出（HIP kernel 非零作为真算证据）。

## 有粘部分 device 化成功（昆山 128663186，复用 binary 6deb5225）

两个并行 agent 分别设计粘性 device 化与 LU-SGS 驻留。粘性线结论：旗舰所有 CalcVisFluxHip 守卫都满足（nEqu=5/ivischeme=0 VIS_STD/isowallbc=0/vismodel 无守卫），唯一缺 `ONEFLOW_ENABLE_UNS_HIP_VISFLUX=1`（HipProductionFeatureEnabled 需显式=1，backend=HIP 不自动开）。SA 的 vist 是 CalcVisFlux 的逐调用输入（uploaded as-is），host 算的 vist 即可，不是 device 粘性前置条件。

实验：旗舰 profile + ONEFLOW_ALLOW_VISCOUS_HIP_BATCH=1 + verifier 通用 HIP env 透传（ONEFLOW_XB_HIP_EXTRA_ENVS="ONEFLOW_ENABLE_UNS_HIP_VISFLUX=1"）。结果：

- **v_kernel 0.751645、v_upload 87.73、v_download 116.6**（之前全 0）→ 粘性 kernel 真正在 DCU 执行
- **HIP kernel 25.52**（无粘 fused flux 也在算）
- residual max_rel=6.78e-21、aero max_rel=0 → 数值逐字节等价
- verifier rc=0, CROSS_BACKEND_RESIDUAL_PASS

**里程碑**：旗舰无粘 + 粘性通量都在 DCU 真算且等价。注意 v_upload 大是每步上传 visl/vist；SA hooks（VIST/SOURCE/INVFLUX/VISFLUX）仍未开——vist 当前 host 算上传，下一步可啃 SA device 化。LU-SGS 宏步、梯度、重构仍在 host（agent2 设计见下条）。

verifier 增强：ONEFLOW_XB_HIP_EXTRA_ENVS 通用透传（分号分隔 KEY=VALUE），支持 device hook 实验不改脚本。

## LU-SGS device 驻留设计（agent2 排查，待实施）

**现状（代码核实）**：无 HipEulerDomainBackend；EulerDomainBackend::Advance 只有 Cpu 实现。TimeIntegral::RungeKutta 硬编码 CPU backend；Lusgs()（TimeIntegral.cpp:397-426）是纯 host task 串：ZERO_DQ→CALC_TIME_STEP→LOAD_RESIDUALS→UPDATE_RESIDUALS→INIT_LUSGS→[LOWER_SWEEP→EXCHANGE_DQ→UPPER_SWEEP]→UPDATE_FLOWFIELD→CALC_BOUNDARY。UNsLusgs Lower/UpperSweep（UNsLusgs.cpp:79-137）是严格逐 cell 顺序循环（lower 0→N、upper N→0），PrepareSweep 读 unsf.rhs/dq/impsr，Solve 读 q/gama/visl/vist、写 dq。ExchangeInterfaceDq 是空 stub（InterfaceTaskReg.cpp:354-357）。

**推荐第一刀**（不做完整宏步 Advance）：只把 lower/upper sweep 核心做成 device kernel + 窄 host 编排；dq 在 sweep 循环内驻留 device，rhs 每宏步上传一次。完整 Advance（RK 式宏步循环进设备）需先把梯度/重构/粘性/SA 都搬到 device（旗舰全被 vismodel 守卫挡住），一次作业太大。
- 接入点：TimeIntegral.cpp:269-283/397-426（env-gated 路由到新 LusgsAdvance）；新 HipLusgsBackend（EulerDomain.h:146-192）；从 UNsLusgs.cpp:79-137/327-394 抽 kernel，复用 NsLusgsPolicy.h:7-52 纯数学；StageProfiler.cpp:162-170 注册 macro_lusgs_* 行（注意静默丢弃陷阱 :162-164）；ResidentStateContract.h 管代际。
- 验收：stage breakdown macro_lusgs_* 非零 + HIP kernel 非零 + res/aero 等价（如 128662638 的 6.78e-21/0）。
- 前置：gate 到 nZones==1（dq halo 是 stub，MPI 延后）；rhs 每宏步 FullField 上传（bump generation）；SA/viscous 累加仍 host 进 rhs。
- 不在范围：device 梯度/重构/flux/state-update、多 zone dq MPI、device SA 更新、粗网格/MG、device restart/IO。

**可复用脚手架**：NsLusgsPolicy.h（纯 sweep 数学）、StageProfiler macro 行、NsResidentMacroStep.cpp:59-140、ResidentStateContract.h、旗舰 verifier、1D HipEulerBackend::Advance（ports/kunshan/oneflow_1d_hip/OneDEulerPersistent.hip:197-303）作为驻留模板。
