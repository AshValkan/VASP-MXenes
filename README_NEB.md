# CI-NEB 离子迁移势垒计算：方法学总结与立项

> 日期：2026-10-11 ｜ 状态：**方法论已复核通过（VASP Wiki SPRING/NFREE 页逐条核对，见 §3 复核记录），先导计算运行中**
> 先导案例：Li 在 Ti₃C₂O₂（y=0）表面相邻 O-top 位点间的表面迁移势垒

---

## 1. 方法原理

CI-NEB（Climbing-Image Nudged Elastic Band，Henkelman & Jónsson 2000）在初、末态之间构造一串由弹簧连接的"镜像"（images），迭代优化使整条带收敛到**最小能量路径（MEP）**；攀爬镜像（climbing image）取消最高能量镜像的弹簧力、使其沿梯度自由爬升至**一阶鞍点**，其能量与初态能量之差即为迁移势垒：

$$E_a = E_{\text{saddle}} - E_{\text{initial}}$$

对吸附离子表面迁移这一场景，CI-NEB 的输入只有两样：**充分弛豫的初态与末态**（VASP Wiki 明确强调端点质量决定路径质量）。

## 2. 标准流程（六步）

1. **初态弛豫**：Q3 已完成（ads_Li 4/4 收敛，CONTCAR 即初态）
2. **末态构建与弛豫**：Li 平移至相邻 O-top 位点（本项目 3.000 Å）→ 独立弛豫
3. **镜像构建**：两端点间线性插值（简单位点跳跃适用；复杂路径用 IDPP）
4. **CI-NEB 运行**：3 个中间镜像 + 2 个固定端点
5. **收敛判据**：所有镜像最大原子力 < |EDIFFG|
6. **势垒提取**：$E_a = \max_k(E_k) - E_{00}$；可选鞍点虚频验证

## 3. 参数表（检索推荐值 → 本项目取值）

| 参数 | 文献/VASP Wiki 推荐 | 本项目取值 | 依据 |
|---|---|---|---|
| IMAGES（中间镜像数） | 越少越好，1 起步逐级加密；3–7 常用 | **3** | 对称位点跳跃、路径简单；3 镜像 × np5 恰好占满 3 个核块 |
| SPRING | 负值，−5 推荐，勿用大值 | **−5** | Wiki 原文推荐 |
| LCLIMB | .TRUE.（CI-NEB 灵魂） | **.TRUE.** | 直接获得鞍点，无需两阶段 |
| IBRION | 1（RMM-DIIS）或 3（quick-min）；CG 易在非线性约束下发散 | **1 + NFREE=2** | Wiki：低维参数起步稳健；且本体系 IBRION=1 有 ZBRENT 教训后的成功记录 |
| POTIM | IBRION=1 时 0.1–0.3 | **0.1** | 保守起步 |
| EDIFF | 力需可信，1e-5 足够（Henkelman：1e-4 即安全） | **1e-5** | 与端点弛豫口径一致，保证能量可比 |
| EDIFFG | −0.03~−0.05 定性可用；发表级 −0.02 | **−0.03** | 平衡机时与精度；必要时事后精修至 −0.02 |
| NSW | 数百步 | **200** | — |
| ISYM | 关闭（ISYM=0） | **0** | 防止对称化强行扭曲路径 |
| KPOINTS/泛函/+U/SIGMA | 与端点弛豫**完全一致** | 3×3×1、GGA+U、SIGMA=0.1 | 能量口径可比性 |
| LCHARG/LWAVE | 关闭省盘 | **.FALSE.** | 镜像多、WAVECAR 巨大 |
| 电荷混合 | — | **AMIX=0.1 预置** | 本体系 SCF 震荡三次教训 |

### 复核记录（2026-10-11，VASP Wiki 原文核对）

- **SPRING**（[SPRING 页](https://www.vasp.at/wiki/index.php/SPRING)）：*"The nudged elastic band method is applied when SPRING is set to a negative value e.g. SPRING=-5. This is also the recommended setting."*；同页明确 *"we recommended to use the RMM-DIIS algorithm (IBRION=1) or the quick-min algorithm (IBRION=3)"*（CG 在非线性约束下常发散）→ 与本项目 SPRING=−5、IBRION=1 取值一致。
- **NFREE**（[NFREE 页](https://www.vasp.at/wiki/index.php/NFREE)）：IBRION=1 下 NFREE 限制 DIIS 迭代历史秩（近似 Hessian 秩 ≤ NFREE）；SPRING 页对 NEB 起步阶段 *"a very low dimensionality parameter (IBRION=1, NFREE=2) should be applied in the first few steps"* → **IBRION=1 + NFREE=2 为 Wiki 逐字推荐组合**，本项目取值正确。

## 4. 并行方案（本机 20 核）

- VASP 原生 NEB：**单个 mpirun 内部按镜像分组**——`mpirun -np 15` + `IMAGES=3` → 3 镜像 × 5 ranks
- 核块绑定：`taskset -c 0-29`（15 ranks 均布于 3 个 10-CPU 块），**必须 `--bind-to none`**（PRRTE 覆盖 taskset 事故教训）
- 端点不参与运行（00/04 固定），无单独计算量

## 5. 风险提示

| 风险 | 触发症状 | 处置 |
|---|---|---|
| 端点未充分弛豫 | 路径端点能量不齐、势垒失真 | 末态弛豫收敛判据与初态一致（EDIFFG=−0.02） |
| 线性插值路径差 | 镜像原子间距过近、收敛到错误通道 | 先导路径短（3 Å）且无避让动作，风险低；异常时改 IDPP（需装 pymatgen） |
| SCF 震荡 | OSZICAR 步数不动、dE 平台 | AMIX=0.1 已预置；仍震荡则 CONTCAR 续跑 |
| NEB 收敛慢/力震荡 | Fmax 长期 >0.1 | 切 IBRION=3（SMASS=2）预收敛后再回 IBRION=1 |
| ZBRENT 类崩溃 | IBRION=1 报错 | CONTCAR→POSCAR 续跑（成熟流程） |
| 多通道竞争 | 势垒异常低 | 与空位/hollow 通道对比（后继立项） |
| 虚频验证缺失 | 鞍点存疑 | 表面迁移文献常省略；发表级需对鞍点做频率分析（恰 1 虚频） |

## 6. 立项范围与机时预算

**先导（已启动）**：Li @ y=0（Ti₃C₂O₂），无 V 磁性复杂度，验证全流程。

| 阶段 | 内容 | 状态 | 预计机时 |
|---|---|---|---|
| 终态弛豫 | ads_Li_fB（Li 移位 3.000 Å，IBRION=1，AMIX=0.1） | ✅ 14:58 已启动 | 0.5–2 天 |
| CI-NEB | neb_Li 3 镜像 × np5，脚本已预部署 | 待终态收敛触发 | 2–4 天 |

**扩展矩阵（先导验证后逐批）**：{Li, Na, K, Zn} × {y=0, 0.5, 1.0, 1.5} 共 16 条路径；科学主线：ε_d 上移 → 吸附增强（Q3/Q4 已证）vs 迁移势垒 E_a 变化 → "容量-倍率"权衡叙事（仓库 neb.py 注释中预设的故事线）。

## 7. 产出物

- 迁移势垒 E_a(Li, y=0) + 能量路径曲线（E vs 反应坐标）
- 鞍点构型（VESTA 可视化过渡态）
- 扩展后：E_a vs y 曲线（与 E_ads vs y 并列，论文第 4 章）

## 8. 服务器端已部署脚本（已 commit 至仓库根目录）

| 脚本 | 作用 |
|---|---|
| `make_final_state.py` | 从 CONTCAR 构建终态 POSCAR（PBC 感知邻位搜索） |
| `launch_final_state.sh` | 终态目录搭建 + INCAR 修正 + 启动（已执行） |
| `verify_final.py` | 核对终态 POSCAR 的逐 Li 位移（防止 CRLF 显示伪影误判） |
| `check_neb_prereq.sh` | NEB 前置核查（端点 CONTCAR、INCAR 关键标签、二进制 NEB 字符串） |
| `make_neb_images.py` | 线性插值生成 00–04 镜像 + 根目录 INCAR/KPOINTS/POTCAR（10-11 修正：源 INCAR 剔除 EDIFF/LCHARG/LWAVE 再追加，消除重复标签） |
| `run_neb.sh` | 15-rank NEB 启动器（参数化 neb_root/np/cpu_block，taskset + --bind-to none） |

## 9. 参考

- VASP Wiki: Nudged elastic bands（SPRING=−5、IBRION=1/3、端点优化、镜像数最小化）
- Henkelman & Jónsson, JCP 113, 9901 (2000)（CI-NEB 原始文献）
- Henkelman et al., JCP 113, 9978 (2000)
- VASP Forum #5546：EDIFFG=−0.01~−0.05 对势垒影响 <0.05 eV（graeme/Henkelman 回复）
- Mills & Jónsson, Surf. Sci. 324, 305 (1995)
