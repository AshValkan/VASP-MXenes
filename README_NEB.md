# CI-NEB 离子迁移势垒计算：方法与参数

> 适用于吸附离子在 MXene 表面相邻位点间的最小能量路径（MEP）与迁移势垒计算。
> 工具脚本见仓库根目录（§6），体系与计算参数以 `alloy_config.yaml` 为唯一来源。

## 1. 方法原理

CI-NEB（Climbing-Image Nudged Elastic Band，Henkelman & Jónsson 2000）在初、末态之间构造一串由弹簧连接的"镜像"（images），迭代优化使整条带收敛到**最小能量路径（MEP）**；攀爬镜像（climbing image）取消最高能量镜像的弹簧力、使其沿梯度自由爬升至**一阶鞍点**，其能量与初态能量之差即为迁移势垒：

$$E_a = E_{\text{saddle}} - E_{\text{initial}}$$

对吸附离子表面迁移场景，CI-NEB 的输入只有两样：**充分弛豫的初态与末态**（VASP Wiki 明确强调端点质量决定路径质量）。

## 2. 计算流程（六步）

1. **初态弛豫**：吸附构型弛豫至 EDIFFG=−0.02，CONTCAR 即初态
2. **末态构建与弛豫**：吸附离子平移至相邻 O-top 位点（本项目 3.000 Å）→ 独立弛豫（收敛判据与初态一致）
3. **镜像构建**：两端点间线性插值（简单位点跳跃适用；复杂路径用 IDPP）
4. **CI-NEB 运行**：3 个中间镜像 + 2 个固定端点
5. **收敛判据**：所有镜像最大原子力 < |EDIFFG|
6. **势垒提取**：$E_a = \max_k(E_k) - E_{00}$；可选鞍点虚频验证

## 3. 参数表（检索推荐值 → 本项目取值）

| 参数 | 文献/VASP Wiki 推荐 | 本项目取值 | 依据 |
|---|---|---|---|
| IMAGES（中间镜像数） | 越少越好，1 起步逐级加密；3–7 常用 | **3** | 对称位点跳跃、路径简单；np = 3 × 5 恰好按镜像分组 |
| SPRING | 负值启用 NEB，−5 推荐，勿用大值 | **−5** | Wiki 原文推荐 |
| LCLIMB | .TRUE.（CI-NEB 核心） | **.TRUE.** | 直接获得鞍点，无需两阶段 |
| IBRION | 1（RMM-DIIS）或 3（quick-min）；CG 在非线性约束下易发散 | **1 + NFREE=2** | Wiki 对 NEB 起步的低维参数推荐组合 |
| POTIM | IBRION=1 时 0.1–0.3 | **0.1** | 保守起步 |
| EDIFF | 力需可信，1e-5 足够（Henkelman：1e-4 即安全） | **1e-5** | 与端点弛豫口径一致，保证能量可比 |
| EDIFFG | −0.03~−0.05 定性可用；发表级 −0.02 | **−0.03** | 平衡精度与算力；必要时事后精修至 −0.02 |
| NSW | 数百步 | **200** | — |
| ISYM | 关闭 | **0** | 防止对称化强行扭曲路径 |
| KPOINTS/泛函/+U/SIGMA | 与端点弛豫**完全一致** | 3×3×1、GGA+U、SIGMA=0.1 | 能量口径可比性 |
| LCHARG/LWAVE | 关闭省盘 | **.FALSE.** | 镜像多、WAVECAR 巨大 |
| 电荷混合 | — | **AMIX=0.1** | 电荷混合稳定性 |

### 关键参数出处（VASP Wiki 原文）

- SPRING：*"The nudged elastic band method is applied when SPRING is set to a negative value e.g. SPRING=-5. This is also the recommended setting."*
- IBRION：*"we recommended to use the RMM-DIIS algorithm (IBRION=1) or the quick-min algorithm (IBRION=3)"*（非线性约束下 CG 常发散）
- NFREE：IBRION=1 下限制 DIIS 迭代历史秩（近似 Hessian 秩 ≤ NFREE）；NEB 起步阶段 *"a very low dimensionality parameter (IBRION=1, NFREE=2) should be applied in the first few steps"*

## 4. 并行方案

- VASP 原生 NEB：**单个 mpirun 内部按镜像分组**——总 rank 数 = IMAGES × 每镜像 rank 数（如 3 × 5 = 15）
- CPU 绑定：`taskset -c <CPU 块>` + **必须 `--bind-to none`**（PRRTE 默认把 rank 绑到机器前几个核，会覆盖 taskset 掩码）
- 端点目录 00 / N+1 固定，不参与运行，无单独计算量

## 5. 风险提示

| 风险 | 触发症状 | 处置 |
|---|---|---|
| 端点未充分弛豫 | 路径端点能量不齐、势垒失真 | 末态弛豫收敛判据与初态一致（EDIFFG=−0.02） |
| 线性插值路径差 | 镜像原子间距过近、收敛到错误通道 | 短路径（~3 Å）且无避让动作时风险低；异常时改 IDPP（需 pymatgen） |
| SCF 震荡 | OSZICAR 步数不动、dE 平台震荡无下降趋势 | AMIX=0.1；仍震荡则 CONTCAR→POSCAR 续跑 |
| NEB 收敛慢/力震荡 | Fmax 长期 >0.1 | 切 IBRION=3（SMASS=2）预收敛后再回 IBRION=1 |
| ZBRENT 类崩溃 | IBRION=1 报 ZBRENT 错误 | CONTCAR→POSCAR 续跑（IBRION=1 近极小区稳健） |
| 多通道竞争 | 势垒异常低 | 与空位/hollow 通道对比 |
| 虚频验证缺失 | 鞍点存疑 | 表面迁移文献常省略；发表级需对鞍点做频率分析（恰 1 虚频） |

## 6. 工具脚本（仓库根目录）

| 脚本 | 作用 |
|---|---|
| `make_final_state.py` | 从 CONTCAR 构建终态 POSCAR（PBC 感知邻位搜索，断言 hop 距离合理区间） |
| `launch_final_state.sh` | 终态目录搭建 + INCAR 修正（IBRION=1、AMIX=0.1）+ 启动 |
| `verify_final.py` | 核对终态 POSCAR 的逐离子位移（防 CRLF 显示伪影误判） |
| `check_neb_prereq.sh` | NEB 前置核查（端点 CONTCAR、INCAR 关键标签、二进制 NEB 支持字符串） |
| `make_neb_images.py` | 线性插值生成 00–04 镜像 + 根目录 INCAR/KPOINTS/POTCAR（分数坐标 nearest-image 对齐） |
| `run_neb.sh` | NEB 启动器（参数化 neb_root / np / cpu_block，taskset + --bind-to none） |

## 7. 参考

- VASP Wiki: Nudged elastic bands（SPRING=−5、IBRION=1/3、端点优化、镜像数最小化）
- VASP Wiki: SPRING / NFREE / LCLIMB 标签页
- Henkelman & Jónsson, JCP 113, 9901 (2000)（CI-NEB 原始文献）
- Henkelman et al., JCP 113, 9978 (2000)
- VASP Forum #5546：EDIFFG=−0.01~−0.05 对势垒影响 <0.05 eV（graeme/Henkelman 回复）
- Mills & Jónsson, Surf. Sci. 324, 305 (1995)
