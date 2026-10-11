# VASP-MXenes — 合金化 MXene 电子结构计算工具集

>
> 本仓库是 [diegonti/VASP-MXenes](https://github.com/diegonti/VASP-MXenes) 的**分支**
> （branch `alloy-mxene-tiv`），在其"单金属 MXene → 批量生成 / 优化 / 计算 / 分析"的骨架上，
> 扩展为可处理 **Ti<sub>3−y</sub>V<sub>y</sub>C<sub>2</sub>O<sub>2</sub>** 双金属固溶体系，
> 并把 **d 带中心**作为核心描述符。
>
> 上游工具集作者：**Diego Ontiveros**（原文文档见 [`README_upstream.md`](README_upstream.md)）。
> 分支扩展：alloyed-MXene 工作流。

---

## 一、项目简介与核心功能

### 要回答的三个科学问题

| | 问题 | 本仓库提供的手段 |
|---|---|---|
| Q1 | V 原子在 M 位**占据外层还是内层**？ | `alloy.py` 生成 `outer`/`inner`/`sqs`/… 显式构型；`energetics.site_preference()` 给出能量差 |
| Q2 | V 的引入是否使 **d 带中心向费米能级移动**？ | `dband.py` 计算分元素、分原子层的 ε_d，并给出线性拟合的 **R²** |
| Q3 | 电子结构变化如何**传递**到容量、倍率与循环？ | `adsorption.py`（吸附能）+ `neb.py`（扩散势垒）+ `analyze_alloy.py`（与实测数据交叉作图） |

### 核心功能

**结构建模与输入生成**
- 二元 M 位固溶 MXene，任意组分 y 与超胞尺寸（自动校验公度性）
- V 占位模式：`outer` / `inner` / `random` / `random-outer` / `random-inner` / `sqs` / `cluster`
- 自动生成 POSCAR / POTCAR / KPOINTS / INCAR，含 **DFT+U**（U(Ti)=U(V)=2.5 eV）、**DFT-D3**、slab 偶极校正
- k 网格按超胞自动缩减，避免约 4 倍算力浪费

**后处理与分析**
- d 带中心 ε_d、带宽、偏度；支持**按原子层**（outer_bottom / inner / outer_top）分辨
- H / Li / Na / K 吸附建模，可按"表面 O 下方金属"选位，默认双面对称吸附
- CI-NEB 建模（线性 / IDPP）与能垒提取
- 形成能、混合能（凸包）、构型排序
- 一条命令输出论文图表：ε_d vs y、ε_d + E_ads + 实测 C_g 三联图

**工程化**
- `alloy_config.yaml` 是体系与计算参数的**唯一来源**
- `tests/selftest.py`：**68 项自检，不需要 VASP**
- `tests/make_fake_results.py`：伪造完整结果树，在真 VASP 跑通前即可预演分析流程

---

## 二、技术栈与依赖说明

| 类别 | 名称 | 版本 | 必需 | 说明 |
|---|---|:---:|:---:|---|
| 语言 | Python | ≥ 3.9 | ✅ | 仅标准库即可建模 |
| 数值 | numpy | 2.5.x | ✅ | 兼容 numpy 2.x |
| 绘图 | matplotlib | 3.11.x | ✅ | Agg 后端，无需显示器 |
| 配置 | PyYAML | 6.x | ✅ | 读取 `alloy_config.yaml` |
| 结构 | pymatgen | — | ⬜ | 仅 NEB 的 IDPP 初始猜测需要；缺失时自动降级为线性插值 |
| 求解器 | **VASP** | 5.4.1+ / 6.x | ✅ | PAW-PBE；NEB 用 VASP 原生 LCLIMB（无需 VTST） |
| 辅助 | VASPKIT / Multiwfn / VESTA | — | ⬜ | 后处理与可视化 |

> **VASP 与 POTCAR 均需授权许可，不随本仓库分发。**

---

## 三、目录结构与模块职责

```
VASP-MXenes/
├── .gitignore
├── README.md                  本文件（项目主文档）
├── README_upstream.md         上游原作者文档（原文保留）
├── README_ALLOY.md            分支改动说明
├── alloy_config.yaml          体系定义与计算参数（唯一来源）
│
├────── 合金化扩展（本分支新增）──────
├── run_alloy.py               配置驱动的统一入口（生成 / 分析）
├── alloy.py                   AlloyMX 类：合金建模、占位控制、VASP 输入生成
├── dband.py                   PDOS 解析、d 带中心、能带矩、组分趋势作图
├── adsorption.py              吸附位点建模与吸附能
├── neb.py                     CI-NEB 建模与能垒提取
├── energetics.py              形成能 / 混合能 / 占位偏好能 / 构型排序
├── analyze_alloy.py           一键产出论文图表与数据表
│
├────── 上游核心（diegonti）──────
├── VASP.py                    MX() 类与输入文件树生成
├── VASPread.py                OUTCAR 解析（能量、优化状态）
├── structure.py               POSCAR / CONTCAR 编辑（加真空、加终止、平移）
├── searcher.py                结果树检索与批量文件搬运
├── optimizer.py               自动优化机器人（isif2 ↔ isif7 迭代）
├── opt.py                     批量派发 optimizer.py
├── calculate.py               批量提交 DOS / BS / WF 计算
├── analyzer.py                批量 DOSCAR 分析与作图
├── database.py                汇总结构与电子性质数据库
├── DOS.py                     DOSCAR 读取、DOS 作图、带隙
├── LOCPOT.py                  LOCPOT → 静电势 → 功函数
│
├── tests/
│   ├── selftest.py            68 项自检（无需 VASP）
│   ├── make_fake_results.py   伪造完整结果树
│   ├── measured.json          示例：实测比电容 C_g(y)
│   └── adsorption.json        示例：吸附能 E_ads(y)
│
├── car/                       集群脚本与输入模板
│   ├── INCAR  KPOINTS  KPOINTS2  KPOINTS3  control.in  bader  OUTCAR
│   ├── script  scriptBSC  scriptIQTC  scriptIQTC9
│   ├── scriptBSC_aims  scriptIQTC_aims
│   └── graphical_abstract (v4).png
├── PP/                        赝势库（约 343 个元素目录，见注意事项）
├── species_defaults/          FHI-aims 默认基组
│   └── light  tight  really_tight  non-standard  README
│
├── MXenes_alloy/              〔生成〕VASP 输入树与计算结果
└── results/                   〔生成〕分析产物：CSV + PNG
```

**结果树形态**（`MXenes_alloy/`，自动生成）

```
MXenes_alloy/Ti2V1C2O2_s2x2_outer/ABC/HM/
├── CONTCAR          优化后的结构（dband.py 会向上回溯找到它）
├── opt/             POSCAR POTCAR KPOINTS INCAR     结构优化
├── DOS/             …  LORBIT=11, NEDOS=5000        → dband.py
├── BS/              …  ICHARG=11                    能带
├── WF/              …  LVTOT, LVHAR, 偶极校正       → 功函数
└── bader/           …  LAECHG                       Bader 电荷
```

---

## 四、环境要求与安装配置

### 1. 获取代码

```bash
#fork VASP-MXenes
git clone https://github.com/AshValkan/VASP-MXenes.git
cd VASP-MXenes
git checkout alloy-mxene-tiv        # 若 fork 未带该分支：git fetch origin alloy-mxene-tiv

```

### 2. Python 环境

```bash
python -m venv .venv
source .venv/bin/activate           # Windows: .venv\Scripts\activate
pip install numpy matplotlib pyyaml
pip install pymatgen                # 可选：NEB 的 IDPP
```

### 3. 赝势（POTCAR）

`PP/<元素>/POTCAR` 需自备。默认使用增强赝势，映射见 `alloy.PP_ALIAS`：

| 元素 | 使用的 POTCAR | 理由 |
|---|---|---|
| Ti / V / Nb / Mo / Cr / Zr … | `Ti_pv` `V_pv` `Nb_pv` … | 半芯层 p 电子参与成键 |
| Sc / Y | `Sc_sv` `Y_sv` | s 半芯层 |
| C / N / O / F / H | 标准版本 | — |

> 分支已修复上游的一个静默缺陷：`VASP.MX.POTCAR` 只对第一个元素做赝势映射，
> 第二个金属会退回裸 `V` 而非 `V_pv`。

### 4. 集群脚本

把你的作业调度脚本放到 `car/script`（上游示例为 IQTC / MN4 的 SGE 脚本），
VASP 二进制也在此处指定。

---

## 五、启动 / 运行与构建命令

```bash
# ① 自检（无需 VASP，约 3 秒）
python tests/selftest.py                      # 68 项
python tests/selftest.py -v                  # 保留工作目录便于排查

# ② 生成全部 VASP 输入
python run_alloy.py --config alloy_config.yaml
python run_alloy.py --dry-run                # 只列出将创建的算例

# ③ 在集群上优化与计算（上游流程）
python opt.py  -n 2 -T O2                    # 批量结构优化
python calculate.py -n 2 -T O2               # 批量 DOS / BS
python calculate.py -n 2 -T O2 -WF           # 功函数

# ④ 分析出图（回到本地即可跑）
python dband.py --scan MXenes_alloy --plot
python run_alloy.py --analyze \
       --exp tests/measured.json --ads tests/adsorption.json \
       --host-end <E(Ti3C2O2)> --dopant-end <E(V3C2O2)>
```

**无 VASP 时的完整预演**（验证第 ④ 步）：

```bash
python tests/make_fake_results.py            # 伪造 POSCAR / OUTCAR / DOSCAR
python run_alloy.py --analyze \
       --exp tests/measured.json --ads tests/adsorption.json \
       --host-end -600.0 --dopant-end -607.5
# → results/
#     dband_table.csv    dband_table.json    energetics.csv
#     eps_d_Ti_vs_y.png  eps_d_V_vs_y.png    combined_trend.png
```

> `dband.py --scan` 的趋势图默认写在当前目录（可用 `-o` 指定）；
> `run_alloy.py --analyze` 的产物统一写入 `results/`。

---

## 六、可用脚本与配置项

### 脚本一览

| 脚本 | 主要参数 | 用途 |
|---|---|---|
| `run_alloy.py` | `--config` `--analyze` `--dry-run` `--exp` `--ads` `--host-end` `--dopant-end` | 统一入口 |
| `alloy.py` | `-y` `-n` `-X` `-T` `-s`(超胞) `-p`(模式) `--seed` `--stacking` `--hollows` `-v`(真空) `--encut` `--vdw` | 生成合金输入 |
| `dband.py` | `<DOSCAR>` `--poscar` `--scan ROOT` `--plot` `--emin` `--emax` `-o` | d 带中心与趋势 |
| `adsorption.py` | `<POSCAR>` `--inventory` `--adsorbate` `--site` `--metal-below` `--height` `--one-side` `-o` | 吸附位点与结构 |
| `neb.py` | `<初态>` `<终态>` `-o` `-n`(镜像数) `--method` `--potcar` `--analyze` | CI-NEB |
| `energetics.py` | `--root` `--host-end` `--dopant-end` `-o` | 形成能 / 混合能 / 排序 |
| `analyze_alloy.py` | `--root` `--exp` `--ads` `--pattern` `-o` | 一键出图 |
| `tests/selftest.py` | `-v` | 自检 |
| `tests/make_fake_results.py` | `--root` `--ys` `--patterns` | 伪造结果树 |

### `alloy_config.yaml` 主要配置项

| 配置块 | 关键字段 | 说明 |
|---|---|---|
| `system` | `host` / `dopant` | 主体与客体 M 元素 |
| | `n` | MXene 宽度（2 → M₃X₂） |
| | `X` / `T` | X 元素 / 终止基团（`""` / `O2` / `F2` / `O2H2`） |
| | `stacking` / `hollows` | `ABC`/`ABA`；`HM`(FCC) / `H`(HCP) / `HX` / `HMX` |
| | `ys` | V 含量列表，如 `[0.0, 0.5, 1.0, 1.5]` |
| | `supercell` | `[2, 2, 1]`，须与 `ys` 公度 |
| | `patterns` | 占位模式列表 |
| | `seed` | 随机模式种子（保证可复现） |
| `geometry` | `vacuum` / `interlayer` / `lattice_scale` | 真空层、初始层间距、初始晶格 |
| `dft` | `encut`(550) / `vdw`(d3) / `spin` / `ediff` / `ediffg` / `sigma` | 基础参数 |
| | `lorbit`(**11**) / `nedos` / `dipole` | 电子性质（`lorbit=11` 是 dband 的硬前提） |
| | `hubbard_u: {Ti: 2.5, V: 2.5}` | U_eff（eV，Dudarev 形式） |
| | `kpoints_primitive: [7, 7, 1]` | 原胞 k 网格，超胞自动除 |
| `dband` | `emin` / `emax` | ε_d 积分区间（相对 E_F） |
| `adsorption` | `species` / `site` / `both_sides` / `references` | 参考态能量**须自算** |
| `neb` | `images` / `method` / `spring` / `climb` / `optimizer` | NEB 参数 |
| `experimental` | `capacitance_file` | 实测 C_g(y) 的 JSON / CSV |
| `energetics` | `host_end` / `dopant_end` | 两个端点成员总能量（eV） |
| `output` | `root` / `results` / `dpi` | 输出位置 |

---

## 七、基本使用示例

### 示例 1：生成 Ti–V 系列输入

```bash
python run_alloy.py --config alloy_config.yaml
# 11 个算例写入 MXenes_alloy/（y=1.5 的 inner 模式因内层位点不足被自动跳过并提示）
```

### 示例 2：单个 DOSCAR 的 d 带中心

```bash
python dband.py MXenes_alloy/Ti2V1C2O2_s2x2_outer/ABC/HM/DOS/DOSCAR --plot
```

```python
from dband import PDOS
d = PDOS(".../DOS/DOSCAR")            # 自动向上查找 POSCAR / CONTCAR
d.dband_center("Ti")                  # → -1.70  (eV, 相对 E_F)
d.dband_center("V", emin=-10, emax=6) # 全带一阶矩约定
d.d_moments("Ti")                     # {'center', 'width', 'skewness', 'filling'}
d.layers(n=2)                         # outer_bottom / inner / outer_top
d.summary()                           # 论文表格所需的全部数字
```

### 示例 3：吸附位点与结构

```bash
python adsorption.py <slab POSCAR> --inventory       # 列出表面 O 及其下方金属
python adsorption.py <slab POSCAR> --adsorbate H --site Otop --metal-below V
```

### 示例 4：扩散势垒

```bash
python neb.py ini.vasp fin.vasp -o NEB -n 5 --potcar <POTCAR>
# 用 vasp_vtst 运行后
python neb.py --analyze -o NEB        # → E_a_forward / E_a_reverse / TS_image
```

### 示例 5：一键出图

```bash
python run_alloy.py --analyze --exp tests/measured.json --ads tests/adsorption.json \
       --host-end -600.0 --dopant-end -607.5
```

输出 `combined_trend.png`：ε_d、E_ads(H)、实测 C_g 三条曲线归一化后同轴展示，
并给出 ε_d–y 的 **R²**。没有这个 R²，"V 使 ε_d 上移 ⇒ 容量提升"只是一张图，不是结论。

---

## 八、注意事项

**许可与合规**
- VASP 与 POTCAR 均需授权，**不**随仓库公开分发。在你提交前请确认你的 `PP/` 已被 `.gitignore` 排除。
- 上游原作者文档位于 `README_upstream.md`。

**计算参数**
- `LORBIT = 11` 是 `dband.py` 的硬前提；缺少它 DOSCAR 无 lm 分解，ε_d 无法计算。
- ε_d 的积分区间一旦选定，**全部组分必须一致**，并在论文中写明用的是"占据态"还是"全带"约定。
- U_eff 对 ε_d 有实质影响，建议做 U = 0 / 2 / 3 eV 的敏感性分析并说明取值依据。

**物理与化学**
- **H_M（FCC）终止下，表面 O 正对中间层金属**。因此 `pattern=outer` 时吸附位下方是 Ti，
  `pattern=inner` 时下方是 V——两者通过**完全不同**的机制影响吸附，不可当作两个随机猜测。
- 分数组分需要**公度超胞**：`y=0.5` / `1.5` 要求超胞面积为偶数（2×2 可行，3×3 不可行）。
  `suggest_supercell()` 会给出最小可行胞（**面积优先**，方正次之——先小后方正可省下约 20 倍算力）。
- `V₃C₂Tₓ 仅作为 E_mix 的计算端点。
- 吸附参考态（½E(H₂)、体相 Li/Na/K）**必须用你自己的赝势和泛函计算**，不可引用文献数值。

**软件与环境**
- 兼容 numpy 2.x：`ndarray.ptp()` 与 `np.nanptp` 在 2.0 后不可用，勿改回。
- POTCAR 查找路径是相对当前工作目录的 `./PP/...`，**请在仓库根目录运行脚本**。
- 涉及到使用wsl环境的，Windows 下 git-bash 的 `/tmp` 与 Python 不互通，临时文件请放在项目内。
- 输出文件名按 key 全量 sanitize，`eps_d(Ti)` / `eps_d(V)` 分别输出，不会互相覆盖。

**数据处理**
- 加入端点成员（y=0 无 V、y=3 无 Ti）时，作图逻辑按"至少两个有效点"判断，
  不会因个别组分缺少某元素而静默丢图。
- 三电极比电容按**单电极质量**计，两电极按**双电极总质量**计，二者理论相差 4 倍，
  论文中必须标明所用定义。
- 若图中 H 吸附能（E_ads）随 V 含量变得更负，按原始值绘制会向下走；
  论文中更可读的画法是绘制 −E_ads 使其与 ε_d 同步向上。脚本保留你的符号约定权。

---
