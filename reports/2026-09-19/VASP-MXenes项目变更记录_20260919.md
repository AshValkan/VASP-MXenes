# VASP-MXenes 项目变更记录

**项目**：Ti₃₋ᵧVᵧC₂O₂ MXene 合金化体系 DFT 计算（硕士论文 Q2/Q3 支撑数据）
**整理日期**：2026-09-19
**仓库**：fork AshValkan/VASP-MXenes，分支 alloy-mxene-tiv
**服务器**：<user>@<server> → `~/calc/VASP-MXenes-alloy`（本地工作目录 `C:\Users\Atlas\WorkBuddy\2026-09-13-17-06-37`）

---

## 1. 计算环境要点

| 项 | 内容 |
|---|---|
| VASP | vasp.6.4.2 std，二进制 `~/calc/vasp_std` |
| MPI 栈 | OpenMPI 5.0.10（/usr/local）+ PMIx **5.0.1**（`~/calc/pmix-5.0.1`；用户自装 6.1.0 与 PRRTE 不兼容，勿用） |
| MPI 启动 | `LD_LIBRARY_PATH=/home/atlas/calc/pmix-5.0.1/lib:/usr/local/lib` + `/usr/local/bin/mpirun` |
| 硬件 | Xeon E5-2630v4 ×2，WSL2 虚拟机（内核 microsoft-standard-WSL2）：物理双路被抹平为单 socket 视图，20 虚拟核 ×2 SMT = 40 逻辑 CPU；numactl 不可用，taskset 可用 |
| 核块划分 | 0-9 / 10-19 / 20-29 / 30-39（每块 = 5 物理核 = 10 SMT 线程），4 个 np5 任务并行 |

## 2. 计算体系与基准参数

- **体系**：Ti₃₋ᵧVᵧC₂O₂（y = 0 / 0.5 / 1.0 / 1.5），s2×2 超胞、outer 层 V 占位、O 双面封端
- **路径**：`MXenes_alloy/{Ti3V0C2O2,Ti2.5V0.5C2O2,Ti2V1C2O2,Ti1.5V1.5C2O2}_s2x2_outer/ABC/HM/{opt,DOS,ads_H,ads_Li,ads_Na,ads_K}`
- **INCAR 基准**（opt）：PBE/PAW、ENCUT=550 eV、GGA+U（U=2.5 eV 于 Ti/V d，LMAXMIX=4）、ISPIN=2、IVDW=11、SIGMA=0.10、EDIFF=1e-5、EDIFFG=-0.02、AMIX=0.2/BMIX=0.0001（含 _MAG 变体）、ALGO=Normal、NELM=200、LREAL=AUTO、k=3×3×1；**无 LDIPOL/IDIPOL/DIPOL**（对称 O 封端无需偶极校正，模板已修）
- **DOS**：k=5×5×1、ISMEAR=-5、LORBIT=11、NEDOS=5000

## 3. 服务器系统级修复（时间线）

### 09-15 loopback 策略路由劫持（致命）
- 现象：所有本机 TCP 与 MPI 全断。
- 根因：策略路由 `1: ipproto tcp lookup 127` 把 127.0.0.1 劫持到坏 relay（loopback0/169.254.73.152）。
- 修复：`sudo ip rule del priority 1 ipproto tcp lookup 127`；已加 @reboot cron 自愈。
- 恢复命令（如复发）：`sudo ip rule add priority 1 ipproto tcp lookup 127` 是肇事规则，删除即可。

### 09-18 PRRTE rank 堆叠（性能，第一形态）
- 现象：4 个 np5 任务并发时每进程仅 ~60% CPU，其余核闲置。
- 根因：PRRTE 默认 bind-to-core，每个独立启动的任务都从机器核 0-4 开始绑，多任务叠在同一批核上。
- 修复：运行中任务在线 `taskset -pc <核块>` 重绑（实测 ~60%→~110%）；启动器加第 4 参数 CPUS。注意：**WSL2 内无物理 NUMA 视图，核块绑定是等效替代方案**。

### 09-19 PRRTE 覆盖 taskset 事故（致命，已根治）
- 现象：队列 01:09–01:39 把 11 个碱金属任务全部启动，55 个 rank 堆在 CPU 0-9，每 rank 仅 ~16% CPU 爬行 15 小时；队列日志反复显示 "launching … on cpus 30-39"。
- 根因：`taskset -c $CPUS mpirun` 只约束 wrapper 本身；**PRRTE 启动 rank 时以自身 bind-to-core 策略重新 sched_setaffinity，完全无视继承的 taskset 掩码**（rank 掩码呈 0,1 / 2,3 / 4,5 … 序列即此特征）。block_busy() 扫描掩码时 30-39 永远"空闲"，队列无限补位。
- 修复（三件套）：
  1. `run_direct_mpi.sh` CPUS 分支改为 `taskset -c $CPUS mpirun --bind-to none …`（旧版备份 `run_direct_mpi.sh.bak_0919`）；探针实测 `--bind-to none` 下 rank 掩码继承 30-39 ✓
  2. 保留 4 个进度最多任务在线重绑到互不重叠核块（重绑后每 rank ~115% 满速）
  3. 其余 7 个堆叠任务终止；旧输出移入各自 `restart_bak_0919/`；Ti2V1/ads_K 的 CONTCAR 回存 POSCAR 保住 1 步弛豫；队列文件重写、队列重启待命
- 教训：**taskset 包 mpirun 对 PRRTE 无效，必须 `--bind-to none`**；绑定正确性要用真实探针读 `/proc/<pid>/status` 验证，不能只看启动日志。

## 4. 仓库代码变更（相对上游 fork）

| 文件 | 变更 |
|---|---|
| `alloy.py` | INCAR 模板 **DIPOL 动态化**：从 POSCAR 解析 slab 原子 z 范围、wrap 处理后计算板中心写入 DIPOL（原硬编码 0.5 0.5 0.5 落在 21 Å 真空层中央，引发偶极反馈振荡，最终删除偶极校正但保留动态计算能力）；selftest 68/68 通过 |
| `neb.py` | 删除硬编码 DIPOL |
| `dband.py` | 适配 VASP 6.4.2 真实 DOSCAR：① NEDOS 头行 4→4-5 列；② 数据起始偏移修正；③ 自旋判据（ncol≥5）；④ 块分隔符改为识别重复头行（6.4.2 不用空行）。selftest 68/68 |
| `run_direct_mpi.sh`（新增） | MPI 启动器：`run_direct_mpi.sh <CASE> [NP=5] [SUB=opt] [CPUS]`；opt 收敛跳过；CPUS 分支 taskset + **--bind-to none**（09-19 修正） |
| `run_ads_queue.sh`（新增） | 核块感知自适应队列：扫 `/proc/*/status` Cpus_allowed_list，整块空闲才启动下一个任务；依赖 `.refs_done` 哨兵 |
| `start_queue.sh`（新增） | 前台启动器（内部 setsid nohup），规避 ssh 含 `&` 命令挂起问题 |
| `setup_adsorption.py` / `setup_alkali.py`（新增） | H / 碱金属吸附任务批量部署（INCAR 沿用 opt + MAGMOT/LDAU 块按物种数补齐） |
| `run_bulk_refs.sh`（新增） | Li/Na/K bcc 参考态（ISIF=3、13×13×13、ISMEAR=0/SIGMA=0.05） |
| `analyze_adsorption.py` / `analyze_alkali.py`（新增） | 吸附能分析 + matplotlib 绘图（E_ads = [E(slab+2M) − E(slab) − 2E_ref(M)]/2）；H 版含出版级样式（sci 计数轴、加粗、1:1 方图） |
| `analyze_alloy.py`（上游已有） | d 带中心/宽度提取 |

## 5. 计算配置坑（重要经验）

1. **LMAXMIX 必须 =4**（Ti/V d + U 体系）；=2 会导致磁性体系 SCF 不收敛。
2. **DIPOL 参考点错位**：slab 在 z∈[0.037,0.259] 而 DIPOL=0.5 0.5 0.5 落在真空层 → 偶极锯齿势失真 → ±100 eV 震荡。判断特征：能量下探到物理值（~−1480 eV）后弹回 0 附近震荡；对称封端却报 ~9 e·Å 假偶极。
3. **OSZICAR 有 Fortran 4K 缓冲滞后**：监控文件进度看 vasp.out/OUTCAR 的 mtime，别信 OSZICAR 实时性。
4. **LDAU 维度错**：向已有体系加新物种（H/碱金属）时，LDAUU/LDAUJ/LDAUL 条目数必须按新物种数补齐（H/K 位填 -1/0），否则 VASP 报维度错。
5. **K_sv 赝势**：仓库 PP 无普通 K，用 K_sv（ZVAL=9，含 3s3p 半芯）；ENCUT=550 覆盖其 ENMAX=259。
6. **SSH 工具链**：Windows 本地无 sshpass，用 managed venv 的 paramiko（工作区 `ssh_exec.py` / `ssh_put.py` / `ssh_get.py`）；ssh_exec 对含 `&` 的远程命令会挂起，一律走"脚本内 setsid + 前台调用"。

## 6. 计算结果汇总（截至 2026-09-19 18:00）

### 6.1 弛豫 opt（4/4 收敛）
| y | E0 (eV) | mag (μB) |
|---|---|---|
| 0 | −236.320 | 0.00 |
| 0.5 | −234.983 | 0.008 |
| 1.0 | −233.109 | 2.45 |
| 1.5 | −231.306 | 4.42 |

每 V 磁矩随 y 递增至 ~0.74 μB/V。

### 6.2 DOS 电子结构（4/4 完成）
- E_fermi：−1.83 / −1.72 / −1.77 / −2.04 eV
- ε_d(Ti)：−3.572 / −3.599 / −3.494 / −3.296 eV（线性拟合 0.187·y − 3.630，R²=0.775；y=0.5 为稀掺杂反常点）
- ε_d(V)（y≥0.5）：−3.367 / −3.368 / −3.148 eV，恒高于 Ti 位
- N(E_F)：0.55 / 7.90 / 9.97 / 9.88 states/eV（较基线 ×18，金属性显著增强）

### 6.3 H 吸附 Q3（4/4 收敛，正式）
| y | E_ads (eV/H) |
|---|---|
| 0 | −0.427 |
| 0.5 | −0.397 |
| 1.0 | −0.616 |
| 1.5 | −0.734 |

E_ref(H)=−3.38623 eV（E_H2=−6.77245）；趋势与 ε_d(Ti) 上移一致（Hammer–Nørskov）。附产物：H 吸附使磁矩量子化重整为 ~2 μB/V。

### 6.4 碱金属吸附 Q3（进行中）
- bulk 参考态 3/3 收敛：Li −2.0027、Na −1.4648、K_sv −1.2245 eV/atom（残余压力 ≤0.4 kB）
- ads_Li(y=0) 已收敛（16 步，E0=−242.769）
- 09-19 16:30 修复后运行中：Ti3V0/ads_Na（步 4+）、Ti2V1/ads_Li、Ti2V1/ads_Na、Ti3V0/ads_K（各步 1+），核块 0-9/10-19/20-29/30-39 隔离满速
- 队列待启动 7 个：Ti2.5V0.5×{Li,Na,K}、Ti1.5V1.5×{Li,Na,K}、Ti2V1/ads_K

## 7. 本地工作目录清单（`C:\Users\Atlas\WorkBuddy\2026-09-13-17-06-37`）

- `计算进展汇报/`：汇报 PPTX（9 页）+ slides + STORY.md + results/（opt_dos_merged.csv 含注释与文献讨论列、eps_d_Ti_vs_y_barline.svg、adsorption_H.csv、E_ads_vs_y.png）
- 远程工具：`ssh_exec.py` / `ssh_put.py` / `ssh_get.py`（paramiko）
- 部署/分析脚本：`setup_adsorption.py`、`analyze_adsorption*.py`、`setup_alkali.py`、`analyze_alkali.py`、`run_bulk_refs.sh`、`run_ads_queue.sh`、`start_queue.sh`、`run_direct_mpi*.sh`
- 数据/图表生成：`make_merged_and_chart.py`（纯标准库 SVG 柱线图）
- 事故处置：`remediation_0919.sh`、`verify_fix_0919.sh`、`run_direct_mpi_fixed.sh`
- PPT 工具：`run_slidep.py`、`rebuild_pptx.py`
- 项目记忆：`.workbuddy/memory/`（MEMORY.md + 每日日志 09-13 ~ 09-19）

## 8. 备份内容说明（本次 Z 盘备份）

```
Z:\WorkBuddy_Backup\2026-09-19\
├── workspace\                                  # 本地工作目录全量（含 .workbuddy/memory 项目记忆与本文档）
├── user_memory_cache\                          # 云端用户画像本地缓存（~/.workbuddy/memory；用户级 MEMORY.md 从未创建，无此文件）
└── server_snapshot\mxene_snapshot_20260919.tar.gz   # 服务器仓库快照（46.6 MB，含 .git 历史与全部输入/日志，排除 POTCAR/OUTCAR/CHGCAR/WAVECAR/vasprun.xml/DOSCAR 等大体积二进制输出与可重建的 .venv）
```

**服务器快照恢复**：解包至 `~/calc/VASP-MXenes-alloy`，POTCAR 需从仓库 PP/ 目录按 POSCAR 物种顺序重新拼接；大体积输出（OUTCAR/CHGCAR 等）不在快照内，如需完整输出请另行从服务器按 case 拷贝。

## 9. 后续待办

1. 碱金属 12 任务全部收敛后跑 `analyze_alkali.py` → `adsorption_alkali.csv` + `E_ads_alkali_vs_y.png`
2. E_ads(H) vs ε_d(Ti) 关联图（论文 Q2/Q3 交叉证据）
3. 候选：CI-NEB 离子迁移势垒、形成能交叉验证、图表插入 PPT/论文
