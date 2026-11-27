# 统计分析计划（Statistical Analysis Plan, SAP）v1.0

- 项目：ICU 糖尿病患者新发房颤（new-onset atrial fibrillation, NOAF）与已有房颤、无房颤的预后比较，以及糖尿病的效应修饰
- 版本：**v1.0（锁定）**
- 锁定日期：2026-11-20（Pacific/Apia，UTC+13）
- 批准：史教授批准按心电图验证后的修正定义撰写本 SAP 并开始主分析
- 报告规范：STROBE（队列）
- 改动规则：锁定后的任何修改记入 `docs/SAP_changelog.md`，并升版本（v1.1…），不覆盖本文件

本 SAP 在跑主分析前锁定。可行性与心电图验证见  
the local feasibility project (`projects/af_dm_noaf_feasibility/REPORT.md` and `REPORT_ECG_VALIDATION.md`; path is environment-specific and is not shipped with this repository).

---

## 0. 合规

- 原始区只读；原 `mimic.duckdb` / `.wal` 仅 `ATTACH … (READ_ONLY)`，不写。
- 个体数据不出本机；脚本只用 `Rscript`；不打开、不打印、不把个体行写到项目目录。中间分析只在内存或系统临时空间，导出仅汇总。
- 对外表计数或率出现 ≤9 时掩码为「≤9」，精确值只放 `output/INTERNAL_ONLY/`。
- 不新装 R 包。本机已有：R 4.5.1，duckdb 1.5.6，survival 3.8.3，survminer 0.5.1，ggplot2 4.0.0，cmprsk 2.2.12，tidycmprsk 1.1.0，mice 3.18.0，rms 8.0.0，tableone 0.13.2，forestplot 3.1.7。**未安装 EValue**，E-value 按 VanderWeele & Ding (2017) 公式手算，不新装包。
- 数据使用权：史教授已于 2026-11-20 确认（见 `metadata/licenses_TODO.md`）。MIMIC-IV / eICU 为 credentialed DUA；MIMIC-IV-ECG v1.0 为 ODbL 开放，但与临床表连接仍须遵守 MIMIC-IV DUA。发表时按 PhysioNet 要求写数据来源与致谢。
- 随机种子：`20261120`。软件与种子写入结果附录。

---

## 1. PECO 与库角色

**P（人群）**：MIMIC-IV 3.1 成年（`anchor_age≥18`）首次 ICU 入住者。主分析限制 `anchor_year_group` 为 2008–2019（心电图模块有覆盖的年份）。糖尿病为关注人群；同期建立非糖尿病对照，用于修饰分析。

**E（暴露）**：入 ICU 时/后的房颤状态，三组——已有房颤（pre-existing AF）、新发房颤（NOAF）、无房颤（no AF）。定义见第 3 节（已纳入入院前心电图机器判读）。

**C（对照）**：主比较为 NOAF vs 无房颤；次比较为已有房颤 vs 无房颤，以及 NOAF vs 已有房颤。非糖尿病队列用于检验糖尿病是否修饰上述关联。

**O（结局）**：主要——住院死亡、28 天全因死亡、1 年全因死亡。次要——ICU 住院时间（ICU LOS）、本次住院卒中诊断、出院时（出 ICU 时）仍为房颤节律。

**库角色**

| 库 | 角色 |
|---|---|
| MIMIC-IV 3.1 临床表 + mimic-code 派生表 | 发现队列：暴露、协变量、住院与随访死亡 |
| MIMIC-IV-ECG v1.0 元数据（机器报告，无波形） | 修正「已有房颤」：入院前诊断心电图机器判读房颤；敏感性：仅 ECG 确认的 NOAF |
| eICU-CRD 2.0 | 仅住院死亡的**粗**外部验证（口径不可与 MIMIC 对齐，见第 9 节） |

**与已发表研究的区别（基于前期查新，含 Alomari 2025, PMID 40496587；文献 NOAF 发生率约 4.6%–24%）**：多数 ICU NOAF 研究以监护心律和/或 ICD 定义，较少用诊断心电图把「入院前已有房颤」从新发中剔除，也较少把 NOAF 当作入 ICU 后的时变暴露以处理 immortal time，并在同一首次 ICU 队列中同时报告糖尿病与非糖尿病的乘法与加法（RERI）交互，以及血糖变异度/低血糖分层。本研究在修正定义的基础上做这些分析。若后续查新有高度重叠工作，再在 changelog 中注明。

---

## 2. 人群与流程图

**纳入**

1. MIMIC-IV `icustays` 与 `patients` 连接，`anchor_age≥18`。
2. 每名 `subject_id` 仅保留时间上第一次 ICU（`intime`、`stay_id` 排序，`ROW_NUMBER=1`）。
3. 能与 `admissions` 连接（同一 `subject_id, hadm_id`）。

**主分析年份**：`anchor_year_group ∈ {2008 - 2010, 2011 - 2013, 2014 - 2016, 2017 - 2019}`。  
**敏感性**：加入 `2020 - 2022`（该年组心电图覆盖为 0，定义退化为无 ECG 信息）。

**糖尿病（主定义）**：本次住院 `diagnoses_icd` 中 ICD-9 `250*` 或 ICD-10 `E10–E14`（去点后正则 `^E1[0-4]`；不用 `SIMILAR TO` 的 `%` 写法）。此为**主分析糖尿病队列**。

**糖尿病敏感性扩展**（不改变主定义）：在主定义阴性者中，若满足以下任一，则纳入「扩展糖尿病」队列并重复主模型：

- 入院前 90 天至入院后 24 小时内，距入院最近一次 HbA1c（itemid 50852）≥6.5%；或
- 医院入院至入 ICU 后 24 小时内医嘱中出现降糖药（胰岛素、二甲双胍、磺脲类、SGLT2 抑制剂、GLP-1RA 等，药名对照见 codebook）。

**重要限制（预声明）**：当前仓库的 MIMIC-IV 3.1 包**没有 ED / medrecon**，无法可靠识别真正的院前用药。上述药物标准是住院早期医嘱的**代理**，不是门诊处方。可行性研究已表明 SGLT2i/GLP-1RA 暴露极稀。SAP 不把「院前降糖药」写成已实现的金标准。

**非糖尿病对照**：主定义糖尿病为阴性者（扩展定义仅用于敏感性，不改变对照的主口径）。

**排除出三组主分析、但在流程图中计数者**：无 Heart Rhythm（itemid 220048）记录者；以及「无法判定 / indeterminate」（见第 3 节）。灰区单列，不进入三组主模型（敏感性再归组）。

记录每步人数，用于纳入流程图（STROBE）。

---

## 3. 暴露定义（锁定的修正定义）

数据来源：监护仪心律 `chartevents` itemid **220048**（Heart Rhythm）；既往与本次 ICD；MIMIC-IV-ECG 机器报告文本。心律 AF/AFL 规则与可行性研究相同：`atrial fib|a-fib|afib|(^|[^a-z])af([^a-z]|$)|a flut|flutter`。窦律家族：取值含 `sinus`（SR/ST/SB/SA）。心电图机器房颤正则与 `af_dm_noaf_feasibility/output/05_ecg_regex.csv` 一致（否定词剥离后匹配 atrial fibrillation / flutter / a-fib / afib 等，不用单独 `af` 词；主分析用当前房颤或房扑，含可疑；确定性作为敏感性）。

时间窗：入院前 ECG 为 `ecg_time < intime`（任意时点，不限于 1 年）；入 ICU 2 小时窗为 `first_af_time ≤ intime + 2 hours`。`ecg_time` 为 TIMESTAMP（非纯日期）。

### 3.1 已有房颤（pre-existing AF）

满足任一：

- 本次住院之前任一住院的 AF ICD（ICD-9 `42731` 或 ICD-10 `I48*`）；或
- 入 ICU 前任意时点诊断心电图机器判读为当前房颤/房扑；或
- 入 ICU 2 小时内监护已记录 AF/AFL。

### 3.2 NOAF（主定义，推荐）

同时满足：

- 不具备 3.1 中任一条；
- 有 220048 记录，且**首条**心律为窦律家族，且首条不是 AF/AFL；
- 此后（同一 ICU 住院，至 `outtime`）出现至少一次 AF/AFL。

首次 AF 时刻记为 `t_af`（用于时变暴露）。

### 3.3 无房颤（no AF）

不具备 3.1；有 220048；整个 ICU 期间从未记录 AF/AFL；且**不是**灰区（即本次住院无 AF ICD）。

### 3.4 灰区（gray zone，单列）

本次住院有 AF ICD，但：无既往住院 AF ICD、无入院前 ECG 机器房颤、ICU 内监护从未记录 AF/AFL（须有心律记录）。主分析不并入三组。敏感性归组：（A）归入已有房颤；（B）归入无房颤；（C）保持排除。

### 3.5 无法判定（indeterminate）

例如：无心律记录；或无 3.1，且首条非窦律、但后来出现 AF（扩大定义下可为 NOAF，主定义排除）；或其他不能归入 3.1–3.4 者。流程图报告人数，主模型不纳入。

### 3.6 预指定的暴露敏感性

- 扩大 NOAF（rec_2h 风格）：无 3.1，首条为任意非 AF 节律，其后出现 AF/AFL。
- ≥2 次连续 AF 记录。
- 连续 AF 段持续时间 ≥1 小时。
- 仅 ECG 确认的 NOAF：主定义 NOAF 且 ICU 窗或「入住至出院」窗内至少一份机器报告当前 AF/AFL（无窗口内 ECG 者在此敏感性中排除，不把缺失当成阴性）。
- 不加入入院前 ECG 标准（回到仅 ICD+2h 监护的旧已有房颤定义），以量化心电图修正的影响。

---

## 4. 结局

时间零点：ICU `intime`（landmark 分析除外，见第 5 节）。

**住院死亡**：`hospital_expire_flag=1`。主分析同时报告：(i) 粗率与 Logistic 回归（OR）；(ii) 以出院为删失、死亡为事件的 Cox（住院期间）。存活出院者在 `dischtime` 删失。

**28 天 / 1 年全因死亡**：事件 = 住院死亡，或 `dod` 非缺失且 `dod ≤ intime + 28/365 天`。存活者在窗口结束时删失。若无 `dod` 且存活出院，视为该窗口内删失（即当作仍存活），**会低估死亡率**（MIMIC `dod` 不完整）。此限制在所有相关结果中写明，1 年死亡为次要/敏感性结局。

**ICU LOS**：`outtime - intime`（小时或天）。描述用中位数（IQR）；比较用分位数回归或 Wilcoxon（完整描述），并以出 ICU 存活为事件、ICU 内死亡为竞争风险的 Fine-Gray 模型（亚分布风险）。

**卒中**：本次住院 ICD 为脑卒中（缺血或出血：ICD-9 430–434、436；ICD-10 I60–I64，去点后匹配；详表见 codebook）。**无发病时间、无 POA 标志**，不能区分入院即存在与院内新发。因此卒中做住院期间二分类（Logistic），**不做**以臆造时间为基础的 Fine-Gray。此为方法学限制，不是分析结果。

**出院时 AF**：出 ICU 前 6 小时内最后一次 220048 为 AF/AFL。分母：存活出 ICU 且该 6 小时内有心律记录者。死亡者排除出该二分类（另报告死亡人数）。敏感性：本次出院诊断含 AF ICD。

---

## 5. 时间偏倚与生存模型（主分析必须处理 immortal time）

NOAF 在入 ICU **之后**才发生，从 `intime` 到 `t_af` 的存活时间属于「尚未暴露」。若把 NOAF 当作基线固定暴露，会引入 immortal time bias，高估或扭曲 NOAF 与死亡的关联。

**主分析：时变 Cox（time-varying Cox, counting-process）**

- 已有房颤：从 `intime` 起即为「已有房颤」，不随时间改变。
- 无房颤：全程为「无房颤」。
- NOAF：`[intime, t_af)` 编码为「无房颤」（或单独的 at-risk 状态，主模型中与无房颤相同的参照编码，即暴露=0）；`[t_af, 随访结束]` 编码为 NOAF。若事件发生在 `t_af` 之前，该例对 NOAF 暴露不贡献事件（避免 immortal time）。
- 模型：`coxph(Surv(tstart, tstop, status) ~ exposure_tvc + covariates)`。暴露为分类（时变 NOAF + 基线已有房颤，参照无房颤）。
- 分别拟合住院死亡（出院删失）、28 天死亡、1 年死亡。报告 HR 与 95% CI。
- 比例风险：Schoenfeld 残差；若 NOAF 明显违例，加报告按时间分段的 HR（如 0–7 天、7–28 天）。

**预指定替代：48 小时 landmark**

- 分析集：`intime` 后 48 小时仍存活者（landmark 前死亡排除，并报告人数）。
- 在 48 小时处分类：已有房颤（定义同 3.1）；NOAF（`t_af ≤ 48h` 且符合 NOAF）；无房颤（48h 内无 AF）。48h 后才出现的 AF 在该模型中仍算无房颤（限制在结果中说明）。
- 从 48 小时起随访 28 天 / 1 年（及住院死亡，若仍在院）。

**竞争风险**：非死亡结局中，ICU 存活出院用 Fine-Gray（ICU 内死亡为竞争）。卒中因无时间，不用 Fine-Gray（第 4 节）。

**描述**：按组报告粗死亡率、随访人时发生率（28 天 / 1 年），以及 Kaplan–Meier 曲线（固定基线分类的 KM **仅用于描述**，不作为主因果估计；图注写明未校正 immortal time）。主推断以时变 Cox 为准。

---

## 6. 协变量与缺失

均为**基线**（入 ICU 后 24 小时内，或本次住院入 ICU 前已确定的特征），避免把暴露后的中介因子当基线，除非另有说明。

| 协变量 | 来源（官方表/字段） | 说明 |
|---|---|---|
| 年龄、性别 | `patients.anchor_age`, `gender` | |
| 种族 | `admissions.race` | 合并为 White / Black / Hispanic / Asian / Other or Unknown |
| 入院类型 | `admissions.admission_type` | 急诊（EW EMER.、DIRECT EMER.、URGENT 等）、择期（ELECTIVE、SURGICAL SAME DAY ADMISSION）、其他 |
| 首日 SOFA | `mimiciv_derived.first_day_sofa.sofa` | 只读挂载原 duckdb |
| SAPS II | `mimiciv_derived.sapsii.sapsii` | |
| Charlson | 由本次住院 `diagnoses_icd` 按 Quan 等 ICD-9/10 权重计算 | 仓库无现成 charlson 派生表 |
| 脓毒症 | `mimiciv_derived.sepsis3` 在该 `stay_id` 有记录 | Sepsis-3 派生定义 |
| 心衰、高血压、CKD | 本次住院 ICD（codebook） | |
| 心脏手术 | 本次住院服务含 `CSURG`，或 CABG/瓣膜等操作 ICD（codebook） | 排除心脏手术为敏感性 |
| 机械通气 | `mimiciv_derived.ventilation` 在入 ICU 24h 内 `ventilation_status='InvasiveVent'` | |
| 血管活性药 | 入 ICU 24h 内：派生表去甲肾上腺素/肾上腺素/多巴胺/多巴酚丁胺，或 `inputevents` itemid 221906、221289、221662、221653、221749（去氧肾上腺素）、222315（血管加压素） | |
| 钾 | `first_day_lab.potassium_min` / `_max`（取首次 24h 最低与最高，模型用最低或均值，表 1 两者都列） | |
| 镁 | `labevents` itemid 50960，入 ICU 24h 内 | 派生 first_day_lab 无镁 |
| 首日血糖 | 化验 50931、50809、52569、52027 与 chart 220621、226537、225664、228388，入 ICU 24h，按分钟去重 | 均值用于协变量；CV、低血糖用于分层 |
| 年份组 | `anchor_year_group` | 主分析已限制 2008–2019，模型仍可调整年份组 |

**缺失**

- 先出缺失表（按三组 × 糖尿病）。
- **主分析用多重插补**（`mice`，已安装）：只插补基线协变量，不插补暴露与结局；`m=10`，`seed=20261120`，方法以 predictive mean matching（数值）和 logistic/多项式（分类）为主。各套插补分别拟合时变 Cox，再按 Rubin 规则合并。
- 敏感性：完整病例分析（complete-case）。
- HbA1c / SHR：已知缺失率高且随房颤分类不同（已有房颤组更高），属非随机缺失风险。SHR **仅在有 HbA1c 的子集做探索**，不插补 HbA1c 后当全队列协变量，并在文中说明选择偏倚。

---

## 7. 修饰（interaction / effect modification）

1. **糖尿病 vs 非糖尿病**（主修饰）  
   - 乘法：在合并队列的时变 Cox 中纳入 `NOAF（时变）× has_dm`、`已有房颤 × has_dm`，报告交互项 P 与各层 HR。  
   - 加法：以「无糖尿病且无房颤」为参照，估计糖尿病、NOAF 及其并存的 HR，计算 RERI、attributable proportion（AP）、synergy index（S）及 95% CI（delta 方法或 bootstrap，B=500，种子同上）。主 RERI 针对 NOAF 与糖尿病（已有房颤另报）。  
   - 若交互与多重插补同时进行，在各套插补中计算后合并，或在完整病例中做加法交互作为敏感性；SAP 主报告：插补后的乘法交互 + 完整病例或插补合并的 RERI（实施时在结果中写明实际采用的合并方式）。

2. **糖尿病内部地层**（仅糖尿病队列，时变 Cox，与主模型协变量相同，不另做三重交互，除非样本足够）：  
   - 入 ICU 24h 血糖变异度（CV = SD/mean，按三分位）；  
   - 24h 内是否出现低血糖（任一血糖 &lt;70 mg/dL）；  
   - 首日平均血糖三分位。  
   报告各层 HR 与层内人数；若某层事件数过少（对外 ≤9），只在 INTERNAL_ONLY 写精确数，对外不估计或不展示该层模型。

3. **SHR**（血糖/糖化校正）：仅 HbA1c 非缺失子集，探索性，明确非随机缺失。

---

## 8. 敏感性（预指定）

1. 第 3.6 节各 NOAF / 已有房颤定义。  
2. 仅 ECG 确认的 NOAF 子集。  
3. 排除心脏手术。  
4. 加入 2020–2022 年组。  
5. 灰区归组 A/B/C。  
6. 完整病例 vs 多重插补。  
7. 48h landmark。  
8. E-value（主模型 NOAF vs 无房颤的 28 天及住院死亡 HR；点估计与 CI 靠近无效值的一侧）。公式（HR&gt;1）：\(E=\mathrm{HR}+\sqrt{\mathrm{HR}(\mathrm{HR}-1)}\)；HR&lt;1 先取倒数。不新装 EValue 包。  
9. 糖尿病扩展定义（第 2 节）。

---

## 9. eICU 粗验证（预声明限制）

**不能**复制 MIMIC 的心律时序 + 心电图定义。可行口径（与可行性报告一致）：

- 成年首次 ICU；  
- 可能 NOAF：无 `pastHistory` 房颤，且诊断表中最早 AF 的 `diagnosisoffset>120` 分钟；  
- 已有房颤：`pastHistory` 或 offset≤120 的 AF；  
- 无房颤：无 AF 诊断且无既往史；  
- 糖尿病：诊断或既往史文本/ICD（可行性脚本口径）。  

仅比较**住院死亡**的粗率与（尽量）校正 Logistic/Cox。明确这是近似口径，`diagnosisoffset` 更可能是录入时间而非起病时间，且 eICU **不能**评估 28 天/1 年死亡，也无 HbA1c。不把 eICU 结果解释为对主定义的严格外部验证。若方向一致，写「粗验证方向一致」；不一致则讨论口径差异，不做 meta 强行合并。

---

## 10. 统计与图表

- 检验水准双侧 α=0.05，不因多重比较而改主结论的声称，但敏感性与分层解释为支持性。  
- Table 1：三组 × 糖尿病（主分析年份、进入三组者）。连续变量中位数（IQR）或均数（SD），分类变量 n（%）。组间比较： Kruskal–Wallis / χ² 或 Fisher（小格）。对外掩码 ≤9。灰区另表。  
- 发生率：NOAF 在「有风险」（非已有房颤）中的累积发生率；按糖尿病分层。  
- KM：描述性，PDF + PNG，300 dpi，期刊风格（无灰底、色盲友好配色、风险表、95% CI 带）。  
- 森林图：主模型与关键敏感性的 HR（95% CI），同样 PDF+PNG、300 dpi。  
- 三线表导出为 csv（对外 + INTERNAL_ONLY）；不把个体数据或中间文件放进对外目录。

---

## 11. 软件与复现

- 入口脚本在 `scripts/`，`nohup` 可续跑；已完成的汇总输出若完整则跳过。  
- 可复用 `af_dm_noaf_feasibility` 中的队列/心律/心电图 SQL，但暴露必须以本 SAP 第 3 节为准。  
- 不把个体分析数据集写进项目目录。

---

## 12. 已预见到的偏离（实施前写明，不算事后改口）

1. 真正「院前降糖药」在本仓库无法实现，改为住院早期医嘱代理，并可能因暴露过稀而不建模、只做描述。  
2. 卒中无时间戳，不用 Fine-Gray。  
3. E-value 用手算公式。  
4. `dod` 不全导致 28 天/1 年死亡低估。  
5. 2020–2022 无 ECG，加入该年组的敏感性会使「已有房颤」定义不完整。  
6. 仓库无现成 Charlson 表，用 Quan 编码自行计算，编码表写入 codebook，供核对。

锁定后若再改定义、窗口或主模型，必须升版本并在 changelog 写明理由与对已跑结果的影响。
