# Primer Panel Pipeline

从人类 hg38 基因符号生成 PCR 引物面板（Stage 1: Target Generation → Stage 2: Primer Design → Stage 3: Specificity）。

## 核心概念

本流水线严格分离三个阶段的职责：

| 阶段 | 职责 | 输入 | 输出 |
|------|------|------|------|
| **Stage 1** | 决定"要扩增哪个区域" | 基因符号 | SEQUENCE_TEMPLATE + SEQUENCE_TARGET |
| **Stage 2** | Primer3Plus-like 引物设计 | SEQUENCE_TEMPLATE + SEQUENCE_TARGET | primer pair + Primer3 指标 |
| **Stage 3** | in-silico PCR 特异性检查 | primer pair | genomic product 坐标 + 特异性状态 |

### 坐标层级

| 概念 | 定义 | 说明 |
|------|------|------|
| **required_region** | raw CDS exons | 必须被 PCR 产物覆盖的区域；`--cds-buffer` 已废弃并被忽略 |
| **extended_target** | required_region 向基因内部扩展到 `--target-size` 下限 | Primer3 的 SEQUENCE_TARGET 区域 |
| **design_template** | extended_target + `primer-flank` 两侧 | Primer3 的 SEQUENCE_TEMPLATE |
| **SEQUENCE_TARGET** | extended_target 在 design_template 内的相对坐标 | 传给 Primer3 的 0-based 坐标 |

- required_region 的长度**不要求均一**，也不强制扩展到 target-size 中点。
- 短的 required_region 完全合法——先向基因内部扩展成 extended_target，再加 flank 生成 design_template。
- **genomic product 坐标只由 Stage 3 产生**，Stage 1/2 不输出基因组产物坐标。

### 坐标层级示例（HFE_cds1）

```
required_region:   chr6:26087410-26087546   136 bp   (CDS exon 1)
extended_target:   chr6:26087410-26090110   2700 bp  (向基因内部扩展到 target-size 下限)
design_template:   chr6:26087110-26090410   3300 bp  (extended_target + primer-flank 两侧)
SEQUENCE_TARGET:   300,2700 (0-based)  /  301,2700 (1-based, Primer3Plus 显示)
```

## 功能概述

### Stage 1 — Target Generation

- 输入基因符号列表（如 `HFE HJV TFR2`）
- 通过 Ensembl REST API 获取转录本信息
- 优先选择 MANE Select → MANE Plus Clinical → canonical → 最长 protein-coding 转录本
- **覆盖 CDS（编码序列），而非完整 exon**；UTR 区域不作为必须覆盖区域
- 使用原始 CDS 坐标作为必须覆盖区间（required interval）；`--cds-buffer` 保留为兼容参数但不生效
- 合并相邻 CDS 区间，按 `--target-size` 约束分组生成 target
- 输出 SEQUENCE_TEMPLATE（FASTA）和 SEQUENCE_TARGET（相对坐标）
- Stage 1 QC 检查区域构建，不涉及引物

### Stage 2 — Primer3Plus-like Primer Design（需 `--design-primers`）

- **模拟 Primer3Plus 默认流程**（https://www.primer3plus.com）
- 输入 Stage 1 的 SEQUENCE_TEMPLATE + SEQUENCE_TARGET
- 使用 **Primer3Plus 默认 Product Size Ranges**（501-600, 601-700, ..., 10001-20000）
- **不使用 `--target-size` 作为 Primer3 产物大小硬过滤条件**
- Primer3 参数：`PRIMER_TASK=generic`, `PRIMER_PICK_INTERNAL_OLIGO=0`
- QC 检查：验证 primer product 是否覆盖 SEQUENCE_TARGET（仅 warning，非 filter）
- 输出 `primers.tsv`（Primer3 设计指标，无基因组坐标）

### Stage 3 — In-silico PCR Specificity（需 `--check-specificity`）

- 输入 Stage 2 的 primer pair
- 使用 UCSC isPcr 做 genome-wide in-silico PCR
- **Stage 3 不使用 `--target-size` 过滤 PCR 产物**，保留所有 isPcr 命中
- **Stage 3 是唯一产生 genomic product 坐标的阶段**
- 分类逻辑：
  - `unique_pass`：单个命中，且 chrom/start/end 与预期匹配（tolerance 内）
  - `unique_off_target`：单个命中，但 chrom/start/end 不匹配预期
  - `multi_hit`：多个命中（不论其中一个是否匹配预期）
  - `no_hit`：无命中
- 输出特异性状态和 `specificity_explain`（人类可读的分类原因）

## 环境准备（WSL / Ubuntu）

```bash
# 推荐：用 micromamba 建独立环境
micromamba create -n primer_panel -c conda-forge -c bioconda \
  python=3.11 requests openpyxl pandas pyfaidx primer3 -y

micromamba activate primer_panel

# 验证
python -c "import requests, openpyxl, pandas; print('OK')"
primer3_core --version
```

依赖说明：
- `requests`（必填）：Ensembl API 通信
- `openpyxl`（可选）：XLSX 输出
- `pyfaidx`（可选）：从 genome FASTA 提取真实序列
- `primer3`（Stage 2 必填）：引物设计引擎

## 运行示例

### Stage 1 only（不需要 Primer3）

```bash
python -m primer_panel --genes HFE HJV TFR2 SLC40A1 HAMP FTH1 \
    --target-size 2700-3300 \
    --cds-buffer 30 \
    --primer-flank 300 \
    --output-dir outputs/test_targets
```

### Stage 1 + Stage 2 + Stage 3（完整流程）

```bash
python -m primer_panel --genes HFE HJV TFR2 SLC40A1 HAMP FTH1 \
    --target-size 2700-3300 \
    --cds-buffer 30 \
    --primer-flank 300 \
    --genome-fasta /path/to/hg38.fa \
    --design-primers \
    --check-specificity \
    --output-dir outputs/hcc6_primers
```

## 参数说明

### Stage 1 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--genes` | (必填) | 基因符号列表 |
| `--target-size MIN-MAX` | 2700-3300 | Stage 1 target 合并/扩展尺度（bp）。**仅用于 CDS 合并、拆分、扩展，不约束 Primer3 产物大小。** 兼容 `--product-size`。 |
| `--cds-buffer` | 0 | 已废弃并忽略；required_region 使用原始 CDS 坐标 |
| `--primer-flank` | 300 | extended_target 两侧各加的 flank（bp），用于 Primer3 搜索 |
| `--output-dir` | outputs | 输出目录 |
| `--genome-fasta` | 无 | bgzipped+索引的 hg38 FASTA |
| `-v / --verbose` | off | 调试日志 |

### Stage 2 参数（需 `--design-primers`）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--design-primers` | off | 启用 Primer3 引物设计 |
| `--primer-num-return` | 10 | 每个 target 返回的 primer pair 数量 |
| `--primer-opt-size` | 20 | 引物最优长度 |
| `--primer-min-size` | 18 | 引物最小长度 |
| `--primer-max-size` | 25 | 引物最大长度 |
| `--primer-opt-tm` | 60.0 | 引物最优 Tm |
| `--primer-min-tm` | 58.0 | 引物最小 Tm |
| `--primer-max-tm` | 64.0 | 引物最大 Tm |
| `--primer-max-tm-diff` | 2.0 | F/R 引物最大 Tm 差 |
| `--primer-min-gc` | 35.0 | 最小 GC% |
| `--primer-max-gc` | 65.0 | 最大 GC% |
| `--primer3-bin` | primer3_core | primer3_core 二进制路径 |

**注意：Stage 2 使用 Primer3Plus 默认 Product Size Ranges，不使用 `--target-size` 约束产物大小。**

### Stage 3 参数（需 `--check-specificity`）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--check-specificity` | off | 启用 in-silico PCR 特异性检查 |
| `--is-pcr-bin` | isPcr | isPcr 二进制路径 |
| `--pcr-tolerance` | 10 | 坐标匹配容差 (bp) |

## Stage 3 单独运行（从已有 primers.tsv）

如果已有 Stage 2 的 `primers.tsv` 和 Stage 1 的 `target_summary.tsv`，可以直接运行 Stage 3：

```bash
python -m primer_panel.insilico_pcr_cli \
    --primers-tsv outputs/hcc6_primers/primers.tsv \
    --targets-tsv outputs/hcc6_primers/target_summary.tsv \
    --genome-fasta /path/to/hg38.fa \
    --is-pcr-bin /path/to/isPcr \
    --output-dir outputs/stage3_from_existing
```

输出文件：`primers_specificity.tsv`、`primers_unique.tsv`、`stage3_summary.txt`。

## 输出文件

### Stage 1 输出

| 文件 | 说明 |
|------|------|
| `targets.bed` | BED6 格式（design_template 坐标） |
| `required_regions.bed` | BED6 格式（required_region 坐标，即 CDS+buffer） |
| `target_summary.tsv` | 详细 TSV 摘要（含 SEQUENCE_TARGET 坐标） |
| `target_summary.xlsx` | 同上，Excel 格式（需 openpyxl） |
| `targets.fa` | design_template FASTA 序列（SEQUENCE_TEMPLATE） |
| `failed_targets.tsv` | 处理失败的基因（仅在有失败时生成） |
| `run_summary.txt` | QC 汇总 |

### Stage 2 输出（需 `--design-primers`）

| 文件 | 说明 |
|------|------|
| `primers.tsv` | 引物设计结果 TSV（Primer3 设计指标，无基因组坐标） |
| `primers.xlsx` | Excel 格式，含 Targets + Primers 两个 sheet |

### Stage 3 输出（需 `--check-specificity`）

| 文件 | 说明 |
|------|------|
| `primers_specificity.tsv` | 含 in-silico PCR 特异性检查结果和基因组坐标 |
| `primers_unique.tsv` | 仅 unique_pass 的引物 |
| `stage3_summary.txt` | Stage 3 QC 汇总 |

### target_summary.tsv / xlsx 列说明

| 列 | 说明 |
|----|------|
| gene | 基因符号 |
| transcript_id | 选中的转录本 ID |
| selection_reason | 转录本选择原因 |
| target_id | 靶标 ID，如 `HFE_cds1_3` 表示覆盖 CDS exon 1-3 |
| required_chrom / required_start / required_end / required_length | required_region 坐标 |
| extended_chrom / extended_start / extended_end / extended_length | extended_target 坐标 |
| template_chrom / template_start / template_end / template_length | design_template 坐标 |
| sequence_target_start_0based | SEQUENCE_TARGET 在 template 内的 0-based 起始 |
| sequence_target_length | SEQUENCE_TARGET 长度 |
| sequence_target_for_primer3plus_1based | Primer3Plus 1-based 显示格式 |
| strand | 链方向 (+/-) |
| product_min / product_max | Stage 1 target 扩展尺度约束 |
| cds_exon_numbers | 覆盖的 CDS exon 编号（如 `1,2,3`） |
| cds_exon_ids | 对应的 Ensembl exon ID |
| covered_cds_count | 覆盖的 CDS exon 数量 |
| cds_exon_coords | CDS exon 基因组坐标 |
| status | 靶标状态 |
| needs_review | 是否需要人工审查 |
| sequence_status | 序列来源（placeholder / real） |
| target_qc_status | Stage 1 QC 状态（ok 或问题描述） |

### primers.tsv 列说明（Stage 2 输出）

| 列 | 说明 |
|----|------|
| target_id | 靶标 ID |
| primer_rank | 引物对排名（1=最佳） |
| forward_primer / reverse_primer | 引物序列 |
| forward_tm / reverse_tm / tm_diff | Tm 值 |
| forward_gc / reverse_gc | GC% |
| primer_pair_penalty | Primer3 penalty 值 |
| primer_left_start / primer_left_len | 正向引物在 template 上的位置（0-based） |
| primer_right_start / primer_right_len | 反向引物在 template 上的位置（0-based） |
| primer3_product_size | Primer3 报告的产物大小（template-relative） |
| primer3_status | ok / no_primer |
| primer3_explain | Primer3 解释或 QC 警告 |
| sequence_target_start_0based | SEQUENCE_TARGET 起始 |
| sequence_target_length | SEQUENCE_TARGET 长度 |

**注意：primers.tsv 不包含基因组产物坐标。基因组坐标由 Stage 3 产生。**

### primers_specificity.tsv 列说明（Stage 3 输出）

在 Stage 2 基础上增加：

| 列 | 说明 |
|----|------|
| insilico_status | unique_pass / no_hit / multi_hit / unique_off_target / pcr_error |
| insilico_hit_count | isPcr 命中数 |
| insilico_hits | 所有命中坐标 |
| insilico_best_chrom / insilico_best_start / insilico_best_end / insilico_best_size | 最佳命中基因组坐标 |
| specificity_pass | 是否通过特异性检查 |
| expected_target_chrom / expected_target_start / expected_target_end | 预期靶标区域（来自 Stage 1） |
| specificity_explain | 特异性检查说明 |

### 靶标生成规则

1. 从 Ensembl 提取 CDS/Translation 区域，使用原始 CDS 坐标作为 required_region（`--cds-buffer` 已废弃并忽略）
2. 相邻/重叠的 required_region 合并
3. 合并后区间按 `--target-size` 约束分组：
   - 如果多个 CDS 区间合并后 required_span ≤ `product_max`，则合并以减少扩增次数
   - 如果合并后超过 `product_max`，在 CDS 间最大 gap 处拆分
   - 如果单个区间超过 `product_max`，自动 tile 拆分（target 间重叠 200 bp）
4. 短的 required_region 向基因内部扩展到 `--target-size` 下限 → extended_target
5. extended_target 两侧加 `primer-flank` → design_template
6. 计算 SEQUENCE_TARGET = extended_target 在 design_template 内的相对坐标

## 当前阶段限制

- **Ensembl API 依赖**: 需要网络连接；离线模式不可用
- **占位序列**: 默认输出 N 填充序列；需 `--genome-fasta` + pyfaidx 获取真实序列
- **FASTA 链方向**: 当前按 hg38 正链输出；负链靶标未做 reverse-complement
- **Stage 2 需要真实序列**: N 占位序列会阻止 Primer3 运行，需提供 `--genome-fasta`
- **单物种**: 仅支持 homo_sapiens (hg38)

## 下一阶段计划

- Stage 4: 多基因批量优化与报告

## 项目结构

```
primer_panel/
├── __init__.py          # 版本标记
├── __main__.py          # python -m 入口
├── main.py              # CLI 主流程（Stage 1 + 2 + 3）
├── config.py            # 参数配置
├── ensembl_client.py    # Ensembl REST API 客户端 + CDS 提取
├── cds_handler.py       # CDS 必须覆盖区间处理
├── exon_handler.py      # 外显子区间处理（遗留）
├── target_planner_adapter.py  # primer-target-planner 适配层
├── primer3_runner.py    # Primer3 Boulder IO 运行器（Stage 2, Primer3Plus-like）
├── insilico_pcr.py      # UCSC isPcr 封装（Stage 3）
├── insilico_pcr_cli.py  # Stage 3 独立 CLI（从已有 primers.tsv 运行）
├── writers.py           # 输出格式（BED/TSV/XLSX/FASTA/primers/specificity）
└── utils.py             # 工具函数
```
