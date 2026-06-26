# Target Planning Audit: v1 vs v2 Grouper

## 1. Coordinate Summary (v2)

### HFE (chr6, +strand)

| target_id | required | extended | cds_exons |
|-----------|----------|----------|-----------|
| HFE_cds1 | 26087410-26087546 (136bp) | 26087410-26090110 (2700bp) | 1 |
| HFE_cds2_5 | 26090810-26093262 (2452bp) | 26090562-26093262 (2700bp) | 2,3,4,5 |
| HFE_cds4_6 | 26092654-26094256 (1602bp) | 26091556-26094256 (2700bp) | 4,5,6 |

### TFR2 (chr7, −strand)

| target_id | required | extended | cds_exons |
|-----------|----------|----------|-----------|
| TFR2_cds18 | 100620826-100621156 (330bp) | 100620826-100623526 (2700bp) | 18 |
| TFR2_cds10_17 | 100626732-100629402 (2670bp) | 100626717-100629417 (2700bp) | 10-17 |
| TFR2_cds5_9 | 100630858-100633370 (2512bp) | 100630764-100633464 (2700bp) | 5,6,7,8,9 |
| TFR2_cds4 | 100633385-100633586 (201bp) | 100632136-100634836 (2700bp) | 4 |
| TFR2_cds1_3 | 100640655-100641539 (884bp) | 100638839-100641539 (2700bp) | 1,2,3 |

### FTH1 (chr11, −strand)

| target_id | required | extended | cds_exons |
|-----------|----------|----------|-----------|
| FTH1_cds2_4 | 61964696-61965545 (849bp) | 61964696-61967396 (2700bp) | 2,3,4 |
| FTH1_cds1_3 | 61964956-61967455 (2499bp) | 61964856-61967556 (2700bp) | 1,2,3 |

### SLC40A1 / HAMP / HJV — identical in v1 and v2

---

## 2. HFE: v1 cds6 → v2 cds4_6

### v1
- HFE_cds6: required=26094155-26094256 (101bp), extended=26091556-26094256 (2700bp), cds=6

### v2
- HFE_cds4_6: required=26092654-26094256 (1602bp), extended=26091556-26094256 (2700bp), cds=4,5,6

### Analysis

CDS exon 4 (26092684-26092960) and exon 5 (26093118-26093232) are already covered by HFE_cds2_5 (extended to 26093262). In v1, these exons were covered by cds2_5, and cds6 was a separate target.

In v2, the backward sweep from cds6 (the last uncovered region) extends leftward by product_min=2700 → absorbs cds4 and cds5. This is the "terminal sweep" strategy working correctly.

**Result**: cds4 and cds5 appear in BOTH HFE_cds2_5 and HFE_cds4_6. This is redundant coverage but not harmful — each target's extended_target independently covers its required_regions.

**Verdict**: Acceptable. The v2 strategy naturally produces overlapping targets at gene interior via the terminal sweep. v1's split (cds2_5 + cds6) was slightly more efficient (no overlap), but v2's split is also valid. **No algorithm change needed.**

---

## 3. TFR2: v1 cds4_9 → v2 cds5_9 + cds4

### v1
- TFR2_cds4_9: required=100630858-100633586 (2728bp), extended=100630858-100633586 (2728bp), cds=4-9

### v2
- TFR2_cds5_9: required=100630858-100633370 (2512bp), extended=100630764-100633464 (2700bp), cds=5-9
- TFR2_cds4: required=100633385-100633586 (201bp), extended=100632136-100634836 (2700bp), cds=4

### Root Cause: Forward Window Anchor Bug

The v2 grouper's forward window anchors from `regions[start_idx].start` instead of `regions[last_grouped].end`. For TFR2:

The cds_handler produces 18 separate required intervals (one per CDS exon with buffer, gaps > 60bp). The v2 grouper starts forward sweep from cds9 (start=100630858):

- product_min window_end = 100630858 + 2700 = 100633558
- Covers cds8 (end=100631975 ✓), cds7 (end=100632228 ✓), cds6 (end=100633153 ✓), cds5 (end=100633370 ✓)
- cds4.end = 100633586 > 100633558 → **NOT covered** (missed by 28bp!)
- product_max window_end = 100630858 + 3300 = 100634158
- cds4.end = 100633586 ≤ 100634158 → covered by product_max

So the grouper uses product_max and should group cds5-9 + cds4. But the actual output shows them as separate targets. The issue is that `_count_covered_forward` iterates from `start_idx` in order and breaks at the first non-matching region. cds5 (region[13]) has end=100633370 ≤ 100633558, but cds4 (region[14]) has end=100633586 > 100633558. The function returns count_min=4 (cds8-5), then checks product_max which covers cds4 → count_max=5.

**However**, the actual pipeline output shows cds5_9 and cds4 as separate targets. This means the grouper is NOT using product_max to group them. The likely reason: the grouper processes cds9-5 as one group (product_min covers them), then cds4 becomes the next `i`. Since cds4 is the last region, it goes through the backward sweep path, not the forward path.

**The fix**: anchor the forward window from `regions[last_grouped_idx].end` instead of `regions[start_idx].start`. After grouping cds9-5, the window should start from cds5.end (100633370), giving:
- product_min window_end = 100633370 + 2700 = 100636070
- cds4.end = 100633586 ≤ 100636070 → **covered!**

### Specificity Impact

| Target | v1 unique_pass | v2 unique_pass |
|--------|---------------|----------------|
| TFR2_cds4_9 | 9/10 | — |
| TFR2_cds5_9 | — | 9/10 |
| TFR2_cds4 | — | 8/10 |

Both v2 sub-targets have good unique primer coverage. **No functional problem, but unnecessary split.**

### Verdict

The split is a **bug consequence**, not a strategy choice. The 15bp gap between cds5.end and cds4.start means they should be grouped. **Recommend fixing the grouper's window anchor.**

---

## 4. FTH1: v1 cds1_4 → v2 cds2_4 + cds1_3

### v1
- FTH1_cds1_4: required=61964696-61967455 (2759bp), extended=61964696-61967455 (2759bp), cds=1,2,3,4

### v2
- FTH1_cds2_4: required=61964696-61965545 (849bp), extended=61964696-61967396 (2700bp), cds=2,3,4
- FTH1_cds1_3: required=61964956-61967455 (2499bp), extended=61964856-61967556 (2700bp), cds=1,2,3

### Analysis

FTH1 CDS exons (−strand, numbered from 3' end):
- cds4: 61964726-61964891 (genomically leftmost)
- cds3: 61964986-61965112
- cds2: 61965368-61965515
- cds1: 61967311-61967425 (genomically rightmost)

Required intervals (after +30bp buffer, cds_handler produces separate intervals since gaps > 60bp):
- cds4+buffer: 61964696-61964921
- cds3+buffer: 61964956-61965142 (gap=35bp from cds4 → merged!)
- cds2+buffer: 61965338-61965545
- cds1+buffer: 61967281-61967455

Wait — cds4.end=61964921, cds3.start=61964956. Gap=35bp < 60bp → **merged** into 61964696-61965142.

So cds_handler produces: [61964696-61965142, 61965338-61965545, 61967281-61967455]

v2 grouper forward sweep from first interval (61964696):
- product_min window_end = 61964696 + 2700 = 61967396
- Covers second interval (end=61965545 ✓)
- Third interval end=61967455 > 61967396 → NOT covered by product_min
- product_max window_end = 61964696 + 3300 = 61967996
- Third interval end=61967455 ≤ 61967996 → covered by product_max

So product_min covers cds4+cds3 + cds2, but not cds1. Product_max covers all three. Since product_min covers at least one more (cds2), the grouper uses product_min → groups cds4+cds3 + cds2 → target cds2_4.

Then cds1 is the last region → backward sweep from cds1.end (61967455) - 2700 = 61964755. This covers cds2.start=61965338 ✓ and cds3.start=61964956 ✓, but NOT cds4.start=61964696 (61964696 < 61964755 by 59bp). So backward absorbs cds2 and cds3 but not cds4 → target cds1_3.

**Result**: cds2 and cds3 appear in BOTH targets (redundant coverage). cds4 is only in cds2_4. cds1 is only in cds1_3.

### Specificity Impact

| Target | unique_pass/10 | multi_hit/10 |
|--------|---------------|-------------|
| v1 FTH1_cds1_4 | 4 | 6 |
| v2 FTH1_cds2_4 | 6 | 4 |
| v2 FTH1_cds1_3 | 6 | 4 |

v2 sub-targets have **better unique primer ratios** (6/10 vs 4/10). Splitting FTH1 improves primer design options — shorter templates give Primer3 more flexibility to pick primers in unique flanking regions.

### Verdict

**v2 split is beneficial.** Each sub-target has more unique_pass primers than the single v1 target. The redundant coverage of cds2,3 is acceptable — it provides primer design flexibility. **Recommend keeping v2's 2-target split for FTH1.**

---

## 5. CDS Coverage Completeness

| Gene | CDS Exons | v1 Coverage | v2 Coverage |
|------|-----------|-------------|-------------|
| HFE | 1,2,3,4,5,6 | all ✓ | all ✓ |
| TFR2 | 1-18 | all ✓ | all ✓ |
| FTH1 | 1,2,3,4 | all ✓ | all ✓ |
| SLC40A1 | 1-8 | all ✓ | all ✓ |
| HAMP | 1,2,3 | all ✓ | all ✓ |
| HJV | 1,2,3 | all ✓ | all ✓ |

**All CDS exons are fully covered in both v1 and v2.**

---

## 6. Unnecessary Redundancy Assessment

| Overlap | v1 | v2 | Harmful? |
|---------|----|----|----------|
| HFE cds4,5 in both cds2_5 and cds4_6 | No overlap | Yes, ~600bp overlap | No — provides primer flexibility |
| FTH1 cds2,3 in both cds2_4 and cds1_3 | No overlap | Yes, ~600bp overlap | No — improves unique primer count |

**No harmful redundancy.**

---

## 7. Strategy Compliance

| Rule | Compliance | Notes |
|------|-----------|-------|
| Forward: try product_min first | ✓ | All targets |
| Forward: fall back to product_max | Partial | TFR2 cds4 should be absorbed but isn't (anchor bug) |
| Terminal: backward sweep for last region | ✓ | FTH1_cds1_3, TFR2_cds4 |
| Single region > product_max → tile | ✓ | No instances in HCC6 |
| chr0 clamp → compensate rightward | ✓ | No instances in HCC6 |

---

## 8. Recommendations

1. **Fix the forward window anchor bug.** Anchor from `regions[last_grouped_idx].end` instead of `regions[start_idx].start`. This will correctly group TFR2 cds5_9 + cds4.

2. **Keep v2 as default** after the anchor fix. The sliding-window strategy is fundamentally correct.

3. **Keep v1 fallback** for now. Can be removed after the anchor fix is validated.

4. **FTH1**: Keep 2 targets (cds2_4 + cds1_3). The split improves unique primer ratios.

5. **No other algorithm changes needed.** The min-first/max-rescue/terminal-sweep strategy is sound.
