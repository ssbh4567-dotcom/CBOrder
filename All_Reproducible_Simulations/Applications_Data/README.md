# Applications Data

This folder contains the datasets used in the two real-data applications reported in the paper. Both are fully reported in their original published sources; the CSV files here transcribe those published values exactly, so the scripts in `applications/` (or `R/`) can be run end-to-end without needing to consult the original sources.

---

## `Application1_toothgrowth.csv`

**Source:** Crampton, E. W. (1947), *Journal of Nutrition*, 33(5), 491–504; data later analyzed by Bliss, C. I. (1952).

**Description:** Odontoblast (tooth growth) lengths in guinea pigs given Vitamin C via two delivery methods — ascorbic acid (VC) and orange juice (OJ) — at three dose levels (0.5, 1.0, 2.0 mg/day). Six groups, n = 10 each.

**Columns:**
- `VC_0.5mg`, `OJ_0.5mg` — control groups (partition **P1**, 0.5 mg dose)
- `VC_1.0mg`, `VC_2.0mg`, `OJ_1.0mg`, `OJ_2.0mg` — treatment groups (partition **P2**, higher doses)

**Note:** Values are close to but not identical to R's built-in `ToothGrowth` dataset — use this file, not `data(ToothGrowth)`, to reproduce the paper's results.

---

## `Application2_cholesterol_reduction.csv`

**Source:** Braat, S. et al. (2008); data reproduced from Table 1 (example dataset, Westfall et al., 1999).

**Description:** Cholesterol reduction under three test-drug dose regimens and two active comparator drugs. Five groups, n = 10 each.

**Columns:**
- `1-time`, `2-times`, `4-times` — test drug regimens, 20 mg once/10 mg twice/5 mg four-times daily (partition **P1**)
- `drug D`, `drug E` — active comparators (partition **P2**)

---

## Reproducibility note

Values in both files are identical to those hardcoded in the corresponding analysis scripts, so the data and code stay in sync as a single, citable source of the exact inputs used to generate the paper's results.
