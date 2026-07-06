# Data

## `toothgrowth.csv` — Application 1
Vitamin C tooth growth data (Crampton, 1947; Bliss, 1952). Guinea pig odontoblast lengths, 6 groups (n=10 each): 2 delivery methods (VC, OJ) × 3 doses (0.5, 1.0, 2.0 mg).

Columns: `VC_0.5mg`, `OJ_0.5mg` (controls, P1); `VC_1.0mg`, `VC_2.0mg`, `OJ_1.0mg`, `OJ_2.0mg` (treatments, P2).

Note: similar to R's built-in `ToothGrowth` but not identical — use this file to reproduce the paper's results.

Used by: `applications/analysis1_toothgrowth.R`

## `cholesterol_reduction.csv` — Application 2
Cholesterol reduction data (Braat et al., 2008; Table 1, Westfall et al., 1999). 5 groups (n=10 each): 3 test-drug regimens + 2 active comparators.

Columns: `1-time`, `2-times`, `4-times` (test drug doses, P1); `drug D`, `drug E` (comparators, P2).

Used by: `applications/analysis2_cholesterol.R`

---
Both datasets are transcribed exactly from the published sources and match the values hardcoded in their respective analysis scripts.
