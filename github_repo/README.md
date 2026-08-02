# Network-CMIN: Biomarker Network Dysregulation and Cardiovascular Mortality

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)

Analysis code and figures for:

> **Mahalanobis distance-based biomarker network dysregulation and residual cardiovascular mortality: a binational cohort study**
>
> Gong Q, Fan Y, Shao S, Li X, Du Y, Guo Q. *BMC Medicine* (submitted).

## Overview

Network-CMIN is a Mahalanobis distance-based score that quantifies multi-system biomarker network dysregulation from 7–10 routine clinical biomarkers (CRP, glucose, HbA1c, total cholesterol, HDL-C, LDL-C, triglycerides, SBP, DBP, BMI). Higher scores indicate greater deviation of an individual's biomarker profile from the tightly coordinated configuration observed in healthy middle-aged adults.

## Key Findings

| Analysis | CHARLS (N=12,436) | NHANES (N=17,804) | NHANES Wave I (N=2,349) |
|----------|:---:|:---:|:---:|
| All-cause mortality HR (per SD) | 1.14 (1.10–1.18) | 1.13 (1.09–1.18) | 1.27 (1.04–1.54) |
| CVD mortality HR (per SD) | — | 1.25 (1.17–1.34) | 1.57 (1.03–2.41) |
| Discordant vs Concordant-Low HR | 1.86 (1.35–2.57) | — | 1.71 (0.88–3.32)\* |
| NRI (category-free) | 0.195 (P < 0.001) | — | — |

\*Underpowered (52% power for HR=1.86; 42 combined events). P = 0.11.

## Repository Structure

```
├── manuscript_bmc.tex          # LaTeX manuscript source
├── references_verified.bib     # Verified reference list (26 entries)
├── cover_letter.tex            # Submission cover letter
├── requirements.txt            # Python dependencies
├── README.md
├── scripts/
│   ├── phase0_load.py          # Data loading: CHARLS 2015 + NHANES 1999–2016
│   ├── phase1_7_charls_analysis.py  # CHARLS: DM → Cox → NRI → RCS → DCA
│   ├── nhanes_validation.py    # NHANES external validation
│   ├── nhanes_framingham_h2.py # NHANES H2 with Framingham 2008 risk categories
│   ├── nhanes_waveI_full.py    # NHANES Wave I CRP-inclusive sub-analysis
│   ├── plot_figures.py         # Figure generation (matplotlib)
│   └── gen_fig1.py             # STROBE flow diagram (SVG)
├── figures/
│   ├── fig1_strobe.svg         # Study flow diagram
│   ├── fig2_km_discordance.svg # KM curves by risk group
│   ├── fig3_forest.svg         # Forest plot across cohorts
│   ├── fig4_doseresponse.svg   # Dose-response quartiles
│   ├── efigS1_nhanes_km.svg    # NHANES KM curves (supplementary)
│   └── efigS2_dm_distribution.svg  # DM Z-score distributions
└── output/
    └── table1.csv              # Baseline characteristics by DM tertile
```

## Data Availability

This study uses publicly available, de-identified secondary data:

- **CHARLS 2015**: [http://charls.pku.edu.cn/](http://charls.pku.edu.cn/) (registration required)
- **NHANES 1999–2016**: [https://wwwn.cdc.gov/nchs/nhanes/](https://wwwn.cdc.gov/nchs/nhanes/)
- **NHANES Linked Mortality Files**: [https://www.cdc.gov/nchs/data-linkage/mortality-public.htm](https://www.cdc.gov/nchs/data-linkage/mortality-public.htm)

Raw data files are NOT included in this repository. Users must download them from the original sources. The `phase0_load.py` script automatically reads files from the expected directory structure.

## Reproducing the Analysis

### Requirements

Python 3.12 with the following packages:

```
pandas>=3.0
numpy>=2.5
scipy>=1.14
lifelines>=0.30
matplotlib>=3.9
```

### Data Setup

1. Download CHARLS 2015 data to `D:/CHARLS/CHARLS_2015/`
2. Download NHANES 1999–2016 `.xpt` files to `D:/NHANES/`
3. Download NHANES linked mortality `.dat` files to `D:/NHANES_mortality/`
4. Update paths in `scripts/phase0_load.py` if needed

### Run

```bash
# Step 0: Load and merge data
python scripts/phase0_load.py

# Step 1-7: Full CHARLS analysis
python scripts/phase1_7_charls_analysis.py

# NHANES validation
python scripts/nhanes_validation.py
python scripts/nhanes_framingham_h2.py
python scripts/nhanes_waveI_full.py

# Figures
python scripts/plot_figures.py
python scripts/gen_fig1.py
```

Output files are saved to `./output/`.

## License

MIT License. See [LICENSE](LICENSE) file.

## Citation

If you use this code or data, please cite:

```
Gong Q, Fan Y, Shao S, Li X, Du Y, Guo Q. Network-based multi-system
biomarker dysregulation identifies residual cardiovascular mortality
risk beyond traditional scores: a binational cohort study. BMC Medicine.
2026.
```
