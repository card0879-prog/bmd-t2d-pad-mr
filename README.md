# Bone mineral density, type 2 diabetes and peripheral artery disease — Mendelian randomization

Reproducible analysis code for:

> **Type 2 diabetes rather than bone mineral density underlies the osteoporosis–peripheral artery disease association: a multivariable Mendelian randomization study.**
> Zhang H†, Wang Z†, Guo W, Ma Y, Tao J, Wang L, Wang J, Luo C, Xu K, Ning K, Wang L\*, Yang Y\*. *Cardiovascular Diabetology* (under review), 2026.
> († joint first authors; \* joint corresponding authors)

This repository reproduces every Mendelian randomization (MR) result in the paper:
the primary two-sample MR (eBMD → PAD), multivariable MR (MVMR), phenotype-specificity
controls, DXA site-specific sensitivity, and the FinnGen R12 endpoint-definition
sensitivity analysis.

---

## Data availability — fully public

**This study uses only publicly available GWAS summary statistics — no individual-level
data are used or required.** All exposures and outcomes are pulled programmatically from
the [OpenGWAS database](https://gwas.mrcieu.ac.uk/) (except the FinnGen R12 sensitivity
endpoint, which is downloaded from FinnGen). The analysis is therefore reproducible by
anyone with a free OpenGWAS token.

| Trait | Role | Source | Accession |
|---|---|---|---|
| Heel eBMD | exposure (primary) | UK Biobank (Morris 2019) | `ebi-a-GCST006979` |
| PAD | outcome (primary) | FinnGen R5 | `finn-b-I9_PAD` |
| Peripheral atherosclerosis | outcome (sensitivity) | FinnGen R12 | `finngen_R12_I9_ATHSCLE` * |
| Body mass index | MVMR exposure | GIANT (Locke 2015) | `ieu-a-2` |
| Type 2 diabetes | MVMR exposure | Xue 2018 | `ebi-a-GCST006867` |
| Smoking initiation | MVMR exposure | GSCAN (Liu 2019) | `ieu-b-4877` |
| Coronary artery disease | control outcome | CARDIoGRAMplusC4D (Nikpay 2015) | `ieu-a-7` |
| Ischemic stroke | control outcome | MEGASTROKE (Malik 2018) | `ebi-a-GCST005843` |
| Femoral-neck / lumbar / forearm BMD (DXA) | sensitivity exposure | GEFOS (Zheng 2015) | `ieu-a-980` / `982` / `977` |

\* FinnGen R12 is not on OpenGWAS — download `I9_ATHSCLE` (release R12) from
https://www.finngen.fi/en/access_results and place it at `data_raw/finngen_R12_I9_ATHSCLE.gz`.

> The companion **observational** analysis (referenced in the paper) uses individual-level
> UK Biobank data, which **cannot** be shared and is available only via UK Biobank
> application (https://www.ukbiobank.ac.uk). That data and code are not part of this repository.

---

## Repository structure

```
.
├── R/
│   ├── 00_config.R                  # packages, OpenGWAS token note, accessions, settings
│   ├── 01_main_mr_eBMD_PAD.R        # primary two-sample MR + sensitivity (Q, Egger, MR-PRESSO)
│   ├── 02_mvmr.R                    # MVMR (eBMD + T2D + BMI + smoking) + conditional F
│   ├── 03_phenotype_specificity.R  # eBMD -> CAD / ischemic stroke control outcomes
│   ├── 04_dxa_sensitivity.R        # femoral-neck / lumbar / forearm DXA BMD -> PAD
│   ├── 05_finngen_r12_sensitivity.R# eBMD -> FinnGen R12 peripheral atherosclerosis
│   └── 06_figures.R                # scatter / funnel / leave-one-out / forest plots (R)
├── run_all.R                        # runs the pipeline in order
├── results/                         # created on run: estimates + per-SNP CSVs + plots
├── data_raw/                        # (you create) holds the FinnGen R12 download only
├── README.md  ·  LICENSE  ·  CITATION.cff  ·  .gitignore
```

---

## How to reproduce

**1. Install R (≥ 4.3; developed on 4.6.0) and the packages:**

```r
install.packages(c("remotes", "data.table", "ggplot2"))
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/ieugwasr")
remotes::install_github("WSpiller/MVMR")
remotes::install_github("rondolab/MR-PRESSO")
```

**2. Add your OpenGWAS token** (free, from https://api.opengwas.io) to `~/.Renviron`:

```
OPENGWAS_JWT=your_token_here
```

Restart R. **Do not put the token in any script.**

**3. Run** (from the repository root, so that `R/` and `results/` resolve correctly):

```r
source("run_all.R")
```

Script `05` additionally needs the FinnGen R12 file in `data_raw/` (see above); it is
commented out in `run_all.R` and can be run on its own once the file is in place.

Results (effect estimates, per-SNP tables, and figures) are written to `results/`.
Exact numbers may differ marginally from the paper if OpenGWAS updates a dataset version.

---

## Citation

Please cite the paper above. An archived, citable version of this repository is on
Zenodo: **DOI: pending** (a DOI will be added here after the Zenodo release is created).

## License

Code released under the MIT License (`LICENSE`). The license covers the code only; the
underlying GWAS summary statistics remain under their original providers' terms.
