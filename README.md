# Bone mineral density, type 2 diabetes, and peripheral artery disease: a Mendelian randomization study

Analysis code for a two-sample and multivariable Mendelian randomization (MR) study testing whether genetically predicted bone mineral density (BMD) causally affects peripheral artery disease (PAD), and whether the observational osteoporosis–PAD association reflects shared cardiometabolic risk (notably type 2 diabetes).

**Status:** Submitted to *BMC Medical Genomics* (2026).

## Authors
Wenjuan Guo, Li Wang, Zhijun Wang, Yinxiang Ma, Zhengyi Cao, Jie Wang, Cheng Luo, Kailang Xu, Kai Ning, Haoshaqiang Zhang.
Corresponding author: Haoshaqiang Zhang (zhanghaoshaqiang@126.com).

## Summary
- **Exposure:** heel estimated BMD (eBMD; UK Biobank, Morris 2019).
- **Primary outcome:** PAD (FinnGen Data Freeze 5); sensitivity endpoint: peripheral atherosclerosis (FinnGen Data Freeze 12).
- **Methods:** random-effects IVW (primary estimator), MR-Egger, weighted median and mode, MR-PRESSO, leave-one-out; multivariable MR (BMI, type 2 diabetes, smoking initiation); phenotype-specificity controls (coronary artery disease, ischemic stroke); site-specific DXA-derived BMD; and an instrument-level screen for cardiometabolic confounders.
- **Main result:** no evidence of a causal effect of genetically predicted eBMD on PAD (IVW OR 1.03, 95% CI 0.96–1.11); type 2 diabetes was the strongest identifiable causal signal (OR 1.26, p = 7.4×10⁻⁹).

## Data sources
All analyses use publicly available GWAS summary statistics obtained through the OpenGWAS database (https://gwas.mrcieu.ac.uk/) and FinnGen (https://www.finngen.fi/en/access_results). No individual-level data were used. GWAS accessions and sample sizes are listed in the manuscript (Table 1 and Additional file 1 / Table S1).

## Repository contents
- `R/` — analysis scripts (`00_config.R` … `07_supplementary_analyses.R`), including `PLOSONE_supplementary.R`, which runs the statistical-precision and cardiometabolic-confounder-screen analyses.
- `results/harmonised_eBMD_PAD_R5.csv` — harmonised per-variant dataset for the primary eBMD→PAD analysis (493 instruments: rsID, effect/other alleles, betas, effect-allele frequencies, and the mr_keep flag).
- `results/PLOSONE_supplementary_summary.txt` — output of the confounder screen and precision analysis: the 52 instruments associated with a cardiometabolic confounder (p < 5×10⁻⁸) and the IVW estimate before (OR 1.03) and after (OR 1.02) their exclusion.

## Requirements
R (≥ 4.0) with `TwoSampleMR`, `ieugwasr`, `MVMR`, `MRPRESSO`, `phenoscanner`, and `data.table`. An OpenGWAS API token is required for live data retrieval (obtain one at https://api.opengwas.io).

## Reproducing the analyses
1. Set your OpenGWAS token as `OPENGWAS_JWT` in `~/.Renviron` (or at the top of the script).
2. Run the scripts in `R/` in order; outputs are written to `results/`.

## License
Released for academic, non-commercial use.
