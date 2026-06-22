## =====================================================================
## 00_config.R  —  shared settings for all analysis scripts
## Study: Type 2 diabetes rather than bone mineral density underlies the
##        osteoporosis–PAD association: a multivariable MR study
## =====================================================================

## ---- packages ----
## install once (see README):
##   remotes::install_github("MRCIEU/TwoSampleMR")
##   remotes::install_github("MRCIEU/ieugwasr")
##   remotes::install_github("WSpiller/MVMR")
##   remotes::install_github("rondolab/MR-PRESSO")
suppressMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
})

## ---- OpenGWAS access token ----
## Get a free token at https://api.opengwas.io and put it in ~/.Renviron as:
##   OPENGWAS_JWT=your_token_here
## NEVER hard-code the token in these scripts (it is a secret).
## ieugwasr/TwoSampleMR pick it up automatically from the environment.

## ---- GWAS data sources (OpenGWAS IDs) ----
## These match Supplementary Table S1. Verify/update at https://gwas.mrcieu.ac.uk
EBMD    <- "ebi-a-GCST006979"   # Heel eBMD (exposure), Morris 2019, UK Biobank
PAD_R5  <- "finn-b-I9_PAD"      # PAD (primary outcome), FinnGen R5
T2D     <- "ebi-a-GCST006867"   # Type 2 diabetes, Xue 2018
BMI     <- "ieu-a-2"            # Body mass index, GIANT (Locke 2015)
SMOKING <- "ieu-b-4877"         # Smoking initiation, GSCAN (Liu 2019)
CAD     <- "ieu-a-7"            # Coronary artery disease, CARDIoGRAMplusC4D (Nikpay 2015)
STROKE  <- "ebi-a-GCST005843"   # Ischemic stroke, MEGASTROKE (Malik 2018)
DXA_FN  <- "ieu-a-980"          # Femoral-neck BMD (DXA), GEFOS (Zheng 2015)
DXA_LS  <- "ieu-a-982"          # Lumbar-spine BMD (DXA), GEFOS
DXA_FA  <- "ieu-a-977"          # Forearm BMD (DXA), GEFOS

## FinnGen R12 peripheral atherosclerosis (sensitivity) is NOT on OpenGWAS;
## download from https://www.finngen.fi/en/access_results and place the file in
## data_raw/  (see 05_finngen_r12_sensitivity.R).

## ---- clumping / instrument settings (as reported in the paper) ----
P_THRESH  <- 5e-8     # genome-wide significance
CLUMP_R2  <- 0.001    # LD clumping r^2
CLUMP_KB  <- 10000    # LD clumping window

## ---- output dir ----
OUTDIR <- "results"
dir.create(OUTDIR, showWarnings = FALSE)

cat("Config loaded. Output ->", normalizePath(OUTDIR), "\n")
