## Run the full Paper-2 MR analysis pipeline in order.
## Requires: OPENGWAS_JWT in ~/.Renviron, and (for script 05) the FinnGen R12
## summary file in data_raw/. See README.md.
source("R/00_config.R")
source("R/01_main_mr_eBMD_PAD.R")
source("R/02_mvmr.R")
source("R/03_phenotype_specificity.R")
source("R/04_dxa_sensitivity.R")
## source("R/05_finngen_r12_sensitivity.R")   # needs the FinnGen R12 download
source("R/06_figures.R")
cat("\nAll done. Outputs are in results/\n")
