## =====================================================================
## 02_mvmr.R
## Multivariable MR: direct effects of eBMD, T2D, BMI, smoking on PAD
## Reproduces: T2D OR ~1.26 (1.16-1.36), p~7e-9; eBMD null (OR ~1.03);
##             BMI OR ~1.30 (p~0.06); smoking OR ~1.25 (p~0.12);
##             conditional F: eBMD ~133, T2D ~16 (gencov = 0 caveat)
## =====================================================================
source("R/00_config.R")

## 1. Assemble instruments across the four exposures
mv_exp <- mv_extract_exposures(c(EBMD, T2D, BMI, SMOKING))

## 2. Outcome (PAD R5) for the combined instrument set
mv_out <- extract_outcome_data(snps = mv_exp$SNP, outcomes = PAD_R5)

## 3. Harmonise
mvdat <- mv_harmonise_data(mv_exp, mv_out)

## 4. Multivariable IVW (direct effects)
mv_res <- mv_multiple(mvdat)
mv_or  <- generate_odds_ratios(mv_res$result)
print(mv_or[, c("exposure", "nsnp", "or", "or_lci95", "or_uci95", "pval")])
write.csv(mv_or, file.path(OUTDIR, "mvmr_results.csv"), row.names = FALSE)

## 5. Conditional instrument strength + MVMR heterogeneity (MVMR package)
## NOTE: gencov = 0 assumes no covariance between exposure effects; because the
## exposure GWAS share UK Biobank participants, conditional F may be modestly
## optimistic (discussed in the paper's limitations).
if (requireNamespace("MVMR", quietly = TRUE)) {
  bx <- mvdat$exposure_beta            # SNP x exposure matrix
  sx <- mvdat$exposure_se
  by <- mvdat$outcome_beta
  sy <- mvdat$outcome_se
  mvmr_in <- MVMR::format_mvmr(BXGs = bx, BYG = by,
                               seBXGs = sx, seBYG = sy,
                               RSID = rownames(bx))
  cat("\nConditional F-statistics (gencov = 0):\n")
  print(MVMR::strength_mvmr(mvmr_in, gencov = 0))
  cat("\nMVMR Cochran's Q (horizontal pleiotropy):\n")
  print(MVMR::pleiotropy_mvmr(mvmr_in, gencov = 0))
}
cat("02_mvmr done.\n")
