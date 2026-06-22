## =====================================================================
## 01_main_mr_eBMD_PAD.R
## Primary analysis: genetically predicted heel eBMD -> PAD (FinnGen R5)
## Reproduces: IVW OR ~1.03 (0.96-1.11), p~0.43; 493 variants; mean F~167;
##             Egger intercept p~0.21; Cochran's Q~608; MR-PRESSO corrected p~0.48
## =====================================================================
source("R/00_config.R")

## 1. Instruments for eBMD (genome-wide significant, LD-clumped)
exp <- extract_instruments(EBMD, p1 = P_THRESH, clump = TRUE,
                           r2 = CLUMP_R2, kb = CLUMP_KB)

## 2. Outcome: PAD (FinnGen R5)
out <- extract_outcome_data(snps = exp$SNP, outcomes = PAD_R5)

## 3. Harmonise exposure and outcome to the same effect allele
dat <- harmonise_data(exp, out, action = 2)

## Instrument strength (per-SNP F = beta^2 / se^2)
dat$Fstat <- dat$beta.exposure^2 / dat$se.exposure^2
cat(sprintf("Variants used: %d | mean F = %.0f\n",
            sum(dat$mr_keep), mean(dat$Fstat[dat$mr_keep])))

## 4. Main MR estimates
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
          "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
res_or <- generate_odds_ratios(res)
print(res_or[, c("method", "nsnp", "or", "or_lci95", "or_uci95", "pval")])

## 5. Sensitivity / pleiotropy
egger_intercept <- mr_pleiotropy_test(dat)   # MR-Egger intercept (directional pleiotropy)
heterogeneity   <- mr_heterogeneity(dat)     # Cochran's Q
print(egger_intercept)
print(heterogeneity)

## MR-PRESSO (global test + outlier-corrected estimate)
if (requireNamespace("MRPRESSO", quietly = TRUE)) {
  set.seed(1)
  presso <- MRPRESSO::mr_presso(
    BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome", SdExposure = "se.exposure",
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
    data = dat[dat$mr_keep, ], NbDistribution = 1000, SignifThreshold = 0.05)
  print(presso$`Main MR results`)
  print(presso$`MR-PRESSO results`$`Global Test`)
}

## 6. Per-SNP outputs (used to draw scatter / funnel / leave-one-out figures)
single <- mr_singlesnp(dat)
loo    <- mr_leaveoneout(dat)

write.csv(dat,    file.path(OUTDIR, "harmonised_eBMD_PAD_R5.csv"), row.names = FALSE)
write.csv(res_or, file.path(OUTDIR, "univariable_results_R5.csv"), row.names = FALSE)
write.csv(single, file.path(OUTDIR, "singlesnp_R5.csv"),           row.names = FALSE)
write.csv(loo,    file.path(OUTDIR, "loo_R5.csv"),                 row.names = FALSE)
write.csv(egger_intercept, file.path(OUTDIR, "egger_intercept_R5.csv"), row.names = FALSE)
write.csv(heterogeneity,   file.path(OUTDIR, "heterogeneity_R5.csv"),   row.names = FALSE)
cat("01_main_mr done.\n")
