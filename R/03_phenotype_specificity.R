## =====================================================================
## 03_phenotype_specificity.R
## Control outcomes: genetically predicted eBMD -> CAD and ischemic stroke
## Reproduces: CAD OR ~1.05 (1.005-1.106), p~0.031 (fails Bonferroni 0.017);
##             ischemic stroke OR ~1.00, p~0.97
## Bonferroni threshold = 0.05 / 3 vascular outcomes (PAD, CAD, stroke) = 0.017
## =====================================================================
source("R/00_config.R")

exp <- extract_instruments(EBMD, p1 = P_THRESH, clump = TRUE,
                           r2 = CLUMP_R2, kb = CLUMP_KB)

controls <- c("Coronary artery disease" = CAD,
              "Ischemic stroke"          = STROKE)

out_all <- list()
for (nm in names(controls)) {
  o   <- extract_outcome_data(snps = exp$SNP, outcomes = controls[[nm]])
  d   <- harmonise_data(exp, o, action = 2)
  r   <- generate_odds_ratios(mr(d, method_list = "mr_ivw"))
  r$outcome_label <- nm
  out_all[[nm]] <- r
  cat(sprintf("%-26s OR = %.3f (%.3f-%.3f), p = %.3g\n",
              nm, r$or, r$or_lci95, r$or_uci95, r$pval))
}
res <- do.call(rbind, out_all)
res$bonferroni_0.017 <- ifelse(res$pval < 0.05/3, "significant", "ns")
write.csv(res, file.path(OUTDIR, "phenotype_specificity.csv"), row.names = FALSE)
cat("03_phenotype_specificity done.\n")
