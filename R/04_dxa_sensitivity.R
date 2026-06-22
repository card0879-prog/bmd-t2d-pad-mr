## =====================================================================
## 04_dxa_sensitivity.R
## Site-specific DXA BMD (femoral neck, lumbar spine, forearm) -> PAD
## Reproduces: nominal femoral-neck IVW signal invalidated by Egger intercept
##             (directional pleiotropy); not robust -> not a primary result.
## =====================================================================
source("R/00_config.R")

dxa <- c("Femoral neck" = DXA_FN,
         "Lumbar spine" = DXA_LS,
         "Forearm"      = DXA_FA)

rows <- list()
for (nm in names(dxa)) {
  e <- tryCatch(extract_instruments(dxa[[nm]], p1 = P_THRESH, clump = TRUE,
                                    r2 = CLUMP_R2, kb = CLUMP_KB),
                error = function(z) NULL)
  if (is.null(e) || nrow(e) < 3) {
    cat(sprintf("%-14s : too few instruments (uninformative)\n", nm)); next
  }
  o <- extract_outcome_data(snps = e$SNP, outcomes = PAD_R5)
  d <- harmonise_data(e, o, action = 2)
  r <- generate_odds_ratios(mr(d, method_list = c("mr_ivw", "mr_egger_regression",
                                                  "mr_weighted_median")))
  plei <- mr_pleiotropy_test(d)   # Egger intercept p < 0.05 => pleiotropy flag
  cat(sprintf("\n== %s (%d SNPs) ==\n", nm, sum(d$mr_keep)))
  print(r[, c("method", "or", "or_lci95", "or_uci95", "pval")])
  cat(sprintf("Egger intercept = %.4f, p = %.3g\n",
              plei$egger_intercept, plei$pval))
  r$site <- nm; r$egger_intercept_p <- plei$pval
  rows[[nm]] <- r
}
if (length(rows)) write.csv(do.call(rbind, rows),
                            file.path(OUTDIR, "dxa_sensitivity.csv"), row.names = FALSE)
cat("04_dxa_sensitivity done.\n")
