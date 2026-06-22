## =====================================================================
## 05_finngen_r12_sensitivity.R
## Endpoint-definition sensitivity: eBMD -> FinnGen R12 peripheral
## atherosclerosis (I9_ATHSCLE). Reproduces: IVW OR ~1.04 (0.99-1.09),
## consistent with the null primary finding in a larger, broader endpoint.
##
## NOTE: FinnGen R12 is NOT on OpenGWAS. Download the summary file from
##   https://www.finngen.fi/en/access_results
## (endpoint I9_ATHSCLE, release R12) and place it at:
##   data_raw/finngen_R12_I9_ATHSCLE.gz
## =====================================================================
source("R/00_config.R")

r12_file <- "data_raw/finngen_R12_I9_ATHSCLE.gz"
if (!file.exists(r12_file))
  stop("Place the FinnGen R12 I9_ATHSCLE summary file at ", r12_file)

exp <- extract_instruments(EBMD, p1 = P_THRESH, clump = TRUE,
                           r2 = CLUMP_R2, kb = CLUMP_KB)

## FinnGen summary columns: #chrom pos ref alt rsids pval mlogp beta sebeta af_alt ...
r12 <- data.table::fread(r12_file)
out <- format_data(as.data.frame(r12), type = "outcome",
                   snp_col = "rsids", beta_col = "beta", se_col = "sebeta",
                   effect_allele_col = "alt", other_allele_col = "ref",
                   eaf_col = "af_alt", pval_col = "pval")
out <- out[out$SNP %in% exp$SNP, ]

dat <- harmonise_data(exp, out, action = 2)
res <- generate_odds_ratios(mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                                    "mr_weighted_median")))
print(res[, c("method", "nsnp", "or", "or_lci95", "or_uci95", "pval")])

## per-SNP outputs for Supplementary Figure S1 (scatter / funnel / leave-one-out)
write.csv(dat,              file.path(OUTDIR, "harmonised_eBMD_ATHSCLE_R12.csv"), row.names = FALSE)
write.csv(res,             file.path(OUTDIR, "results_R12.csv"),  row.names = FALSE)
write.csv(mr_singlesnp(dat), file.path(OUTDIR, "singlesnp.csv"),  row.names = FALSE)
write.csv(mr_leaveoneout(dat), file.path(OUTDIR, "loo.csv"),      row.names = FALSE)
cat("05_finngen_r12_sensitivity done.\n")
