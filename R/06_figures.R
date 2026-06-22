## =====================================================================
## 06_figures.R
## Regenerate the primary diagnostic plots (scatter, funnel, leave-one-out,
## forest) in pure R using TwoSampleMR's built-in plotting functions.
## The publication figures were styled separately, but these reproduce the
## same underlying content from the analysis above.
## =====================================================================
source("R/00_config.R")
suppressMessages(library(ggplot2))

exp <- extract_instruments(EBMD, p1 = P_THRESH, clump = TRUE,
                           r2 = CLUMP_R2, kb = CLUMP_KB)
out <- extract_outcome_data(snps = exp$SNP, outcomes = PAD_R5)
dat <- harmonise_data(exp, out, action = 2)
res <- mr(dat)

p_scatter <- mr_scatter_plot(res, dat)
p_funnel  <- mr_funnel_plot(mr_singlesnp(dat))
p_loo     <- mr_leaveoneout_plot(mr_leaveoneout(dat))
p_forest  <- mr_forest_plot(mr_singlesnp(dat))

ggsave(file.path(OUTDIR, "fig_scatter_R5.pdf"), p_scatter[[1]], width = 6, height = 6)
ggsave(file.path(OUTDIR, "fig_funnel_R5.pdf"),  p_funnel[[1]],  width = 6, height = 6)
ggsave(file.path(OUTDIR, "fig_leaveoneout_R5.pdf"), p_loo[[1]], width = 6, height = 8)
ggsave(file.path(OUTDIR, "fig_forest_R5.pdf"),  p_forest[[1]], width = 6, height = 8)
cat("06_figures done -> see", OUTDIR, "\n")
