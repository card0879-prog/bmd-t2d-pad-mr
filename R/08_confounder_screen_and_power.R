############################################################
## Paper 2 — PLOS ONE 补充分析
## PhenoScanner 混杂筛查 + 剔除混杂 SNP 后的 IVW 再分析
## 在 RStudio 里整段运行即可。只需改第 10 行的项目路径。
## 跑完把 results/confounder_screen_summary.txt 发给我。
############################################################

## ---- 0. 设置 ---------------------------------------------------------------
PROJECT_DIR <- "~/Desktop/ukbank延续/bmd-t2d-pad-mr"   # <-- 改成你的项目根目录
setwd(PROJECT_DIR)

EBMD_ID  <- "ebi-a-GCST006979"   # heel eBMD (UKB, Morris 2019)
PAD_ID   <- "finn-b-I9_PAD"      # PAD, FinnGen R5
P_THRESH <- 5e-8
CLUMP_R2 <- 0.001
CLUMP_KB <- 10000

## OpenGWAS token：优先从 .Renviron 读 OPENGWAS_JWT；读不到就把 token 粘到引号里
TOKEN <- Sys.getenv("OPENGWAS_JWT")
if (TOKEN == "") TOKEN <- ""     # <-- 若上面为空，把你的 token 粘进这对引号
if (TOKEN != "") Sys.setenv(OPENGWAS_JWT = TOKEN)

## ---- 1. 依赖包（缺了自动装）------------------------------------------------
need <- function(pkg, gh = NULL){
  if (!requireNamespace(pkg, quietly = TRUE)){
    message("installing ", pkg, " ...")
    if (is.null(gh)) install.packages(pkg) else {
      if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
      remotes::install_github(gh, upgrade = "never")
    }
  }
}
need("data.table")
need("TwoSampleMR",  gh = "MRCIEU/TwoSampleMR")
need("ieugwasr",     gh = "MRCIEU/ieugwasr")
need("phenoscanner", gh = "phenoscanner/phenoscanner")
suppressMessages({library(TwoSampleMR); library(data.table)})

## ---- 2. 重建 493 个 SNP 的 harmonised 数据 --------------------------------
csv <- "results/harmonised_eBMD_PAD_R5.csv"
if (file.exists(csv)) {
  dat <- as.data.frame(fread(csv))
  message("已读取现成 harmonised 数据：", nrow(dat), " SNPs")
} else {
  message("未找到现成数据，从 OpenGWAS 重新提取 ...")
  exp <- extract_instruments(EBMD_ID, p1 = P_THRESH, clump = TRUE,
                             r2 = CLUMP_R2, kb = CLUMP_KB)
  out <- extract_outcome_data(snps = exp$SNP, outcomes = PAD_ID)
  dat <- harmonise_data(exp, out, action = 2)
  dir.create("results", showWarnings = FALSE)
  fwrite(dat, csv)
  message("重建完成：", nrow(dat), " SNPs")
}
dat   <- dat[dat$mr_keep %in% c(TRUE, NA) | is.na(dat$mr_keep), ]
if ("mr_keep" %in% names(dat)) dat <- dat[dat$mr_keep == TRUE, ]
rsids <- unique(dat$SNP)
cat("\n纳入分析的工具变量数：", length(rsids), "\n")

## ---- 3. PhenoScanner 混杂筛查 ---------------------------------------------
library(phenoscanner)
chunks <- split(rsids, ceiling(seq_along(rsids) / 100))   # 每批 100 个
ph <- do.call(rbind, lapply(seq_along(chunks), function(i){
  Sys.sleep(1)
  message("  PhenoScanner 批次 ", i, "/", length(chunks))
  phenoscanner(snpquery = chunks[[i]], pvalue = 5e-8)$results
}))

# 针对 eBMD->PAD 问题的心血管代谢混杂关键词
conf <- paste0("(?i)(smok|cigarett|type 2 diabet|t2d|glucose|hba1c|",
               "ldl|hdl|cholesterol|triglycerid|lipid|",
               "blood pressure|hypertens|systolic|diastolic|",
               "body mass|bmi|obesit|adiposit|waist|c-reactive|crp)")
ph$trait  <- as.character(ph$trait)
hits      <- ph[grepl(conf, ph$trait, perl = TRUE), c("snp","trait","beta","p")]
conf_snps <- sort(unique(hits$snp))
cat("\n与心血管代谢混杂关联的 SNP (p<5e-8)：",
    length(conf_snps), " / ", length(rsids), "\n", sep = "")

## ---- 4. 主 IVW vs 剔除混杂 SNP 后的 IVW -----------------------------------
res_all  <- mr(dat, method_list = "mr_ivw")
dat_excl <- dat[!dat$SNP %in% conf_snps, ]
res_excl <- mr(dat_excl, method_list = "mr_ivw")
show <- function(r) sprintf("OR %.3f (95%% CI %.3f-%.3f), p=%.3g, nSNP=%d",
  exp(r$b), exp(r$b - 1.96*r$se), exp(r$b + 1.96*r$se), r$pval, r$nsnp)
cat("\nIVW 全部工具变量 ：", show(res_all),  "\n")
cat(  "IVW 剔除混杂后   ：", show(res_excl), "\n")

## ---- 5. 保存汇总 -----------------------------------------------------------
dir.create("results", showWarnings = FALSE)
sink("results/confounder_screen_summary.txt")
cat("Paper 2 — PLOS ONE supplementary analyses\n")
cat("Date:", as.character(Sys.Date()), "\n\n")
cat("Instruments in primary analysis:", length(rsids), "\n")
cat("SNPs linked to a cardiometabolic confounder (PhenoScanner, p<5e-8):",
    length(conf_snps), "\n\n")
cat("Confounder-associated SNPs:\n"); print(conf_snps)
cat("\nTrait hits:\n"); print(hits)
cat("\nIVW all instruments :", show(res_all),  "\n")
cat(  "IVW excl. confounded:", show(res_excl), "\n")
sink()
cat("\n已保存：results/confounder_screen_summary.txt\n")
cat(">>> 把这个 txt 发我即可 <<<\n")