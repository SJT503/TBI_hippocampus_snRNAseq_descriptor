# 30_geo_gene_celltype_lookup.R — 基因×细胞类型速查表: pct_detected + mean counts per nucleus
# 服务 Q2.2 (CSV/Excel 可及性): 让"基因X在哪些细胞类型表达"一行查表即答, 零代码。
# 数据源与 29 号完全相同: 已发布 GEO filtered_feature_bc_matrix + per_nucleus_metadata.csv
suppressPackageStartupMessages({library(Seurat); library(Matrix)})

BASE <- "SciData_descriptor"
libs <- c(Sham="Sham", CCI_Scr="CCI_Scr", CCI_TGFbetaR1_KO="CCI_TGFbetaR1_KO")
meta <- read.csv(file.path(BASE,"geo_submission/processed/per_nucleus_metadata.csv"))
stopifnot(nrow(meta)==35047)
stopifnot(all(c("barcode","library","novogene_celltype") %in% names(meta)))

cat("cell types:", paste(sort(unique(meta$novogene_celltype)), collapse=" | "), "\n")
cat("n per type:\n"); print(table(meta$novogene_celltype))

stripbc <- function(x) sub("-\\d+$","",x)
det <- list(); cnt <- list(); ncell <- list(); ref_genes <- NULL
for (lib in libs) {
  m <- Read10X(file.path(BASE,"geo_submission/processed",lib,"filtered_feature_bc_matrix"))
  if (is.null(ref_genes)) ref_genes <- rownames(m)
  stopifnot(identical(rownames(m), ref_genes))          # 三库基因集一致
  ml <- meta[meta$library==lib,]; ml$bc <- stripbc(ml$barcode)
  bc <- match(stripbc(colnames(m)), ml$bc)
  ct <- ml$novogene_celltype[bc]
  cat(lib, "unmatched barcodes:", sum(is.na(ct)), "of", ncol(m), "\n")
  for (g in sort(unique(ct[!is.na(ct)]))) {
    idx <- which(ct==g)
    sub_m <- m[,idx,drop=FALSE]
    d <- Matrix::rowSums(sub_m > 0)
    c_ <- Matrix::rowSums(sub_m)
    if (is.null(det[[g]])) { det[[g]] <- d; cnt[[g]] <- c_; ncell[[g]] <- ncol(sub_m) }
    else { det[[g]] <- det[[g]] + d; cnt[[g]] <- cnt[[g]] + c_; ncell[[g]] <- ncell[[g]] + ncol(sub_m) }
  }
}
types <- names(det)
stopifnot(sum(unlist(ncell)) == 35047)                   # 全部 35,047 核均入表
stopifnot(length(ref_genes) == 32285)

out <- data.frame(gene=ref_genes, check.names=FALSE)
for (g in types) {
  out[[paste0(g,"_pct_detected")]] <- round(100 * det[[g]] / ncell[[g]], 2)
  out[[paste0(g,"_mean_counts")]]  <- round(cnt[[g]] / ncell[[g]], 3)
}
outdir <- file.path(BASE,"geo_submission/update_20260824")
write.csv(out, file.path(outdir,"gene_by_celltype_detection.csv"), row.names=FALSE)

# 生物学 sanity(打印给 Response 引用): 每基因取 pct_detected 前三的细胞类型
for (gn in c("Gfap","Aqp4","Snap25","Syt1","Plp1","Mbp","P2ry12","Hexb")) {
  r <- out[out$gene==gn,]
  pc <- unlist(r[grep("_pct_detected", names(r))])
  top <- head(sort(pc, decreasing=TRUE), 3)
  cat(gn, " top3 pct_detected:", paste(sprintf("%s=%.2f%%", names(top), top), collapse=", "), "\n")
}

# README_UPDATE 同步(4 个文件全列, 保持与 29 号的 [PENDING-IMAGES] 注记一致)
writeLines(c(
 "GSE333879 supplementary update (2026-08-24 revision)",
 "1. per_nucleus_metadata_v2.csv  — adds hippocampal-subregion annotation for excitatory nuclei",
 "   (subregion label, z_DG/z_CA3/z_CA2/z_CA1 marker scores, z_margin; Unassigned where below thresholds;",
 "   standalone copy of the subregion columns: per_nucleus_subregion.csv).",
 "2. pseudobulk_by_condition.csv — per-gene summed counts per library (3 columns), CSV for direct use.",
 "3. pseudobulk_by_condition_celltype.csv — per-gene summed counts per library x cell type (21 columns).",
 "4. gene_by_celltype_detection.csv — per-gene lookup across 7 cell types pooled over all libraries",
 "   (32,285 genes x 7 types: percent of nuclei detecting the gene, and mean UMI counts per nucleus).",
 "Validation-image files (in-vitro sgRNA screening Western blots; in-vivo HA-tag and TGF-betaRI IF) will be",
 "added in the same update once final high-resolution versions are supplied. [PENDING-IMAGES]"
), file.path(outdir,"README_UPDATE.txt"))
cat("\nDONE ->", normalizePath(file.path(outdir,"gene_by_celltype_detection.csv")), "\n")
