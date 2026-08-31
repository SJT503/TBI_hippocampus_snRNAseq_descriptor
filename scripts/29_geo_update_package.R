# 29_geo_update_package.R — R2 数据可及性更新包:合并metadata(含亚区) + 假_bulk计数CSV
suppressPackageStartupMessages({library(Seurat); library(Matrix)})

BASE <- "SciData_descriptor"
libs <- c(Sham="Sham", CCI_Scr="CCI_Scr", CCI_TGFbetaR1_KO="CCI_TGFbetaR1_KO")
meta <- read.csv(file.path(BASE,"geo_submission/processed/per_nucleus_metadata.csv"))
subr <- read.csv(file.path(BASE,"results/revision_analysis/subregion/per_nucleus_subregion.csv"))
stopifnot(nrow(meta)==35047, nrow(subr)==18131)

# 亚区并入 metadata(barcode 直接同源,无需剥后缀——都来自 atlas)
meta2 <- merge(meta, subr[,c("barcode","subregion","z_DG","z_CA3","z_CA2","z_CA1","z_margin")], by="barcode", all.x=TRUE)
stopifnot(nrow(meta2)==35047)
meta2$subregion[is.na(meta2$subregion)] <- ""
outdir <- file.path(BASE,"geo_submission/update_20260824")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
write.csv(meta2, file.path(outdir,"per_nucleus_metadata_v2.csv"), row.names=FALSE)

# 假_bulk:逐库读已发布 filtered 矩阵,求列和;矩阵 barcode 需剥 -1 后与 metadata 对齐
stripbc <- function(x) sub("-\\d+$","",x)
pb_all <- NULL; pb_ct <- NULL
for (lib in libs) {
  m <- Read10X(file.path(BASE,"geo_submission/processed",lib,"filtered_feature_bc_matrix"))
  genes <- rownames(m)
  cs_all <- Matrix::rowSums(m)
  pb_all <- cbind(pb_all, cs_all)
  colnames(pb_all)[ncol(pb_all)] <- lib
  ml <- meta[meta$library==lib,]
  ml$bc <- stripbc(ml$barcode)   # 两侧都剥 GEM 后缀
  bc <- match(stripbc(colnames(m)), ml$bc)
  ct <- ml$novogene_celltype[bc]
  for (g in sort(unique(ct[!is.na(ct)]))) {
    idx <- which(ct==g)
    v <- Matrix::rowSums(m[,idx,drop=FALSE])
    nm <- paste0(lib,"|",g)
    pb_ct <- if (is.null(pb_ct)) { m2 <- matrix(v, ncol=1, dimnames=list(names(v), nm)); m2 }
             else cbind(pb_ct, v)
    colnames(pb_ct)[ncol(pb_ct)] <- nm
  }
  cat(lib, "unmatched barcodes:", sum(is.na(bc)), "of", ncol(m), "\n")
  cat(lib, "done:", ncol(m), "barcodes\n")
}
stopifnot(identical(rownames(pb_all), rownames(pb_ct)))
write.csv(data.frame(gene=rownames(pb_all), as.matrix(pb_all), check.names=FALSE),
          file.path(outdir,"pseudobulk_by_condition.csv"), row.names=FALSE)
write.csv(data.frame(gene=rownames(pb_ct), as.matrix(pb_ct), check.names=FALSE),
          file.path(outdir,"pseudobulk_by_condition_celltype.csv"), row.names=FALSE)

# 释放注记
writeLines(c(
 "GSE333879 supplementary update (2026-08-24 revision)",
 "1. per_nucleus_metadata_v2.csv  — adds hippocampal-subregion annotation for excitatory nuclei",
 "   (subregion label, z_DG/z_CA3/z_CA2/z_CA1 marker scores, z_margin; Unassigned where below thresholds).",
 "2. pseudobulk_by_condition.csv — per-gene summed counts per library (3 columns), CSV for direct use.",
 "3. pseudobulk_by_condition_celltype.csv — per-gene summed counts per library x cell type (21 columns).",
 "Validation-image files (in-vitro sgRNA screening Western blots; in-vivo HA-tag and TGF-betaRI IF) will be",
 "added in the same update once final high-resolution versions are supplied. [PENDING-IMAGES]"
), file.path(outdir,"README_UPDATE.txt"))
cat("DONE ->", normalizePath(outdir), "\n")
