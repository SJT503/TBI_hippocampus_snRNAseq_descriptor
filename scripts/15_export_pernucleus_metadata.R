# 15_export_pernucleus_metadata.R — 导出每核 metadata CSV，含本 descriptor 增值列
# 子项目 2 · 2026-05-21 · GEO 投稿包 supplementary file
# 经门禁读 atlas + B2/B3 csv 合并

suppressPackageStartupMessages(library(Seurat))
source("00_gates.R")
PROJ <- normalizePath("..")
TAB  <- file.path(PROJ, "results", "qc_tables")
OUT  <- file.path(PROJ, "geo_submission", "processed", "per_nucleus_metadata.csv")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

seu <- load_definitive_rds("atlas")
DefaultAssay(seu) <- "RNA"

mt <- grep("^mt-",   rownames(seu), value = TRUE, ignore.case = TRUE)
rb <- grep("^Rp[sl]", rownames(seu), value = TRUE, ignore.case = TRUE)
seu[["percent.mt"]]   <- PercentageFeatureSet(seu, features = mt)
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, features = rb)

b23 <- read.csv(file.path(TAB, "B2_B3_per_cell.csv"), check.names = FALSE)
rownames(b23) <- b23$barcode
b23 <- b23[colnames(seu), ]
stopifnot(all(b23$barcode == colnames(seu)))

umap <- as.data.frame(Embeddings(seu, "umap")); colnames(umap) <- c("UMAP_1","UMAP_2")
tsne <- as.data.frame(Embeddings(seu, "tsne")); colnames(tsne) <- c("tSNE_1","tSNE_2")

out <- data.frame(
  barcode               = colnames(seu),
  library               = as.character(seu$orig.ident),
  nCount_RNA            = round(seu$nCount_RNA),
  nFeature_RNA          = seu$nFeature_RNA,
  percent.mt            = round(seu$percent.mt, 4),
  percent.ribo          = round(seu$percent.ribo, 4),
  seurat_clusters       = as.character(seu$seurat_clusters),
  novogene_celltype     = as.character(seu$novogene_celltype),
  novo_broad            = as.character(seu$novo_broad),
  SingleR_label         = as.character(seu$SingleR_label),
  SingleR_broad         = as.character(seu$SingleR_broad),
  scDblFinder_class     = b23$scDblFinder_class,
  scDblFinder_score     = round(b23$scDblFinder_score, 4),
  decontX_contamination = round(b23$decontX_contamination, 4),
  UMAP_1                = round(umap$UMAP_1, 4),
  UMAP_2                = round(umap$UMAP_2, 4),
  tSNE_1                = round(tsne$tSNE_1, 4),
  tSNE_2                = round(tsne$tSNE_2, 4),
  stringsAsFactors      = FALSE
)
write.csv(out, OUT, row.names = FALSE)
cat(sprintf("Wrote %d rows × %d cols to:\n  %s\n", nrow(out), ncol(out), OUT))
cat("File size:", round(file.info(OUT)$size/1e6, 2), "MB\n")
