# 10_inspect_astrocyte.R — 探查 astrocyte_v2_unified.rds 结构
# 子项目 2 · 2026-05-18 · A1 supplementary figure 前置探查

suppressPackageStartupMessages(library(Seurat))
rds <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/seurat/astrocyte_v2_unified.rds"
seu <- readRDS(rds)

cat("class:", class(seu), "\n")
cat("ncol:", ncol(seu), "  nrow:", nrow(seu), "\n")
cat("Seurat ver:", as.character(seu@version), "\n")
cat("Assays:", paste(Assays(seu), collapse=", "), "\n")
cat("DefaultAssay:", DefaultAssay(seu), "\n")
cat("Reductions:", paste(Reductions(seu), collapse=", "), "\n")
cat("\nmeta.data columns:\n"); print(colnames(seu@meta.data))
cat("\n各列样例:\n")
for (cn in colnames(seu@meta.data)) {
  v <- seu@meta.data[[cn]]
  if (is.numeric(v)) {
    cat(sprintf("  %-30s [num] range %.2g~%.2g\n", cn, min(v,na.rm=T), max(v,na.rm=T)))
  } else {
    u <- unique(as.character(v))
    cat(sprintf("  %-30s [%s] %d uniq: %s\n", cn, class(v)[1], length(u),
                paste(head(u,8), collapse=" | ")))
  }
}
cat("\n候选分组/亚群列 table:\n")
for (cn in intersect(c("orig.ident","stim","seurat_clusters","integrated_snn_res.0.6",
                        "integrated_snn_res.0.4","integrated_snn_res.0.8",
                        "subcluster","subtype","astrocyte_subtype","ast_subtype"),
                     colnames(seu@meta.data))) {
  cat("\n>>", cn, "\n"); print(table(seu@meta.data[[cn]]))
}
