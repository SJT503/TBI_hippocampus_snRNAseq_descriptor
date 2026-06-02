# 01_inspect_atlas.R — 探查 novogene_atlas_annotated.rds 真实结构
# 一次性探查脚本，用于校准 00_gates.R 与撰写 B1/B6。直接 readRDS（建门禁前的例外）。
# 子项目 2 · 2026-05-18

suppressPackageStartupMessages(library(Seurat))

rds <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/seurat/novogene_atlas_annotated.rds"
cat("== 读取:", rds, "==\n")
seu <- readRDS(rds)

cat("\n## 基本信息\n")
cat("class      :", class(seu), "\n")
cat("ncol(细胞) :", ncol(seu), "\n")
cat("nrow(基因) :", nrow(seu), "\n")
cat("Seurat ver :", as.character(seu@version), "\n")

cat("\n## Assays / layers\n")
cat("assays     :", paste(Assays(seu), collapse = ", "), "\n")
cat("DefaultAssay:", DefaultAssay(seu), "\n")
for (a in Assays(seu)) {
  ly <- tryCatch(Layers(seu[[a]]), error = function(e) NA)
  cat(sprintf("  [%s] layers: %s\n", a, paste(ly, collapse = ", ")))
}

cat("\n## meta.data 列名 (", ncol(seu@meta.data), "列)\n")
print(colnames(seu@meta.data))

cat("\n## meta.data 各列类型与样例\n")
for (cn in colnames(seu@meta.data)) {
  v <- seu@meta.data[[cn]]
  if (is.numeric(v)) {
    cat(sprintf("  %-22s [num]  range %.3g ~ %.3g  median %.3g\n",
                cn, min(v, na.rm = TRUE), max(v, na.rm = TRUE), median(v, na.rm = TRUE)))
  } else {
    u <- unique(as.character(v))
    cat(sprintf("  %-22s [%s] %d 唯一值: %s\n",
                cn, class(v)[1], length(u),
                paste(head(u, 12), collapse = " | ")))
  }
}

cat("\n## 候选分组列 table()\n")
for (cn in intersect(c("group","Group","orig.ident","condition","Condition",
                        "sample","Sample","stim","treatment"), colnames(seu@meta.data))) {
  cat(">>", cn, "\n"); print(table(seu@meta.data[[cn]])); cat("\n")
}

cat("\n## 候选细胞类型列 table()\n")
for (cn in intersect(c("cell_type","celltype","CellType","cell.type","annotation",
                        "Annotation","cellType","SingleR","singleR","celltype_singleR",
                        "predicted.id"), colnames(seu@meta.data))) {
  cat(">>", cn, "\n"); print(table(seu@meta.data[[cn]])); cat("\n")
}

cat("\n## QC 列是否预存\n")
for (cn in c("nFeature_RNA","nCount_RNA","percent.mt","percent.mito","percent.ribo",
             "percent_mt","pct_counts_mt")) {
  cat(sprintf("  %-16s : %s\n", cn,
              if (cn %in% colnames(seu@meta.data)) "存在" else "无"))
}

cat("\n## meta.data head(3)\n")
print(head(seu@meta.data, 3))

cat("\n== 探查完成 ==\n")
