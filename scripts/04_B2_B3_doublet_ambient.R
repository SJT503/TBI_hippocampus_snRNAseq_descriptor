# 04_B2_B3_doublet_ambient.R — B2 doublet 检测 + B3 ambient RNA 评估
# 子项目 2 · Technical Validation · 2026-05-18
# 经 00_gates.R 门禁读 atlas；scDblFinder 估 doublet，decontX 估 ambient 污染
# 均按 orig.ident 分库运行（doublet/ambient 是库内现象）；纯描述性，不做组间定量声明

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(decontX)
  library(ggplot2)
})
set.seed(20260518)

source("00_gates.R")
PROJ <- normalizePath("..")
TAB  <- file.path(PROJ, "results", "qc_tables")
FIG  <- file.path(PROJ, "results", "figures")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

seu <- load_definitive_rds("atlas")
DefaultAssay(seu) <- "RNA"
GROUP_COL <- "orig.ident"
CT_COL    <- "novogene_celltype"
GROUP_LV  <- c("Sham", "CCI_Scr", "CCI_TGFbetaR1_KO")

# ---- 转 SCE（仅 RNA counts）----
sce <- SingleCellExperiment(
  assays  = list(counts = GetAssayData(seu, layer = "counts", assay = "RNA")),
  colData = seu@meta.data)

# ============ B2 · scDblFinder doublet 检测（按库）============
cat("\n[B2] scDblFinder 运行中（按 orig.ident 分库）...\n")
sce <- scDblFinder(sce, samples = sce[[GROUP_COL]])
seu$scDblFinder_class <- sce$scDblFinder.class
seu$scDblFinder_score <- sce$scDblFinder.score

b2 <- do.call(rbind, lapply(c(GROUP_LV, "All"), function(g) {
  cl <- if (g == "All") seu$scDblFinder_class
        else seu$scDblFinder_class[seu@meta.data[[GROUP_COL]] == g]
  data.frame(Group = g, n = length(cl),
             n_doublet = sum(cl == "doublet"),
             doublet_rate_pct = round(100 * mean(cl == "doublet"), 2))
}))
write.csv(b2, file.path(TAB, "B2_doublet_summary.csv"), row.names = FALSE)
cat("\n## B2 doublet 汇总\n"); print(b2, row.names = FALSE)

# ============ B3 · decontX ambient RNA 评估（按库）============
cat("\n[B3] decontX 运行中（按 orig.ident 分库）...\n")
sce <- decontX(sce, batch = sce[[GROUP_COL]])
seu$decontX_contamination <- sce$decontX_contamination

b3 <- do.call(rbind, lapply(c(GROUP_LV, "All"), function(g) {
  ct <- if (g == "All") seu$decontX_contamination
        else seu$decontX_contamination[seu@meta.data[[GROUP_COL]] == g]
  data.frame(Group = g, n = length(ct),
             mean_contam_pct   = round(100 * mean(ct), 2),
             median_contam_pct = round(100 * median(ct), 2),
             pct_cells_contam_gt20 = round(100 * mean(ct > 0.2), 2))
}))
write.csv(b3, file.path(TAB, "B3_ambient_summary.csv"), row.names = FALSE)
cat("\n## B3 ambient 污染汇总\n"); print(b3, row.names = FALSE)

# ---- 每细胞结果导出（供 descriptor processed metadata / Usage Notes）----
cell_tab <- data.frame(
  barcode       = colnames(seu),
  group         = seu@meta.data[[GROUP_COL]],
  cell_type     = seu@meta.data[[CT_COL]],
  scDblFinder_class = seu$scDblFinder_class,
  scDblFinder_score = round(seu$scDblFinder_score, 4),
  decontX_contamination = round(seu$decontX_contamination, 4))
write.csv(cell_tab, file.path(TAB, "B2_B3_per_cell.csv"), row.names = FALSE)

# ============ 图 ============
seu@meta.data[[GROUP_COL]] <- factor(seu@meta.data[[GROUP_COL]], levels = GROUP_LV)

# B2 doublet UMAP + 每组分面
p_db <- DimPlot(seu, reduction = "umap", group.by = "scDblFinder_class",
                cols = c(singlet = "grey80", doublet = "#C44E52"),
                order = "doublet") +
  ggtitle("scDblFinder doublet classification") + theme(aspect.ratio = 1)
ggsave(file.path(FIG, "B2_doublet_umap.pdf"), p_db, width = 7, height = 6)

# B3 ambient 污染 UMAP + 小提琴
p_c1 <- FeaturePlot(seu, "decontX_contamination", reduction = "umap", order = TRUE) +
  scale_color_viridis_c() + ggtitle("decontX ambient contamination") +
  theme(aspect.ratio = 1)
p_c2 <- VlnPlot(seu, "decontX_contamination", group.by = GROUP_COL, pt.size = 0) +
  theme(axis.title.x = element_blank()) + NoLegend()
ggsave(file.path(FIG, "B3_ambient_umap_violin.pdf"),
       patchwork::wrap_plots(p_c1, p_c2, ncol = 2), width = 12, height = 5)

cat("\n== B2+B3 完成 ==\n")
cat("表:", file.path(TAB, c("B2_doublet_summary.csv","B3_ambient_summary.csv",
                            "B2_B3_per_cell.csv")), sep = "\n  ")
cat("\n图:", file.path(FIG, c("B2_doublet_umap.pdf","B3_ambient_umap_violin.pdf")),
    sep = "\n  ")
cat("\n")
