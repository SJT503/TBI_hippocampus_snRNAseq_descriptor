# 03_B5_marker_validation.R — B5 细胞类型 marker 发表级验证
# 子项目 2 · Technical Validation · 2026-05-18
# 经 00_gates.R 门禁读 atlas；产出注释 UMAP + marker 气泡图（证明 7 类注释可靠）
# 纯描述性，不做组间定量声明（venue 红线）

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

source("00_gates.R")
PROJ <- normalizePath("..")
FIG  <- file.path(PROJ, "results", "figures")
TAB  <- file.path(PROJ, "results", "qc_tables")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

seu <- load_definitive_rds("atlas")
DefaultAssay(seu) <- "RNA"                 # marker 表达用 RNA assay

CT_COL    <- "novogene_celltype"
GROUP_COL <- "orig.ident"
# 细胞类型展示顺序：神经元 → 大胶质 → 小胶质 → 血管
CT_ORDER  <- c("Excitatory", "Inhibitory", "Astrocytes",
               "Oligodendrocytes", "OPCs", "Microglia", "Pericyte")
stopifnot(setequal(CT_ORDER, unique(seu@meta.data[[CT_COL]])))
seu@meta.data[[CT_COL]] <- factor(seu@meta.data[[CT_COL]], levels = CT_ORDER)
seu@meta.data[[GROUP_COL]] <- factor(seu@meta.data[[GROUP_COL]],
                                     levels = c("Sham", "CCI_Scr", "CCI_TGFbetaR1_KO"))

# ---- 经典 marker（按细胞类型顺序；均已核实存在于对象）----
markers <- list(
  Excitatory       = c("Slc17a7", "Satb2"),
  Inhibitory       = c("Gad1", "Gad2"),
  Astrocytes       = c("Aqp4", "Gfap", "Slc1a3"),
  Oligodendrocytes = c("Plp1", "Mbp", "Mog"),
  OPCs             = c("Pdgfra", "Cspg4", "Olig1"),
  Microglia        = c("Csf1r", "Hexb", "P2ry12", "C1qa"),
  Pericyte         = c("Pdgfrb", "Rgs5", "Vtn")
)
mk_flat <- unlist(markers, use.names = FALSE)
miss <- setdiff(mk_flat, rownames(seu))
if (length(miss)) stop("MARKER FAIL: 缺失 ", paste(miss, collapse = ", "))

ct_pal <- setNames(
  c("#4C72B0","#5B9BD5","#DD8452","#55A868","#8CD17D","#C44E52","#937860"),
  CT_ORDER)

# ============ 图1 · 注释 UMAP（7 类）============
p_umap <- DimPlot(seu, reduction = "umap", group.by = CT_COL,
                  cols = ct_pal, label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("Cell type annotation (n = 35,047 nuclei)") +
  theme(aspect.ratio = 1)
ggsave(file.path(FIG, "B5_umap_celltype.pdf"), p_umap, width = 7, height = 6)

# ============ 图2 · UMAP 按组分面 ============
p_grp <- DimPlot(seu, reduction = "umap", group.by = CT_COL, split.by = GROUP_COL,
                 cols = ct_pal, ncol = 3) +
  ggtitle("Cell type annotation by group") +
  theme(aspect.ratio = 1)
ggsave(file.path(FIG, "B5_umap_by_group.pdf"), p_grp, width = 14, height = 5)

# ============ 图3 · marker 气泡图 ============
p_dot <- DotPlot(seu, features = markers, group.by = CT_COL,
                 cols = c("lightgrey", "#C44E52"), dot.scale = 6) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8),
        strip.text  = element_text(size = 8)) +
  labs(x = NULL, y = NULL, title = "Canonical marker expression by cell type")
ggsave(file.path(FIG, "B5_dotplot_markers.pdf"), p_dot, width = 12, height = 4.5)

# ============ 图4 · 每类 1 个代表 marker FeaturePlot ============
rep_mk <- c(Excitatory="Slc17a7", Inhibitory="Gad1", Astrocytes="Aqp4",
            Oligodendrocytes="Plp1", OPCs="Pdgfra", Microglia="Csf1r",
            Pericyte="Rgs5")
p_feat <- FeaturePlot(seu, features = unname(rep_mk), reduction = "umap",
                      ncol = 4, order = TRUE) &
  theme(aspect.ratio = 1)
ggsave(file.path(FIG, "B5_featureplot_markers.pdf"), p_feat, width = 14, height = 7)

# ============ 表 · marker 平均表达 + 阳性细胞比例（供文中引用）============
expr <- GetAssayData(seu, layer = "data", assay = "RNA")
ct   <- seu@meta.data[[CT_COL]]
mk_tab <- do.call(rbind, lapply(names(markers), function(cn) {
  do.call(rbind, lapply(markers[[cn]], function(g) {
    in_ct  <- ct == cn
    data.frame(
      CellType        = cn,
      Marker          = g,
      MeanExpr_inType = round(mean(expr[g, in_ct]), 3),
      MeanExpr_other  = round(mean(expr[g, !in_ct]), 3),
      PctPos_inType   = round(100 * mean(expr[g, in_ct]  > 0), 1),
      PctPos_other    = round(100 * mean(expr[g, !in_ct] > 0), 1)
    )
  }))
}))
write.csv(mk_tab, file.path(TAB, "B5_marker_specificity.csv"), row.names = FALSE)
cat("\n## B5 marker 特异性（组内 vs 组外）\n"); print(mk_tab, row.names = FALSE)

cat("\n== B5 完成 ==\n")
cat("图:", file.path(FIG, c("B5_umap_celltype.pdf","B5_umap_by_group.pdf",
                            "B5_dotplot_markers.pdf","B5_featureplot_markers.pdf")),
    sep = "\n  ")
cat("\n表:", file.path(TAB, "B5_marker_specificity.csv"), "\n")
