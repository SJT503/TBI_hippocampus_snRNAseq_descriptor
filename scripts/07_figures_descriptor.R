# 07_figures_descriptor.R — Data Descriptor 发表级 Technical Validation 图集
# 子项目 2 · 2026-05-18 (rev2) · Scientific Data 风格（Okabe-Ito 色盲安全配色，矢量 PDF）
# 经 00_gates.R 门禁读 atlas；doublet/ambient 复用 B2_B3_per_cell.csv（不重算）
# 4 张主图 + 2 补充图。纯技术验证，不做组间定量生物学声明（venue 红线）
# rev2 修复：Fig1c 分面标签截断、Fig2b FeaturePlot 对比度与面板标号、Fig4a 着色 bug

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})

source("00_gates.R")
PROJ <- normalizePath("..")
TAB  <- file.path(PROJ, "results", "qc_tables")
FIG  <- file.path(PROJ, "results", "figures_publication")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ============ 配色与主题 ============
CT  <- c("Excitatory","Inhibitory","Astrocytes","Oligodendrocytes","OPCs","Microglia","Pericyte")
CTP <- setNames(c("#0072B2","#56B4E9","#D55E00","#009E73","#E69F00","#CC79A7","#999999"), CT)
GR  <- c("Sham","CCI_Scr","CCI_TGFbetaR1_KO")
GRP <- setNames(c("#5B8FA8","#E69F00","#BB3E03"), GR)

theme_sd <- function(bs = 7) {
  theme_bw(base_size = bs, base_family = "sans") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = bs - 1, colour = "black"),
          plot.title = element_text(size = bs + 1, hjust = 0),
          legend.text = element_text(size = bs - 1),
          legend.title = element_text(size = bs),
          legend.key.size = unit(3, "mm"),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = bs))
}
theme_set(theme_sd())
umap_blank <- theme(aspect.ratio = 1, axis.title = element_blank(),
                    axis.text = element_blank(), axis.ticks = element_blank())

# ============ 数据加载 ============
seu <- load_definitive_rds("atlas")
DefaultAssay(seu) <- "RNA"
seu$novogene_celltype <- factor(seu$novogene_celltype, levels = CT)
seu$orig.ident        <- factor(seu$orig.ident, levels = GR)

mt <- grep("^mt-", rownames(seu), value = TRUE, ignore.case = TRUE)
rb <- grep("^Rp[sl]", rownames(seu), value = TRUE, ignore.case = TRUE)
seu[["percent.mt"]]   <- PercentageFeatureSet(seu, features = mt)
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, features = rb)

b23 <- read.csv(file.path(TAB, "B2_B3_per_cell.csv"), check.names = FALSE)
rownames(b23) <- b23$barcode
b23 <- b23[colnames(seu), ]
stopifnot(all(b23$barcode == colnames(seu)))            # 对齐校验
seu$scDblFinder_class     <- b23$scDblFinder_class
seu$scDblFinder_score     <- b23$scDblFinder_score
seu$decontX_contamination <- b23$decontX_contamination

b1 <- read.csv(file.path(TAB, "B1_per_group_qc.csv"), check.names = FALSE)
b4 <- read.csv(file.path(TAB, "B4_sequencing_metrics.csv"), check.names = FALSE)
b6 <- read.csv(file.path(TAB, "B6_correlations.csv"), check.names = FALSE)
md <- seu@meta.data

# UMAP 坐标（手动绘图用，着色完全可控）
umap <- as.data.frame(Embeddings(seu, "umap")); colnames(umap) <- c("U1","U2")
umap <- cbind(umap, md[, c("novogene_celltype","orig.ident","scDblFinder_class",
                           "scDblFinder_score","decontX_contamination")])

# ================================================================
# Figure 1 — Dataset overview & study design
# ================================================================
stages <- data.frame(x = 1:4,
  lab = c("AAV injection\n(hippocampal CA3)","CCI or sham\nsurgery",
          "Injured-side\nhippocampus","snRNA-seq\n(10x 3' v3)"))
links <- data.frame(x = c(1.34,2.34,3.34), xe = c(1.66,2.66,3.66),
                    lab = c("3 weeks","30 days",""))
f1a <- ggplot() +
  geom_segment(data = links, aes(x, 1, xend = xe, yend = 1),
               arrow = arrow(length = unit(1.6,"mm"), type = "closed"),
               linewidth = .4, colour = "grey40") +
  geom_text(data = links, aes((x+xe)/2, 1.17, label = lab), size = 2.2,
            colour = "grey25") +
  geom_tile(data = stages, aes(x, 1), width = .64, height = .44, fill = "#EAF1F5",
            colour = "#5B8FA8", linewidth = .4) +
  geom_text(data = stages, aes(x, 1, label = lab), size = 2.3, lineheight = .92) +
  annotate("text", x = 2.5, y = .66,
           label = "3 conditions (Sham / CCI_Scr / CCI_TGFbetaR1_KO), 8 mice pooled per condition",
           size = 2.4, colour = "grey25") +
  scale_y_continuous(limits = c(.5,1.32)) + scale_x_continuous(limits = c(.55,4.45)) +
  labs(title = "Study design") +
  theme_void(base_size = 7) + theme(plot.title = element_text(size = 8, hjust = 0))

f1b <- ggplot(umap[sample(nrow(umap)), ], aes(U1, U2, colour = novogene_celltype)) +
  geom_point(size = .22, stroke = 0) +
  scale_colour_manual(values = CTP, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.8))) +
  labs(title = "Cell-type annotation (35,047 nuclei)", x = "UMAP 1", y = "UMAP 2") +
  theme_sd() + theme(aspect.ratio = 1)

comp <- as.data.frame(table(md$orig.ident, md$novogene_celltype))
names(comp) <- c("Group","CellType","n")
f1d <- ggplot(comp, aes(Group, n, fill = CellType)) +
  geom_col(position = "fill", width = .7) +
  scale_fill_manual(values = CTP, name = NULL) +
  scale_y_continuous("Proportion", expand = c(0,0)) +
  labs(title = "Cell-type composition", x = NULL) +
  theme_sd() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

f1c <- ggplot(umap[sample(nrow(umap)), ], aes(U1, U2, colour = novogene_celltype)) +
  geom_point(size = .2, stroke = 0) +
  scale_colour_manual(values = CTP, name = NULL) +
  facet_wrap(~ orig.ident, nrow = 1) +
  guides(colour = guide_legend(override.aes = list(size = 1.8))) +
  labs(title = "Cell-type annotation by condition", x = "UMAP 1", y = "UMAP 2") +
  theme_sd() + theme(aspect.ratio = 1)

fig1 <- (f1a / (f1b | f1d) / f1c) +
  plot_layout(heights = c(.42, 1, .72)) +
  plot_annotation(tag_levels = "a", title = "Figure 1. Dataset overview",
                  theme = theme(plot.title = element_text(size = 9, face = "bold")))
ggsave(file.path(FIG, "Figure1_overview.pdf"), fig1, width = 175, height = 190,
       units = "mm")

# ================================================================
# Figure 2 — Cell-type identity validation
# ================================================================
markers <- list(
  Excitatory = c("Slc17a7","Satb2"), Inhibitory = c("Gad1","Gad2"),
  Astrocytes = c("Aqp4","Gfap","Slc1a3"),
  Oligodendrocytes = c("Plp1","Mbp","Mog"),
  OPCs = c("Pdgfra","Cspg4","Olig1"),
  Microglia = c("Csf1r","Hexb","P2ry12","C1qa"),
  Pericyte = c("Pdgfrb","Rgs5","Vtn"))
mk_flat <- unlist(markers, use.names = FALSE)

f2a <- DotPlot(seu, features = markers, group.by = "novogene_celltype",
               cols = c("grey90","#BB3E03"), dot.scale = 4.5) +
  RotatedAxis() +
  labs(title = "Canonical marker expression", x = NULL, y = NULL) +
  theme_sd() + theme(axis.text.x = element_text(size = 5.5, angle = 90, hjust = 1, vjust = .5),
                     strip.text = element_text(size = 5.6), legend.position = "bottom")

rep_mk <- c("Slc17a7","Gad1","Aqp4","Plp1","Pdgfra","Csf1r","Rgs5")
# 改纯矢量 (raster=FALSE)：35K 点完全能矢量，Sci Data 投稿要求 >=300 DPI，
# 矢量永远清晰，不受 DPI 限制
f2b <- FeaturePlot(seu, rep_mk, reduction = "umap", raster = FALSE,
                   pt.size = 0.3, ncol = 7, order = TRUE,
                   cols = c("grey88","#08306B")) &
  theme_sd() & umap_blank &
  theme(plot.title = element_text(size = 7, face = "italic"),
        legend.position = "none")
f2b <- wrap_elements(full = f2b)        # 整行作为单一面板

expr <- GetAssayData(seu, layer = "data", assay = "RNA")
ct   <- seu$novogene_celltype
avg  <- sapply(CT, function(g) Matrix::rowMeans(expr[mk_flat, ct == g, drop = FALSE]))
zavg <- t(scale(t(avg)))
hm <- reshape2::melt(zavg, varnames = c("Marker","CellType"), value.name = "z")
hm$Marker   <- factor(hm$Marker, levels = rev(mk_flat))
hm$CellType <- factor(hm$CellType, levels = CT)
f2c <- ggplot(hm, aes(CellType, Marker, fill = z)) +
  geom_tile(colour = "white", linewidth = .25) +
  scale_fill_gradient2("Row\nz-score", low = "#2166AC", mid = "white",
                       high = "#B2182B", midpoint = 0) +
  labs(title = "Marker specificity (row-scaled mean expression)", x = NULL, y = NULL) +
  theme_sd() + theme(axis.text.x = element_text(angle = 40, hjust = 1),
                     axis.text.y = element_text(size = 5.5, face = "italic"))

cm  <- table(seu$novogene_celltype, seu$SingleR_broad)
cm  <- cm[, colSums(cm) > 0, drop = FALSE]
cmp <- sweep(cm, 1, rowSums(cm), "/") * 100
cmd <- reshape2::melt(cmp, varnames = c("novogene","SingleR"), value.name = "pct")
cmd$novogene <- factor(cmd$novogene, levels = rev(CT))
f2d <- ggplot(cmd, aes(SingleR, novogene, fill = pct)) +
  geom_tile(colour = "white", linewidth = .25) +
  geom_text(data = subset(cmd, pct >= 5), aes(label = round(pct)), size = 1.9,
            colour = ifelse(subset(cmd, pct >= 5)$pct > 50, "white", "black")) +
  scale_fill_viridis_c("% of row", option = "mako", direction = -1) +
  labs(title = "Concordance with independent SingleR annotation",
       x = "SingleR label", y = "Provider annotation") +
  theme_sd() + theme(axis.text.x = element_text(angle = 40, hjust = 1))

fig2 <- (f2a / f2b / (f2c | f2d)) +
  plot_layout(heights = c(1, .5, 1)) +
  plot_annotation(tag_levels = "a", title = "Figure 2. Cell-type identity validation",
                  theme = theme(plot.title = element_text(size = 9, face = "bold")))
ggsave(file.path(FIG, "Figure2_celltype_validation.pdf"), fig2,
       width = 175, height = 200, units = "mm")

# ================================================================
# Figure 3 — Sequencing and data quality
# ================================================================
qcm <- c("nFeature_RNA","nCount_RNA","percent.mt")
qcl <- c("Genes per nucleus","UMIs per nucleus","Mitochondrial %")
f3a <- wrap_plots(lapply(seq_along(qcm), function(i)
  ggplot(md, aes(orig.ident, .data[[qcm[i]]], fill = orig.ident)) +
    geom_violin(scale = "width", linewidth = .3, colour = "grey30") +
    geom_boxplot(width = .12, outlier.shape = NA, linewidth = .3, fill = "white") +
    scale_fill_manual(values = GRP, guide = "none") +
    labs(title = qcl[i], x = NULL, y = NULL) +
    theme_sd() + theme(axis.text.x = element_text(angle = 25, hjust = 1))), nrow = 1)

dep <- data.frame(Group = factor(GR, levels = GR),
                  MeanReads = b4$Mean_reads_per_cell[match(GR, b4$Group)],
                  MedianUMI = b1$Median_UMI[match(GR, b1$Group)])
mkbar <- function(y, ttl) ggplot(dep, aes(Group, .data[[y]], fill = Group)) +
  geom_col(width = .65) +
  geom_text(aes(label = format(.data[[y]], big.mark = ",")), vjust = -.4, size = 2) +
  scale_fill_manual(values = GRP, guide = "none") +
  scale_y_continuous(expand = expansion(c(0,.15))) +
  labs(title = ttl, x = NULL, y = NULL) +
  theme_sd() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
f3b <- mkbar("MeanReads","Mean reads per cell") |
       mkbar("MedianUMI","Median UMIs per nucleus (2.6-fold range)")

f3c <- ggplot(md, aes(nCount_RNA, nFeature_RNA)) +
  geom_point(size = .15, alpha = .25, colour = "#0072B2", stroke = 0) +
  geom_text(data = b6, aes(x = Inf, y = -Inf,
            label = paste0("r=", sprintf("%.2f", r_nFeature_nCount))),
            hjust = 1.1, vjust = -1, size = 1.9, inherit.aes = FALSE) +
  facet_wrap(~ orig.ident) + scale_x_log10() +
  labs(title = "Genes vs UMIs per nucleus (Pearson r)",
       x = "UMIs per nucleus (log10)", y = "Genes per nucleus") +
  theme_sd()

qm <- data.frame(Group = factor(rep(GR, 4), levels = GR),
  Metric = rep(c("Q30 barcode","Q30 RNA read","Reads mapped to genome",
                 "Fraction reads in cells"), each = 3),
  Value = c(b4$Q30_barcode_pct, b4$Q30_RNA_read_pct,
            b4$Reads_mapped_genome_pct, b4$Fraction_reads_in_cells_pct))
f3d <- ggplot(qm, aes(Metric, Value, fill = Group)) +
  geom_col(position = position_dodge(.8), width = .72) +
  scale_fill_manual(values = GRP) +
  coord_cartesian(ylim = c(80,100)) +
  labs(title = "Sequencing and alignment quality", x = NULL, y = "%") +
  theme_sd() + theme(axis.text.x = element_text(angle = 20, hjust = 1),
                     legend.position = "bottom")

fig3 <- (f3a / f3b / f3c / f3d) +
  plot_layout(heights = c(1, .9, .95, 1)) +
  plot_annotation(tag_levels = "a", title = "Figure 3. Sequencing and data quality",
                  theme = theme(plot.title = element_text(size = 9, face = "bold")))
ggsave(file.path(FIG, "Figure3_data_quality.pdf"), fig3,
       width = 175, height = 215, units = "mm")

# ================================================================
# Figure 4 — Doublet and ambient RNA assessment（手动绘图，着色可控）
# ================================================================
ud <- umap[order(umap$scDblFinder_class == "doublet"), ]   # singlet 先、doublet 后(上层)
f4a <- ggplot(ud, aes(U1, U2, colour = scDblFinder_class)) +
  geom_point(size = .22, stroke = 0) +
  scale_colour_manual(values = c(singlet = "grey82", doublet = "#BB3E03"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.8))) +
  labs(title = "scDblFinder doublet classification", x = "UMAP 1", y = "UMAP 2") +
  theme_sd() + theme(aspect.ratio = 1)

ua <- umap[order(umap$decontX_contamination), ]
f4b <- ggplot(ua, aes(U1, U2, colour = decontX_contamination)) +
  geom_point(size = .22, stroke = 0) +
  scale_colour_viridis_c("Ambient\nfraction", option = "viridis") +
  labs(title = "decontX ambient RNA contamination", x = "UMAP 1", y = "UMAP 2") +
  theme_sd() + theme(aspect.ratio = 1)

dbg <- data.frame(Group = factor(GR, levels = GR),
  rate = sapply(GR, function(g) 100*mean(md$scDblFinder_class[md$orig.ident==g]=="doublet")))
f4c1 <- ggplot(dbg, aes(Group, rate, fill = Group)) +
  geom_col(width = .65) + geom_text(aes(label = sprintf("%.1f%%", rate)),
            vjust = -.4, size = 2) +
  scale_fill_manual(values = GRP, guide = "none") +
  scale_y_continuous(expand = expansion(c(0,.18))) +
  labs(title = "Doublet rate per library", x = NULL, y = "%") +
  theme_sd() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
dbc <- data.frame(CellType = factor(CT, levels = CT),
  rate = sapply(CT, function(g) 100*mean(md$scDblFinder_class[md$novogene_celltype==g]=="doublet")))
f4c2 <- ggplot(dbc, aes(CellType, rate, fill = CellType)) +
  geom_col(width = .7) +
  scale_fill_manual(values = CTP, guide = "none") +
  scale_y_continuous(expand = expansion(c(0,.1))) +
  labs(title = "Doublet rate by cell type", x = NULL, y = "%") +
  theme_sd() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
f4c <- f4c1 | f4c2

f4d <- ggplot(md, aes(novogene_celltype, decontX_contamination,
                      fill = novogene_celltype)) +
  geom_violin(scale = "width", linewidth = .3, colour = "grey30") +
  geom_boxplot(width = .12, outlier.shape = NA, linewidth = .3, fill = "white") +
  scale_fill_manual(values = CTP, guide = "none") +
  labs(title = "Ambient contamination by cell type", x = NULL,
       y = "decontX contamination") +
  theme_sd() + theme(axis.text.x = element_text(angle = 35, hjust = 1))

fig4 <- ((f4a | f4b) / f4c / f4d) +
  plot_layout(heights = c(1, .72, .72)) +
  plot_annotation(tag_levels = "a",
                  title = "Figure 4. Doublet and ambient RNA assessment",
                  theme = theme(plot.title = element_text(size = 9, face = "bold")))
ggsave(file.path(FIG, "Figure4_doublet_ambient.pdf"), fig4,
       width = 175, height = 190, units = "mm")

# ================================================================
# Supplementary
# ================================================================
ts <- as.data.frame(Embeddings(seu, "tsne")); colnames(ts) <- c("T1","T2")
ts$ct <- seu$novogene_celltype
fs1 <- ggplot(ts[sample(nrow(ts)), ], aes(T1, T2, colour = ct)) +
  geom_point(size = .22, stroke = 0) +
  scale_colour_manual(values = CTP, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.8))) +
  labs(title = "Supplementary Figure S1. t-SNE embedding (cell types)",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_sd() + theme(aspect.ratio = 1)
ggsave(file.path(FIG, "FigureS1_tsne.pdf"), fs1, width = 120, height = 105,
       units = "mm")

fs2 <- FeaturePlot(seu, mk_flat, reduction = "umap", raster = FALSE,
                   pt.size = 0.25, ncol = 5, order = TRUE,
                   cols = c("grey88","#08306B")) &
  theme_sd() & umap_blank &
  theme(plot.title = element_text(size = 7, face = "italic"), legend.position = "none")
fs2 <- fs2 + plot_annotation(title = "Supplementary Figure S2. All canonical markers",
                             theme = theme(plot.title = element_text(size = 9, face = "bold")))
ggsave(file.path(FIG, "FigureS2_all_markers.pdf"), fs2, width = 175, height = 165,
       units = "mm")

cat("\n== 发表级图集完成 (rev2) ==\n")
for (f in list.files(FIG, pattern = "pdf$", full.names = TRUE)) cat("  ", f, "\n")
