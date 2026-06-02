# 02_B1_B6_qc_metrics.R — B1 每组 QC 总表 + B6 相关性质量指标
# 子项目 2 · Technical Validation · 2026-05-18
# 经 00_gates.R 门禁读 novogene_atlas_annotated.rds（post-QC，35,047 核）
# 产出纯描述性 QC 指标，不做组间定量生物学声明（venue 红线）

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# ---- 路径（脚本从 scripts/ 运行）----
source("00_gates.R")
PROJ <- normalizePath("..")                       # SciData_descriptor/
TAB  <- file.path(PROJ, "results", "qc_tables")
FIG  <- file.path(PROJ, "results", "figures")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ---- 门禁加载 ----
seu <- load_definitive_rds("atlas")
DefaultAssay(seu) <- "RNA"

GROUP_COL <- "orig.ident"            # 已核实 (CONSTITUTION 2.1)
CT_COL    <- "novogene_celltype"
GROUP_LV  <- c("Sham", "CCI_Scr", "CCI_TGFbetaR1_KO")   # 展示顺序
seu@meta.data[[GROUP_COL]] <- factor(seu@meta.data[[GROUP_COL]], levels = GROUP_LV)

# ---- 现算线粒体 / 核糖体比例（percent.mt 未预存）----
mt_genes   <- grep("^mt-",   rownames(seu), value = TRUE, ignore.case = TRUE)
ribo_genes <- grep("^Rp[sl]", rownames(seu), value = TRUE, ignore.case = TRUE)
cat(sprintf("\n[B1] 线粒体基因 %d 个: %s\n", length(mt_genes),
            paste(head(mt_genes, 20), collapse = ", ")))
cat(sprintf("[B1] 核糖体基因 %d 个\n", length(ribo_genes)))
if (length(mt_genes) == 0)
  stop("MT FAIL: 未匹配到 ^mt- 线粒体基因——确认基因命名（atlas 是否已剔除线粒体基因？）")

seu[["percent.mt"]]   <- PercentageFeatureSet(seu, features = mt_genes)
seu[["percent.ribo"]] <- if (length(ribo_genes))
  PercentageFeatureSet(seu, features = ribo_genes) else NA_real_

md <- seu@meta.data

# ============ B1 · 每组 QC 总表 ============
b1 <- do.call(rbind, lapply(c(GROUP_LV, "All"), function(g) {
  sub <- if (g == "All") md else md[md[[GROUP_COL]] == g, ]
  data.frame(
    Group                = g,
    Nuclei_preQC         = NA_integer_,                    # 见诺禾结题报告，待填
    Nuclei_postQC        = nrow(sub),
    Median_genes         = round(median(sub$nFeature_RNA)),
    Mean_genes           = round(mean(sub$nFeature_RNA)),
    Median_UMI           = round(median(sub$nCount_RNA)),
    Mean_UMI             = round(mean(sub$nCount_RNA)),
    Median_pct_mt        = round(median(sub$percent.mt), 3),
    Mean_pct_mt          = round(mean(sub$percent.mt), 3),
    Pct_nuclei_mt_gt5    = round(100 * mean(sub$percent.mt > 5), 2),
    Pct_nuclei_mt_gt10   = round(100 * mean(sub$percent.mt > 10), 2),
    Median_pct_ribo      = round(median(sub$percent.ribo), 3),
    check.names = FALSE
  )
}))
write.csv(b1, file.path(TAB, "B1_per_group_qc.csv"), row.names = FALSE)
cat("\n## B1 每组 QC 总表\n"); print(b1, row.names = FALSE)

# ---- B1 附：组 × 细胞类型 核数交叉表（供 Data Records）----
comp <- as.data.frame.matrix(table(md[[GROUP_COL]], md[[CT_COL]]))
comp <- cbind(Group = rownames(comp), comp, Total = rowSums(comp))
write.csv(comp, file.path(TAB, "B1_celltype_composition.csv"), row.names = FALSE)
cat("\n## B1 附 · 组 × 细胞类型核数\n"); print(comp, row.names = FALSE)

# ============ B6 · 相关性质量指标 ============
b6 <- do.call(rbind, lapply(GROUP_LV, function(g) {
  sub <- md[md[[GROUP_COL]] == g, ]
  data.frame(
    Group              = g,
    n                  = nrow(sub),
    r_nFeature_nCount  = round(cor(sub$nFeature_RNA, sub$nCount_RNA), 3),
    r_pctmt_nCount     = round(cor(sub$percent.mt,  sub$nCount_RNA), 3),
    r_pctmt_nFeature   = round(cor(sub$percent.mt,  sub$nFeature_RNA), 3)
  )
}))
write.csv(b6, file.path(TAB, "B6_correlations.csv"), row.names = FALSE)
cat("\n## B6 相关性指标（Pearson）\n"); print(b6, row.names = FALSE)

# ============ 图 · B1 QC 小提琴 ============
qc_vln <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  group.by = GROUP_COL, pt.size = 0, ncol = 3) &
  theme(axis.title.x = element_blank())
ggsave(file.path(FIG, "B1_qc_violin.pdf"), qc_vln, width = 11, height = 4)

# ============ 图 · B6 散点 ============
p1 <- ggplot(md, aes(nCount_RNA, nFeature_RNA)) +
  geom_point(size = .2, alpha = .3) +
  facet_wrap(~ get(GROUP_COL)) + scale_x_log10() +
  labs(x = "UMI / nucleus (log10)", y = "Genes / nucleus",
       title = "B6 · nFeature vs nCount") + theme_bw()
p2 <- ggplot(md, aes(nCount_RNA, percent.mt)) +
  geom_point(size = .2, alpha = .3) +
  facet_wrap(~ get(GROUP_COL)) + scale_x_log10() +
  labs(x = "UMI / nucleus (log10)", y = "Mitochondrial %",
       title = "B6 · percent.mt vs nCount") + theme_bw()
ggsave(file.path(FIG, "B6_scatter_qc.pdf"),
       patchwork::wrap_plots(p1, p2, ncol = 1), width = 10, height = 8)

cat("\n== B1+B6 完成 ==\n")
cat("表:", file.path(TAB, c("B1_per_group_qc.csv",
                            "B1_celltype_composition.csv",
                            "B6_correlations.csv")), sep = "\n  ")
cat("\n图:", file.path(FIG, c("B1_qc_violin.pdf", "B6_scatter_qc.pdf")), sep = "\n  ")
cat("\n")
