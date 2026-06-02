# 08_critical_checks.R — 批判性自审两项真实分析
# 子项目 2 · 2026-05-18 · 经 00_gates.R 门禁读 atlas
# 项1：KO 臂星形胶质 Tgfbr1 表达验证（敲除是否在数据中可见）
# 项4：Pericyte 簇纯度复查（是否血管混合体）
# 纯描述性，不做组间定量统计声明（venue 红线；n=1 仅描述）

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})
source("00_gates.R")
PROJ <- normalizePath("..")
TAB  <- file.path(PROJ, "results", "qc_tables")
FIG  <- file.path(PROJ, "results", "figures_publication")

seu <- load_definitive_rds("atlas")
DefaultAssay(seu) <- "RNA"
CT <- c("Excitatory","Inhibitory","Astrocytes","Oligodendrocytes","OPCs","Microglia","Pericyte")
GR <- c("Sham","CCI_Scr","CCI_TGFbetaR1_KO")
seu$novogene_celltype <- factor(seu$novogene_celltype, levels = CT)
seu$orig.ident        <- factor(seu$orig.ident, levels = GR)
expr <- GetAssayData(seu, layer = "data", assay = "RNA")
md   <- seu@meta.data

cat("\n###################  项1：Tgfbr1 敲除验证  ###################\n")
tg <- grep("^Tgfbr1$", rownames(seu), value = TRUE, ignore.case = TRUE)
if (length(tg) == 0) stop("Tgfbr1 不在对象中——核对基因名")
cat("基因名:", tg, "\n")
tgv <- expr[tg, ]

# 1a. 各细胞类型 × 组：Tgfbr1 均值 + 阳性率
t1 <- do.call(rbind, lapply(CT, function(ct) {
  do.call(rbind, lapply(GR, function(g) {
    idx <- md$novogene_celltype == ct & md$orig.ident == g
    data.frame(CellType = ct, Group = g, n = sum(idx),
               mean_Tgfbr1 = round(mean(tgv[idx]), 4),
               pct_pos = round(100 * mean(tgv[idx] > 0), 2))
  }))
}))
write.csv(t1, file.path(TAB, "C1_Tgfbr1_by_celltype_group.csv"), row.names = FALSE)
cat("\n## Tgfbr1 表达（细胞类型 × 组）\n"); print(t1, row.names = FALSE)

# 1b. 星形胶质聚焦：KO vs Scr vs Sham
cat("\n## 星形胶质 Tgfbr1（敲除靶细胞）\n")
ast <- t1[t1$CellType == "Astrocytes", ]
print(ast, row.names = FALSE)
cat(sprintf("\n解读：KO 星形 mean=%.4f vs Scr=%.4f vs Sham=%.4f；",
            ast$mean_Tgfbr1[ast$Group=="CCI_TGFbetaR1_KO"],
            ast$mean_Tgfbr1[ast$Group=="CCI_Scr"],
            ast$mean_Tgfbr1[ast$Group=="Sham"]))
cat(" 阳性率 KO=%.1f%% Scr=%.1f%%\n")

# 1c. 特异性：Tgfbr1 在非星形细胞是否随组变化（GFAP 启动子应只敲星形）
cat("\n## astrocyte-specific 检查：各细胞类型 KO/Scr Tgfbr1 均值比\n")
spec <- do.call(rbind, lapply(CT, function(ct) {
  ko  <- t1$mean_Tgfbr1[t1$CellType==ct & t1$Group=="CCI_TGFbetaR1_KO"]
  scr <- t1$mean_Tgfbr1[t1$CellType==ct & t1$Group=="CCI_Scr"]
  data.frame(CellType = ct, KO = ko, Scr = scr,
             KO_vs_Scr_ratio = round(ifelse(scr>0, ko/scr, NA), 3))
}))
print(spec, row.names = FALSE)

# 1d. 图：星形胶质 Tgfbr1 小提琴（三组）
ast_cells <- colnames(seu)[md$novogene_celltype == "Astrocytes"]
dfa <- data.frame(group = md[ast_cells, "orig.ident"], Tgfbr1 = tgv[ast_cells])
pC1 <- ggplot(dfa, aes(group, Tgfbr1, fill = group)) +
  geom_violin(scale = "width", linewidth = .3, colour = "grey30") +
  geom_boxplot(width = .12, outlier.shape = NA, linewidth = .3, fill = "white") +
  scale_fill_manual(values = c(Sham="#5B8FA8",CCI_Scr="#E69F00",
                               CCI_TGFbetaR1_KO="#BB3E03"), guide = "none") +
  labs(title = expression(italic("Tgfbr1")~"in astrocytes by condition"),
       x = NULL, y = "Expression (log-normalised)") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(FIG, "C1_Tgfbr1_astrocyte_validation.pdf"), pC1,
       width = 80, height = 70, units = "mm")

cat("\n\n###################  项4：Pericyte 簇纯度复查  ###################\n")
peri <- colnames(seu)[md$novogene_celltype == "Pericyte"]
cat("Pericyte 细胞数:", length(peri), "\n")
vasc <- list(
  Pericyte    = c("Pdgfrb","Rgs5","Vtn","Kcnj8","Notch3","Anpep"),
  Endothelial = c("Cldn5","Pecam1","Flt1","Slco1a4","Ly6c1"),
  VLMC_Fibro  = c("Dcn","Col1a1","Col1a2","Lum","Pdgfra"),
  SMC         = c("Acta2","Myh11","Tagln"))
cat("\n## Pericyte 簇内各血管标志物表达（mean / %阳性）\n")
p4 <- do.call(rbind, lapply(names(vasc), function(cat0) {
  gs <- intersect(vasc[[cat0]], rownames(expr))
  do.call(rbind, lapply(gs, function(g) {
    v <- expr[g, peri]
    data.frame(MarkerClass = cat0, Gene = g,
               mean_inPericyte = round(mean(v), 4),
               pct_pos = round(100 * mean(v > 0), 1))
  }))
}))
print(p4, row.names = FALSE)
write.csv(p4, file.path(TAB, "C4_pericyte_purity.csv"), row.names = FALSE)

# 4b. 每个 Pericyte 细胞主导身份（按四类标志物平均表达 argmax）
score <- sapply(names(vasc), function(cat0) {
  gs <- intersect(vasc[[cat0]], rownames(expr))
  colMeans(expr[gs, peri, drop = FALSE])
})
dom <- colnames(score)[max.col(score, ties.method = "first")]
cat("\n## Pericyte 簇内细胞按主导血管标志物分类\n")
print(table(dom))
cat(sprintf("→ 纯 Pericyte 主导占比: %.1f%%\n", 100*mean(dom=="Pericyte")))

cat("\n== 08 完成 ==\n")
cat("表:", file.path(TAB, c("C1_Tgfbr1_by_celltype_group.csv","C4_pericyte_purity.csv")),
    sep="\n  ")
cat("\n图:", file.path(FIG, "C1_Tgfbr1_astrocyte_validation.pdf"), "\n")
