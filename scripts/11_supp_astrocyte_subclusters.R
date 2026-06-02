# 11_supp_astrocyte_subclusters.R — Supplementary Figure S4：星形胶质亚型
# 子项目 2 · 2026-05-18 · A1 升级（对标 PMID 41667502 Fig 3 结构）
# 经门禁读 atlas，atlas 内 subcluster 1038 星形胶质；与 astrocyte_v2_unified 精炼注释交叉验证
# 纯描述性，不做组间定量声明（venue 红线，n=1）

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})
set.seed(20260518)
source("00_gates.R")
PROJ <- normalizePath("..")
FIG  <- file.path(PROJ, "results", "figures_publication")
TAB  <- file.path(PROJ, "results", "qc_tables")

# ---- 主对象：atlas 内 subcluster ----
seu <- load_definitive_rds("atlas")
ast <- subset(seu, subset = novogene_celltype == "Astrocytes")
cat("atlas 内星形胶质 n =", ncol(ast), "\n")
DefaultAssay(ast) <- "integrated"
ast <- ScaleData(ast, verbose = FALSE)
ast <- RunPCA(ast, npcs = 20, verbose = FALSE)
ast <- FindNeighbors(ast, dims = 1:15, verbose = FALSE)
ast <- FindClusters(ast, resolution = 0.2, verbose = FALSE)
ast <- RunUMAP(ast, dims = 1:15, verbose = FALSE)
DefaultAssay(ast) <- "RNA"
cat("subcluster 数:", length(levels(ast$seurat_clusters)), "\n")
print(table(ast$seurat_clusters, ast$orig.ident))

# ---- 拉精炼注释（astrocyte_v2_unified）做交叉验证 ----
av2 <- readRDS(file.path("..","..","wet_lab_data","results","seurat","astrocyte_v2_unified.rds"))
cat("astrocyte_v2_unified n =", ncol(av2), "\n")
av2_labels <- setNames(as.character(av2$cluster_label), colnames(av2))
cl <- unname(av2_labels[colnames(ast)])
cl[is.na(cl)] <- "not_in_v2_object"
ast <- AddMetaData(ast, metadata = cl, col.name = "curated_label")
cat("\n## atlas subcluster vs astrocyte_v2_unified 精炼注释交叉表\n")
xt <- table(atlas_subcluster = ast$seurat_clusters, curated = ast$curated_label)
print(xt)
write.csv(as.data.frame.matrix(xt),
          file.path(TAB, "C_astrocyte_subcluster_vs_curated.csv"))

# ---- 标志物：每亚群 top markers ----
mk_state <- c("Apoe","Id3","Gfap","Slc1a3","Aqp4","Vim","Cd44","Hes5","Sparc",
              "Aldoc","Mt1","Mt2","Igfbp5","Cst3","S100a6")
mk_state <- intersect(mk_state, rownames(ast))
mks <- FindAllMarkers(ast, only.pos = TRUE, min.pct = .25, logfc.threshold = .25,
                      verbose = FALSE)
top3 <- by(mks, mks$cluster, function(d) head(d[order(-d$avg_log2FC), "gene"], 3))
cat("\n## 每亚群 top3 markers (data-driven)\n")
print(top3)

# ---- 构图 ----
nC <- length(levels(ast$seurat_clusters))
ast_pal <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00")[seq_len(nC)]
names(ast_pal) <- levels(ast$seurat_clusters)

# a. 数据驱动 subcluster UMAP
emb <- as.data.frame(Embeddings(ast, "umap")); colnames(emb) <- c("U1","U2")
emb$sub <- ast$seurat_clusters; emb$grp <- ast$orig.ident
emb$cur <- ast$curated_label

pA <- ggplot(emb, aes(U1, U2, colour = sub)) +
  geom_point(size = .5, stroke = 0) +
  scale_colour_manual(values = ast_pal, name = "Subcluster") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(title = sprintf("a  Astrocyte subclusters (n = %d)", nrow(emb)),
       x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank(), aspect.ratio = 1)

# b. 精炼注释叠在同一 UMAP 上
cur_levels <- c("Homeostatic-A","Homeostatic-B(Apoe+)","Reactive(Id3+)",
                "not_in_v2_object")
cur_pal <- c("#0072B2","#56B4E9","#BB3E03","grey75")
names(cur_pal) <- cur_levels
emb$cur <- factor(emb$cur, levels = cur_levels)
pB <- ggplot(emb, aes(U1, U2, colour = cur)) +
  geom_point(size = .5, stroke = 0) +
  scale_colour_manual(values = cur_pal, name = "Curated label\n(astrocyte_v2_unified)") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(title = "b  Curated state annotation (cross-reference)",
       x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank(), aspect.ratio = 1)

# c. marker 气泡图
pC <- DotPlot(ast, features = mk_state, group.by = "seurat_clusters",
              cols = c("grey90","#BB3E03"), dot.scale = 4) +
  RotatedAxis() +
  labs(title = "c  Canonical astrocyte state markers", x = NULL, y = "Subcluster") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 6))

# d. subcluster × group 组成（仅描述性核数）
comp <- as.data.frame(table(emb$grp, emb$sub))
names(comp) <- c("Group","Subcluster","n")
comp$Group <- factor(comp$Group, levels = c("Sham","CCI_Scr","CCI_TGFbetaR1_KO"))
pD <- ggplot(comp, aes(Group, n, fill = Subcluster)) +
  geom_col(position = "fill", width = .7) +
  scale_fill_manual(values = ast_pal) +
  scale_y_continuous("Proportion", expand = c(0,0)) +
  labs(title = "d  Subcluster composition by condition (descriptive)", x = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Supplementary Figure S4. Astrocyte subclusters in the released atlas",
    caption = paste0("Subclusters resolved by re-clustering the 1,038 atlas astrocytes ",
                     "(integrated assay, resolution 0.2); curated labels in (b) are from ",
                     "the more refined astrocyte_v2_unified.rds (",ncol(av2)," nuclei); ",
                     "(d) shows descriptive composition (n = 1 library per condition)."),
    theme = theme(plot.title   = element_text(size = 9, face = "bold"),
                  plot.caption = element_text(size = 6.5, hjust = 0)))
ggsave(file.path(FIG, "FigureS4_astrocyte_subclusters.pdf"), fig,
       width = 175, height = 165, units = "mm")
cat("\n== S4 完成 ==\n")
cat("图:", file.path(FIG, "FigureS4_astrocyte_subclusters.pdf"), "\n")
