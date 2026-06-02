# 09_figS3_Tgfbr1.R — Supplementary Figure S3：Tgfbr1 跨细胞类型×组
# 子项目 2 · 2026-05-18 · 据 C1_Tgfbr1_by_celltype_group.csv（08 脚本产出）
# 诚实呈现：KO 库 Tgfbr1 在所有细胞类型都偏低 → 与深度差混杂

suppressPackageStartupMessages(library(ggplot2))
PROJ <- normalizePath("..")
d <- read.csv(file.path(PROJ, "results/qc_tables/C1_Tgfbr1_by_celltype_group.csv"))
CT <- c("Excitatory","Inhibitory","Astrocytes","Oligodendrocytes","OPCs","Microglia","Pericyte")
GR <- c("Sham","CCI_Scr","CCI_TGFbetaR1_KO")
d$CellType <- factor(d$CellType, levels = CT)
d$Group    <- factor(d$Group, levels = GR)

p <- ggplot(d, aes(CellType, mean_Tgfbr1, fill = Group)) +
  geom_col(position = position_dodge(.8), width = .72) +
  scale_fill_manual(values = c(Sham = "#5B8FA8", CCI_Scr = "#E69F00",
                               CCI_TGFbetaR1_KO = "#BB3E03"), name = NULL) +
  labs(title = "Supplementary Figure S3. Tgfbr1 expression by cell type and condition",
       x = NULL, y = "Mean expression (log-normalised)") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(PROJ, "results/figures_publication/FigureS3_Tgfbr1_by_celltype.pdf"),
       p, width = 130, height = 85, units = "mm")
cat("FigS3 saved OK\n")
