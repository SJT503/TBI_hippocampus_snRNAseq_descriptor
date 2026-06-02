# 05_B4_seq_metrics.R — B4 测序与比对质量指标（取自诺禾结题报告）
# 子项目 2 · Technical Validation · 2026-05-18
# 数据源：诺禾 Cell Ranger metrics_summary（all_metrics_summary.txt），非重算
# 同时回填 B1_per_group_qc.csv 的 pre-QC 核数

PROJ <- normalizePath("..")
TAB  <- file.path(PROJ, "results", "qc_tables")

REPORT_TXT <- file.path(
  "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/novogene_data",
  "单细胞核测序/数据",
  "Report-X101SC24097123-Z01-J001-B1-1_20241019005301 (1)",
  "Report-X101SC24097123-Z01-J001-B1-1/src/file/all_metrics_summary.txt")
if (!file.exists(REPORT_TXT)) stop("诺禾报告不存在: ", REPORT_TXT)

raw <- read.delim(REPORT_TXT, check.names = FALSE, stringsAsFactors = FALSE)
raw <- raw[raw$Sample_Name != "" & !is.na(raw$Sample_Name), ]
cat("## 诺禾 Cell Ranger 原始指标\n"); print(raw)

# ---- 数值清洗（去逗号/百分号）----
num <- function(x) as.numeric(gsub(",", "", gsub("%", "", x)))

b4 <- data.frame(
  Group                       = raw$Sample_Name,
  Estimated_cells_preQC       = as.integer(num(raw$`Estimated.Number.of.Cells`)),
  Number_of_reads             = as.integer(num(raw$`Number.of.Reads`)),
  Mean_reads_per_cell         = as.integer(num(raw$`Mean.Reads.per.Cell`)),
  Valid_barcodes_pct          = num(raw$`Valid.Barcodes`),
  Sequencing_saturation_pct   = num(raw$`Sequencing.Saturation`),
  Q30_barcode_pct             = num(raw$`Q30.Bases.in.Barcode`),
  Q30_RNA_read_pct            = num(raw$`Q30.Bases.in.RNA.Read`),
  Q30_UMI_pct                 = num(raw$`Q30.Bases.in.UMI`),
  Reads_mapped_genome_pct     = num(raw$`Reads.Mapped.to.Genome`),
  Reads_conf_genome_pct       = num(raw$`Reads.Mapped Confidently.to.Genome`),
  Reads_conf_intronic_pct     = num(raw$`Reads.Mapped.Confidently.to.Intronic.Regions`),
  Reads_conf_exonic_pct       = num(raw$`Reads.Mapped.Confidently.to.Exonic.Regions`),
  Reads_conf_transcriptome_pct= num(raw$`Reads.Mapped.Confidently.to.Transcriptome`),
  Fraction_reads_in_cells_pct = num(raw$`Fraction.Reads.in.Cells`),
  Total_genes_detected        = as.integer(num(raw$`Total Genes.Detected`)),
  Median_genes_per_cell_CR    = as.integer(num(raw$`Median.Genes.per Cell`)),
  Median_UMI_per_cell_CR      = as.integer(num(raw$`Median.UMI.Counts.per.Cell`)),
  check.names = FALSE)
write.csv(b4, file.path(TAB, "B4_sequencing_metrics.csv"), row.names = FALSE)
cat("\n## B4 测序/比对质量指标表\n"); print(b4, row.names = FALSE)

# ---- 回填 B1 pre-QC 核数 + QC 保留率 ----
b1f <- file.path(TAB, "B1_per_group_qc.csv")
if (file.exists(b1f)) {
  b1 <- read.csv(b1f, check.names = FALSE)
  pre <- setNames(b4$Estimated_cells_preQC, b4$Group)
  pre["All"] <- sum(b4$Estimated_cells_preQC)
  b1$Nuclei_preQC <- as.integer(pre[b1$Group])
  b1$QC_retention_pct <- round(100 * b1$Nuclei_postQC / b1$Nuclei_preQC, 2)
  write.csv(b1, b1f, row.names = FALSE)
  cat("\n## B1 已回填 pre-QC 核数 + QC 保留率\n")
  print(b1[, c("Group","Nuclei_preQC","Nuclei_postQC","QC_retention_pct")],
        row.names = FALSE)
}

cat("\n== B4 完成 ==\n")
cat("表:", file.path(TAB, "B4_sequencing_metrics.csv"), "\n")
cat("    ", b1f, "(已回填 pre-QC)\n")
