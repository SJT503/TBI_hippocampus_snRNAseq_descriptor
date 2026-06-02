# 00_gates.R — 子项目 2 门禁函数库
# 生信分析门禁铁律落地。下游脚本禁止直接 readRDS()，必须调 load_definitive_rds()。
# 建立 2026-05-17。

suppressPackageStartupMessages({
  library(Seurat)
})

# 数据物理路径前缀（从 SciData_descriptor/ 出发；脚本以本目录为工作目录运行）
.SEURAT_DIR <- normalizePath(file.path("..", "..", "wet_lab_data", "results", "seurat"),
                             mustWork = FALSE)

# 唯一 RDS 对照表（PROJECT_CONSTITUTION.md 第一章）
.DEFINITIVE_RDS <- list(
  atlas       = "novogene_atlas_annotated.rds",
  astrocyte   = "astrocyte_v2_unified.rds",
  excitatory  = "excitatory_v2_unified.rds",
  inhibitory  = "inhibitory_v2_unified.rds",
  microglia   = "microglia_v2_unified.rds",
  dg_granule  = "DG_granule_subclustered.rds"
)

# 焊死参数（PROJECT_CONSTITUTION.md 第二章）
.EXPECTED <- list(
  atlas_total   = 35047,
  n_groups      = 3,
  group_names   = c("Sham", "CCI_Scr", "CCI_TGFbetaR1_KO"),  # 已核实 (01_inspect_atlas.R)
  n_celltypes   = 7
)

# CLAUDE.md 标志物对照表（纯度门禁）
.PURITY_MARKERS <- list(
  astrocyte  = list(pos = c("Aqp4", "Gfap"), neg = "Snap25"),
  microglia  = list(pos = c("Hexb", "Csf1r"), neg = "Snap25"),
  excitatory = list(pos = c("Snap25", "Syt1"), neg = "Aqp4"),
  inhibitory = list(pos = c("Gad1", "Gad2"),  neg = "Aqp4")
)

# 纯度检验：pos 标志物均值须 >> neg 标志物均值（比值 >= 3）
.check_purity <- function(seu, celltype) {
  mk <- .PURITY_MARKERS[[celltype]]
  if (is.null(mk)) return(invisible(TRUE))  # atlas 等无单一纯度标准，跳过
  expr <- GetAssayData(seu, layer = "data", assay = "RNA")
  have <- intersect(c(mk$pos, mk$neg), rownames(expr))
  if (!all(c(mk$pos, mk$neg) %in% have))
    stop(sprintf("PURITY FAIL [%s]: 标志物缺失 %s",
                 celltype, paste(setdiff(c(mk$pos, mk$neg), have), collapse = ",")))
  pos_mean <- mean(colMeans(expr[mk$pos, , drop = FALSE]))
  neg_mean <- mean(expr[mk$neg, ])
  if (pos_mean < neg_mean * 3)
    stop(sprintf("PURITY FAIL [%s]: pos(%.3f) 未 >> neg(%.3f)。是否读错 RDS？",
                 celltype, pos_mean, neg_mean))
  message(sprintf("  ✓ 纯度 [%s]: pos=%.3f neg=%.3f (%.1fx)",
                  celltype, pos_mean, neg_mean, pos_mean / neg_mean))
  invisible(TRUE)
}

#' 加载定稿 RDS（唯一允许的数据入口）
#' @param key   .DEFINITIVE_RDS 的键：atlas/astrocyte/excitatory/inhibitory/microglia/dg_granule
#' @param group_col  分组列名（默认自动探测）
#' @return Seurat 对象
load_definitive_rds <- function(key, group_col = NULL) {
  if (!key %in% names(.DEFINITIVE_RDS))
    stop(sprintf("非定稿 RDS 键: '%s'。允许: %s",
                 key, paste(names(.DEFINITIVE_RDS), collapse = ", ")))
  fpath <- file.path(.SEURAT_DIR, .DEFINITIVE_RDS[[key]])
  if (!file.exists(fpath))
    stop(sprintf("RDS 不存在: %s", fpath))

  # 检验 5 — 数据来源：路径必须在 wet_lab_data/results/seurat 下
  if (!grepl("wet_lab_data[/\\\\]results[/\\\\]seurat", fpath))
    stop(sprintf("SOURCE FAIL: 路径不在 wet_lab_data/results/seurat 下: %s", fpath))

  message(sprintf("[gate] 读取 %s ...", .DEFINITIVE_RDS[[key]]))
  seu <- readRDS(fpath)

  # 检验 1 — 细胞数
  n <- ncol(seu)
  message(sprintf("  ✓ 细胞数: %d", n))
  if (key == "atlas" && n != .EXPECTED$atlas_total)
    stop(sprintf("CELLNUM FAIL: atlas 核数 %d != 焊死值 %d", n, .EXPECTED$atlas_total))

  # 检验 2 — 分组完整
  if (is.null(group_col)) {
    cand <- intersect(c("group", "Group", "orig.ident", "condition", "sample"),
                      colnames(seu@meta.data))
    group_col <- if (length(cand)) cand[1] else NA
  }
  if (!is.na(group_col)) {
    grps <- unique(as.character(seu@meta.data[[group_col]]))
    message(sprintf("  ✓ 分组 [%s]: %s", group_col, paste(grps, collapse = ", ")))
    if (length(grps) != .EXPECTED$n_groups)
      stop(sprintf("GROUP FAIL: 分组数 %d != %d", length(grps), .EXPECTED$n_groups))
  } else {
    warning("未找到分组列，跳过分组检验——请人工确认 metadata 列名")
  }

  # 检验 3 — 注释完整（atlas 须有细胞类型列且 7 类）
  ann_cand <- intersect(c("novogene_celltype", "cell_type", "celltype",
                          "CellType", "annotation"),
                        colnames(seu@meta.data))
  if (key == "atlas") {
    if (!length(ann_cand)) stop("ANNOTATION FAIL: atlas 未找到细胞类型注释列")
    n_ct <- length(unique(seu@meta.data[[ann_cand[1]]]))
    message(sprintf("  ✓ 细胞类型 [%s]: %d 类", ann_cand[1], n_ct))
    if (n_ct != .EXPECTED$n_celltypes)
      stop(sprintf("ANNOTATION FAIL: 细胞类型数 %d != %d", n_ct, .EXPECTED$n_celltypes))
  }

  # 检验 4 — 纯度（单类型对象）
  .check_purity(seu, key)

  message("[gate] 五检验通过。")
  seu
}
