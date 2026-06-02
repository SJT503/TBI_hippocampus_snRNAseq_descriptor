# 17_md_to_docx.R — 把投稿干净版 MD 转 Sci Data 标准 .docx
# 子项目 2 · 2026-05-25 · 用 rmarkdown 自带 pandoc

suppressPackageStartupMessages(library(rmarkdown))
MD  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/SciData_descriptor/manuscript/data_descriptor_submission.md"
DOCX<- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/SciData_descriptor/manuscript/data_descriptor_submission.docx"

rmarkdown::pandoc_convert(
  input     = MD,
  to        = "docx",
  output    = DOCX,
  options   = c(
    "--from=markdown+pipe_tables+fenced_code_blocks",
    "--reference-doc=" # let pandoc use default
  )
)
cat("Output:", DOCX, "\n")
cat("Size:", file.info(DOCX)$size, "bytes\n")
