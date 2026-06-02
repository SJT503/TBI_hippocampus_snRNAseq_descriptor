# TBI Hippocampus snRNA-seq Data Descriptor — Analysis Code

> Companion code repository for: *A single-nucleus RNA sequencing dataset of the
> mouse hippocampus after controlled cortical impact with astrocyte-specific
> TGF-βRI knockout* (Scientific Data, submitted).
>
> Data: GEO accession **GSE333879**
> (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333879).
> Private until 31 May 2027; a reviewer access token is provided in the
> submission cover letter.
>
> Co-corresponding authors: Kangsheng Li (ksli2013@yeah.net); Weiqiang Chen
> (wqchen@stu.edu.cn).
> Data submitter: Jiangtao Sheng (jtsheng@stu.edu.cn).

## What this repository contains

All analysis scripts used to produce the Technical Validation figures and tables
of the Data Descriptor, plus utilities for building the GEO submission package.

The scripts read the released GEO dataset (filtered Cell Ranger output and the
integrated Seurat object) and reproduce every per-library QC table, every
descriptor figure, and the per-nucleus metadata file released with the data.

## Repository layout

```
code_repo/
├── README.md             # this file
├── LICENSE               # MIT
├── scripts/
│   ├── 00_gates.R                       # gate library — definitive RDS loader + 5 entry checks
│   ├── 01_inspect_atlas.R               # one-off atlas structure probe
│   ├── 02_B1_B6_qc_metrics.R            # per-library QC summary + correlation metrics
│   ├── 03_B5_marker_validation.R        # canonical-marker validation of 7 cell types
│   ├── 04_B2_B3_doublet_ambient.R       # scDblFinder + decontX re-assessment
│   ├── 05_B4_seq_metrics.R              # parse Cell Ranger metrics; backfill QC table
│   ├── 07_figures_descriptor.R          # main Figures 1–4 + Suppl S1–S2 (publication-grade)
│   ├── 08_critical_checks.R             # Tgfbr1-by-celltype + Pericyte purity audit
│   ├── 09_figS3_Tgfbr1.R                # Suppl Figure S4 (Tgfbr1 by celltype)
│   ├── 10_inspect_astrocyte.R           # astrocyte_v2_unified structure probe
│   ├── 11_supp_astrocyte_subclusters.R  # Suppl Figure S3 (astrocyte subclusters)
│   ├── 12_crop_fig4_4_BC.py             # (utility) render thesis page
│   ├── 13_crop_BC_panels.py             # (utility) crop image
│   ├── 14_extract_geo_processed.py      # build the GEO submission package from Novogene Result.zip
│   ├── 15_export_pernucleus_metadata.R  # export per_nucleus_metadata.csv for GEO release
│   ├── 16_check_novogene_software.py    # (utility) inspect Novogene software_list
│   └── 17_md_to_docx.R                  # (utility) MD → DOCX for journal submission
```

## Software environment

All analyses ran on Windows 11, R 4.5.3 (Git Bash with Rscript), Python 3.14.

### R packages (Bioconductor 3.21)

| Package      | Version  | Use                                       |
|--------------|----------|-------------------------------------------|
| Seurat       | 4.1.3*   | atlas object (Novogene-produced)          |
| Seurat       | 5.5.0    | re-analysis (load + auxiliary computation) |
| scDblFinder  | 1.24.10  | doublet re-assessment                     |
| decontX      | 1.8.0    | ambient RNA estimation                    |
| SingleR      | 2.12.0   | independent cell-type cross-check         |
| celldex      | 1.20.0   | MouseRNAseqData reference for SingleR     |
| ggplot2      | 3.5+     | figures                                   |
| patchwork    | 1.x      | multi-panel figure composition            |
| reshape2     | 1.4+     | matrix → long format                      |

\* The released Seurat object was created with Seurat v4.1.3 (Novogene
processing); our re-analysis loads it under Seurat v5.5.0 which preserves the
v4-style assays/slots.

### Python packages

| Package        | Version | Use                                   |
|----------------|---------|---------------------------------------|
| PyMuPDF (fitz) | 1.27+   | PDF rendering for figure provenance   |
| python-docx    | 1.2     | docx utilities                        |
| pypandoc       | 1.17    | markdown → docx conversion            |
| Pillow         | 10+     | image cropping                        |

## Reproduce the analysis

1. **Get the data from GEO** (GSE333879). Download the three per-library
   `filtered_feature_bc_matrix/` triplets and the integrated Seurat object
   (released as supplementary file with the Series record). Place them so that
   `scripts/00_gates.R`'s `.SEURAT_DIR` resolves to the directory containing
   the released Seurat object.
2. **Set working directory** to `scripts/`.
3. **Run scripts in numerical order**:
   - `01_inspect_atlas.R` — confirms the gate's expected structure
   - `02_B1_B6_qc_metrics.R` — Tables / Figure 3 inputs
   - `03_B5_marker_validation.R` — Figure 2 marker panels
   - `04_B2_B3_doublet_ambient.R` — Figure 4 inputs
   - `05_B4_seq_metrics.R` — Table 2 (sequencing metrics from Novogene)
   - `07_figures_descriptor.R` — composes Figures 1–4 + Suppl S1–S2
   - `08_critical_checks.R` — Tgfbr1 / Pericyte audits
   - `09_figS3_Tgfbr1.R` — Supplementary Figure S4
   - `11_supp_astrocyte_subclusters.R` — Supplementary Figure S3
   - `15_export_pernucleus_metadata.R` — re-creates `per_nucleus_metadata.csv`

The numbered prefixes are the order in which the scripts were developed; gates
in `00_gates.R` enforce the same five entry checks (cell count, group
completeness, annotation completeness, purity, source-path) for every script
that loads a definitive RDS.

## License

MIT — see `LICENSE`.

## Citation

If you use this code, please cite the Data Descriptor (Sheng *et al.*,
*Scientific Data*, submitted; GEO **GSE333879**) and the companion mechanistic
study (Liu *et al.*, manuscript in preparation).

## Contact

Co-corresponding authors:
- Kangsheng Li — ksli2013@yeah.net
- Weiqiang Chen — wqchen@stu.edu.cn

Data submitter / first author:
- Jiangtao Sheng — jtsheng@stu.edu.cn
