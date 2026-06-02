# 14_extract_geo_processed.py — 从 Result.zip 抽取 GEO 投稿需要的 processed 文件
# 子项目 2 · 2026-05-21
# 抽取每组的 filtered_feature_bc_matrix（barcodes/features/matrix.mtx.gz）+ .h5 + metrics_summary.csv + web_summary.html

import zipfile, os, shutil, sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from pathlib import Path

ZIP = r"F:\单细胞核测序 原始数据（朱丽红）\X101SC24097123-Z01-J001-B1-1_10X_release_20241019\Result-X101SC24097123-Z01-J001-B1-1.20241019.zip"
OUT = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\SciData_descriptor\geo_submission\processed"
SAMPLES = ["Sham", "CCI_Scr", "CCI_TGFbetaR1_KO"]
PREFIX = "Result-X101SC24097123-Z01-J001-B1-1"

# 每样本需要的 4 个相对路径
def members_for(sample):
    base = f"{PREFIX}/2.Summary/{sample}"
    return [
        (f"{base}/filtered_feature_bc_matrix/barcodes.tsv.gz", f"{sample}/filtered_feature_bc_matrix/barcodes.tsv.gz"),
        (f"{base}/filtered_feature_bc_matrix/features.tsv.gz", f"{sample}/filtered_feature_bc_matrix/features.tsv.gz"),
        (f"{base}/filtered_feature_bc_matrix/matrix.mtx.gz",   f"{sample}/filtered_feature_bc_matrix/matrix.mtx.gz"),
        (f"{base}/filtered_feature_bc_matrix.h5",              f"{sample}/{sample}_filtered_feature_bc_matrix.h5"),
        (f"{base}/metrics_summary.csv",                        f"{sample}/{sample}_metrics_summary.csv"),
        (f"{base}/{sample}_web_summary.html",                  f"{sample}/{sample}_web_summary.html"),
    ]

print(f"Opening zip: {ZIP}")
with zipfile.ZipFile(ZIP, 'r') as z:
    all_names = set(z.namelist())
    for s in SAMPLES:
        print(f"\n=== {s} ===")
        for src, dst in members_for(s):
            if src not in all_names:
                print(f"  [MISSING] {src}")
                continue
            target = os.path.join(OUT, dst)
            os.makedirs(os.path.dirname(target), exist_ok=True)
            with z.open(src) as fin, open(target, 'wb') as fout:
                shutil.copyfileobj(fin, fout, length=1024*1024)
            size = os.path.getsize(target)
            print(f"  [OK] {dst}  ({size/1e6:.1f} MB)")

print("\n=== 完成 ===")
print(f"输出目录: {OUT}")
for s in SAMPLES:
    sample_dir = os.path.join(OUT, s)
    total = sum(os.path.getsize(os.path.join(r,f)) for r,_,fs in os.walk(sample_dir) for f in fs)
    print(f"  {s}: {total/1e6:.1f} MB")
