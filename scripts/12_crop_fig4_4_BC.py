# 12_crop_fig4_4_BC.py — 渲染朱丽红论文 page 75 + 裁剪图 4-4 (B+C) 子图
# 用作 descriptor Supplementary Figure S5（占位，待朱丽红原图替换）

import fitz  # PyMuPDF
from PIL import Image
import os

PDF = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\wet_lab_data\朱丽红 毕业论文定稿.pdf"
OUT = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\SciData_descriptor\results\figures_publication"
os.makedirs(OUT, exist_ok=True)

DPI = 350
PAGE_INDEX = 74  # 0-indexed; PDF 第 75 页

doc = fitz.open(PDF)
page = doc[PAGE_INDEX]
mat = fitz.Matrix(DPI/72, DPI/72)
pix = page.get_pixmap(matrix=mat)
img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
print(f"Full page size: {img.size}")

# 先存全页供裁剪坐标确认
img.save(os.path.join(OUT, "_thesis_page075_full.png"))
print(f"Saved full page: _thesis_page075_full.png ({img.size[0]}x{img.size[1]})")
