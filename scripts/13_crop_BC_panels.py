# 13_crop_BC_panels.py — 裁剪图 4-4 B+C 子图（剔除 A 质粒图）
# 输入：_thesis_page075_full.png（2894x4093, 350 DPI）
# 输出：Supplementary Figure S5 占位

from PIL import Image
import os

IN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\SciData_descriptor\results\figures_publication\_thesis_page075_full.png"
OUT = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\SciData_descriptor\results\figures_publication"

img = Image.open(IN)
W, H = img.size
print(f"Source: {W}x{H}")

# 估计裁剪框（基于 2894x4093 全页布局）
# 图 4-4 在页面下半，A 左 / B 右上 / C 右下
# B+C 占右侧约 x = 0.50W ~ 0.93W，y = 0.66H ~ 0.87H
left   = int(W * 0.49)   # 紧贴 B 标签左侧，排除质粒图右缘
top    = int(H * 0.625)  # 上扩抓 "B" 面板标签
right  = int(W * 0.94)
bottom = int(H * 0.855)  # 上移排除中文 caption
print(f"Crop box: left={left} top={top} right={right} bottom={bottom}")

crop = img.crop((left, top, right, bottom))
print(f"Cropped: {crop.size}")

# 高分辨率 PNG（占位用）+ PDF（投稿用）
png_path = os.path.join(OUT, "FigureS5_invitro_TGFbRI_KO_validation.png")
pdf_path = os.path.join(OUT, "FigureS5_invitro_TGFbRI_KO_validation.pdf")
crop.save(png_path, optimize=True, dpi=(350, 350))
crop.save(pdf_path, "PDF", resolution=350.0)
print(f"Saved: {png_path}")
print(f"Saved: {pdf_path}")

# 清理 v1 占位
import os as _os
old = os.path.join(OUT, "_supp_S5_crop_v1.png")
if _os.path.exists(old): _os.remove(old); print("Removed v1 placeholder")
