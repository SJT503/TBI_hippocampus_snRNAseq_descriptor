# 16_check_novogene_software.py — 探查诺禾报告里的软件版本
# 子项目 2 · 2026-05-21

import sys, os, re, zipfile, shutil
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# 1. software_list.xls (实为 HTML 伪装)
xls = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\novogene_data\单细胞核测序\数据\Report-X101SC24097123-Z01-J001-B1-1_20241019005301 (1)\Report-X101SC24097123-Z01-J001-B1-1\src\images\software_list.xls"
print("=== software_list.xls 原始内容 ===")
try:
    raw = open(xls, 'rb').read()
    txt = raw.decode('utf-8', errors='replace')
    print(txt[:1500])
except Exception as e:
    print("read fail:", e)

# 2. Readme.pdf in Result.zip
print("\n=== 抽 Result.zip 中的 Readme.pdf ===")
ZIP = r"F:\单细胞核测序 原始数据（朱丽红）\X101SC24097123-Z01-J001-B1-1_10X_release_20241019\Result-X101SC24097123-Z01-J001-B1-1.20241019.zip"
OUT_DIR = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\SciData_descriptor\refs"
os.makedirs(OUT_DIR, exist_ok=True)
target = os.path.join(OUT_DIR, "novogene_Readme.pdf")
with zipfile.ZipFile(ZIP, 'r') as z:
    members = [m for m in z.namelist() if m.endswith('Readme.pdf')]
    print("Readme.pdf in zip:", members)
    if members:
        with z.open(members[0]) as fin, open(target, 'wb') as fout:
            shutil.copyfileobj(fin, fout, length=1024*1024)
        print(f"Saved: {target} ({os.path.getsize(target)} bytes)")
