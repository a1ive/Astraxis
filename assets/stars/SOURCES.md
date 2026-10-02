# 星表来源

`bsc5.csv` 来自 Yale Bright Star Catalogue 第 5 修订版（Hoffleit & Warren 1991），取自 CDS 目录 V/50（<https://cdsarc.cds.unistra.fr/ftp/V/50/>，文件 `ReadMe` 和 `catalog.gz`），由 `tools/stars/convert_bsc5.py` 转换而来。

- 共 9096 颗星（原表中没有 J2000 坐标的条目已跳过）
- 坐标为 J2000，历元 2000.0，未做自行改正
- 列：`hr,ra_deg,dec_deg,vmag,bv,name`
- 该星表是 NASA/CDS 公开分发的天文数据
