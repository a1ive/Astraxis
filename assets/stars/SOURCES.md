# 星表来源

## 亮星表

`bsc5.csv` 来自 Yale Bright Star Catalogue 第 5 修订版（Hoffleit & Warren 1991），取自 CDS 目录 V/50（<https://cdsarc.cds.unistra.fr/ftp/V/50/>，文件 `ReadMe` 和 `catalog.gz`），由 `tools/stars/convert_bsc5.py` 转换而来。

- 共 9096 颗星（原表中没有 J2000 坐标的条目已跳过）。
- 坐标为 J2000，历元 2000.0，未做自行改正。
- 列：`hr,ra_deg,dec_deg,vmag,bv,name`
- 该星表是 NASA/CDS 公开分发的天文数据。

## 球状星团成员星

`m4.csv` 保存球状星团 M4（NGC 6121）的成员星，供 PSR B1620-26 场景绘制天空。`tools/stars/make_m4.py` 联网查询 VizieR 和 Pecaut & Mamajek 的表，生成该文件。

- 来源：Gaia DR3（ESA/Gaia/DPAC），CDS 目录 I/355/gaiadr3，经 VizieR 查询。
- 许可证：CC BY-NC 3.0 IGO（<https://www.cosmos.esa.int/web/gaia-users/license>），署名 ESA/Gaia/DPAC，使用范围为非商业用途，限制比仓库中其他公有领域数据更严格。
- 选星：距星团中心 20′ 以内、G < 19、自行与星团相差 1.5 mas/yr 以内、视差与 1.85 kpc 相差 3σ 以内，共 13516 颗。
- 处理：按 A_V = 1.31（Sigurdsson et al. 2003）和 Wang & Chen 2019 的 Gaia 波段系数扣除消光。根据 Bp−Rp，使用 Pecaut & Mamajek 2013 的矮星表换算 V 和 B−V；巨星也使用该表，因此结果是近似值。绝对星等按距离 1.85 kpc（Baumgardt & Vasiliev 2021）计算。
- 天球位置取自 Gaia 实测数据。视线方向的深度未测得，按 King（1962）密度分布随机抽取，参数取自 Harris 2010：c = 1.65，r_c = 1.16′。随机种子固定，三维结构表示统计分布。
- 列：`east_pc,north_pc,away_pc,abs_vmag,bv`（相对星团中心，中心处的切平面坐标）。
- 缺失：G > 19（约 0.65 M☉ 以下）的暗星，以及核心里因为拥挤而 Gaia 没测到的星。
