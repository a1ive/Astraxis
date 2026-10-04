# 星表来源

`bsc5.csv` 来自 Yale Bright Star Catalogue 第 5 修订版（Hoffleit & Warren 1991），取自 CDS 目录 V/50（<https://cdsarc.cds.unistra.fr/ftp/V/50/>，文件 `ReadMe` 和 `catalog.gz`），由 `tools/stars/convert_bsc5.py` 转换而来。

- 共 9096 颗星（原表中没有 J2000 坐标的条目已跳过）
- 坐标为 J2000，历元 2000.0，未做自行改正
- 列：`hr,ra_deg,dec_deg,vmag,bv,name`
- 该星表是 NASA/CDS 公开分发的天文数据

`m4.csv`：球状星团 M4（NGC 6121）的成员星，供 PSR B1620-26 场景的天空使用，由 `tools/stars/make_m4.py` 生成（会联网访问 VizieR 和 Pecaut & Mamajek 的表）。

- 来源：Gaia DR3（ESA/Gaia/DPAC），CDS 目录 I/355/gaiadr3，经 VizieR 查询
- **许可证：CC BY-NC 3.0 IGO**（<https://www.cosmos.esa.int/web/gaia-users/license>），署名 ESA/Gaia/DPAC，仅限非商业用途。这比仓库里其他数据（公有领域）严格
- 选星：距星团中心 20′ 以内、G < 19、自行与星团相差 1.5 mas/yr 以内、视差与 1.85 kpc 相差 3σ 以内，共 13516 颗
- 处理：消光按 A_V = 1.31（Sigurdsson et al. 2003）和 Wang & Chen 2019 的 Gaia 波段系数扣除；V 和 B−V 由 Bp−Rp 按 Pecaut & Mamajek 2013 的矮星表换算（巨星也用这个表，是近似）；绝对星等用距离 1.85 kpc（Baumgardt & Vasiliev 2021）
- 天球上的位置是 Gaia 实测的，视线方向的深度没有测量，按 King（1962）密度分布（Harris 2010：c = 1.65，r_c = 1.16′）随机抽取（固定种子），所以三维结构只是统计上的
- 列：`east_pc,north_pc,away_pc,abs_vmag,bv`（相对星团中心，中心处的切平面坐标）
- 缺失：G > 19（约 0.65 M☉ 以下）的暗星，以及核心里因为拥挤而 Gaia 没测到的星
