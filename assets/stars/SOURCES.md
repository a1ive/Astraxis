# 星表来源

## 星表（带距离）

`hyg.csv` 来自 HYG 数据库 v4.4（astronexus，<https://codeberg.org/astronexus/hyg>，文件 `data/hyg/CURRENT/hyg_v44.csv.gz`，13,636,362 字节，SHA-256 `00b349893b9a53106dd488d8371e8d2fa586043e500bb3cdb8bff3931682197d`，2026-10-07 下载），使用 `tools/stars/convert_hyg.py` 转换。原文件 34 MB（解压后），不入库。

- 来源：HYG 合并了 Hipparcos、Yale Bright Star Catalogue 和 Gliese 近星表，距离来自 Hipparcos 视差。
- 许可证：CC BY-SA 4.0（<https://creativecommons.org/licenses/by-sa/4.0/>），署名 astronexus / HYG database。`hyg.csv` 是它的衍生数据，同样按 CC BY-SA 4.0 分发。
- 选星：观察者在离太阳 15 pc 以内的任何位置时，可能亮于 V = 6.5 的星都保留（距离 d 的星离这样的观察者至少 d − 15 pc）。没有可用视差的星（HYG 中 dist ≥ 100000）按地球上看到的星等保留，当作无穷远。共 13821 颗。从地球上看亮于 6.5 等的约 8900 颗，与 BSC5 相当；其余是近处的暗星，例如 HIP 112325（地球上看 V = 9.4），它离 TRAPPIST-1 只有 1.3 pc，在那里是 4.5 等。
- 去掉了太阳：需要太阳的场景都把它作为天体来画。
- 坐标 J2000，历元 2000.0，未做自行改正。星等是 Hipparcos 的 V 星等，颜色是 B−V。
- 列：`hip,hr,ra_deg,dec_deg,dist_pc,vmag,bv,name`。`dist_pc` 或 `bv` 为空表示未知；`name` 依次取专名、拜耳/弗兰斯蒂德名、Gliese 编号。文件头的 `# max_viewer_pc = 15` 由程序读取：场景的观察者离太阳更远时拒绝加载。
- 程序用法（`scene/star_catalog.cpp` 中的 `catalog_sky`）：把每颗星放到它的三维位置，从观察者处重新计算方向和星等，只画亮于 6.5 等的星，并去掉离观察者 0.5 pc 以内的星（属于场景自己的系统，例如 α Cen A、B 和比邻星）。

## 球状星团成员星

`m4.csv` 保存球状星团 M4（NGC 6121）的成员星，供 PSR B1620-26 场景绘制天空。`tools/stars/make_m4.py` 联网查询 VizieR 和 Pecaut & Mamajek 的表，生成该文件。

- 来源：Gaia DR3（ESA/Gaia/DPAC），CDS 目录 I/355/gaiadr3，经 VizieR 查询。
- 许可证：CC BY-NC 3.0 IGO（<https://www.cosmos.esa.int/web/gaia-users/license>），署名 ESA/Gaia/DPAC，限非商业用途。仓库中的公有领域数据没有这一使用限制。
- 选星：距星团中心 20′ 以内、G < 19、自行与星团相差 1.5 mas/yr 以内、视差与 1.85 kpc 相差 3σ 以内，共 13516 颗。
- 处理：按 A_V = 1.31（Sigurdsson et al. 2003）和 Wang & Chen 2019 的 Gaia 波段系数扣除消光。根据 Bp−Rp，使用 Pecaut & Mamajek 2013 的矮星表换算 V 和 B−V；巨星也使用该表，因此结果是近似值。绝对星等按距离 1.85 kpc（Baumgardt & Vasiliev 2021）计算。
- 天球位置取自 Gaia 实测数据。视线方向的深度未测得，按 King（1962）密度分布随机抽取，参数取自 Harris 2010：c = 1.65，r_c = 1.16′。随机种子固定，三维结构表示统计分布。
- 列：`east_pc,north_pc,away_pc,abs_vmag,bv`（相对星团中心，中心处的切平面坐标）。
- 缺失：G > 19（约 0.65 M☉ 以下）的暗星，以及星团核心中因拥挤而未被 Gaia 检出的星。
