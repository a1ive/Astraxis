# 星表来源

## 星表（带距离）

`hyg.csv` 来自 HYG 数据库 v4.4（astronexus，<https://codeberg.org/astronexus/hyg>，文件 `data/hyg/CURRENT/hyg_v44.csv.gz`，13,636,362 字节，SHA-256 `00b349893b9a53106dd488d8371e8d2fa586043e500bb3cdb8bff3931682197d`，2026-10-07 下载），使用 `tools/stars/convert_hyg.py` 转换。原文件 34 MB（解压后），不入库。

- 来源：HYG 合并了 Hipparcos、Yale Bright Star Catalogue 和 Gliese 近星表，距离来自 Hipparcos 视差。
- 许可证：CC BY-SA 4.0（<https://creativecommons.org/licenses/by-sa/4.0/>），署名 astronexus / HYG database。`hyg.csv` 是它的衍生数据，同样按 CC BY-SA 4.0 分发。
- 选星：观察者在离太阳 30 pc 以内的任何位置时，可能亮于 V = 6.5 的星都保留（距离 d 的星离这样的观察者至少 d − 30 pc）。没有可用视差的星（HYG 中 dist ≥ 100000）按地球上看到的星等保留，当作无穷远。共 21430 颗（2026-10-20 为 Algol 从 15 pc 放宽到 30 pc；15 pc 时 13821 颗）。离太阳 30 pc 处，Hipparcos 以外更暗的星要非常近才能亮到 6.5 等，漏掉的很少。从地球上看亮于 6.5 等的约 8900 颗，与 BSC5 相当；其余是近处的暗星，例如 HIP 112325（地球上看 V = 9.4），它离 TRAPPIST-1 只有 1.3 pc，在那里是 4.5 等。
- 去掉了太阳：需要太阳的场景都把它作为天体来画。
- 坐标 J2000，历元 2000.0，未做自行改正。星等是 Hipparcos 的 V 星等，颜色是 B−V。
- 列：`hip,hr,ra_deg,dec_deg,dist_pc,vmag,bv,name`。`dist_pc` 或 `bv` 为空表示未知；`name` 依次取专名、拜耳/弗兰斯蒂德名、Gliese 编号。文件头的 `# max_viewer_pc = 30` 由程序读取：场景的观察者离太阳更远时拒绝加载。
- 星表里的距离是 Hipparcos 视差，可能和场景采用的距离不同：Algol 在星表里是 27.6 pc（36.27 mas），场景用轨道视差 28.8 pc，星表里的 Algol 离观察者 1.25 pc，超出 `kOwnSystemPc`，会成为一颗 −4.6 等的星。场景用 `[sky] hide_hip = [14576]` 把它去掉。
- 程序用法（`scene/star_catalog.cpp` 中的 `catalog_sky`）：把每颗星放到它的三维位置，从观察者处重新计算方向和星等，只画亮于 6.5 等的星，并去掉离观察者 0.5 pc 以内的星（属于场景自己的系统，例如 α Cen A、B 和比邻星）。

## 造父变星的脉动表

`delta_cep_spips.csv` 是 δ Cephei 一个脉动周期里的半径、有效温度和光度，供 `delta_cephei.toml` 使用。`tools/cepheid/spips_table.py` 只用标准库，从 Mérand et al. 2015（A&A 584, A80）的 SPIPS 模型生成该文件。

- 来源：CDS 目录 J/A+A/584/A80 的 `delcep2.dat`（<https://cdsarc.cds.unistra.fr/ftp/J/A+A/584/A80/>，173,826 字节，SHA-256 `ee9ca65ff91da83767e91c5ece6745450470b56a36ea91b4df6b275924d4ab76`，2026-10-20 下载）。这是论文拟合出的模型（不是观测），按 999 个相位列出 Rosseland 角直径、T_eff 和光度。相位 0 是热光度极大，V 极大在相位 0.995
- 许可证：CDS 分发的已发表研究数据，使用时引用 Mérand et al. 2015 和 VizieR
- 处理：每 5 个相位取 1 个，共 200 个。半径是角直径的一半乘以 274 pc（论文固定的距离，即 Benedict et al. 2002 的 HST 视差），按 IAU 名义太阳半径 695,700 km 换成 R☉。表里的光度和 R²T⁴ 相差约 0.5%（模型的太阳常数不同），场景直接用表里的光度
- 列：`phase,radius_rsun,teff_k,luminosity_lsun`

`polaris_pulsation.csv` 是北极星（Polaris Aa）一个脉动周期里的半径、有效温度和光度，供 `polaris.toml` 使用。北极星没有 SPIPS 那样按相位列出的模型，`tools/cepheid/polaris_table.py`（只用标准库，用到的数都写在脚本里，不需要下载）由论文里的测量拼出：

- 半径：Anderson 2019（A&A 623, A146）表 1 的 Hermes 脉动视向速度（两阶傅里叶级数，2011–2018 年振幅稳定，全振幅 3.526 km/s），乘投影因子 p = 1.29（Mérand 2015 给 δ Cep 的值，北极星没有测过，是假设）后对时间积分，加到 Evans et al. 2024（ApJ 971, 190）用 CHARA 测得的平均半径 46.27 R☉（136.90 pc）上。半振幅 0.39%，和 Evans 引的“约 0.4%”一致。表 1 的图注把 sin 和 cos 写反了：按图注的写法，速度极小落在相位 0.25；CDS 的 Hermes 数据（J/A+A/623/A146）在参考历元 t0 处正好是速度极小，只有 aⱼ 配 cos、bⱼ 配 sin 才对得上（对参考期的数据重新拟合得到 −1.760、−0.093、0.022、0.047，与表 1 相符）
- 有效温度：Usenko et al. 2016（arXiv:1610.03813）表 1 的 21 个谱线深度比温度（2015 年 9–12 月），按 Anderson 的历表（速度极小 BJD 2456553.62553，P = 3.97198 d）定相位，按误差加权拟合正弦：6025 ± 42 K（误差约 8 K），速度极小前 0.076 周最热。他们的周期 3.97987 d 只用了三个月的数据，不用
- 光度：(R/R☉)²(T/5772 K)⁴
- 相位 0 是速度极小。550 nm 黑体的光在相位 0.9635 最亮，即速度极小前 0.145 d（Torres 2023 用 0.21 d），V 振幅 0.061 等
- 脚本还用 Evans 2024 表 2 的 6 个光极大时刻（2016–2021 年，他们按 Torres 2023 图 6 的 O−C 样条算的）拟合二次历表，打印场景的 `epoch_jd_tdb`、`period_days`、`period_change_s_per_year`（−5.12 s/年：周期在 2010 年前后由变长转为变短）
- 许可证：论文发表的数值，引用 Anderson 2019、Evans et al. 2024、Usenko et al. 2016
- 列：`phase,radius_rsun,teff_k,luminosity_lsun`

## 球状星团成员星

`m4.csv` 保存球状星团 M4（NGC 6121）的成员星，供 PSR B1620-26 场景绘制天空。`tools/stars/make_m4.py` 联网查询 VizieR 和 Pecaut & Mamajek 的表，生成该文件。

- 来源：Gaia DR3（ESA/Gaia/DPAC），CDS 目录 I/355/gaiadr3，经 VizieR 查询。
- 许可证：CC BY-NC 3.0 IGO（<https://www.cosmos.esa.int/web/gaia-users/license>），署名 ESA/Gaia/DPAC，限非商业用途。仓库中的公有领域数据没有这一使用限制。
- 选星：距星团中心 20′ 以内、G < 19、自行与星团相差 1.5 mas/yr 以内、视差与 1.85 kpc 相差 3σ 以内，共 13516 颗。
- 处理：按 A_V = 1.31（Sigurdsson et al. 2003）和 Wang & Chen 2019 的 Gaia 波段系数扣除消光。根据 Bp−Rp，使用 Pecaut & Mamajek 2013 的矮星表换算 V 和 B−V；巨星也使用该表，因此结果是近似值。绝对星等按距离 1.85 kpc（Baumgardt & Vasiliev 2021）计算。
- 天球位置取自 Gaia 实测数据。视线方向的深度未测得，按 King（1962）密度分布随机抽取，参数取自 Harris 2010：c = 1.65，r_c = 1.16′。随机种子固定，三维结构表示统计分布。
- 列：`east_pc,north_pc,away_pc,abs_vmag,bv`（相对星团中心，中心处的切平面坐标）。
- 缺失：G > 19（约 0.65 M☉ 以下）的暗星，以及星团核心中因拥挤而未被 Gaia 检出的星。
