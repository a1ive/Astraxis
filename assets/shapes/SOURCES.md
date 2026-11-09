# 形状模型来源

场景中设置 `shape = "shapes/<名称>.mesh"` 后，程序用三角网格绘制不规则天体。网格采用天体固连坐标系，x 指向本初子午线，z 沿自转轴正极，单位为 km。每个顶点保存相对反照率，渲染时与场景的 `color` 相乘。

| 文件 | 来源 | 许可 |
|---|---|---|
| `67p.mesh` | SHAP5 SPC 形状模型（R. Gaskell、L. Jorda；Jorda et al. 2016, Icarus 277, 257），PDS SBN `RO-C-MULTI-5-67P-SHAPE-V2.0`，<https://pdssbn.astro.umd.edu/holdings/ro-c-multi-5-67p-shape-v2.0/>，文件 `data/triplate/spc_lam_psi/shap5/cg_spc_shap5_050k_cart.wrl`（SHA-256 `3304263aaae23de3a82cdf13f754b3c075a57ccd1bdbe4c80117e01a02859a68`） | PDS 存档数据（ESA Rosetta / OSIRIS），公开使用，请引用数据集与作者 |
| `steins.mesh` | Jorda et al. 2012, Icarus 221, 1089（OSIRIS），PDS SBN `RO-A-OSINAC_OSIWAC-5-STEINS-SHAPE-V1.0`，文件 `data/steins_cart.wrl`（SHA-256 `4ef16c6d805336dc4dd878835e900b62c8ba4f0f8056da1ae18189cb6dc67aff`） | 同上 |
| `lutetia.mesh` | Sierks et al. 2011, Science 334, 487；Jorda & Vincent（LAM，OSIRIS），PDS SBN `RO-A-OSINAC_OSIWAC-5-LUTETIA-SHAPE-V1.0`，文件 `data/lutetia_048k_cart.wrl`（SHA-256 `b7185983fba3aa67f342894c7d4e5f1edfed8b2fae74b3e4a7efc6d595016337`） | 同上 |
| `amalthea.mesh`、`thebe.mesh` | P. Stooke，“Stooke Small Body Shape Models” V2.0（PDS SBN，`small_bodies.stooke.shape-models`；EAR-A-5-DDR-STOOKE-SHAPE-MODELS-V2.0），<https://sbnarchive.psi.edu/pds4/non_mission/small_bodies.stooke.shape-models/>，文件 `data/j5amalthea.tab`、`data/j14thebe.tab`；Amalthea：Stooke 1994, EMP 64, 87 | PDS 存档数据，公有领域（请引用：Stooke, P., Stooke Small Body Shape Models V2.0, NASA PDS, 2016） |
| `halley.mesh` | 同上 Stooke 数据集，文件 `data/1682q1halley.tab`（SHA-256 `e9437ac850060a4c527d1045e028ac28f58bfe0386fc1da0bf50405b8128b824`）；1P/Halley，Giotto HMC 与 Vega 1/2 TVS 图像的边缘和明暗界线拟合，指向由 A. Abergel 计算；自转状态按 Belton et al. 1991, Icarus 93, 183 | PDS 存档数据，公有领域（标签要求同时署名 Abergel 与 Stooke） |
| `hyperion.mesh` | P. C. Thomas，“Saturn Small Moon Shape Models” V1.0（PDS SBN，`saturn_satellite_shape_models`，卡西尼 ISS），<https://sbnarchive.psi.edu/pds4/cassini/saturn_satellite_shape_models_V1_0/>，文件 `data/hyperion_30k_plt.tab`；Thomas et al. 2007, Nature 448, 50 | PDS 存档数据，公有领域（请引用数据集与作者） |
| `janus.mesh`、`epimetheus.mesh` | 同上 Thomas 数据集，文件 `data/janus_30k_plt.tab`、`data/epimetheus_30k_plt.tab`（标签所列 PCK：`janus_mst2013.bpc` 等）；Thomas et al. 2013, Icarus 226, 999 | 同上 |
| `phoebe.mesh` | R. W. Gaskell，“Gaskell Phoebe Shape Model” V1.0（PDS4 `urn:nasa:pds:gaskell.phoebe.shape-model::1.0`，由 PDS3 `CO-SA-ISSNA-5-PHOEBESHAPE-V2.0` 迁移；卡西尼 ISS 的 208 张图像，立体光度法），doi:[10.26033/ehkj-xj95](https://doi.org/10.26033/ehkj-xj95)，<https://sbnarchive.psi.edu/pds4/non_mission/gaskell.phoebe.shape-model/>，文件 `data/phoebe_quad64q.tab` | PDS 存档数据，公有领域（请引用：Gaskell, R.W., Gaskell Phoebe Shape Model V1.0, NASA PDS, 2020） |
| `pallas.mesh` | VLT/SPHERE 大型观测计划的 MPCD 形状模型（Marsset et al. 2020, Nature Astronomy 4, 569；Vernazza et al. 2021, A&A 654, A56），Sketchfab 用户 “marin”（@marin14）上传的 “(2) Pallas - MPCD model”，<https://sketchfab.com/3d-models/2-pallas-mpcd-model-59c1ec418b2f4f13bc182633a2cddc66>，文件 `2_Pallas_mpcd.obj`（2022-06-13） | CC BY 4.0（Sketchfab 页面标注；署名：marin (@marin14) / Marsset et al. 2020, Vernazza et al. 2021） |
| `vesta.mesh` | USGS Astrogeology / DLR，Vesta Dawn FC HAMO Global DTM 93m（2013-12-10 发布，黎明号分幅相机 HAMO 立体图像），<https://astrogeology.usgs.gov/search/map/vesta_dawn_fc_hamo_global_dtm_93m>，文件 `Vesta_Dawn_HAMO_DTM_DLR_Global_48ppd.tif` | 公有领域（请引用作者） |
| `phobos.mesh`、`deimos.mesh` | P. C. Thomas，“Small Body Optical Shape Models” V1.0（PDS SBN，`ast-sat.thomas.shape-models`，由 PDS3 数据集 EAR-A-5-DDR-SHAPE-MODELS-V2.1 迁移），<https://sbnarchive.psi.edu/pds4/non_mission/ast-sat.thomas.shape-models_V1_0/>，文件 `data/m1phobos.tab`、`data/m2deimos.tab`；Thomas 1993, Icarus 105, 326 | PDS 存档数据，公有领域（请引用数据集与作者） |
| `endurance_hull.mesh`、`endurance_ring.mesh` | 本项目原创的简化造型，致敬电影《星际穿越》（2014）里的飞船 Endurance，由 `tools/shapes/make_endurance.py` 生成（只用标准库）；不是电影的模型或制作数据 | 随本项目（GPL-3.0-or-later） |
| `voyager.mesh`、`jwst.mesh`、`parker.mesh`（和 `../textures/` 里同名的图集） | NASA 3D Resources，<https://github.com/nasa/NASA-3D-Resources>，见下面的“探测器模型” | NASA 素材，美国不受版权保护（NASA 媒体使用指南），徽标除外 |
| `didymos.mesh`、`dimorphos.mesh` | DART Shapemodel Archive Bundle（Daly, T.; Barnouin, O.; Ernst, C.; Nair, H.; Espiritu, R.; Waller, D.，PDS SBN `urn:nasa:pds:dart_shapemodel::1.0`，doi:[10.26007/96fn-p578](https://doi.org/10.26007/96fn-p578)），<https://pdssbn.astro.umd.edu/holdings/pds4-dart_shapemodel-v1.0/>：Didymos v003 的 9.3 m 版本 `data_derived_didymos_model_v003/didymos_g_9309mm_spc_obj_0000n00000_v003.obj`（SHA-256 `b5349374186d254ca21b7b00ddd9f6dea7e2762ca6bbff13f43151ee448e3f20`），Dimorphos v004 的 1.9 m 版本 `data_derived_dimorphos_model_v004/dimorphos_g_1940mm_spc_obj_0000n00000_v004.obj`（SHA-256 `833eaeac8c49338abe64012a679b1b9695f774f1b5e4842a8e5b9f0af1d28664`）；Daly et al. 2023, Nature 616, 443；Daly et al. 2024, PSJ 5, 24 | NASA DART 任务的 PDS 数据，公有领域（请引用数据集与作者） |
| `arrokoth.mesh` | Porter, S. et al. 2024, “New Horizons Porter (2024) Arrokoth Shape Model Collection”，PDS Small Bodies Node，doi:[10.26007/97r3-1e19](https://doi.org/10.26007/97r3-1e19)，<https://pdssbn.astro.umd.edu/holdings/pds4-nh_derived-v3.0/arrokoth_shapemodel_porter2024/>；方法见同目录的 `porteretal2024b.pdf`（Porter et al. 2024, “The Shape of (486958) Arrokoth”） | NASA 新视野号项目的 PDS 数据，公有领域（请引用作者） |
| `homunculus.mesh` | η Carinae 的 Homunculus 星云三维模型（Steffen et al. 2014, MNRAS 442, 3316，arXiv:1407.4096；SHAPE 软件拟合 VLT/X-shooter 的 H₂ 2.12 μm 速度图），ESO “A 3D Model of the Homunculus Nebula around the star Eta Carinae”（3dmodel_004，2014-07-28），<https://www.eso.org/public/products/models3d/3dmodel_004/>，文件 `model_3d_obj/3dmodel_004.zip` 里的 `Eta_Car_Homunuculus_model.obj` | ESO 公开素材，CC BY 4.0（模型页面的使用条款指向 ESO 的 CC BY 4.0 声明）；署名：NASA's Goddard Space Flight Center/Ed Campion & UNAM/W. Steffen |

原始文件（同一目录下）：

- `arrokoth_porter_2024_v01.obj`（3.7 MB）：20,484 个顶点、40,960 个三角面，原点为质心，轴为惯量主轴。
- `albedo_arrokoth4_fp36h2_masked1.png`（0.2 MB）：反照率图，OBJ 的纹理坐标指向它。

`tools/shapes/make_arrokoth.py` 将原始模型转换为 800 KB 的网格文件，依赖 Pillow。

文件格式见 `tools/shapes/axmesh.py`（小端）：`char[8] "AXMESH2\0"`、`uint32` 顶点数、`uint32` 三角形数、`uint32` 标志（第 0 位：带地图坐标；第 1 位：硬边；第 2 位：纹理图集），之后是每个顶点的 `float32` x、y、z（km），每个顶点的 `float32` 相对反照率，（有标志时）每个顶点的 `float32` 地图 u（东经 / 360°，接缝处的顶点重复，取 0 和 1），（有第 2 位时）每个顶点的 `float32` v，每个三角形 3 个 `uint32` 顶点索引（从外面看逆时针）。带纹理图集的网格（探测器）里 u、v 是普通的纹理坐标，“反照率”存的是太阳在 +z 时到达这个顶点的阳光（0 或 1）。法线在加载时按面积加权平均算出，位置完全相同的顶点（接缝、极点）共用一个法线。带硬边标志的网格（人造的平面体，如探测器）不做这种合并，只有共用同一个顶点索引的三角形才平滑过渡：每个平面单独用自己的顶点，棱角才是硬的。带地图坐标的网格可以在场景里同时写 `texture`：u 取自文件，纬度取自顶点方向（行星中心纬度）。

## Arrokoth

- 形状：两瓣接触双星，大瓣 Wenu 在 +x 方向，小瓣 Weeyo 在 −x 方向。整体 34.546 × 19.838 × 13.822 km，与论文表 2 一致。网格体积 4123.8 km³，换算成等体积球直径 19.896 km，也与表 2 一致（`tests/shape_tests.cpp` 的 `test_arrokoth_shape`）。模型标签注明两瓣在连接处略有重叠。
- 极轴与本初子午线：用新视野号项目的 PCK `nh_arrokoth_002.tpc`（<https://naif.jpl.nasa.gov/pub/naif/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/pck/>）：极轴 RA 319.3685820°、Dec −25.5875063°，W = 71.9541475° + 559.3936173°/天（自 J2000 TDB）。文件注释注明 2024 版形状模型使用这组数值。极轴与论文表 2（319.37°、−25.588°）一致。
  - 形状模型标签和 `arrokoth_geophysics` 概述里写的极轴（317.488°、−24.888°）和周期（0.6632553 天 = 15.92 h）是 2021 年的参数；本模型的朝向按上述 PCK 参数计算。
  - 由此推出的周期是 15.45 h，比 Spencer et al. 2020 的 15.92 h 短 3%。New Horizons 在飞掠前后两天拍到了可分辨的图像，观测时段较短，对周期的约束较弱。飞掠时刻的朝向按拟合时使用的这组参数计算。
  - 交叉检查：极轴与 Arrokoth 日心轨道法向的夹角算得 100.39°（表 2：100.39°），与新视野号来向的夹角 41.04°（表 2：41.096°）。
- 反照率：PNG 是 FITS 数组上下翻转后的图像，按显示方向存储。OBJ 纹理坐标 v = 0 对应图像最下一行，符合 OBJ 的惯例。按这个方向采样，亮的“颈部”落在两瓣连接处。像素换算为 albedo = 像素 / 1361975.975 + 0.03188（标签）。新视野号拍到的部分约占 40% 的顶点，主要在南半球。其余顶点对应固定的掩膜值，转换工具从相邻已知顶点沿网格扩散填充，远处渐变到平均值。因此背面显示为平淡的表面。文件存储反照率与已知部分平均值（0.0438）的比值，范围为 0.74 至 1.83。
- 颜色：场景使用示意颜色 `#a8604a`。Arrokoth 偏红，LORRI 是全色相机，图像本身无法给出这里使用的颜色。
- 光照采用 Lambert 漫反射，自阴影尚未实现。太阳很低时，颈部的凹处也会被照亮。

## Phobos 与 Deimos（Thomas 形状模型）

- 形状：Thomas 的数值形状模型，只用海盗号轨道器的图像，按行星中心纬度 / 经度的网格给出半径（0° 和 360° 两列都有）。`tools/shapes/make_grid_table_shape.py`（只用标准库）把网格原样转成网格文件，去掉极点处退化的三角形。
  - Phobos：2° 网格（91 × 181 个点），16,471 个顶点、32,040 个三角形，697 KB。
  - Deimos：5° 网格（37 × 73 个点），2,701 个顶点、5,040 个三角形，111 KB。标签注明 200°至355° 经度（西经）一带的误差约 400 m。
- 经度方向：表里的经度是西经。依据是 Phobos 表：减去拟合的三轴椭球后，49° 处有清楚的陨石坑特征（坑内低约 0.7 km、外圈高起），311° 处没有；Stickney 坑在 49°W。网格里换成东经（u = 东经 / 360°）。形状里 Stickney 凹陷的中心在 2°S 50°W；网格里该处 12° 以内的平均半径比外圈低 1.12 km，镜像位置 50°E 没有（`test_mars_moon_shapes`）。
- 体积与 JPL SSD 卫星物理参数表的平均半径（Archinal et al. 2018）对照（`test_mars_moon_shapes`）：
  - Phobos：5748.7 km³，等体积半径 11.11 km，表中为 11.08 ± 0.04 km。
  - Deimos：1013.5 km³，等体积半径 6.23 km，表中为 6.2 ± 0.24 km。
- 两颗卫星的最长轴都指向火星（本初子午线），各轴的范围接近场景里 PCK 的三轴椭球（Deimos 沿 x / y / z 离中心最远约 8.4 / 6.7 / 5.8 km，PCK 为 7.8 × 6.0 × 5.1 km）。场景的 `radii_km` 使用 PCK 值（影子、相机距离）。
- 自转参数：Deimos 的标签给出极轴 RA 316.65°、Dec 53.53°，W = 79.41° + 285.1618970°/天，与场景所用 PCK 的值一致，两者采用相同的本初子午线。
- 纹理是同一数据集的拼接图 `m1phobosm.fit`、`m2deimosm.fit`（见 `assets/textures/SOURCES.md`），拼接时的位置就是用这些形状模型控制的。备选的 Stooke 拼接图（USGS）基于 DLR 的另一套控制网，地貌比形状模型偏西约 5°（Limtoc 在 12°S 59°W，Thomas 图和 IAU 地名为 11°S 54°W 左右），偏差各处不同，贴在形状上会错开，因此使用 Thomas 拼接图。
- Deimos 的 5° 网格很粗，近看外形是多面体状的；两颗卫星的光照都没有自身阴影。

## Vesta（DLR HAMO DTM）

- 原始文件：<https://asc-pds-services.s3.us-west-2.amazonaws.com/mosaic/Vesta_Dawn_HAMO_DTM_DLR_Global_48ppd.tif>（597 MB，17280 × 8640，32 位浮点，未压缩、每行一个条带）。像素值是从 Vesta 中心起算的半径（米），USGS 页面给出的平均交会误差为 ±8 m。简单圆柱投影、行星中心纬度、东经向右、左边缘 −180°E；经度是 Claudia″ 系统（Claudia 坑在 146°E），与纹理 `vesta.jpg` 和场景所用 PCK（W0 = 285.39°）相同。DTM 与纹理同属 DLR 2013-12-10 的发布（ProductId 20131210），控制网一致。
- USGS 页面说 HAMO 立体图像覆盖约 95% 的表面，但文件里没有无数据像素（半径 211.9至293.0 km），未覆盖的部分在产品里已经补上了。
- `tools/shapes/make_dtm_shape.py`（只用 Pillow 读 TIFF 的条带位置，数据用标准库按行读取）按 1.5° 的格子求平均，每个网格顶点取周围四格的平均，极点取整圈格子的平均：121 × 241 个点，29,161 个顶点、57,120 个三角形，1.2 MB。1.5° 约合 6.9 km，Rheasilvia 盆地（约 500 km）和它的中央峰都分辨得出，小坑的细节靠纹理。
- 交叉检查：体积 7.495 × 10⁷ km³，等体积直径 523.11 km，SBDB 为 522.77 ± 0.1 km（Park et al. 2025，Nature Astronomy）；差 0.34 km，来自 2013 年的 DTM 和网格的平滑（`test_vesta_shape`）。SBDB 的外形尺寸 569.24 × 554.48 × 452.66 km 应是拟合椭球的轴长，网格的包围盒是 574.5 × 554.0 × 467.5 km（包含地形起伏），不直接比较。
- 场景的 `radii_km` 使用 PCK 的三轴椭球（影子、相机距离）；光照没有自身阴影。

## Amalthea 与 Thebe（Stooke 形状模型）

- 表的列是经度、纬度、半径（与 Thomas 表的顺序不同，`make_grid_table_shape.py --lon-first`），5° 网格，行星中心坐标；标签写明卫星的经度向西增加，本初子午线正对木星。转换后各 2,701 个顶点、5,040 个三角形，111 KB。
- Amalthea 基于旅行者 1、2 号图像（修正了 Stooke 1994 里经度 315° 附近的一个“鼓包”，尚未用伽利略号的结果）；Thebe 基于伽利略号的低分辨率图像，标签说只是初步模型。坐标原点不一定是图形中心（标签说明），所以 Amalthea 的 x 范围是 −147..120 km。
- 交叉检查（`test_jupiter_saturn_small_moon_shapes`）：等体积半径 Amalthea 81.7 km、Thebe 45.5 km，JPL SSD 平均半径为 83.5 ± 3.0、49.3 ± 4.0 km；两者的最长轴都沿 x（指向木星），与 PCK 的三轴椭球（125 × 73 × 64、58 × 49 × 42 km）一致。场景的 `radii_km` 使用 PCK 三轴值（影子、相机距离）。
- Amalthea 贴 Stooke 的晕渲图（见 `assets/textures/SOURCES.md`），Thebe 没有图，按场景颜色均匀着色。

## Halley（Stooke / Abergel 形状模型）

- 表的列是经度、纬度、半径，5° 网格，行星中心坐标。和卫星不同，彗星的经度向东增加（V2.0 按 IAU 2009 的约定把小行星和彗星的经度反过来了），所以转换时加 `--lon-first --east`。表里只有 360° 列，没有 0° 列，转换工具按 360° 取模后当作 0°。转换后 2,701 个顶点、5,040 个三角形，111 KB。
- 坐标系（标签）：参考轴是长轴，北极指向“大头”；经度零点的定义是 Vega 2 第一张高分辨率图像（2:00:30）拍摄时星下点经度为 270°。所以网格的 z 轴就是彗核长轴，x、y 在短轴平面内。
  - +z 确实是粗的一头：离长轴中点 2至6 km 处，+z 一侧的最大截面半径是 4.6至3.3 km，−z 一侧是 3.6至2.6 km。
  - 270° 是 V1.0 的西经：V1.0 的表（PDS3 `EAR_A_5_DDR_STOOKE_SHAPE_MODELS_V1_0`）与 V2.0 逐点比较，V2.0 的经度恰好是 360° 减 V1.0 的经度，而标签里这句话一字未改。换成网格的东经就是 90°。
  - 2:00:30 不是 UT（Vega 2 最近点在 07:20:00 UT），应是电视系统开机后的星上时间。Vega 1 的图像头（PDS SBN `VEGA1-C-TVS-3-RDR-HALLEY-PROCESSED-V1.0`）有星上计数 `TONBHRDW`，按数据集说明 `tvs_proc.doc` 的 4.766485 计数/秒和 5.5 s 读出修正，最近点前 6.6 s 那张（`t01135`）是开机后 2:00:32，开机约在最近点前 2 小时。Vega 2 的图像头没有这个计数，这里假定程序相同，取最近点前 1.5 s 那张（`VEGA2-C-TVS-3-RDR-HALLEY-PROCESSED-V1.0` 的 `t01190`，头文件 `TIME-ENC = -1.5`）。Belton et al. 1991 图 8 对照的也是这张（“encounter − 1.5 sec”）。前后一张图差约 100 s，星下点经度在最近点后变化约 10°，最近点前变化快得多。
- 误差（标签）：图像质量较差，边缘和明暗界线的位置本身就不确定，每个点的绝对误差估计为 0.5至1 km，相邻点之间约 100 m；凸包的可靠程度估计与模型本身相当。Stooke 自己的比较认为他的模型偏多棱角、凹陷偏深。
- 交叉检查：Szegő 1991（Comets in the Post-Halley Era 2, 723，第 2 节）的 Vega/Giotto 独立模型装在 15.3 × 7.2 × 7.22 km 的盒子里，体积 365 km³。本网格沿 z 长 15.14 km，x、y 方向 7.53、7.59 km，体积 402.2 km³（等体积半径 4.58 km），差约 10%（`test_halley_shape`）。场景的 `radii_km` 用 Szegő 的盒子尺寸的一半（3.61、3.6、7.65 km，长轴沿 z），只用于影子和相机距离。
- 自转：Belton et al. 1991（Icarus 93, 183）的“Slow LAM”（第 3.1 节），对称陀螺的自由进动：长轴与角动量 M 夹角 66.0°，绕 M 进动周期 P_ψ = 3.69 ± 0.05 天，绕长轴自转周期 P_φ = 7.1 ± 0.1 天，合成的总自转周期 2.84 天，总角速度与 M 夹角 21.4°（这两个数可以由前三个算出，`test_halley_shape` 检查）；M 指向 B1950 赤经 6.2°、赤纬 −60.7°，用 SPICE 的 FK4 → J2000 旋转换算为 J2000 的 6.785°、−60.423°。场景用 `[bodies.free_precession]` 按 z-x-z 欧拉角计算朝向（`iau_pole_frame(M) · Rz(φ) · Rx(θ) · Rz(ψ)`，φ 是进动角、ψ 是绕长轴的转角），历元取 Vega 2 最近点 1986-03-09 07:20:00 UTC。
  - 进动角：正文第 2.2 节定义了彗核坐标系，在 Vega 2 交会时刻（JD 2446498.80556）长轴 P 穿过“大头”，指向 B1950 赤经 313.20°、赤纬 −7.52°（J2000 313.867°、−7.328°，与 M 夹角 66.04°，和 66.0° 一致）。Stooke 的北极也指向大头，所以网格 +z 就是 P。由此 φ₀ = −120.02°。
  - 绕长轴的转角：Belton 的经度零点在 M × P 与彗核赤道的交点，正好是上式里 ψ = 0 时的 x 轴（节点线）。Stooke 的经度零点另有定义（见上），需要把 Vega 2 在 `t01190` 时的方向换到 Belton 的坐标系：用场景里的 Vega 2 轨迹，星下点在 Belton 经度 194.84°（纬度 −37.8°），网格里应是东经 90°，所以 ψ₀ = 104.84°。
  - 场景的 Vega 2 轨迹在交会附近是按脱靶距离和 Vega 1 的方向重建的（`assets/ephem/SOURCES.md`），这里用 Vega 2 图像头（IKI 1986 年 3 月的初步轨道）检查：11 张图的距离差 ≤ 70 km，相位角差 ≤ 0.6°（`t01190`：8,032 km、28.7°，场景 8,030 km、28.9°），太阳在图像里的方位角（`SUNANG`）差 ≤ 2°，所以这一几何关系是可靠的。
  - 交叉检查：Samarasinha & A'Hearn 1991（Icarus 93, 194）表 II 列出三次交会时大头的黄道（B1950）经纬度，每个约 ±15°。Belton 的假设（“SV1G2”）对应其中组合 (2,1,2)：Vega 1 反向 259°、−15°，Vega 2 310°、9°，Giotto 镜像 235°、−32°。只用 Vega 2 的 P 定相位后，模型在三次最近点时的大头方向与它们分别差 1.6°、3.5°、2.0°（`test_halley_shape`）。
  - 误差：Belton 拟合三次交会的长轴方向，均方根 2.5°；周期误差使两个相位每天各漂移约 1.3° 和 0.7°，所以 1986 年 3 月初的几次交会前后朝向可信到几度，离交会越远越示意。Samarasinha & A'Hearn 认为 M 的具体坐标没有很好确定；真实彗核也不是对称陀螺（Belton 等人自己说这只是一级近似）。
- 纹理：缺少 Halley 的全球反照率图；Giotto HMC 只拍到一侧，那张著名的彗核合成照片由 ESA 和 HMC 团队发布，许可证不明确，不当作公有领域。按场景颜色（示意性的暗灰色）均匀着色。

## Hyperion（Thomas 卡西尼形状模型）

- 板块模型：首行是顶点数和板块数，之后每行一个顶点（km），再之后每行一个板块的三个顶点索引（从 0 开始，从外面看逆时针）。`tools/shapes/make_plate_shape.py`（只用标准库）原样转换：14,636 个顶点、29,268 个三角形，571 KB。
- 坐标系：Hyperion 是混沌自转，没有 IAU 自转模型。这个模型的经纬网以卡西尼第 15 圈早段观测到的自转矢量为极（UTC 2005 年第 268 天 04:29–17:20，即 2005-09-25，最近一次飞掠的前一天；Thomas et al. 2007 图 1d），z 轴就是当时的自转轴，保留了 Davies et al. 1983 的经度参考（Bahloo 坑在 196°W）；各次飞掠测得的自转轴相差可达几十度。最佳飞掠覆盖的区域相对误差小于 1 km，对面一侧最多 6 km（数据集文档 `hyperion_document.pdf`）。z 方向的范围最大（341 km），但长轴（惯量最小的主轴 A）与 z 轴相差 26°。
- 交叉检查：等体积半径 136.3 km，JPL SSD 平均半径 135 ± 4 km（`test_jupiter_saturn_small_moon_shapes`）。
- 自转：按 Harbison, Thomas & Nicholson 2011（Celest. Mech. Dyn. Astron. 110, 1）的表 1、表 2 设为 2005-09-25 飞掠时的真实朝向。
  - 主轴：按均匀密度算网格的惯量张量，A/C = 0.578、B/C = 0.878，论文第 2 节为 0.58 ± 0.03、0.87 ± 0.03（同样按均匀密度从 Thomas 2007 的模型算出）。网格 z 轴在主轴系里是 (0.896, 0.173, 0.409)，表 1 中 2005-09-25 的自转方向是 (0.902, 0.133, 0.411)，相差 2°。这说明 PDS 的网格与论文用的是同一个模型、同一个坐标系。主轴的正负号按这三个分量都为正来取（A × B = C 自动满足）。
  - 欧拉角：表 1 的 θ、φ、ψ 是从土星坐标系 xyz（z 为土星极轴，x 为土星指向 Hyperion 近土点的方向）到主轴系的 z-x-z 欧拉角：主轴 → xyz = Rz(θ) Rx(φ) Rz(ψ)（式 13–18 的方向余弦由此得出）。用 2005-09-25 的 θ、φ、ψ = 2.989、1.685、1.641 rad 把表 1 的自转方向转到 xyz，得到 (1.154, 2.017, 3.565) n，表 2 为 (1.151, 2.018, 3.565) n。
  - 历元：表 1 的 M = 303° 对应 Horizons（SAT441，GM 取 Horizons 的 37931206.6 km³/s²）的吻切平近点角 2005-09-25 17:33 UTC，比 Thomas 2007 的成像时段（04:29–17:20）略晚；M 只给到整度，对应 ±0.7 h，即本初子午线 ±2°。同一时刻 e = 0.1132（表 1：0.113）。表 1 的 ϖ 是在土星赤道面内、从土星赤道与 J2000 黄道的升交点量起的：三次飞掠算得 103.7°、100.8°、98.1°（表 1：105°、102°、99°），系统性地差约 1°。场景的 x 轴直接用 Horizons 吻切偏心率矢量在土星赤道面上的投影（与用表 1 的 ϖ 差 0.9°）；土星极轴用 PCK。
  - 结果：极轴 RA 333.746°、Dec 53.899°（与表 2 的自转方向差 2.4°，表 2 的 σ 是 10°），W = 50.136° + 72.0797°/天。速率是 4.255 n，n 取论文第 1 节的 16.94°/天：这样换算出的 72.08、75.10°/天与表 2 的 72 ± 1、75 ± 1°/天一致（第 3 节的时间单位 P = 21.43 天是 Horizons 的吻切周期，用它换算是 71.5、74.5°/天，取整对不上表 2）。`test_hyperion_spin_state` 检查主轴、欧拉角约定、历元和场景里的朝向。
  - 局限：绕固定的极轴匀速自转只在 2005-09-25 前后几天接近真实：40 天前（8 月）测得的自转轴与这里差 28°，107 天前（6 月）差 55°。论文的动力学积分给出约 61 天的李雅普诺夫时间，从 6 月的状态积分也对不上 8、9 月的观测（体内经度差 120°–180°），所以不做积分。旅行者 2 号（1981）和卡西尼 2007-02-16 的状态见表 2。
- `radii_km` 用 PCK 的 180.1 × 133.0 × 102.7 km（影子、相机距离）。
- 纹理：Stooke 的 Hyperion 晕渲图基于旅行者 2 号，所用自转轴与卡西尼的相差 100° 以上，经纬度体系对不上这个形状模型；颜色是示意性的。

## Janus 与 Epimetheus（Thomas 卡西尼形状模型）

- 板块模型，格式同 Hyperion，`tools/shapes/make_plate_shape.py` 原样转换：Janus 13,381 个顶点、26,758 个三角形，522 KB；Epimetheus 13,915 个顶点、27,826 个三角形，543 KB。只用卡西尼 ISS 窄角相机的图像，模型中心与形心重合到 45 m、25 m（标签）。
- 坐标系（标签）：+Z 是北极，+X 是本初子午线，“on average points toward Saturn”，+Y 指向轨道运动的反方向。场景里两者设为 `synchronous`（x 始终指向土星中心、z 沿轨道法线），正好是这个坐标系；不画天平动（光学天平动约 2e，以及卡西尼测到的受迫天平动）。PCK 里的自转模型只在 1981 年旅行者号飞掠前后可用（Horizons 的说明），每次交换轨道都会改变平均运动，所以不用。
- 交叉检查：等体积半径 88.98 km、58.59 km，JPL SSD 平均半径 89.2 ± 0.8 km、58.2 ± 1.2 km；范围 206.7 × 189.8 × 153.2 km、146.1 × 112.8 × 107.8 km，x 最长、z 最短（`test_jupiter_saturn_small_moon_shapes`）。`radii_km` 用 PCK 的三轴值 101.7 × 93.0 × 76.3 km、64.9 × 57.3 × 53.0 km。
- 纹理：没有用全球拼接图，颜色是示意性的。

## Phoebe（Gaskell 卡西尼形状模型）

- 格式：ICQ（隐式连接四边形，数据集文档 `document/icqmodel.asc`）。第一行是 Q，之后按面 f = 1..6、行 j = 0..Q、列 i = 0..Q 列出 6(Q+1)² 个顶点（km）；立方体棱和角上的顶点在相邻面里各存一份。`tools/shapes/make_plate_shape.py` 识别首行只有一个数的表为 ICQ，按 Gaskell 文档里的索引表合并重复顶点，每个四边形沿较短的对角线分成两个三角形（与他的转换代码相同），只保留用到的 6Q² + 2 个顶点。
- 分辨率：数据集有 Q = 64、128、256、512 四种（512 是原始模型，157 万个顶点）。用 Q = 64：24,578 个顶点、49,152 个三角形，960 KB，网格间距约 2.6 km；飞掠镜头在几千 km 外，看不出与 Q = 128（9.8 万个顶点，3.8 MB）的差别。两者体积相差 0.02%。
- 坐标系：readme 写明模型用的自转常数是 `BODY609_POLE_RA = 356.90`、`POLE_DEC = +77.88`、`PM = 178.58 + 931.639 d`，与场景所用 pck00011 一致，只有极赤纬是 77.88° 而不是 77.80°（pck00011 注明 IAU 2015 报告未改动这组值；0.08° 在 106 km 的半径上只差 0.15 km）。场景用 PCK 的值。
- 交叉检查：等体积半径 106.56 km，JPL SSD 平均半径 106.50 ± 0.70 km；范围 226.0 × 228.4 × 209.0 km，最短的是 z（自转轴），与 PCK 的三轴 109.4 × 108.5 × 101.8 km 一样（`test_jupiter_saturn_small_moon_shapes`）。`radii_km` 用 PCK 三轴值（影子、相机距离）。
- 误差（readme）：522 个地标，位置均方根误差 147 m。卡西尼只在 2004-06-11 近距离飞掠过一次。
- 纹理：USGS 和 Stooke 的小天体地图里都没有 Phoebe 的全球拼接图，按场景颜色（示意性的暗灰色，几何反照率约 0.08）均匀着色。

## Pallas（VLT/SPHERE MPCD 形状模型）

- 来源：Pallas 尚无航天器探测，模型根据 ESO VLT/SPHERE/ZIMPOL 自适应光学图像（2017 至 2019）重建。重建分两步：先用 ADAM 程序结合光变曲线、掩星和 SPHERE 图像得到形状和自转参数，再以它为输入，用 MPCD（Capanna et al. 2013）从 SPHERE 图像中提取更细的地形，包括“高尔夫球”状的撞击坑（Vernazza et al. 2021 第 3 节）。上传者 “marin” 的 Sketchfab 账号只发布这个计划的 MPCD 模型（Psyche、Kalliope、Kleopatra 等），其身份可能是论文作者之一 Marin Ferrais，页面未注明，尚待确认。DAMIT 和 VizieR（J/A+A/654/A56）都没有收录 MPCD 模型，只有 ADAM 模型。
- 文件：OBJ 由 Fortran 程序 `OBJW_SAVE_RX01` 写出，22,530 个顶点、45,056 个三角面，单位 km，没有纹理坐标。`tools/shapes/make_plate_shape.py`（只用标准库，同时读板块表和 OBJ）原样转换，880 KB。
- 坐标系：MPCD 以 ADAM 模型为输入，两者的体固连坐标系可通过模型比较核对。用 DAMIT 的 ADAM 模型 4395（Marsset et al. 2020，<https://damit.cuni.cz/projects/damit/asteroid_models/view/4395>，DAMIT 内容为 CC BY 4.0）核对：沿 600 个随机方向比较两个模型的半径，绕 z 轴转 0° 时最吻合（相关系数 0.983，均方根差 4.6 km），转 ±10° 后明显变差。MPCD 模型的体积中心偏离原点约 3 km，转换时保持该偏移。
- 自转：DAMIT 模型 4395 的 λ = 42°、β = −15°（J2000 黄道），P = 7.81322 h，在 JD 2433827.77154 时 φ0 = 0。约定为“体固连 → 黄道 = Rz(λ) Ry(90° − β) Rz(φ)”，与 Vernazza et al. 2021 表 A.1（λ 42 ± 3°、β −15 ± 3°、P 7.81321 ± 2 × 10⁻⁵ h）一致。场景换算成 IAU 形式：极轴 RA 44.1136°、Dec 1.1250°，W = 41.608° + 1105.818241°/天。自转速率取 DAMIT 的 IAUspin 文件（它用了 P 的更多位数，7.8132189 h；用四舍五入后的 7.81322 h 推到 J2000 会差约 2.8°）。DAMIT 的 IAUspin 文件给出 W0 = 41.7°，它的极轴取整成 RA 44°、Dec 1°，与这里的结果一致。`test_pallas_shape` 检查体积、周期、W0，并在 t0 时刻比较场景的天体轴与按 DAMIT 约定直接构造的轴（相差小于 0.05°）。
  - PCK（IAU 2015）给出的极轴为 RA 33°、Dec −3°，W = 38° + 1105.8036°/天。这组参数与形状模型的坐标系无关，场景因此使用上述 DAMIT 参数。
- 交叉检查：等体积直径 508.2 km，Vernazza et al. 2021 表 1 为 511 ± 4 km（ADAM 模型为 513 km）；外形尺寸 562 × 529 × 429 km，表 A.1 的三轴椭球为 568 × 530 × 450（± 12）km。`radii_km` 使用 [MARS] 的 568 × 532 × 448 km（影子、相机距离）。
- 局限：SPHERE/ZIMPOL 的角分辨率约 20 mas（600 nm），能探测到的撞击坑直径至少约 25至70 km（随观测时的地心距离而变，Vernazza et al. 2021），模型里只保留了其中最显著的；SPHERE 图像覆盖不到的部分主要由光变曲线约束，较平滑。没有反照率信息，按场景颜色均匀着色。

## 67P、Steins、Lutetia（Rosetta 的 OSIRIS 形状模型）

- **格式**：VRML 2.0 的 IndexedFaceSet：`point [...]` 是顶点（km，天体固连），`coordIndex [a b c -1 ...]` 是三角形（从 0 开始）。`tools/shapes/make_plate_shape.py` 新增了 `.wrl` 读取，原样转换。
- **67P**：
  - 模型：SHAP5 的 5 万面版本，24,877 个顶点、49,748 个三角形，971 KB。
  - 体积 18.74 km³，等体积半径 1.65 km；范围 5.0 × 3.7 × 3.3 km。
  - 坐标系：SHAP5 相对 Cheops 参考系只转了 0.28°、平移了约 17 m（数据集 `USER_GUIDE.ASC`），这里直接当作 Cheops 系使用。
  - 自转：按 Cheops 参考系文档（Scholten et al. 2015，同一数据集 `DOCUMENT/CHEOPS_REF_FRAME_V1.PDF`）：极 RA 69.54°、Dec 64.11°，W = 114.69° + 696.543884683°/天（12.4041 h），只对 2014 年 8–9 月有效。0.14° 的进动没有加；过近日点后自转周期缩短到约 12.06 h，也没有建模。所以 2015 年以后彗核的自转相位不准（轨道不受影响）。
- **Steins**：
  - 模型：10,242 个顶点、20,480 个三角形，400 KB；等体积半径 2.63 km，范围 6.8 × 5.6 × 4.2 km。
  - 自转：数据集说明（`catalog/dataset.cat`）给出极 RA 91°、Dec −62°（±5°），周期 6.04681 h。本初子午线过 Spinel 坑，但没有给出某一时刻的相位，场景里取 0，是任意的。
- **Lutetia**：
  - 模型：48k 版本，23,894 个顶点、47,784 个三角形，933 KB；等体积半径 49.2 km，范围 111.7 × 120.8 × 84.8 km。
  - 自转：数据集文档 `LUTETIA_ROTATION.PDF`（Jorda & Vincent）：极 RA 51.80°、Dec +10.83°（Sierks et al. 2011），W = 289.50° + 1057.751519°/天（周期 8.168270 h，Carry et al. 2010）。零经线过 Lauriacum 坑，模型已按这个定义旋转。

## Didymos 与 Dimorphos（DART 的 SPC 形状模型）

- **格式**：OBJ（`v`、`f` 行，km，三角形按右手定则从外面看逆时针），文件头的 `#` 注释给出面数、体积、惯量张量和范围。`tools/shapes/make_plate_shape.py` 原样转换（文件头有 UTF-8 的弯引号，要 `python -X utf8`），各 24,578 个顶点、49,152 个三角形，960 KB。每个模型都有 4 种分辨率，这里取最粗的（更细的 14–221 MB）
- **Didymos**（v003，DART 的 DRACO 和 LICIACube 图像，Barnouin 等）：体积 0.2034 km³（等体积直径 730 m），范围 818 × 801 × 604 m。Daly et al. 2023 表 1 用的是更早的模型（849 × 851 × 620 m，761 m）；Naidu et al. 2024 按动力学拟合的扁球是 394 × 290 m（体积 0.189 km³，相差 8%，在他们的误差内）。文件头的形心在 (16, 42, −34) m，不在原点（标签说原点是质心），转换时保持原样
- **Dimorphos**（v004，“撞击前的最终模型”，只有 DRACO 图像，阳光相角 60°，背面主要由临边约束）：体积 0.001759 km³（等体积直径 150 m），范围 178 × 169 × 115 m（Daly 2023 表 1：177 × 174 × 116 m、0.00181 km³ ± 10%）。撞击后的形状（Naidu 2024 推断长轴比变成 1.3）不知道，场景始终用这个模型
- **自转**：模型的坐标系就是 DART 的本体系（x 指向零经线，z 沿自转轴）。Didymos 按 Horizons 说明的极（69.32686°, −72.72455°）和 W = 3823.0088°/天（2.2600 h），W0 取 0；Dimorphos 同步自转（x 指向 Didymos，z 沿轨道极），和 Horizons 撞前的自转模型相差不到 0.02°（`test_didymos_scene`），撞后 Horizons 本身就用这个定义。DART 撞在东经 264°（Daly 2023），在同步坐标系里是前导面，测试里 DART 来的方向离它约 15°（撞击点离形心 25 m，DART 在场景里瞄的是中心）

## Endurance（致敬造型）

- 来源：本项目原创。造型只取电影里 Endurance 的大意（中间的船体、一圈 12 个舱段、辐条），尺寸和明暗都是示意的，没有参考或使用电影的模型、图纸和制作数据。单位 m，写入文件时换成 km。
- 两个网格：`endurance_hull.mesh` 是不转的船体（八棱柱主体、两端的环、背面的设备架，设备架不对称，舱段转动时看得出船体的朝向），`endurance_ring.mesh` 是绕 +z 轴转动的环（12 个带倒角的舱段、外侧的设备托架、环管、轮毂和三根辐条；0 号舱段颜色深一些，能看出转动）。场景里船体是天体的 `shape`，环是 `[bodies.rotor]`：同一位置、同一朝向，再绕本体 z 轴转。
- 尺寸：环直径 70 m（托架的角到 70.2 m），舱段在 23–33 m 半径之间，地板（舱段外壁内侧）取 32.5 m；船体宽 25 m、长 20.4 m。场景的 `radii_km` 取 0.036 km，包住整个环。
- 转速：环每 11.44 s 转一圈，32.5 m 处的离心加速度为 1 g（2π√(32.5 m / 9.80665 m/s²)，`test_endurance_shape` 检查）。
- 硬边：文件带硬边标志，每个平面用自己的顶点；圆管的侧壁沿圆周共用顶点，看起来是圆的。
- 反照率：脚本里的明暗值（0.12–0.9）按两个网格合起来的面积加权平均（0.584）归一，文件里是相对值（0.21–1.54），整体颜色由场景的 `color` 给。
- 光照只有吸积盘的漫反射，没有自身阴影：环内侧和环后面的船体也会被照亮。

## 环旁小卫星与特洛伊卫星（Thomas 卡西尼形状模型）

Prometheus、Pandora、Pan、Daphnis、Telesto、Calypso 和 Helene 来自同一套 [Saturn Small Moon Shape Models V1.0](https://sbnarchive.psi.edu/pds4/cassini/saturn_satellite_shape_models_V1_0/)，作者 P. C. Thomas，观测为 Cassini ISS。沿用该 PDS 数据集的公有领域许可，请引用数据集及作者。数据目录下每颗的 `NAME_30k_plt.xml` 给出成像来源、坐标系及版本，板块表为 `NAME_30k_plt.tab`。

下载板块表后运行（例如 Pan）：

```powershell
python tools/shapes/make_plate_shape.py pan_30k_plt.tab assets/shapes/pan.mesh
```

按发布的 km 坐标原样转换，保留全部顶点及三角形，无缩放、补面或纹理。模型 +x 是平均指向土星的本初子午线，+z 是北自转极，+y 朝公转的后随侧。场景按同步自转定向，长轴每时刻指向土星，未加入物理天平动。`radii_km` 使用 IAU 2015 / `pck00011.tpc` 的三轴半径，供阴影和相机距离使用；绘制尺寸由网格本身决定。颜色及均匀反照率是示意值。

| 文件 | 顶点 | 三角形 | 等体积半径（km） | 原始板块表 SHA-256 |
|---|---:|---:|---:|---|
| `prometheus.mesh` | 14183 | 28362 | 42.764 | `78d1d11346d554f969d4c0f0df3969e5c7f50c715274f69a57ad31ec3be5425c` |
| `pandora.mesh` | 13564 | 27124 | 40.035 | `d950c909cf9c991108ef92c07cf63bc6d3bdadc4b19028bb3599aae280f24294` |
| `pan.mesh` | 13736 | 27468 | 13.696 | `18a96fbec0ab0493e78a707644db1a4314db54eec662094ea8d5210a614dbd1b` |
| `daphnis.mesh` | 13465 | 26926 | 3.893 | `8fc85246848ad05de7c68b0686f9fde7d3c23529d02db9c1e7287864b74c34be` |
| `telesto.mesh` | 13804 | 27604 | 12.287 | `46836f5785b12651f4b7ce8db8f14c9023bca897e13f0fe64a62cd05b58bae0f` |
| `calypso.mesh` | 14140 | 28276 | 9.524 | `840fc381659291b29fbaa51bdaf09c870da6bcc96b4ae791709e2a354048fb63` |
| `helene.mesh` | 13495 | 26986 | 18.074 | `d8fbff19525bca341bb3e4bde205ae77d8d05096f5244b77405ab27ebfd83e23` |

Pan 和 Daphnis 的标签引用 `pan_mst2018.bpc`、`daphnis_mst2018.bpc`；这些是晚期卡西尼形状数据。网格上的环脊及大尺度凹凸来自该模型，细小表面纹理和环粒子没有重建。各卫星背面及较少成像区域受原模型的观测覆盖限制。

Thomas 的这个归档没有 Polydeuces。场景用 [NAIF pck00011.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc) 的 `BODY634_RADII = (1.5 1.2 1.0)` km 画光滑三轴椭球；这不是 Thomas 网格，也不代表已重建的表面地形。

## Proteus

`proteus.mesh` 来源为 [Stooke Small Body Shape Models](https://sbn.psi.edu/pds/resource/stkshape.html)，作者 Philip Stooke；论文为 Stooke 1994, *The Surfaces of Larissa and Proteus*, EMP 65, 31–54，doi:10.1007/BF00572198。数据集引用：Stooke, P. 2025, *Stooke small bodies shape models V1.0*, NASA PDS，doi:10.26033/yt84-5y91。

许可：PDS 存档数据，公有领域，使用时请引用作者与数据集（[PDS 数据使用说明](https://pds-geosciences.wustl.edu/faq.html)）。[原始产品标签](https://sbnarchive.psi.edu/pds4/non_mission/small_bodies.stooke.shape-models/data/n8proteus.xml) 将观测来源列为 Voyager 2 ISSN，网格间隔为 5°，半径单位为 km，经度向西增加。

转换命令：

```powershell
python -X utf8 tools/shapes/make_plate_shape.py Proteus_Stooke.obj assets/shapes/proteus.mesh
```

保留 OBJ 的 2,522 个顶点和 5,040 个三角面，坐标单位 km，+z 朝北，+x 指向零经线，+y 指向东经 90°；未缩放或旋转。将 PDS `n8proteus.tab` 的西经转换为右手直角坐标后，OBJ 每个顶点与最近的表格点相差不超过 0.02 m。模型范围约 437.3 × 398.4 × 411.0 km，等体积半径约 201.00 km。场景的 `radii_km` 取 [NAIF pck00011.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc) 的 `BODY808_RADII = (218 208 201)` km，几何绘制使用网格本身。

OBJ 文件头注明由 MeshLab 导出。几何与 PDS 网格的核对支持上述来源；OBJ 没有附带纹理。顶点反照率统一为 1，场景灰色是显示用的示意值。5° 网格只描述大尺度形状，旅行者 2 号未充分成像的区域受原模型约束；同步自转省略物理天平动。

SHA-256：

- 输入 OBJ：`0f28d321156cd74e39167f3f83a0205cf39aab88f8df3e389655bcf04e289faa`
- PDS `n8proteus.tab`：`d079c8eac52d84224d5cb413526ff9415ab33afd48965e7e025544832ddd1868`
- 输出 `proteus.mesh`：`3425f7459d83a0b200f30d34a9f88d4d6d4e1747ebc7467394f786390306fe5a`

## Larissa

`larissa.mesh` 来自 Stooke 用旅行者 2 号图像重建的 Larissa 形状模型。论文与数据集同 Proteus，PDS 产品为 [`n7larissa.tab`](https://sbnarchive.psi.edu/pds4/non_mission/small_bodies.stooke.shape-models/data/)。

```powershell
python -X utf8 tools/shapes/make_plate_shape.py ref/neptune_moons/Larissa_Stooke.obj assets/shapes/larissa.mesh
```

MeshLab 导出的 OBJ 有 2,522 个顶点和 5,040 个三角面。坐标以 km 为单位原样写入。场景按 Proteus 的体轴约定解释 OBJ，+x 为零经线，+z 朝北。网格范围约为 201.1 × 192.6 × 177.7 km，等体积半径为 94.56 km。5° 网格描述大尺度轮廓，表面用场景里的示意灰色均匀着色。

- 输入 OBJ 的 SHA-256 为 `904016d2ad0a8b53e332ff2406e40ca7f2dbb7ae558bef59b4ea6dc4716af5fc`
- 输出网格的 SHA-256 为 `478666ec40174ed77ac25432ed986768b4498f4e8fb5be00a59d24b7c5c483ab`

## Homunculus（η Carinae，Steffen 2014 的 SHAPE 模型）

`homunculus.mesh` 由 ESO 的 3dmodel_004（OBJ）转换，坐标原样保留（单位是模型自己的，不是 km）：

```
python -E tools/shapes/make_plate_shape.py Eta_Car_Homunuculus_model.obj assets/shapes/homunculus.mesh
```

- 8,552 个顶点、17,101 个三角形，一个连通块，拓扑上是球面（χ = 2），有向体积为正（法线朝外）。有小缺陷：2 条边界边、3 条被三个三角形共用的边、1 个退化三角形，不影响按面积撒点
- 极轴就是模型的 z 轴（表面惯性主轴与 z 相差不到 1.5°），两瓣到 z = ±3.8；x 方向 −1.28 .. 2.72，y 方向 −1.80 .. 1.82
- 中心星不在原点：两瓣各截面的中心都在 x ≈ +0.85 的直线上，原点在两瓣连接的“腰”的一侧。网格关于 (0.75–0.85, 0, 0) 点对称得最好（每个顶点镜像后离网格平均 0.13 个单位，约为半径的 3%，Smith 2006 说真实的不对称就是 2–3%；以原点为中心时是 0.63）。Steffen 2014 说极区的沟是关于中心星点对称的，场景取 (0.80, 0, 0)。从这里到两极是 3.87 和 3.88 个单位（以原点为中心时是 4.03 和 3.83）
- 哪端是蓝瓣（朝向我们的东南瓣）：Steffen 2014 图 5 里红瓣上有一片“irregularities”，蓝瓣只有一个孔和一道沟；−z 那瓣凹凸明显多，所以 +z 是蓝瓣
- 两个凸起（表 1 的 Blue/Red Protrusion）：+z 瓣在方位角 48°（从 +x 向 +y 量）、离轴最远处高于腰 0.8 个单位，−z 瓣在 164°；两者相隔 116°（论文约 110°），平分线方位角约 106° 就是远星点方向。按“从地球看逆时针”（绕指向我们的蓝瓣轴）量经度，蓝凸起 −58°、红凸起 +58°，符号与表 1 的 −55°、+55° 一致，所以模型没有被镜像。两极附近凹得最深的地方也大致对得上表 1 的沟（蓝瓣在方位角约 280°，即经度 180°；红瓣在约 100°，即经度 0°），但离轴约 33°，表 1 是 15°，只算旁证
- 在天空中的朝向：论文没有给数字。轴的位置角和倾角用 Smith 2006 4.1 节的 PA 310°、i = 41°（极轴与视线的夹角）；绕轴的转角按论文图 4 左（“从地球看，北上东左”）对出来：远星点方向的投影位置角 330° 对应模型方位角约 88° 朝向我们，两个凸起的投影方向对应约 102°（图上两条虚线画成一条直线，110° 的夹角投影后不可能是 180°，那张图是示意的），取 95° ± 10°。按这个朝向从地球方向渲染，和图 5 上图一致（蓝凸起在两瓣连接处的右下，红凸起在上方）
- 独立核对：这样摆出的极轴和 Madura 2012 的轨道轴差 4.7°，凸起的平分线和轨道的远星点方向差 2.1°（`tests/stellar_system_tests.cpp`）
- 尺度：模型单位不是 10¹⁵ m（那样极半径约 25,400 au，比 Smith 2006 大 11%），看来是任意单位。场景按 Smith 2006 表 1 的极半径（纬度 89.1° 处 21,690 au，2005-03-01 的观测、2350 pc）定为每单位 5,600 au，再从 1847.1 年（Smith 2017）起均匀膨胀。Steffen 的形状在低纬度（20°–35°）比 Smith 2006 的胖 20–60%（他们说外形更接近圆柱），所以只用极半径标定
- 模型描的是 H₂ 层（尘埃在它里面一点），不是散射光的尘埃面；赤道的“裙边”在 H₂ 里描得不完整，没有 1890 年喷发的 Little Homunculus

SHA-256：

- 输入 `3dmodel_004.zip`：`d2dc1cfef4ebccbc8649aac0548a7708288b8e3c9a62e7bff22b66de9ebe717f`
- 输入 OBJ：`51af9b8ce70c71eb57f0118da17288f2ed1b12d9de21689b64d065ff921c64f4`
- 输出 `homunculus.mesh`：`96eaa06eb18132e7a9dd0c9b666d74492a4e47a742195cbe2ef3cf2e6234f433`

## 探测器模型

NASA 3D Resources（<https://github.com/nasa/NASA-3D-Resources>，提交 `11ebb4e`，2025-06-03）和 NASA 科学网站上单独发布的模型（GLB，DART 是 STL），用 `tools/shapes/make_spacecraft.py` 转换。原始文件不入库：仓库克隆到 `ref/NASA-3D-Resources`，单独下载的放在 `ref/probes/`。

| 输出 | 输入文件 | 页面上的署名 | 输入 SHA-256 |
|---|---|---|---|
| `voyager.mesh`、`../textures/voyager.jpg` | `Voyager Probe (B)/Voyager Probe (B).glb` | NASA/Michael D. Carbajal，<https://science.nasa.gov/3d-resources/voyager-probe-b/> | `bd86ded828dd3f459293aee4ffc3cd0998d8db67439317c8299650a1174c3289` |
| `jwst.mesh`、`../textures/jwst.png` | `James Webb Space Telescope (B)/James Webb Space Telescope (B).glb` | NASA/Christopher R. Meaney，<https://science.nasa.gov/3d-resources/james-webb-space-telescope-b/> | `4958e61e5a564f2efcc32cd516b9d36ea0b9bc644f64dbcfe182ce8a422ae526` |
| `parker.mesh`、`../textures/parker.jpg` | `Parker Solar Probe/Parker Solar Probe.glb` | NASA/Matthew J. Garcia，<https://science.nasa.gov/3d-resources/parker-solar-probe/> | `4d26e84efd29a78c5d07e474f6272bb111668b58da5a7555e305a65c72ddb1c5` |
| `cassini.mesh`、`../textures/cassini.jpg` | `Cassini-Huygens (A)/Cassini-Huygens (A) (without Hyugens).glb` | NASA/Brian E. Kumanchik，<https://science.nasa.gov/3d-resources/cassini-huygens-a/> | `87424f2c11b8e0c6a212292221e317396735f5e0876b7bb1204351998dd22f2c` |
| `juno.mesh`、`../textures/juno.jpg` | `Juno (A)/Juno (A).glb` | 页面没写署名，<https://science.nasa.gov/3d-resources/juno-a/> | `13f6bd48ff830a1890b21ef6e0182bb1043926feebb77af9733e7374b050367a` |
| `galileo.mesh`、`../textures/galileo.png` | `Galileo/Galileo.glb` | NASA/JPL/Solar System Simulator，<https://science.nasa.gov/3d-resources/galileo/> | `d192caaa4973c99187aadf9d4758ac286c564e7edb4d15305787579476eba2f7` |
| `pioneer10.mesh`、`../textures/pioneer10.jpg` | `Pioneer 10/Pioneer 10.glb` | NASA/JPL/Eyes on the Solar System，<https://science.nasa.gov/3d-resources/pioneer/> | `81443a7fa01cca0dda878cd27692b3307f9c332c4eb8e55a844499d7b350573c` |
| `chandra.mesh`、`../textures/chandra.jpg` | `Chandra X-ray Observatory/Chandra X-ray Observatory.glb` | NASA/Brian E. Kumanchik，<https://science.nasa.gov/3d-resources/chandra-x-ray-observatory/> | `9383f44fc4c6eba73c40b0c2cfe33c08c90f1fb5230274bfdaa98e0c9f116a34` |
| `roman.mesh`、`../textures/roman.png` | `Nancy Grace Roman Space Telescope (B)/Nancy Grace Roman Space Telescope (B).glb`（比建成的望远镜早的设计；(A) 是更早的 AFTA，30 万个三角形） | 页面没写署名，<https://science.nasa.gov/3d-resources/nancy-grace-roman-space-telescope-b/> | `23b01c53522452bf67608f908240c2e2f909700d2998c04153d025a275e28b71` |
| `spitzer.mesh`、`../textures/spitzer.png` | `Spitzer Space Telescope/Spitzer Space Telescope.glb` | NASA/Janice C. Lee; NASA/Robert L. Hurt，<https://science.nasa.gov/3d-resources/spitzer-space-telescope/> | `fa0b7c89fd2c8632d981d3a0ad471dc70620dda4e663788409b251e0f6e6f558` |
| `kepler.mesh`、`../textures/kepler.jpg` | `Kepler (A)/Kepler (A).glb` | NASA/Brian E. Kumanchik; NASA/Christian A. Lopez，<https://science.nasa.gov/3d-resources/kepler-a/> | `382014d90e2cf3e1f23e931d233547e652fe3629fddfb5acb8a99e2ceb72ff19` |
| `tess.mesh`、`../textures/tess.png` | `Transiting Exoplanet Survey Satellite (TESS) (A)/…(A).glb`（(B) 的机身是假的蓝色） | NASA/Christopher R. Meaney，<https://science.nasa.gov/3d-resources/transiting-exoplanet-survey-satellite-tess-a/> | `2b82a0191af6f06081368bb59d6cbb63b13ad8f019cd7ad08229e911e5b9ee6a` |
| `wmap.mesh`、`../textures/wmap.png` | `Wilkinson Microwave Anisotropy Probe (WMAP)/…(WMAP).glb` | NASA/Christopher R. Meaney，<https://science.nasa.gov/3d-resources/wilkinson-microwave-anisotropy-probe-wmap/> | `e229741cec91e370cbc8be74bf41b66eff4f164ca18b0e08dbd6d6ed4f59aeca` |
| `maven.mesh`、`../textures/maven.jpg` | `Mars Atmosphere and Volatile EvolutioN (MAVEN) (B)/…(B).glb`（(A) 有 175 万个三角形，没有纹理） | NASA/Brian E. Kumanchik，<https://science.nasa.gov/3d-resources/mars-atmosphere-and-volatile-evolution-maven-b/> | `cbb0d721a09bb8bee391d7a9eaa3edf12c7db4cb089465e9061113b9fa49d91c` |
| `mro.mesh`、`../textures/mro.jpg` | `Mars Reconnaissance Orbiter (MRO) (C)/…(C).glb` | NASA/Brian E. Kumanchik; NASA/Christian A. Lopez，<https://science.nasa.gov/3d-resources/mars-reconnaissance-orbiter-mro-c/> | `0f92872363ebe1984738caa4a5739863b125a1bc890b3c77bff23d8daa406747` |
| `odyssey.mesh`、`../textures/odyssey.jpg` | `Mars Odyssey/Mars Odyssey.glb` | NASA/JPL/Eyes on the Solar System，<https://science.nasa.gov/3d-resources/mars-odyssey/> | `ac0faa8dbfb54aabcbe3d06886c3d98dd29d120fd49b1c9fa6abbfbed22b049a` |
| `messenger.mesh`、`../textures/messenger.jpg` | `MErcury Surface, Space ENvironment, GEochemistry, and Ranging (MESSENGER)/…(MESSENGER).glb` | 页面没写署名，<https://science.nasa.gov/3d-resources/mercury-surface-space-environment-geochemistry-and-ranging-messenger/> | `c0268bfed0910e68a54c45925aea8edfc6e479d8131a849de7949252e01b9875` |
| `ladee.mesh`、`../textures/ladee.png` | `Lunar Atmosphere and Dust Environment Explorer (LADEE)/…(LADEE).glb` | NASA/Christopher R. Meaney，<https://science.nasa.gov/3d-resources/lunar-atmosphere-and-dust-environment-explorer-ladee/> | `92d5c7b782ab620ea0f7eb8e26d1905de43daeffdf6bcae8eaac878f123c59a0` |
| `themis.mesh`、`../textures/themis.jpg`（ARTEMIS-P1、P2） | `Time History of Events and Macroscale Interactions During Substorms (THEMIS)/…(THEMIS).glb` | NASA/JPL/Eyes on the Solar System，<https://science.nasa.gov/3d-resources/time-history-of-events-and-macroscale-interactions-during-substorms-themis/> | `2dfe2a2e07de211f6686adbb612c09c5dfac62a397cb88343f119c2124d818e3` |
| `rosetta.mesh`、`../textures/rosetta.jpg` | `Rosetta/Rosetta.glb` | 页面没写署名（ESA 的任务，页面的说明提到 ESA），<https://science.nasa.gov/3d-resources/rosetta/> | `92b3a6a2a24b4fce0ea1ba1ecc9b47d699a396422d096ef9e21b9ba97a5125c4` |
| `hubble.mesh`、`../textures/hubble.jpg` | `Hubble-1.glb`（2023-09-20 的版本，不在仓库里） | NASA Visualization Technology Applications and Development (VTAD)，<https://science.nasa.gov/resource/hubble-space-telescope-3d-model-2/> | `e6a40640d6846e8f5d9438821ebe2d7de730c8dadba1cc9dfdbcc148ff1052cd` |
| `new_horizons.mesh`、`../textures/new_horizons.jpg` | `New_Horizons.glb`（2023-09-09 的版本，不在仓库里） | NASA Visualization Technology Applications and Development (VTAD)，<https://science.nasa.gov/resource/new-horizons-3d-model/> | `cf152d8cea17c4a83711041cb0d29ca1d980bfafa084ada2656569ea92cdd73f` |
| `tgo.mesh`、`../textures/tgo.jpg` | `TGO.glb`（2023-09-09 的版本，不在仓库里） | NASA Visualization Technology Applications and Development (VTAD)，<https://science.nasa.gov/resource/trace-gas-orbiter-3d-model/> | `8d65a22cdbf29b9192506a7c31eafbdc25676db163849251e30683b45aaf536f` |
| `clipper.mesh`、`../textures/clipper.jpg` | `clipper_spacecraft.glb`（2023-09-09 的版本，不在仓库里） | NASA Visualization Technology Applications and Development (VTAD)，<https://science.nasa.gov/resource/europa-clipper-3d-model/> | `a3f5e4533b31bf4bb25f4332f6879055dfd744f3ff9add23bfa9bdfbefd4f7fe` |
| `dart.mesh`、`../textures/dart.png` | `3D Printing/Double Asteroid Redirection Test (DART)/…(DART).stl`（打印用，没有颜色） | 页面没写署名，<https://science.nasa.gov/3d-resources/double-asteroid-redirection-test-dart/> | `324725ab0595f0d69b4921e2d9521c42cceb0655b157a77b53813b9830d3eaf3` |
| `isee3.mesh`、`../textures/isee3.png` | `probes/isee3/` 下打印用的零件（2025-04-25 上传，没有颜色）：`Core 2.stl`、`Exboom1.stl`、`Exboom2.stl`、`Inertia boom.stl`、`Rad antenna.stl` | 页面没写署名，<https://science.nasa.gov/3d-resources/international-sun-earth-explorer-3/> | `04a9f9c0f6ecd9359a4b48640ec1fad03451ad2f265f8d1b28486a1ae1cdb3a6`、`b30c38e63363dd8e5b33cd4df68a4d84bea0244f64258994743783ac3950fa6b`、`4e49b9251c444fa3609acfa34748313c2aeefd9cded371fb430c1094bf3b4e66`、`75e70795cad2c1fa55affa04ed663ab53cc7f793df9355fdcbe804d7d1e2c811`、`2138bbf5ef241014d2a3a345204201476ef121e48788c1fd31897f7a4942fd73` |

- **许可**：仓库 README 写“free and without copyright”，并指向 NASA 媒体使用指南（<https://www.nasa.gov/nasa-brand-center/images-and-media/>）：渲染 3D 模型用的纹理和多边形数据在美国一般不受版权保护，但 NASA 的徽标和标识不在其内，也不能暗示 NASA 背书。这些页面（2026-11-07 查看）都没写单独的许可，也没有列第三方权利人。别的模型入库前要同样看一下页面上的署名：仓库的贡献者名单里有 Ball Aerospace、DigitalSpace 等外部机构。`Hubble Space Telescope (A)` 的页面署名就是 DigitalSpace Corporation（一家公司，版权归属不清楚），所以没有用；`(B)` 没有署名也没有材质，82 万个三角形，还有零件飘在外面。用的是 NASA 网站上 VTAD 署名的那个。
- **徽标**：Voyager 的纹理里有一个 NASA 徽标（像素矩形 828,440–882,480，贴在一块奶油色的板上），工具用周围的颜色涂掉。Hubble 的纹理里有一条横幅，印着 NASA 和 ESA 的徽标，只有 `NASA_logo` 这个贴花用它：工具跳过这个材质（`skip_materials`），并把横幅涂成灰色（`erase_polygon`），徽标不留在发布的图集里。Odyssey 的 `tex_01` 有 HEND 仪器的徽标和标签，`tex_boom` 是 Voyager 那张带 NASA 徽标的贴图，都涂掉。纯色的模型（JWST、Galileo、Roman、Spitzer、TESS、WMAP、LADEE、DART）没有纹理，其余的纹理里没有徽标（MESSENGER 的美国国旗也是真实的，留着；New Horizons 的美国国旗是探测器上真实的，留着；`tex_02` 里还有一小张没用到的探测器结构图，没有文字）。
- **处理**：仓库里的文件用 Draco 压缩，工具用 DracoPy（Apache-2.0，只在这个工具里用，不进程序）解码，Pillow 读图；New Horizons 没有压缩。去掉 Blender 导出时留下的占位方块（`_root`），按配置转到本体系（+z 是指向的轴）。Voyager 按 NASA 的 3.7 m 天线缩放（模型里 3.81 m），New Horizons 按 APL 的 2.1 m 天线缩放（模型里 1.888 m），Galileo 按 NASA 的 5.3 m 总高缩放（模型里 146.7 个单位，单位没写），Hubble 按 NASA 的 4.3 m 宽度缩放（模型里 4.52 m，镜筒也就是 13 m），Roman 按 NASA 的 12.7 m 长度（411.9 个单位），TESS 按 Horizons 说明里 3.9 m 的跨度（37.1 个单位），WMAP 按 Bennett 2003 的 5.1 m 圆盘（278.3 个单位），LADEE 按 NSSDCA 的 1.85 m 直径（机身 28.3 个单位），其余的尺寸不改。为控制文件大小，`merge_m` 把小于这个尺寸、法线相近的顶点合并（JWST、Chandra、Roman、Spitzer、TESS、WMAP、MAVEN、LADEE）。Hubble 的文件把同一张图存了三份，工具只放一份进图集。它的电池片（深色）朝模型的 −z，转到本体系里就是朝太阳的 +z；背面（浅色）是白色的散热器（ESA 和 NASA 都这样写）。New Horizons 的金箔在文件里是金属（近乎纯黄），这里只有漫反射，改成 sRGB (186, 140, 58)。Galileo 文件里的颜色是假的（叫 `bus:black` 的是亮紫色），改成深色和金色，半透明的金网画成不透明。Juno 有几个完全透明的材质，面和可见的面重合，跳过。JWST 用 2 cm 的网格合并顶点（法线相近的才合并），99,239 个三角形减到约 39,000 个（再加上薄片的背面）。

### ESA SciFleet 的模型

模型及原始纹理由欧洲空间局（ESA）提供，来源为 [ESA SciFleet](https://scifleet.esa.int/model-downloads)。采用 [CC BY](https://creativecommons.org/cc-licenses/)。发布的是 `tools/shapes/make_spacecraft.py` 转换后的网格和颜色图集。转换包含重定位、旋转、缩放、顶点合并和图集重采样，材质只保留底色，太阳定向的模型烘焙顶点自阴影。

| 输入文件 | 模型与标定尺寸 | 输入 SHA-256 |
|---|---|---|
| `bepi_mcs.glb` | BepiColombo，[30 m（MTM 翼展）](https://www.esa.int/Science_Exploration/Space_Science/BepiColombo/BepiColombo_factsheet) | `eaa9e016169c2d3a8ac0504b78f5d4d6122133fe81c0522440572295d3d6ab7a` |
| `giotto.glb` | Giotto，[2.848 m（总高）](https://www.esa.int/esapub/achievements/Sc72s4.pdf) | `842624dd18d57a4bd498755229800e67c7f0454642bc4e1d502458a120b2d984` |
| `herschel.glb` | Herschel，[7.5 m（总高）](https://www.esa.int/Science_Exploration/Space_Science/Herschel/Vital_stats) | `17d416a79dbc8fddd495e974cee13b6d9f2babb72c43e9e880a014c2c885c929` |
| `huygens.glb` | Huygens，[2.7 m（防热罩直径）](https://www.esa.int/Science_Exploration/Space_Science/Cassini-Huygens/Huygens_spacecraft) | `913879d43c0743c5becb6bd9150108f9231c81c2b88d9d0b4649b8a62e54e16d` |
| `juice.glb` | Juice，[27.1 m（翼展）](https://www.esa.int/articles/about-us/directorates/science/juice/juice-spacecraft-specs) | `367a7ea534f23d23894c12cd4ec38a71d1f33e205dddaabb683e25f7e9ce24b1` |
| `mars_express.glb` | Mars Express，[12 m（翼展）](https://www.esa.int/esaSC/SEMTV8374OD_0_spk.html) | `8d5901bf2ae337ab9fc55e26f11bd910467210a724d99483c2da39d63cc1c38e` |
| `philae.glb` | Philae，[约 1 m（着陆器机身宽度）](https://www.esa.int/About_Us/ESAC/Philae_found) | `30c4f3c2842a306eb6c4db53d606f0d709643b82bf638a5ac1785046bb3dc40d` |
| `planck.glb` | Planck，[4.22 m（底部宽度）](https://www.esa.int/Science_Exploration/Space_Science/Planck/Vital_stats) | `aabac98019a12c4514b0b843a08b19bbe8eced6ebdf39b23f54b14d39c65af18` |
| `solar_orbiter.glb` | Solar Orbiter，[18 m（翼展）](https://www.esa.int/Science_Exploration/Space_Science/Solar_Orbiter/Solar_Orbiter_factsheet) | `c542189e24a5190ffcc451a589845cb89db4ece55445fff3977c3c082d728bc7` |
| `venus_express.glb` | Venus Express，[8 m（翼展）](https://www.esa.int/Science_Exploration/Space_Science/Venus_Express/The_spacecraft) | `5734b72c9007a8cbfb99838728d5cec12643c70310a11d548adefd00e7bc2f0c` |

运行以下命令可重建这批资源。网格与图集同名，BepiColombo 的巡航组合体和单独的 MPO 共用 `bepi_mcs.png`。

Giotto 的图集去掉了天线反射面和机身上的 ESA 标志及其在小零件上的投影残片。

```powershell
python -E -P tools/shapes/make_spacecraft.py ref bepi_mcs bepi_mpo giotto herschel huygens juice mars_express philae planck solar_orbiter venus_express
```

### Gaia Sky 的模型

Gaia Sky 的探测器仓库（<https://gaia.ari.uni-heidelberg.de/gaiasky/repository/extra/spacecraft/>）里的 OBJ（带 MTL 和纹理），放在 `ref/probes/` 下各自的目录里，也用 `make_spacecraft.py` 转换（`Obj`：v 翻转，因为 OBJ 的纹理坐标从图的底部量起；Kd 按 Blender 导出的线性值读）。许可都按 CC BY，署名照各目录的 `credit.txt` / `credits.txt`。

| 输出 | 输入 | 署名（`credit.txt`） | 输入 SHA-256 |
|---|---|---|---|
| `gaia.mesh`、`../textures/gaia.jpg` | `mission-gaia/gaia.obj` | Toni Sagristà（Gaia Sky），CC BY | `66ddc8ef4ffd34d2751d0763dd9130a8c826099cca76b1f26f6bc60790a8a5d4` |
| `euclid.mesh`、`../textures/euclid.jpg` | `spacecraft-euclid/euclid.obj`、`Euclid2_Base_Color.png` | ESA，CC BY | `c138a35e182e6b70bfcb1906b9c7d7d198411184b1e610f887ece5f75f682526`、`cd4c48630ebd8c256304f68197ee98b71a3796ea9042e6f9a04c13f6b1a09117` |

- **Gaia**：遮阳罩在模型里宽 1 个单位，按 ESA 的 10.2 m（sci.esa.int 的 deployable sunshield 页）缩放，高约 3.5 m。很多纹理在一个三角形上平铺好几次，图集做不到，换成各自纹理的平均色；太阳能电池片按 8 × 8 格先铺好（`tile_images`，板子的 UV 最多跨 7 格）。纹理里没有徽标。
- **Euclid**：Scifleet 的模型，单张 4096 像素的图集，缩到 1/4。从服务舱底到遮阳罩顶 0.869 个单位，按 ESA 概览页的 4.7 m 缩放（下面的天线和推力器不算），宽度因此是 3.5 m，ESA 写 3.7 m（概览页）或收拢时的 3.1 m（fact sheet）。按 2 cm 和图集的 1/256 合并顶点，136,061 个三角形减到约 74,000 个。纹理里没有徽标。

### Orion（本项目自建）

| 输出 | 输入 | 来源 | 输入 SHA-256 |
|---|---|---|---|
| `orion.mesh`、`orion_crew_module.mesh`、`../textures/orion.png`（Artemis I、II） | `ref/artemis_orion/artemis_orion.obj`、`.mtl`（由同目录的 `build_orion.py` 在 Blender 里生成，`ref/` 不入库）；乘员舱的零件名单取自同一 Blender 文件导出的 `artemis_orion.glb`（`2d3dc13368fd7fb67ca813236efb8197bf154fcd381f74f5971a1644e191d5cb`） | 本项目为 Astraxis 生成的几何和材质，GPL-3.0-or-later；没有导入第三方网格、纹理或徽标 | `a93f3d6c703fd494495a54275974510f54b712e3b1f35da3dbb2067a84cb7fb7`、`248b2ca1f9978af691e0cb33d2ae3eaac8621cd2504a381505dc5efaf1d342c0`（`build_orion.py`：`90785e4ce87add1d3ffa9172b59c7fe245933de39823e5561062914c8f67140a`） |

- **外形**：太阳翼展开的在轨组合体：乘员舱、适配器、欧洲服务舱、主发动机和辅助发动机、四片各三块 2 m 面板的太阳翼。尺寸按 NASA “Orion by the Numbers” 2025（乘员舱直径约 5 m、高 3.35 m）、NASA “How to Draw Artemis / Orion”（布局、每片三段）、Airbus 的欧洲服务舱页面（每片约 7 m、展开约 19 m）和 ESA “Orion's wings”（电池片布局），2026-11-07 查阅；服务舱各段长度、喷口、隔热片接缝和电池片数量（每块 24 × 26 片，比实物少）是近似重建，不是工程 CAD。
- **处理**：米制，Z 轴从发动机指向舱顶，就是本体系的 +z；太阳翼平放在垂直于 +z 的平面里，电池片朝 +z，所以场景让 +z 朝太阳并烘焙阳光（乘员舱挡住翼根，翼背面在自己的阴影里）。原点在 7.4 m 的组合体中间。纯色材质（Kd 按 Blender 写的线性值），四种电池片颜色各占一个色板格，没有纹理图。按 2 cm 合并顶点，117,436 个三角形减到约 95,000 个。
- **乘员舱**：`orion_crew_module.mesh` 是和服务舱分离以后只剩的乘员舱，场景用 `[bodies.shape_after]` 在分离时刻换过去。Blender 文件把乘员舱的 237 个零件放在 `01 Return capsule (Crew module)` 集合的父对象 `Orion_Return_Capsule` 下；OBJ 没有层级，所以工具从 GLB 的节点树读出这个父对象下的网格名（`keep_glb_node`，空格换成下划线就是 OBJ 的对象名），网格只留这些对象，图集仍按整个模型做，和 `orion.mesh` 共用 `orion.png`。坐标系和原点不变，分离时乘员舱留在原处。防热罩是平的圆盘（模型在组合体里看不到底面，没有做出曲面），烘焙按同样的“+z 朝太阳”，约 19,000 个三角形。

### LICIACube（本项目自建）

| 输出 | 输入 | 来源 | 输入 SHA-256 |
|---|---|---|---|
| `liciacube.mesh`、`../textures/liciacube.png` | `ref/LICIACube/liciacube.obj`、`.mtl`（由同目录的 `build_liciacube.py` 在 Blender 里生成，`ref/` 不入库） | 本项目为 Astraxis 生成的几何和材质，GPL-3.0-or-later；没有导入第三方网格、纹理或标识 | `a203a3a7a29337726732804019025a62980b59537b8aec02dd9401d9ebce67cf`、`c99b911b86c3204beb0627931a33b3ebf4f6d43caf99be09dc93ee4156d87250`（`build_liciacube.py`：`abbd2046aaf82d14b1b4fc105cc7ef050ebdeacb0efffe08938db260feb04fe1`） |

- **外形**：太阳翼展开构型。机身按 INAF（“LICIACube, nanosatellite tutto italiano”，2023）的 6U、约 10 × 20 × 30 cm；两侧各三块太阳能板、端面的 LEIA 和 LUKE 两台相机、顶部的平板天线按 Argotec 的任务页面（argotecgroup.com/missions/liciacube）和 Poggiali 等 2022 图 1 估计，背面设备是示意的，不是工程尺寸。2026-11-07 查阅。
- **处理**：米制，电池片朝 +z，相机在 Blender 里朝 −y、在 OBJ 里朝 +y（导出时 `forward_axis='NEGATIVE_Y'` 让模型绕 z 转了 180°）。场景让相机对着 Dimorphos（本体 +z = OBJ 的 +y），电池片因此在本体 +y，太阳不在 +z，所以不烘焙阳光。全部是纯色材质，没有纹理图；15,612 个三角形，不合并顶点。

### MOM、Suisei（本项目自建）

| 输出 | 输入 | 来源 | 输入 SHA-256 |
|---|---|---|---|
| `mom.mesh`、`../textures/mom.png` | `ref/mom/mom.obj`、`.mtl`（`build_mom.py` 在 Blender 里生成，`ref/` 不入库） | 本项目为 Astraxis 生成的几何和材质，GPL-3.0-or-later；没有导入第三方网格、纹理或标识 | `cc5037cc59062d23d5effaf694672a977eef631aacdb60c1cb2f46b828159296`、`24ee8ff8191a8824192e3ce512f55b209a8877119a123582e690ff07e063046a`（`build_mom.py`：`664dc9357e620ca59a05012b4e07f099761ec23e35df3c447c6be9e9648308dd`） |
| `suisei.mesh`、`../textures/suisei.png`（Suisei 和 Sakigake） | `ref/suisei/suisei.obj`、`.mtl`（`build_suisei.py` 在 Blender 里生成，`ref/` 不入库） | 同上 | `2e50470552d963968feee2a0937251441e5eced0945e36ea3aa89748ace99032`、`277d9b582ec4f253693662a92e2a1f593717bcfaf1f25d40d81e61fe96194a29`（`build_suisei.py`：`8f39cba9a433253126a1b40c23de1bd14f46cad29424d226527a0d01709a7159`） |

- **MOM**：ISRO 的三块 1.8 × 1.4 m 太阳能板、2.2 m 天线；1.5 m 的机身和其余是估计；参考图片的权利归原作者，没有放进模型。米制，太阳能板在 OBJ 的 +y 一侧、电池片朝 +z，场景让 +z 朝太阳并烘焙阳光。原点在机身中间。按 1 cm 合并顶点，约 45,000 个三角形。
- **Suisei**：ISAS 的 1.4 × 0.7 m 圆筒、标称 0.8 m 的椭圆天线；其余按图片估计；参考图片的权利归原作者，没有放进模型。米制，自旋轴是 +z，天线在 +z 一端。场景让自旋轴在本体 +y（靠近黄道北），圆筒的一个侧面（OBJ 的 +x）朝太阳；它在自转，所以不烘焙阳光。天线的细网（管径 0.9 mm）合并顶点会塌掉，所以不合并，约 80,000 个三角形。
- **Sakigake** 没有单独的模型，用 Suisei 的网格和图集（两者圆筒相同，载荷不同）。

### 希望号（本项目自建）

| 输出 | 输入 | 来源 | 输入 SHA-256 |
|---|---|---|---|
| `hope.mesh`、`../textures/hope.png` | `ref/hope/hope.obj`、`.mtl`（`build_hope.py` 在 Blender 里生成，`ref/` 不入库） | 本项目为 Astraxis 生成的几何和材质，GPL-3.0-or-later；没有导入第三方网格、纹理或标识 | `0dcf89564dde3304d8c628c90cbbb5bbb12ba35f69e5559c9ded4eb6fea5dd1e`、`989f64ab19240b5fcfdf71e75f87e599d96fd4350eb198da469bce9ba260cd90`（`build_hope.py`：`4d2a9112118d44ef304a63ced4ffeaf5eff7c9811327c48176534c08d825de16`） |

- **外形**：太阳翼展开构型。按 Amiri 等 2022（“The Emirates Mars Mission”，Space Sci. Rev. 218, 4）第 6 节：两翼各两块太阳能板、1.85 m 高增益天线、六个主推进器、八个姿控推进器、同一块仪器板上的 EXI、EMIRS、EMUS；外观对照 The National 发表的发射前照片（图片的权利归原作者，没有放进模型）。翼展约 7.9 m、机身各边、仪器外壳和推进器位置是估计的，不是工程尺寸。
- **处理**：米制，机身长轴是 +z；电池片和天线在 Blender 里朝 −y，在 OBJ 里朝 +y（同 LICIACube），场景让它们朝太阳（本体 +z）并烘焙阳光，长轴因此是本体 +y。原点在长轴中间。按 1 cm 合并顶点，约 28,000 个三角形。

### 嫦娥五号 T1、天问一号（本项目自建）

| 输出 | 输入 | 来源 | 输入 SHA-256 |
|---|---|---|---|
| `ce5t1.mesh`、`../textures/ce5t1.png` | `ref/ce5t1/ce5t1.obj`、`.mtl`（`build_ce5t1.py` 在 Blender 里生成，`ref/` 不入库） | 本项目为 Astraxis 生成的几何和材质，GPL-3.0-or-later；没有导入第三方网格、纹理或标识 | `51777e4f170899fd3f4843596061121f01f0b83acaedf90cad9ef394bd0efd4f`、`3ddccb19322a4da59eadb22364ad770936856476822c180672a3f20e13c4753b`（`build_ce5t1.py`：`a23bea23d3f7dd50a6f5ab28006d19d51ae994d345d8eb1de2b6305df2d78ce3`） |
| `tianwen1.mesh`、`tianwen1_orbiter.mesh`、`../textures/tianwen1.png` | `ref/tianwen1/tianwen1.obj`、`.mtl`（`build_tianwen1.py` 在 Blender 里生成，`ref/` 不入库） | 同上 | `523177590663292fd3222ce1da9a4f36001715b39666f6d46e308579dc2d2b10`、`0d49258b8b556c4d0361a3b0b542be6db2c5e0d81fef9fa353b024a856e2f056`（`build_tianwen1.py`：`337846008cd9bab5db6850a6493c0c09831b32802166060fbabee56ee632647c`） |

- **外形**：都按图片建模（嫦娥五号 T1 两张，天问一号一张）。图上没有尺寸，模型的尺寸都是估计的：嫦娥五号 T1 服务舱约 1.68 × 1.48 × 1.52 m、翼展 8.64 m；天问一号机身约 2.65 × 2.25 × 1.35 m、进入舱底部直径约 3.34 m、翼展 14.1 m。没有对照工程资料，不能当实物尺寸引用。
- **嫦娥五号 T1**：场景里的表从返回舱分离之后（2014-10-31 22:19 TDB）开始，所以只画服务舱：工具按 OBJ 的对象名去掉返回舱的零件（`skip_objects`），留下它坐的适配器。米制，+z 从服务舱指向返回舱，电池片朝 +z，场景让 +z 朝太阳并烘焙阳光。按 1 cm 合并顶点，约 15,000 个三角形。
- **天问一号**：两份网格，同一张图集。`tianwen1.mesh` 是带着进入舱的巡航组合体，`tianwen1_orbiter.mesh` 是 2021-05-14 20:00 UTC 放下进入舱之后的轨道器（去掉进入舱、分离环和舱上的小零件），场景用 `[bodies.shape_after]` 在这个时刻换过去。`skip_objects` 只从网格里去掉零件，图集仍按整个模型做，所以两份网格的纹理坐标一样（`texture` 指向同一张图，第二次不再写）。米制，+z 朝进入舱，电池片朝 −z，场景让 −z 朝太阳（本体 +z）并烘焙阳光。Blender 里隔热膜的程序凹凸没有烘焙进 OBJ，看起来是平的。按 1 cm 合并顶点，约 36,000 和 34,000 个三角形。

### LCROSS、CAPSTONE、Akatsuki 与两艘 Pioneer（本项目自建）

运行 `python -E -P tools/shapes/make_spacecraft.py ref lcross capstone akatsuki pioneervo pioneer7`，为每艘探测器生成同名的 `.mesh` 和 `.png` 图集。输入是各自目录中的同名 OBJ、MTL，由 `build_<名称>.py` 在 Blender 中生成。几何和材质采用 GPL-3.0-or-later。各目录的 `SOURCES.md` 记录了建模资料与估算尺寸，`ref/` 由 Git 忽略。

| 名称 | 建模依据 | OBJ、MTL 的 SHA-256 |
|---|---|---|
| `lcross` | NASA 的 LRO/LCROSS 发射资料、飞行运行报告和 PDS 仪器标定图 | `e1e6e64b43c2da7a689d9542c392ea037b280a6bad232aa81275bc8f084015ad`、`f6a0554487ccd0b23fa3aa62d5d602f6538505111bc2d10bc184b7e05ac68350` |
| `capstone` | NASA 的装配资料和洁净室照片、Advanced Space 的飞行运行论文 | `18f0c1b15946e55ab888855bde54e011945fda6d4cb10bb02e362bbe46d3cf3e`、`241d9a20cb0e1c49c858a22c6697e9f20eed55bc7a0b3bbcdf1b8f7b13f01ac3` |
| `akatsuki` | JAXA/ISAS 构型图、NEC Technical Journal 6(1) 第 48 页图 3 | `e75ef4e2ab1a20ebb77fbbcb8d7e07483b55b478775b2072d02cb05844f3e069`、`aa368f1f3fa71bf9b13237b46d692ef000f00753614e8d7f800b78f760ca7159` |
| `pioneervo` | NASA SP-461 第 38–44 页、Ames 装配照和 JPL 42-28 图 1 | `08e9dda984cbe858bad5ed650cc2095d9aa3a5cf7139ffcf25d9b92b5d58414e`、`36afd0f93734f0042a964a6404f6c044f2d8ddbfe03bf570a00138eec8eed59a` |
| `pioneer7` | NASA SP-278 卷一第 50–52、58–61 页和 NASA/TRW 设计报告图 3–4 | `75921f86344b6636532ef134a7d6e9608f7f12a0feed57b3d6411c6866af7c1f`、`79f4a4f16b8f77b169700096da80c855f75110a82bebc5ecc0fa9a84f4b5495f` |

OBJ 使用米制。LCROSS 和 CAPSTONE 的电池面在 OBJ +y，Akatsuki 的在 OBJ +z。转换把这些面转到本体 +z，按场景的朝日方向烘焙自阴影。两艘 Pioneer 的圆筒轴从 OBJ +z 转到本体 +y，和 `point_at` 给出的朝日方向垂直，顶点的自阴影系数统一取 1。

### ISEE-3 的组装

NASA 的 ISEE-3 是给 3D 打印的散件（英寸，鼓 69.6 英寸宽 = 1.77 m），每根吊杆在自己的文件里、从原点沿 +z 伸出，没有装配图。`make_spacecraft.py` 的 `assembly` 按 JHU/APL 的飞行构型图（eoPortal “ISEE-3 / ICE” 图 3）和 NASA 计划文件（spacecollege.org/?p=245，2026-11-07 查看）装起来：鼓的轴（模型 y）是自转轴；四根刚性吊杆和四根导线天线都从鼓中间那条带（y 27.9–30.2 英寸）伸出，实验吊杆（`Exboom1` 是搜索线圈，`Exboom2` 是磁强计，标称 302 cm）装在鼓上有铰链座的两侧（模型 ±z，座的中心偏 7.15 英寸），惯性吊杆（`Inertia boom`，114 cm）装在与之垂直的两侧，导线天线（`Rad antenna` 拉伸到 49 m）在它们之间的 45° 方向，轴向天线（同一个零件拉伸到 701 cm）沿自转轴上下各一根。文件自己的 `Rad antenna` 是 55 m。`Exboom (empty).stl` 是没有探头的空吊杆，`Thruster.stl` 是一个装不到鼓上任何推力器座的散件，都不用。鼓是打印用的空筒，两端补上 16 边形的盖。太阳能电池覆盖侧面，除了中间那条带和上下两道环（按面的三个顶点所在的高度分，大的侧面三角形跨过中间的带，所以带只在少数面上看得到）。颜色都是示意的。它每 3 秒转一圈，所以不烘焙阳光。

- **TGO、Europa Clipper、DART、Huygens**：TGO 和 Clipper 没有压缩，尺寸不改（17.5 m、30.5 m，和 Horizons 的说明一致）；纹理里没有徽标。TGO 的金属材质在这里的漫反射下近乎纯白、纯黄，毯子的纹理又几乎全黑（1024 像素的图上只有一些金点），都换成纯色（颜色是这里选的，见 `make_spacecraft.py` 的 `colors`）。Clipper 的太阳能板每个 0.56 × 0.79 m 的矩形把一张电池片的图平铺 7 × 18 次，图集不能平铺，工具把电池片缩到 24 像素、铺成 7 × 19 格放进图集（`tile_images`），亮度乘 2.6（原图近乎全黑，在天空上看不出来）；雷达天线上有几十个三角形的 UV 也跨过多个平铺，这些被截断，天线太细看不出来。为控制大小，Clipper 按 3 cm 和图集的 1/128 合并顶点（`merge_m`、`merge_uv`），203,662 个三角形减到约 114,000 个。DART 的 STL 有 446,397 个三角形、348 个互不相连的零件，单位不明：跨度 40.34 个单位，按场景的 19 m 缩放后机身宽 1.2 m，和 APL 的 1.2 × 1.3 × 1.3 m 一致；法线按 40° 的折痕角平滑，按 3 cm 合并后约 90,000 个三角形。它没有材质，工具按零件的外接盒归类（电池片、吊杆、散热板、毯子、天线、DRACO、离子发动机、其他），颜色都是示意的（没有找到能据以取色的公开资料）。两翼是同一个形状绕望远镜的轴转 180°，文件看不出电池片在哪一面，这里都放在模型的 +y。Huygens 是进入舱（防热罩和背罩），直径 2.615 m，按 2.7 m 缩放；它朝向 Titan 而不是太阳，所以不烘焙阳光（`bake_sun: False`），顶点通道都是 1。
- **双面**：这些模型有很多不封闭的薄片，程序画图集网格时不剔除背面，从背面看时翻转法线，所以文件里每个面只存一次。
- **图集**：每个探测器一张纹理。模型的图放在上面（周围留 4 像素的重复边），纯色材质各占一个 8×8 的色板格（glTF 的颜色是线性的，转成 sRGB）。伸出 0–1 的 UV 按整数移回来（Parker 的箔材质）。
- **烘焙的阳光**：顶点的“反照率”通道存太阳在 +z 时能否照到（从顶点沿 +z 的光线和整个网格求交，法线背朝太阳的顶点也算，因为薄片的另一面可能朝着太阳），所以只适用于始终用同一面朝太阳的探测器。阳光在顶点之间线性插值，低面数的大三角形要先细分（`subdivide_m`，Galileo 的遮阳板是一个大扇形，中心顶点被挡住时整块都发黑）。

NASA 科学网站上有同一模型的 STL（“Eta Carinae Homunculus Nebula”，835 KB），页面没写许可，没有用。同一网站上 Madura 2015 的“High Mass-Loss Rate”打印模型（碰撞星风在远星点、近星点和相位 1.045 的三个快照，CC BY 4.0）是 Madura 2013 的 SPH 模拟，以后做碰撞星风时可以拿来核对。
