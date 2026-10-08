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
| `pallas.mesh` | VLT/SPHERE 大型观测计划的 MPCD 形状模型（Marsset et al. 2020, Nature Astronomy 4, 569；Vernazza et al. 2021, A&A 654, A56），Sketchfab 用户 “marin”（@marin14）上传的 “(2) Pallas - MPCD model”，<https://sketchfab.com/3d-models/2-pallas-mpcd-model-59c1ec418b2f4f13bc182633a2cddc66>，文件 `2_Pallas_mpcd.obj`（2022-06-13） | CC BY 4.0（Sketchfab 页面标注；署名：marin (@marin14) / Marsset et al. 2020, Vernazza et al. 2021） |
| `vesta.mesh` | USGS Astrogeology / DLR，Vesta Dawn FC HAMO Global DTM 93m（2013-12-10 发布，黎明号分幅相机 HAMO 立体图像），<https://astrogeology.usgs.gov/search/map/vesta_dawn_fc_hamo_global_dtm_93m>，文件 `Vesta_Dawn_HAMO_DTM_DLR_Global_48ppd.tif` | 公有领域（请引用作者） |
| `phobos.mesh`、`deimos.mesh` | P. C. Thomas，“Small Body Optical Shape Models” V1.0（PDS SBN，`ast-sat.thomas.shape-models`，由 PDS3 数据集 EAR-A-5-DDR-SHAPE-MODELS-V2.1 迁移），<https://sbnarchive.psi.edu/pds4/non_mission/ast-sat.thomas.shape-models_V1_0/>，文件 `data/m1phobos.tab`、`data/m2deimos.tab`；Thomas 1993, Icarus 105, 326 | PDS 存档数据，公有领域（请引用数据集与作者） |
| `arrokoth.mesh` | Porter, S. et al. 2024, “New Horizons Porter (2024) Arrokoth Shape Model Collection”，PDS Small Bodies Node，doi:[10.26007/97r3-1e19](https://doi.org/10.26007/97r3-1e19)，<https://pdssbn.astro.umd.edu/holdings/pds4-nh_derived-v3.0/arrokoth_shapemodel_porter2024/>；方法见同目录的 `porteretal2024b.pdf`（Porter et al. 2024, “The Shape of (486958) Arrokoth”） | NASA 新视野号项目的 PDS 数据，公有领域（请引用作者） |

原始文件（同一目录下）：

- `arrokoth_porter_2024_v01.obj`（3.7 MB）：20,484 个顶点、40,960 个三角面，原点为质心，轴为惯量主轴。
- `albedo_arrokoth4_fp36h2_masked1.png`（0.2 MB）：反照率图，OBJ 的纹理坐标指向它。

`tools/shapes/make_arrokoth.py` 将原始模型转换为 800 KB 的网格文件，依赖 Pillow。

文件格式见 `tools/shapes/axmesh.py`（小端）：`char[8] "AXMESH2\0"`、`uint32` 顶点数、`uint32` 三角形数、`uint32` 标志（第 0 位：带地图坐标），之后是每个顶点的 `float32` x、y、z（km），每个顶点的 `float32` 相对反照率，（有标志时）每个顶点的 `float32` 地图 u（东经 / 360°，接缝处的顶点重复，取 0 和 1），每个三角形 3 个 `uint32` 顶点索引（从外面看逆时针）。法线在加载时按面积加权平均算出，位置完全相同的顶点（接缝、极点）共用一个法线。带地图坐标的网格可以在场景里同时写 `texture`：u 取自文件，纬度取自顶点方向（行星中心纬度）。

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
- 误差（标签）：图像质量较差，边缘和明暗界线的位置本身就不确定，每个点的绝对误差估计为 0.5至1 km，相邻点之间约 100 m；凸包的可靠程度估计与模型本身相当。Stooke 自己的比较认为他的模型偏多棱角、凹陷偏深。
- 交叉检查：Szegő 1991（Comets in the Post-Halley Era 2, 723，第 2 节）的 Vega/Giotto 独立模型装在 15.3 × 7.2 × 7.22 km 的盒子里，体积 365 km³。本网格沿 z 长 15.14 km，x、y 方向 7.53、7.59 km，体积 402.2 km³（等体积半径 4.58 km），差约 10%（`test_halley_shape`）。场景的 `radii_km` 用 Szegő 的盒子尺寸的一半（3.61、3.6、7.65 km，长轴沿 z），只用于影子和相机距离。
- 自转：Belton et al. 1991 的摘要给出长轴模式（对称陀螺的自由进动）：长轴与角动量 M 夹角 66.0°，绕 M 进动周期 3.69 天，绕长轴自转周期 7.1 天，合成的总自转周期 2.84 天，总角速度与 M 夹角 21.4°（这两个数可以由前三个算出，`test_halley_shape` 检查）；M 指向 B1950 赤经 6.2°、赤纬 −60.7°，用 SPICE 的 FK4 → J2000 旋转换算为 J2000 的 6.785°、−60.423°。场景用 `[bodies.free_precession]` 按 z-x-z 欧拉角计算朝向。摘要里没有两个角的相位，论文正文需要付费获取，所以相位取 0（历元 1986-03-14 00:00 UTC），任一时刻的具体朝向是示意性的；也不知道 Belton 的长轴方向与 Stooke 的“大头”是否一致，这里假定 +z 与 M 成 66°。
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
