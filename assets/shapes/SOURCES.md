# 形状模型来源

场景中设置 `shape = "shapes/<名称>.mesh"` 后，程序用三角网格绘制不规则天体。网格采用天体固连坐标系，x 指向本初子午线，z 沿自转轴正极，单位为 km。每个顶点保存相对反照率，渲染时与场景的 `color` 相乘。

| 文件 | 来源 | 许可 |
|---|---|---|
| `vesta.mesh` | USGS Astrogeology / DLR，Vesta Dawn FC HAMO Global DTM 93m（2013-12-10 发布，黎明号分幅相机 HAMO 立体图像），<https://astrogeology.usgs.gov/search/map/vesta_dawn_fc_hamo_global_dtm_93m>，文件 `Vesta_Dawn_HAMO_DTM_DLR_Global_48ppd.tif` | 公有领域（请引用作者） |
| `phobos.mesh`、`deimos.mesh` | P. C. Thomas，“Small Body Optical Shape Models” V1.0（PDS SBN，`ast-sat.thomas.shape-models`，由 PDS3 数据集 EAR-A-5-DDR-SHAPE-MODELS-V2.1 迁移），<https://sbnarchive.psi.edu/pds4/non_mission/ast-sat.thomas.shape-models_V1_0/>，文件 `data/m1phobos.tab`、`data/m2deimos.tab`；Thomas 1993, Icarus 105, 326 | PDS 存档数据，公有领域（请引用数据集与作者） |
| `arrokoth.mesh` | Porter, S. et al. 2024, “New Horizons Porter (2024) Arrokoth Shape Model Collection”，PDS Small Bodies Node，doi:[10.26007/97r3-1e19](https://doi.org/10.26007/97r3-1e19)，<https://pdssbn.astro.umd.edu/holdings/pds4-nh_derived-v3.0/arrokoth_shapemodel_porter2024/>；方法见同目录的 `porteretal2024b.pdf`（Porter et al. 2024, “The Shape of (486958) Arrokoth”） | NASA 新视野号项目的 PDS 数据，公有领域（请引用作者） |

原始文件（同一目录下）：

- `arrokoth_porter_2024_v01.obj`（3.7 MB）：20,484 个顶点、40,960 个三角面，原点为质心，轴为惯量主轴。
- `albedo_arrokoth4_fp36h2_masked1.png`（0.2 MB）：反照率图，OBJ 的纹理坐标指向它。

`tools/shapes/make_arrokoth.py` 将原始模型转换为 800 KB 的网格文件，依赖 Pillow。

文件格式见 `tools/shapes/axmesh.py`（小端）：`char[8] "AXMESH2\0"`、`uint32` 顶点数、`uint32` 三角形数、`uint32` 标志（第 0 位：带地图坐标），之后是每个顶点的 `float32` x、y、z（km），每个顶点的 `float32` 相对反照率，（有标志时）每个顶点的 `float32` 地图 u（东经 / 360°，接缝处的顶点重复，取 0 和 1），每个三角形 3 个 `uint32` 顶点索引（从外面看逆时针）。法线在加载时按面积加权平均算出，位置完全相同的顶点（接缝、极点）共用一个法线。带地图坐标的网格可以在场景里同时写 `texture`：u 取自文件，纬度取自顶点方向（行星中心纬度）。

## Arrokoth

- 形状：两瓣接触双星，大瓣 Wenu 在 +x 方向，小瓣 Weeyo 在 −x 方向。整体 34.546 × 19.838 × 13.822 km，与论文表 2 一致。网格体积 4123.8 km³，换算成等体积球直径 19.896 km，也与表 2 一致（`tests/sim_tests.cpp` 的 `test_arrokoth_shape`）。模型标签注明两瓣在连接处略有重叠。
- 极轴与本初子午线：用新视野号项目的 PCK `nh_arrokoth_002.tpc`（<https://naif.jpl.nasa.gov/pub/naif/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/pck/>）：极轴 RA 319.3685820°、Dec −25.5875063°，W = 71.9541475° + 559.3936173°/天（自 J2000 TDB）。文件注释注明 2024 版形状模型使用这组数值。极轴与论文表 2（319.37°、−25.588°）一致。
  - 形状模型标签和 `arrokoth_geophysics` 概述里写的极轴（317.488°、−24.888°）和周期（0.6632553 天 = 15.92 h）是 2021 年的参数；本模型的朝向按上述 PCK 参数计算。
  - 由此推出的周期是 15.45 h，比 Spencer et al. 2020 的 15.92 h 短 3%。New Horizons 在飞掠前后两天拍到了可分辨的图像，观测时段较短，对周期的约束较弱。飞掠时刻的朝向按拟合时使用的这组参数计算。
  - 交叉检查：极轴与 Arrokoth 日心轨道法向的夹角算得 100.39°（表 2：100.39°），与新视野号来向的夹角 41.04°（表 2：41.096°）。
- 反照率：PNG 是 FITS 数组上下翻转后的图像，按显示方向存储。OBJ 纹理坐标 v = 0 对应图像最下一行，符合 OBJ 的惯例。按这个方向采样，亮的“颈部”落在两瓣连接处。像素换算为 albedo = 像素 / 1361975.975 + 0.03188（标签）。新视野号拍到的部分约占 40% 的顶点，主要在南半球。其余顶点对应固定的掩膜值，转换工具从相邻已知顶点沿网格扩散填充，远处渐变到平均值。因此背面显示为平淡的表面。文件存储反照率与已知部分平均值（0.0438）的比值，范围为 0.74 至 1.83。
- 颜色：场景使用示意颜色 `#a8604a`。Arrokoth 偏红，LORRI 是全色相机，图像本身无法给出这里使用的颜色。
- 光照采用 Lambert 漫反射，自阴影尚未实现。太阳很低时，颈部的凹处也会被照亮。

## Phobos 与 Deimos（Thomas 形状模型）

- 形状：Thomas 的数值形状模型，只用海盗号轨道器的图像，按行星中心纬度 / 经度的网格给出半径（0° 和 360° 两列都有）。`tools/shapes/make_thomas_shape.py`（只用标准库）把网格原样转成网格文件，去掉极点处退化的三角形
  - Phobos：2° 网格（91 × 181 个点），16,471 个顶点、32,040 个三角形，697 KB
  - Deimos：5° 网格（37 × 73 个点），2,701 个顶点、5,040 个三角形，111 KB。标签注明 200°–355° 经度（西经）一带的误差约 400 m
- **经度方向**：表里的经度是**西经**。依据是 Phobos 表：减去拟合的三轴椭球后，49° 处有清楚的陨石坑特征（坑内低约 0.7 km、外圈高起），311° 处没有；Stickney 坑在 49°W。网格里换成东经（u = 东经 / 360°）。形状里 Stickney 凹陷的中心在 2°S 50°W；网格里该处 12° 以内的平均半径比外圈低 1.12 km，镜像位置 50°E 没有（`test_mars_moon_shapes`）
- 体积与 JPL SSD 卫星物理参数表的平均半径（Archinal et al. 2018）对照（`test_mars_moon_shapes`）：
  - Phobos：5748.7 km³，等体积半径 11.11 km，表中为 11.08 ± 0.04 km
  - Deimos：1013.5 km³，等体积半径 6.23 km，表中为 6.2 ± 0.24 km
- 两颗卫星的最长轴都指向火星（本初子午线），各轴的范围接近场景里 PCK 的三轴椭球（Deimos 沿 x / y / z 离中心最远约 8.4 / 6.7 / 5.8 km，PCK 为 7.8 × 6.0 × 5.1 km）。场景的 `radii_km` 仍用 PCK 值（影子、相机距离）
- 自转参数：Deimos 的标签给出极轴 RA 316.65°、Dec 53.53°，W = 79.41° + 285.1618970°/天，与场景所用 PCK 的值一致，所以本初子午线相同
- 纹理是同一数据集的拼接图 `m1phobosm.fit`、`m2deimosm.fit`（见 `assets/textures/SOURCES.md`），拼接时的位置就是用这些形状模型控制的。Phobos 原先用的 Stooke 拼接图（USGS）基于 DLR 的另一套控制网，地貌比形状模型偏西约 5°（Limtoc 在 12°S 59°W，Thomas 图和 IAU 地名为 11°S 54°W 左右），偏差各处不同，贴在形状上会错开，所以换掉了
- Deimos 的 5° 网格很粗，近看外形是多面体状的；两颗卫星的光照都没有自身阴影

## Vesta（DLR HAMO DTM）

- 原始文件：<https://asc-pds-services.s3.us-west-2.amazonaws.com/mosaic/Vesta_Dawn_HAMO_DTM_DLR_Global_48ppd.tif>（597 MB，17280 × 8640，32 位浮点，未压缩、每行一个条带）。像素值是从 Vesta 中心起算的半径（米），USGS 页面给出的平均交会误差为 ±8 m。简单圆柱投影、行星中心纬度、东经向右、左边缘 −180°E；经度是 Claudia″ 系统（Claudia 坑在 146°E），与纹理 `vesta.jpg` 和场景所用 PCK（W0 = 285.39°）相同。DTM 与纹理同属 DLR 2013-12-10 的发布（ProductId 20131210），控制网一致
- USGS 页面说 HAMO 立体图像覆盖约 95% 的表面，但文件里没有无数据像素（半径 211.9–293.0 km），未覆盖的部分在产品里已经补上了
- `tools/shapes/make_dtm_shape.py`（只用 Pillow 读 TIFF 的条带位置，数据用标准库按行读取）按 1.5° 的格子求平均，每个网格顶点取周围四格的平均，极点取整圈格子的平均：121 × 241 个点，29,161 个顶点、57,120 个三角形，1.2 MB。1.5° 约合 6.9 km，Rheasilvia 盆地（约 500 km）和它的中央峰都分辨得出，小坑的细节靠纹理
- 交叉检查：体积 7.495 × 10⁷ km³，等体积直径 523.11 km，SBDB 为 522.77 ± 0.1 km（Park et al. 2025，Nature Astronomy）；差 0.34 km，来自 2013 年的 DTM 和网格的平滑（`test_vesta_shape`）。SBDB 的外形尺寸 569.24 × 554.48 × 452.66 km 应是拟合椭球的轴长，网格的包围盒是 574.5 × 554.0 × 467.5 km（包含地形起伏），不直接比较
- 场景的 `radii_km` 仍用 PCK 的三轴椭球（影子、相机距离）；光照没有自身阴影
