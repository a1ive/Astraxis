# 形状模型来源

场景中设置 `shape = "shapes/<名称>.mesh"` 后，程序用三角网格绘制不规则天体。网格采用天体固连坐标系，x 指向本初子午线，z 沿自转轴正极，单位为 km。每个顶点保存相对反照率，渲染时与场景的 `color` 相乘。

| 文件 | 来源 | 许可 |
|---|---|---|
| `arrokoth.mesh` | Porter, S. et al. 2024, “New Horizons Porter (2024) Arrokoth Shape Model Collection”，PDS Small Bodies Node，doi:[10.26007/97r3-1e19](https://doi.org/10.26007/97r3-1e19)，<https://pdssbn.astro.umd.edu/holdings/pds4-nh_derived-v3.0/arrokoth_shapemodel_porter2024/>；方法见同目录的 `porteretal2024b.pdf`（Porter et al. 2024, “The Shape of (486958) Arrokoth”） | NASA 新视野号项目的 PDS 数据，公有领域（请引用作者） |

原始文件（同一目录下）：

- `arrokoth_porter_2024_v01.obj`（3.7 MB）：20,484 个顶点、40,960 个三角面，原点为质心，轴为惯量主轴。
- `albedo_arrokoth4_fp36h2_masked1.png`（0.2 MB）：反照率图，OBJ 的纹理坐标指向它。

`tools/shapes/make_arrokoth.py` 将原始模型转换为 800 KB 的网格文件，依赖 Pillow。文件格式（小端）：`char[8] "AXMESH1\0"`、`uint32` 顶点数、`uint32` 三角形数，之后是每个顶点的 `float32` x、y、z（km），每个顶点的 `float32` 相对反照率，每个三角形 3 个 `uint32` 顶点索引（从外面看逆时针）。法线在加载时按面积加权平均算出。

## Arrokoth

- 形状：两瓣接触双星，大瓣 Wenu 在 +x 方向，小瓣 Weeyo 在 −x 方向。整体 34.546 × 19.838 × 13.822 km，与论文表 2 一致。网格体积 4123.8 km³，换算成等体积球直径 19.896 km，也与表 2 一致（`tests/sim_tests.cpp` 的 `test_arrokoth_shape`）。模型标签注明两瓣在连接处略有重叠。
- 极轴与本初子午线：用新视野号项目的 PCK `nh_arrokoth_002.tpc`（<https://naif.jpl.nasa.gov/pub/naif/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/pck/>）：极轴 RA 319.3685820°、Dec −25.5875063°，W = 71.9541475° + 559.3936173°/天（自 J2000 TDB）。文件注释注明 2024 版形状模型使用这组数值。极轴与论文表 2（319.37°、−25.588°）一致。
  - 形状模型标签和 `arrokoth_geophysics` 概述里写的极轴（317.488°、−24.888°）和周期（0.6632553 天 = 15.92 h）是 2021 年的参数；本模型的朝向按上述 PCK 参数计算。
  - 由此推出的周期是 15.45 h，比 Spencer et al. 2020 的 15.92 h 短 3%。New Horizons 在飞掠前后两天拍到了可分辨的图像，观测时段较短，对周期的约束较弱。飞掠时刻的朝向按拟合时使用的这组参数计算。
  - 交叉检查：极轴与 Arrokoth 日心轨道法向的夹角算得 100.39°（表 2：100.39°），与新视野号来向的夹角 41.04°（表 2：41.096°）。
- 反照率：PNG 是 FITS 数组上下翻转后的图像，按显示方向存储。OBJ 纹理坐标 v = 0 对应图像最下一行，符合 OBJ 的惯例。按这个方向采样，亮的“颈部”落在两瓣连接处。像素换算为 albedo = 像素 / 1361975.975 + 0.03188（标签）。新视野号拍到的部分约占 40% 的顶点，主要在南半球。其余顶点对应固定的掩膜值，转换工具从相邻已知顶点沿网格扩散填充，远处渐变到平均值。因此背面显示为平淡的表面。文件存储反照率与已知部分平均值（0.0438）的比值，范围为 0.74 至 1.83。
- 颜色：场景使用示意颜色 `#a8604a`。Arrokoth 偏红，LORRI 是全色相机，图像本身无法给出这里使用的颜色。
- 光照采用 Lambert 漫反射，自阴影尚未实现。太阳很低时，颈部的凹处也会被照亮。
