# 纹理来源

天体表面纹理采用简单圆柱投影（等经纬度）的全球地图。Eris、Haumea、Makemake 的模型纹理由 `tools/textures/extract_dwarf_textures.py` 转换，做法见下文。其余观测地图由 `tools/textures/prepare_textures.py` 先按整数倍降采样，再用 Lanczos 重采样，将图像缩放到 2048×1024。工具还会填补 USGS 拼接图的极区无数据像素，最后存为 JPEG（质量 90）。原始文件体积较大，仓库中收录的是处理后的纹理，需要时可按下表重新下载。

| 文件 | 来源 | 许可 |
|---|---|---|
| `jupiter.jpg` | NASA/JPL/Space Science Institute，PIA07782 “Cassini's Best Maps of Jupiter: Cylindrical Map”（Cassini 窄角相机，2000-12-11/12），<https://science.nasa.gov/photojournal/cassinis-best-maps-of-jupiter-cylindrical-map> | NASA 图像，公有领域 |
| `io.jpg` | USGS Astrogeology，Io Galileo SSI / Voyager Color Merged Global Mosaic 1km，<https://astrogeology.usgs.gov/search/map/io_galileo_ssi_voyager_color_merged_global_mosaic_1km> | 公有领域（Use Constraints: None） |
| `europa.jpg` | USGS Astrogeology，Europa Voyager - Galileo SSI Global Mosaic 500m，<https://astrogeology.usgs.gov/search/map/europa_voyager_galileo_ssi_global_mosaic_500m> | 公有领域 |
| `ganymede.jpg` | USGS Astrogeology，Ganymede Voyager - Galileo SSI Global Mosaic 1km，<https://astrogeology.usgs.gov/search/map/ganymede_voyager_galileo_ssi_global_mosaic_1km> | 公有领域 |
| `earth.jpg` | NASA Earth Observatory，Blue Marble Next Generation（2004 年 7 月，Reto Stöckli），`world.200407.3x5400x2700.jpg`，<https://eoimages.gsfc.nasa.gov/images/imagerecords/74000/74092/world.200407.3x5400x2700.jpg> | NASA 图像，公有领域（请注明 NASA Earth Observatory） |
| `moon.jpg` | NASA Scientific Visualization Studio，CGI Moon Kit（LRO 数据），`lroc_color_poles_2k.tif`，<https://svs.gsfc.nasa.gov/4720> | NASA 图像，公有领域（请注明 NASA's Scientific Visualization Studio）；页面说明是“为美观而非科学优化” |
| `callisto.jpg` | USGS Astrogeology，Callisto Galileo/Voyager Global Mosaic 1km，<https://astrogeology.usgs.gov/search/map/callisto_galileo_voyager_global_mosaic_1km> | 公有领域 |
| `mimas.jpg` | DLR / USGS Astrogeology，Mimas Cassini ISS 全球底图（2017-06-30，`MI_170630_DLR_basemap_degrees.tif`，在 `Cassini_DLR_Mimas.zip` 里） | 公有领域 |
| `enceladus.jpg` | USGS Astrogeology，Enceladus Cassini ISS Global Mosaic 100m HPF（高通滤波，反照率被压平），<https://astrogeology.usgs.gov/search/map/enceladus_cassini_iss_global_mosaic_hpf_110m> | 公有领域 |
| `tethys.jpg` | USGS Astrogeology，Tethys Cassini Global Mosaic 293m | 公有领域 |
| `dione.jpg` | USGS Astrogeology，Dione Cassini/Voyager Global Mosaic 154m | 公有领域 |
| `rhea.jpg` | USGS Astrogeology，Rhea Cassini/Voyager Global Mosaic 417m | 公有领域 |
| `iapetus.jpg` | USGS Astrogeology，Iapetus Cassini/Voyager Global Mosaic 783m，<https://astrogeology.usgs.gov/search/map/iapetus_cassini_voyager_global_mosaic_803m> | 公有领域 |
| `ariel.jpg` / `umbriel.jpg` / `titania.jpg` / `oberon.jpg` / `miranda.jpg` | USGS 用旅行者 2 号图像拼接的全球图（灰度，1440×720），由 JPL Solar System Simulator 发布：<https://space.jpl.nasa.gov/tmaps/uranus.html>（`ura1vuu2.tif` .. `ura5vuu2.tif`） | USGS 作品，公有领域；页面署名 Caltech/JPL/USGS |
| `triton.jpg` | USGS Astrogeology，Triton Voyager 2 Global Color Mosaic（`Triton_Voyager2_ClrMosaic_GlobalFill_600m.tif`） | 公有领域 |
| `pluto.jpg` | USGS Astrogeology，Pluto New Horizons LORRI/MVIC Global Mosaic 300m（2017-07，`Pluto_NewHorizons_Global_Mosaic_300m_Jul2017_8bit.tif`，灰度） | 公有领域 |
| `charon.jpg` | USGS Astrogeology，Charon New Horizons LORRI/MVIC Global Mosaic 300m（2017-07，`Charon_NewHorizons_Global_Mosaic_300m_Jul2017_8bit.tif`，灰度） | 公有领域 |
| `mercury.jpg` | USGS Astrogeology，Mercury MESSENGER MDIS Global Mosaic 250m（2013-05，750 nm 单色，NAC/WAC，含两极平均拼接图），<https://astrogeology.usgs.gov/search/map/mercury_messenger_mdis_global_mosaic_250m> | 公有领域（请引用作者） |
| `mars.jpg` | USGS Astrogeology，Mars Viking Global Color Mosaic 925m（约 1000 幅海盗号轨道器红、紫滤镜图像，Minnaert 光度归一化，扣除雾霾模型后着色），<https://astrogeology.usgs.gov/search/map/mars_viking_global_color_mosaic_925m> | 公有领域 |
| `phobos.jpg` | P. C. Thomas，Phobos 影像拼接图（海盗号轨道器图像，高通滤波，用形状模型控制位置；每度 16 像素，5760×2880），PDS SBN “Small Body Optical Shape Models” V1.0，`data/m1phobosm.fit`，<https://sbnarchive.psi.edu/pds4/non_mission/ast-sat.thomas.shape-models_V1_0/> | PDS 存档数据，公有领域（请引用数据集与作者） |
| `vesta.jpg` | USGS Astrogeology / DLR，Vesta Dawn FC HAMO Global Mosaic 60m（黎明号分幅相机，约 2500 幅清晰滤镜图像），<https://astrogeology.usgs.gov/search/map/vesta_dawn_fc_hamo_global_mosaic_60m> | 公有领域（请引用作者） |
| `deimos.jpg` | P. C. Thomas，Deimos 影像拼接图（海盗号轨道器图像，高通滤波，用形状模型控制位置；每度 4 像素，1440×720），PDS SBN “Small Body Optical Shape Models” V1.0，`data/m2deimosm.fit`，<https://sbnarchive.psi.edu/pds4/non_mission/ast-sat.thomas.shape-models_V1_0/> | PDS 存档数据，公有领域（请引用数据集与作者） |
| `amalthea.jpg` | P. Stooke，Amalthea 晕渲图（根据旅行者 1、2 号图像手绘的喷笔图，位置控制由 Stooke 完成；简单圆柱投影，每度 10 像素），PDS SBN “Stooke Small Bodies Maps” V3.0（MULTI-SA-MULTI-6-STOOKEMAPS-V3.0），<https://sbnarchive.psi.edu/pds4/non_mission/small_bodies.stooke.maps/>，`miscellaneous/j5amalthea/amalcyl.jpg` | 公有领域（“should not be used without proper credit”：Stooke, P., Stooke Small Bodies Maps V3.0, NASA PDS, 2015） |
| `venus.jpg` | USGS Astrogeology，Venus Magellan Global C3-MDIR Synthetic Color Mosaic 4641m（麦哲伦号合成孔径雷达 C3-MIDR 拼接图的合成彩色版，颜色是模拟的地表色调；关闭大气层时显示），<https://astrogeology.usgs.gov/search/map/venus_magellan_global_c3_mdir_synthetic_color_mosaic_4641m> | 公有领域（Use Constraints: None） |
| `titan.jpg` | USGS Astrogeology / Cassini ISS Team，Titan Cassini ISS Global Mosaic 4005m（2015-06，938 nm 甲烷窗口，数据截至 2014-04 的 T100 飞掠，灰度；关闭大气层时显示；有几块均匀灰色的填补区），<https://astrogeology.usgs.gov/search/map/titan_cassini_iss_global_mosaic_4005m> | 公有领域（请引用作者） |
| `ceres.jpg` | USGS Astrogeology / DLR，Ceres Dawn FC Global Mosaic 400m（2015-10，`Ceres_Dawn_FC_DLR_global_20ppd_Oct2015.tif`），<https://astrogeology.usgs.gov/search/map/ceres_dawn_fc_global_mosaic_400m> | 公有领域（请引用作者） |
| `uranus.jpg` / `neptune.jpg` | 哈勃 OPAL（Outer Planet Atmospheres Legacy，PI A. Simon），MAST 高级科学产品，<https://archive.stsci.edu/hlsp/opal>，DOI 10.17909/T9G593。用的是 2015–2025 年每次自转的 F657N、F547M、F467M 全球图（FITS），不是 TIF 预览图。由 `tools/textures/make_ice_giant_textures.py` 生成，做法见下文 | CC BY 4.0（FITS 头的 `LICENSE`），署名 NASA, ESA, A.A. Simon, M.H. Wong |

原始下载地址：

- <https://assets.science.nasa.gov/content/dam/science/psd/photojournal/pia/pia07/pia07782/PIA07782.jpg>。
- <https://planetarymaps.usgs.gov/mosaic/Io_GalileoSSI-Voyager_Global_Mosaic_ClrMerge_1km.tif>（189 MB）。
- <https://planetarymaps.usgs.gov/mosaic/Europa_Voyager_GalileoSSI_global_mosaic_500m.tif>（184 MB）。
- <https://planetarymaps.usgs.gov/mosaic/Ganymede_Voyager_GalileoSSI_global_mosaic_1km.tif>（131 MB）。
- <https://planetarymaps.usgs.gov/mosaic/Callisto_Voyager_GalileoSSI_global_mosaic_1km.tif>（110 MB）。
- 土星卫星（原下载站 planetarymaps.usgs.gov 跳转到 <https://asc-pds-services.s3.us-west-2.amazonaws.com/mosaic/>，下列路径都相对于它）：
  - `Mimas/Cassini_DLR_Mimas.zip`（40 MB）。
  - `Enceladus/Cassini/Enceladus_Cassini_ISS_Global_Mosaic_100m_HPF.tif`（126 MB）。
  - `Tethys_Cassini_mosaic_global_293m.tif`（66 MB）。
  - `Dione_Cassini_Voyager_mosaic_global_154m.tif`（265 MB）。
  - `Rhea_Cassini_Voyager_mosaic_global_417m.tif`（66 MB）。
  - `Iapetus_Cassini_Voyager_mosaic_global_783m.tif`（17 MB）。
  - `Triton_Voyager2_ClrMosaic_GlobalFill_600m.tif`（286 MB）。
  - `Pluto_NewHorizons_Global_Mosaic_300m_Jul2017_8bit.tif`（295 MB）。
  - `Charon_NewHorizons_Global_Mosaic_300m_Jul2017_8bit.tif`（77 MB）。
  - `Mercury_MESSENGER_mosaic_global_250m_2013.tif`（1.9 GB）。
  - `Mars_Viking_ClrMosaic_global_925m.tif`（798 MB）。
  - `Vesta_Dawn_FC_HAMO_Mosaic_Global_74ppd.tif`（357 MB）。
  - `Ceres_Dawn_FC_DLR_global_20ppd_Oct2015.tif`（27 MB）。
- Venus：<https://planetarymaps.usgs.gov/mosaic/Venus_Magellan_C3-MDIR_Colorized_Global_Mosaic_4641m.tif>（96 MB）。
- Titan：<https://planetarymaps.usgs.gov/mosaic/Titan_ISS_P19658_Mosaic_Global_4km.tif>（8 MB）。
- Phobos、Deimos：<https://sbnarchive.psi.edu/pds4/non_mission/ast-sat.thomas.shape-models_V1_0/data/m1phobosm.fit>（17 MB）、`m2deimosm.fit`（1 MB），8 位 FITS。
- 天王星卫星：<https://space.jpl.nasa.gov/tmaps/pix/ura1vuu2.tif> .. `ura5vuu2.tif`（各 0.2 至 0.4 MB；1 Ariel、2 Umbriel、3 Titania、4 Oberon、5 Miranda）。

## Eris、Haumea 和 Makemake（NASA VTAD）

`eris.jpg`、`haumea.jpg`、`makemake.jpg` 来自 NASA Visualization Technology Applications and Development（VTAD）的三个 3D 模型。入口是 [Pluto and Dwarf Planet Resources](https://science.nasa.gov/solar-system/resources/resource-packages/pluto-dwarf-planet-resources/#3dmodels)，各模型页均署名 NASA VTAD，未标注第三方版权：

| 纹理 | 模型页 | GLB 下载 | 内嵌颜色图 |
|---|---|---|---|
| `eris.jpg` | [Eris 3D Model](https://science.nasa.gov/resource/eris-3d-model/) | [Eris_1_2326.glb](https://assets.science.nasa.gov/content/dam/science/psd/solar/2023/09/e/Eris_1_2326.glb) | `Eris_diff.jpg`，实际为 PNG，2048×1024 |
| `haumea.jpg` | [Haumea 3D Model](https://science.nasa.gov/resource/haumea-3d-model/) | [Haumea_1_1000.glb](https://assets.science.nasa.gov/content/dam/science/psd/solar/2023/09/h/Haumea_1_1000.glb) | `Haumea_diff.jpg`，实际为 PNG，1024×1024 |
| `makemake.jpg` | [Makemake 3D Model](https://science.nasa.gov/resource/makemake-3d-model/) | [Makemake_1_1430.glb](https://assets.science.nasa.gov/content/dam/science/psd/solar/2023/09/m/Makemake_1_1430.glb) | `makemake_diff.jpg`，实际为 PNG，2048×1536，十字形图集 |

使用依据是 [NASA Images and Media Usage Guidelines](https://www.nasa.gov/nasa-brand-center/images-and-media/)：NASA 的三维模型纹理和多边形数据通常在美国不受版权保护，可用于计算机图形模拟等用途，需要注明 NASA 来源；第三方版权内容另行处理。这里署名 NASA/VTAD，遵守该使用规范。核查日期：2026-10-20；本地 `ref/dwarf/` 的三个文件与上述官方下载逐字节一致。SHA-256：

```text
1b97f1b9eea244b047a72c2a3ac0db315715b7d421033d860d669ffd78681d54  Eris_1_2326.glb
7b31737593054bd4ef8f61138a6ad1a2577bec03b84d1c3e76e8ac78632a87f8  Haumea_1_1000.glb
e7d5fa27f50afae10f0e7290cdf7a14ca856c0b9e99b2db02f25b75441fbdbef  Makemake_1_1430.glb
```

复现命令（需要 numpy、Pillow）：

```powershell
python tools/textures/extract_dwarf_textures.py ref/dwarf assets/textures
```

工具沿各经纬度方向与源模型三角形求交，插值原 UV，再采样材质的 `baseColorTexture`。Eris 和 Makemake 的 Y 轴转为项目的北极 Z 轴；Haumea 的短轴本身是 Z 轴。输出从左向右为东经 0°–360°、从上向下为北纬 90°–南纬 90°；经度零点任意。Haumea 按源网格的位置算行星中心纬度，适配项目的三轴椭球贴图方式。在线性 RGB 中双线性采样，以 4096×2048 烘焙后做 2×2 平均，转回 sRGB，存为 2048×1024 JPEG（质量 90，4:4:4）。每张图都检查全像素有网格覆盖。

模型页未提供表面细节的观测地图来源，本项目按示意表面使用；纹理中的沟槽、坑和斑块位置不代表实测地貌。颜色取模型的颜色图，材质中的金属度、粗糙度及 Haumea 的法线图未接入。源模型网格用于离线转换，运行时的半径、自转和 Haumea 三轴尺寸由场景配置决定（源模型短半轴为 569，场景采用 [HAUM] 的 513 km）。

## 天王星和海王星（OPAL）

OPAL 每年观测两次左右完整自转，每次自转每个滤镜一张简单圆柱图（721×361，每度 2 像素，行星图形纬度），已做过 Minnaert 改正（扣除临边昏暗）。每张图只覆盖朝向地球的半球，边缘（出射角很大处）被拉成彩色条纹，单张图做不了纹理。所以：

- FITS 值乘各年 README 里的 I/F 系数得到线性 I/F，F657N、F547M、F467M 当作红、绿、蓝。I/F 是相对阳光的反射率，直接当反照率用，不另调颜色
- 去掉离图像边缘纬度 22° 以内（或离经度空缺 22° 以内）的像素，每次自转取纬向中值，再对所有自转取中值，得到纬度剖面
- 没有观测的纬度（天王星 35.5°S 以南、海王星 48.5°N 以北）取对称纬度的剖面，在边界处乘一个比例接上，比例在 20° 内渐变为 1。天王星南半球因此有和北极一样的亮极冠（旅行者 2 号 1986 年看到的南极区也偏亮），海王星北极同理
- 叠加一张真实图相对它自己纬向中值的起伏，离覆盖边缘 8° 内渐隐：天王星用 2025a，乘 0.75；海王星用 2019a，不缩放（原始对比度），图里有北半球暗斑 NDS-2018（约 20°N，F467M 里暗 15% 左右，真彩色下不明显）和旁边的亮伴云。2018b 也能看到这个暗斑，但它旁边有一条拼接缝，缝的条纹伸进来 30°
- 从行星图形纬度重采样到行星中心纬度（场景按行星中心纬度取纹理），左边缘 0°E、东经向右递增：天王星 2021 年以后的图是“右边缘 0°E，向左递增”，要左右翻转；海王星的图是“左边缘 360°W，向右递减”，不用翻

天王星 2022b 的 F657N 比 2022a 和其他年份亮 1.4 倍（README 里 I/F 系数相同），合成后偏粉，工具里排除了。

下载：每个周期的页面（如 <https://archive.stsci.edu/hlsp/opal/opal-uranus-cycle-33>）列出文件，FITS 和 README 都在 `https://archive.stsci.edu/missions/hlsp/opal/cycle<N>/<planet>/` 下（周期 22–33，约 140 MB）。把各年的 `*_f467m_v1_globalmap.fits`、`*_f547m_*`、`*_f657n_*` 和 `*_all_v*_readme.txt` 放进同一个目录，再运行 `python make_ice_giant_textures.py <目录> <输出目录>`。

## Kerbol 系统（程序生成）

Kerbol 场景 16 个天体的纹理（`kerbin.jpg`、`mun.jpg`、`minmus.jpg`、`moho.jpg`、`eve.jpg`、`gilly.jpg`、`duna.jpg`、`ike.jpg`、`dres.jpg`、`jool.jpg`、`laythe.jpg`、`vall.jpg`、`tylo.jpg`、`bop.jpg`、`pol.jpg`、`eeloo.jpg`）是虚构表面，由 `tools/textures/make_kerbol_textures.py` 生成，不含《坎巴拉太空计划》的任何游戏文件，也不模仿游戏里的具体地形；属于本项目的产物（GPL-3.0-or-later）。每个天体只保留辨识度最高的一两个特征：Kerbin 的海洋（约占 64% 面积）、绿色陆地、沙漠、雪山和冰盖，Mun 的灰色陨石坑，Minmus 的薄荷绿和边缘清晰的浅色平原，Moho 的深棕色坑，Eve 的紫色陆地和暗色海，Duna 的锈红色平原、暗色高地和极冠，Ike 的深灰坑，Dres 的灰色坑和几条暗色峡谷，Jool 的绿色条带和几个风暴卵形斑，Laythe 的海洋（约占 86%）和沙色群岛，Vall 的浅蓝冰面裂纹，Tylo 的大片明暗斑块，Eeloo 的白冰和棕色裂纹，Gilly、Bop、Pol 的棕灰、深棕和黄褐色。

- 裂纹和峡谷是扭曲噪声的零等值线，按 |n| / |∇n| 估计到等值线的距离，线宽在各处一致

- 每个像素对应单位球面上的一点，噪声是在球面上采样的三维梯度噪声，所以 180° 经线处没有接缝，两极也不会被拉扯；倍频程只加到波长约 2 像素为止
- 噪声用自带的整数哈希，不用 numpy 的随机数生成器：种子固定在工具里，换 numpy 版本重新生成结果不变
- 场景的天体着色没有法线贴图，陨石坑只体现在反照率上（坑底略暗、坑缘和溅射物偏亮，年轻的大坑带辐射纹）
- 颜色在线性空间混合，存为 sRGB JPEG（质量 90），2048×1024（Gilly、Bop、Pol 1024×512）；亮度参照地球（海洋很暗）和月球（平均约 150）的纹理

重新生成：`python -E -P tools/textures/make_kerbol_textures.py assets/textures [名称 ...]`（需要 numpy 和 Pillow，全部约 10 分钟）。

## 经度约定

这些 USGS 图像按东经向右递增排列，与元数据中的“positive west”标注相反。左边缘经度如下：

- Io（经度域 −180..180）：左边缘为 180°W，即 −180°E。
- Europa、Ganymede、Callisto（经度域 0..360）：左边缘为 360°W，即 0°E。

经度方向已通过以下地貌核对：

- Callisto：Valhalla 盆地（14°N 56°W）。
- Io：Pele（18.7°S 255°W）、Loki（13°N 309°W）。
- Ganymede：Galileo Regio。

- 土星卫星：全部是东经向右递增（Rhea、Iapetus 的 PositiveWest 标签与图像排列相反）。左边缘：Enceladus 为 0°E，其余（Mimas、Tethys、Dione、Rhea、Iapetus）为 −180°E。用以下地貌验证过：
  - Mimas：Herschel 坑（1.7°N 111.8°W）。
  - Tethys：Odysseus 坑（32.8°N 128.9°W）。
  - Rhea：Tirawa 盆地（34.2°N 151.7°W）。
  - Enceladus：Ali Baba（55.1°N 22.3°W）与 Aladdin（60.7°N 26.7°W）这对坑，Aladdin 在西北。
  - Iapetus：暗色的 Cassini Regio 以前导半球（90°W）为中心。
  - Dione：亮条纹地形在后随半球（270°W）。
  - 这些灰度图以 JPEG 格式保存。
- 地球（Blue Marble）和月球（CGI Moon Kit）都以 0° 为中心、东经递增，左边缘为 −180°E。月球的方向用危海（17°N 59°E）和第谷环形山（43°S 11°W）验证过。

- 天王星卫星（JPL Solar System Simulator 的图没有元数据）：东经向右递增，左边缘为 −180°E。用于核对方向的地貌为：Oberon 的暗底 Hamlet 坑（46°S 44°E）、Titania 的 Ursula 坑（12°S 45°E）、Miranda 的 Arden 冕（29°S 74°E）与 Elsinore 冕（25°S 257°E）（坐标来自 USGS Gazetteer of Planetary Nomenclature，行星中心坐标、东经为正）。若按西经递增解释图像，前两个地貌会落在 315°E。
- Triton、Charon：标签为 PositiveEast、−180..180，左边缘 −180°E；Pluto：PositiveEast、0..360，左边缘 0°E，心形的 Sputnik Planitia（约 175°E）位于图中央。
- Triton 使用旅行者号滤镜合成图像的偏绿色调。

- 类地行星与小天体（以下都是标签为 PositiveEast、行星中心纬度的图，经过地貌核对）：
  - Mercury、Mars、Vesta：经度域 −180..180，左边缘 −180°E。Mercury 用 Caloris 盆地（31.5°N 162.7°E）核对，图里约 14°W 处有一条竖带是 USGS 用高入射角图像补的缺口（阴影明显）。Mars 用 Olympus Mons（18°N 226°E）、Hellas（42°S 70°E）、Syrtis Major 核对。
  - Vesta 的经度是 Claudia″（Claudia Double Prime）系统，Claudia 坑在 146°E（USGS 页面说明）；场景用的 PCK W0 = 285.39° 正是这个系统（Dawn 重力数据集的坐标系文档 `VESTA_COORDINATES_131018` 列出了四种 Vesta 坐标系的 W0）。“雪人”三坑（Marcia、Calpurnia、Minucia）在约 340°E。
  - Ceres：经度域 0..360，左边缘 0°E。Occator 坑（19.8°N 239.3°E，IAU 系统，即 Kait 坑在 0°）落在图中对应位置，与 PCK 的 W0 = 170.65° 一致。
- Phobos、Deimos（Thomas 拼接图，FITS）以 0°N 0°E 为中心，东经向右增加，左边缘 −180°E。这些约定用 Phobos 拼接图核对：Stickney 在中心偏左约 48°（49°W），与形状模型里凹陷的中心（2°S 50°W）重合；Limtoc 在 10°S 54°W（Gazetteer：11°S 54°W）。标签写的是“Line Bottom to Top”，但只有按文件顺序把第一行放在最上（北），Limtoc 坑才落在 Stickney 中心以南（实际为 11°S 与 1°S），凹槽系统的位置也才与 USGS 的 Phobos 图一致。Deimos 用同一批工具制作，采用同样的约定；它的地名（Swift、Voltaire）是 1973 年的坐标，控制网未知，无法用于核对。Phobos 有少量无数据区（一个小方块和零星点），Deimos 约 4% 的无数据区（纯黑，集中在一个经度段和南半球高纬）用 `fill_unimaged` 填成平淡表面。两张图都有来自低分辨率图像的区域，放大后模糊、有马赛克块（Phobos 的 Stickney 内部、Deimos 约一半经度）；高通滤波使对比度偏低。Deimos 原图宽 1440 像素，生成时放大到 2048。
- Amalthea（Stooke 晕渲图）：是带有固定光照阴影的手绘喷笔图（地图索引里的类型 “S”），所以坑的明暗不随太阳方向变化。地图说明写“0 longitude at the center”，没有写经度方向。地名坐标（Pan 55°N 35°W、Gaea 80°S 90°W）在这张模糊的图上无法辨认，所以用 Stooke 自己的形状模型（同一作者）核对：根据形状模型计算不同光照方向的晕渲图，再与原图计算相关系数，“东经向右”的最高相关系数为 0.23（经度偏移 −5°，光从东边来），“西经向右”的最高相关系数为 0.16。相关系数较低，经度对齐的可信度为中等。生成时乘上 Amalthea 的场景颜色 `#9c5a43`（`tint_gray`：灰度换算到线性值后除以全图平均，再乘颜色的线性值），否则这颗很暗、很红的卫星会显示成灰白色。
- Mercury 使用 2013 年的单色底图。备选的彩色拼接图（Global Color Mosaic 665m）将 1000/750/430 nm 映射到 RGB，属于增强假彩色，整体偏蓝，两极还有缺块和杂乱纹理，因此未选用。
- Phobos 使用 Thomas 拼接图，因为它与所用的形状模型共用控制网（见 `assets/shapes/SOURCES.md`）。备选纹理：USGS 的 Stooke 海盗号拼接图（`Phobos_Viking_Mosaic_40ppd_DLRcontrol.tif`）画面更清晰，但基于 DLR 的另一套控制网，地貌偏西约 5° 且各处不一；Mars Express SRC 拼接图（`Phobos_ME_SRC_Mosaic_Global_16ppd.tif`）含有原始图像的阴影和缺块，两侧边缘也无法衔接。
- Vesta 的 HAMO 拼接图中，北纬 60° 以上处于极夜，纹理保留了这些区域的低亮度。Ceres 南纬约 82° 以南缺少数据，按同纬度均值填补。
- Mercury、Phobos、Deimos、Vesta、Ceres 的灰度图存为 JPEG，与其他卫星纹理一样用于反照率。这些图像的亮度经过拉伸，比实际反照率高得多，尤其是表面很暗的 Phobos 和 Ceres。

- Venus：标签为 PositiveEast、−180..180，左边缘 −180°E。用 Maxwell Montes（65°N 3°E，雷达很亮）和它西侧的 Lakshmi Planum 核对过。
- Titan：标签为 PositiveWest、0..360、中心 180°W，但图像按东经向右递增排列（与其他 USGS 图相同），左边缘 0°E。用 Menrva 环形山（19°N 87°W，即 273°E，在图的 0.76 处）、亮的 Xanadu（10°S 100°W）和暗的 Shangri-La（10°S 165°W）核对过；按西经向右解释时 Menrva 会落在 87°E。

场景文件用 `texture_left_lon_deg` 指定图像左边缘的东经。

纹理坐标按行星中心经纬度计算，即从天体中心指向表面的方向。这能让三轴椭球（Miranda、Ariel）上的地貌落在对应经纬度；单位球网格的参数经纬度与行星中心经纬度在 Vesta 这样扁的三轴椭球上最多相差约 7°。有形状模型的天体（Phobos、Deimos、Vesta）按网格顶点的经纬度贴图；Vesta 的形状来自同一 DLR 发布的 DTM。

极区无数据像素的填补适用于 USGS 拼接图。完整地图跳过这一步，因为地球极地海洋等区域本来就暗，按暗像素识别缺失数据会误填这些区域。

天王星卫星、Triton、Pluto、Charon 有大片未拍摄区域。旅行者 2 号在 1986 年拍到的是天王星卫星的南半球；Triton 的北半球、Pluto 和 Charon 30°S 以南当时处于黑夜。这些区域由 `fill_unimaged` 填补：

1. 将纯黑像素标为缺失。天王星卫星图像在拍摄边缘有暗色阴影毛刺，且缺少很黑的地形，因此使用阈值 40。
2. 用闭运算去掉图内零星黑点，让 Pluto 很暗的 Cthulhu 区仍被识别为已知地形。
3. 将已知区域的边界向内腐蚀 15 像素，去掉毛刺。
4. 在缺失区域靠近边界的位置，按附近已知像素的平均值填补。计算使用归一化卷积，高斯半径为图宽的 2%；远处渐变到整幅图的平均亮度。

填补区域的亮度由上述方法估算，显示为缺少地貌细节的表面。当地的实际地貌仍未知。2020 年代，太阳照亮天王星卫星的北半球，画面中朝阳的区域因此多为填补结果。

木星的大气特征会在经度上漂移，所以 `jupiter.jpg` 的经度对齐没有物理意义。

## 天空背景：`milky_way.jpg`

| 文件 | 来源 | 许可 |
|---|---|---|
| `milky_way.jpg` | NASA/Goddard Space Flight Center Scientific Visualization Studio，Deep Star Maps 2020，Milky Way background 图层 `milkyway_2020_4k.exr`（4096×2048，34.7 MB），<https://svs.gsfc.nasa.gov/4851>；数据来自 Gaia DR2（ESA/Gaia/DPAC） | SVS 内容为公有领域（“unless otherwise noted”，见 <https://svs.gsfc.nasa.gov/help/>）。按页面要求署名：“NASA/Goddard Space Flight Center Scientific Visualization Studio. Gaia DR2: ESA/Gaia/DPAC.” |

原始下载地址：<https://svs.gsfc.nasa.gov/vis/a000000/a004800/a004851/milkyway_2020_4k.exr>

- 内容：Gaia DR2 的暗星和弥漫光。该图层排除了 Hipparcos/Tycho 亮星，因此与程序里的星表星点（HYG，到 6.5 等）没有重复。背景缺少 6.5 至 11.5 等的 Tycho 星，对显示影响较小。
- 投影：ICRF（J2000）赤经赤纬的简单圆柱投影。赤经 0h 在图像中央、向左递增（从球内看），上边缘是 +90° 赤纬，即 u = 0.5 − α/360°，v = (90° − δ)/180°。用银心（α = 266.4°，图中最亮的核球）和大麦哲伦云（α ≈ 80°，δ ≈ −69°）的位置验证过。
- 处理：`tools/sky/convert_milkyway.py` 依赖标准库和 Pillow，读取 ZIP 压缩的 half 浮点 EXR。原图存储线性值，最大值为 1.0（1% 分位 7.6e-4，99.99% 分位 0.72）。工具按 sRGB 编码存为 JPEG（质量 92，不做色度抽样），场景中的 `milky_way_brightness` 控制显示亮度。
