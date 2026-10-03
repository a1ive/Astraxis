# 纹理来源

所有纹理都是简单圆柱投影（等经纬度）的全球地图，由 `tools/textures/prepare_textures.py` 缩放到 2048×1024：先按整数倍降采样，再用 Lanczos 重采样，USGS 拼接图还会填补极区无数据像素，最后存为 JPEG（质量 90）。原始文件太大，不放进仓库，需要时按下表重新下载。

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
| `phobos.jpg` | USGS Astrogeology，Phobos Viking Global Mosaic 5m（Phil Stooke，以海盗号高分辨率图像为主，补充其他飞船的图像，DLR 控制网；`Phobos_Viking_Mosaic_40ppd_DLRcontrol.tif`），<https://astrogeology.usgs.gov/search/map/phobos_viking_global_mosaic_5m> | 公有领域（请引用作者） |
| `vesta.jpg` | USGS Astrogeology / DLR，Vesta Dawn FC HAMO Global Mosaic 60m（黎明号分幅相机，约 2500 幅清晰滤镜图像），<https://astrogeology.usgs.gov/search/map/vesta_dawn_fc_hamo_global_mosaic_60m> | 公有领域（请引用作者） |
| `ceres.jpg` | USGS Astrogeology / DLR，Ceres Dawn FC Global Mosaic 400m（2015-10，`Ceres_Dawn_FC_DLR_global_20ppd_Oct2015.tif`），<https://astrogeology.usgs.gov/search/map/ceres_dawn_fc_global_mosaic_400m> | 公有领域（请引用作者） |

原始下载地址：

- <https://assets.science.nasa.gov/content/dam/science/psd/photojournal/pia/pia07/pia07782/PIA07782.jpg>
- <https://planetarymaps.usgs.gov/mosaic/Io_GalileoSSI-Voyager_Global_Mosaic_ClrMerge_1km.tif>（189 MB）
- <https://planetarymaps.usgs.gov/mosaic/Europa_Voyager_GalileoSSI_global_mosaic_500m.tif>（184 MB）
- <https://planetarymaps.usgs.gov/mosaic/Ganymede_Voyager_GalileoSSI_global_mosaic_1km.tif>（131 MB）
- <https://planetarymaps.usgs.gov/mosaic/Callisto_Voyager_GalileoSSI_global_mosaic_1km.tif>（110 MB）
- 土星卫星（planetarymaps.usgs.gov 现在跳转到 <https://asc-pds-services.s3.us-west-2.amazonaws.com/mosaic/>，下列路径都相对于它）：
  - `Mimas/Cassini_DLR_Mimas.zip`（40 MB）
  - `Enceladus/Cassini/Enceladus_Cassini_ISS_Global_Mosaic_100m_HPF.tif`（126 MB）
  - `Tethys_Cassini_mosaic_global_293m.tif`（66 MB）
  - `Dione_Cassini_Voyager_mosaic_global_154m.tif`（265 MB）
  - `Rhea_Cassini_Voyager_mosaic_global_417m.tif`（66 MB）
  - `Iapetus_Cassini_Voyager_mosaic_global_783m.tif`（17 MB）
  - `Triton_Voyager2_ClrMosaic_GlobalFill_600m.tif`（286 MB）
  - `Pluto_NewHorizons_Global_Mosaic_300m_Jul2017_8bit.tif`（295 MB）
  - `Charon_NewHorizons_Global_Mosaic_300m_Jul2017_8bit.tif`（77 MB）
  - `Mercury_MESSENGER_mosaic_global_250m_2013.tif`（1.9 GB）
  - `Mars_Viking_ClrMosaic_global_925m.tif`（798 MB）
  - `Phobos_Viking_Mosaic_40ppd_DLRcontrol.tif`（104 MB）
  - `Vesta_Dawn_FC_HAMO_Mosaic_Global_74ppd.tif`（357 MB）
  - `Ceres_Dawn_FC_DLR_global_20ppd_Oct2015.tif`（27 MB）
- 天王星卫星：<https://space.jpl.nasa.gov/tmaps/pix/ura1vuu2.tif> .. `ura5vuu2.tif`（各 0.2–0.4 MB；1 Ariel、2 Umbriel、3 Titania、4 Oberon、5 Miranda）

## 经度约定

USGS 元数据标注的是“positive west”，但图像的实际排布是东经向右递增。左边缘位置如下：

- Io（经度域 −180..180）：左边缘为 180°W，即 −180°E
- Europa、Ganymede、Callisto（经度域 0..360）：左边缘为 360°W，即 0°E

这一点用已知地貌验证过：
- Callisto：Valhalla 盆地（14°N 56°W）
- Io：Pele（18.7°S 255°W）、Loki（13°N 309°W）
- Ganymede：Galileo Regio

- 土星卫星：全部是东经向右递增（Rhea、Iapetus 的标签写的是 PositiveWest，同样不可信）。左边缘：Enceladus 为 0°E，其余（Mimas、Tethys、Dione、Rhea、Iapetus）为 −180°E。用以下地貌验证过：
  - Mimas：Herschel 坑（1.7°N 111.8°W）
  - Tethys：Odysseus 坑（32.8°N 128.9°W）
  - Rhea：Tirawa 盆地（34.2°N 151.7°W）
  - Enceladus：Ali Baba（55.1°N 22.3°W）与 Aladdin（60.7°N 26.7°W）这对坑，Aladdin 在西北
  - Iapetus：暗色的 Cassini Regio 以前导半球（90°W）为中心
  - Dione：亮条纹地形在后随半球（270°W）
  - 这些都是灰度图，按原样存成 JPEG
- 地球（Blue Marble）和月球（CGI Moon Kit）都以 0° 为中心、东经递增，左边缘为 −180°E。月球的方向用危海（17°N 59°E）和第谷环形山（43°S 11°W）验证过。

- 天王星卫星（JPL Solar System Simulator 的图没有元数据）：东经向右递增，左边缘为 −180°E。三颗卫星的地貌都对得上：Oberon 的暗底 Hamlet 坑（46°S 44°E）、Titania 的 Ursula 坑（12°S 45°E）、Miranda 的 Arden 冕（29°S 74°E）与 Elsinore 冕（25°S 257°E）（坐标来自 USGS Gazetteer of Planetary Nomenclature，行星中心坐标、东经为正）。如果是西经递增，前两个会落在 315°E
- Triton、Charon：标签为 PositiveEast、−180..180，左边缘 −180°E；Pluto：PositiveEast、0..360，左边缘 0°E，心形的 Sputnik Planitia（约 175°E）正好在图中央
- Triton 的颜色来自旅行者号的滤镜合成，偏绿，按原样保留

- 类地行星与小天体（以下都是标签为 PositiveEast、行星中心纬度的图，经过地貌核对）：
  - Mercury、Mars、Phobos、Vesta：经度域 −180..180，左边缘 −180°E。Mercury 用 Caloris 盆地（31.5°N 162.7°E）核对，图里约 14°W 处有一条竖带是 USGS 用高入射角图像补的缺口（阴影明显）；Mars 用 Olympus Mons（18°N 226°E）、Hellas（42°S 70°E）、Syrtis Major 核对；Phobos 的 Stickney 坑在约 55°W（Gazetteer 为 49°W 附近；不规则天体投影到球面，有几度出入）
  - Vesta 的经度是 Claudia″（Claudia Double Prime）系统，Claudia 坑在 146°E（USGS 页面说明）；场景用的 PCK W0 = 285.39° 正是这个系统（Dawn 重力数据集的坐标系文档 `VESTA_COORDINATES_131018` 列出了四种 Vesta 坐标系的 W0）。“雪人”三坑（Marcia、Calpurnia、Minucia）在约 340°E
  - Ceres：经度域 0..360，左边缘 0°E。Occator 坑（19.8°N 239.3°E，IAU 系统，即 Kait 坑在 0°）正好落在图中这个位置，与 PCK 的 W0 = 170.65° 一致
- Mercury 原本想用彩色拼接图（Global Color Mosaic 665m），但它是 1000/750/430 nm 映射到 RGB 的增强假彩色（整体偏蓝），两极还有缺块和乱纹，所以改用 2013 年的单色底图
- Phobos 最初试过 Mars Express SRC 拼接图（`Phobos_ME_SRC_Mosaic_Global_16ppd.tif`），但其中混有原始图像的阴影和缺块，两侧边缘也接不上，所以改用 Stooke 的海盗号拼接图
- Vesta 的 HAMO 拼接图里北纬 60° 以上处于极夜，很暗，按原样保留；Ceres 南纬约 82° 以南没有数据，用同纬度均值填补
- 灰度图（Mercury、Phobos、Vesta、Ceres）按原样存成 JPEG，和其他卫星一样直接作为反照率：它们都经过拉伸，比真实的反照率亮得多（Phobos、Ceres 实际很暗）

场景文件统一写成 `texture_left_lon_deg`（东经）。

纹理坐标按**行星中心**经纬度计算（从天体中心看过去的方向），而不是单位球网格的参数经纬度，所以三轴椭球（Phobos、Vesta、Miranda、Ariel）上的地貌位置是对的；两者在 Vesta 上最多差约 7°。

只有 USGS 拼接图会填补极区无数据像素。完整的地图不填补，否则地球极地海洋这类本来就暗的像素会被误当成缺失数据。

天王星卫星、Triton、Pluto、Charon 有大片从未拍到的区域（旅行者 2 号 1986 年只看到天王星卫星的南半球；Triton 的北半球、Pluto 和 Charon 30°S 以南当时都在黑夜里），用另一种填补（`fill_unimaged`）：纯黑像素算作缺失（天王星卫星的图在拍摄边界有暗色的阴影毛刺，又没有真正很黑的地形，阈值取 40）；先闭运算去掉图内的零星黑点（避免把 Pluto 很暗的 Cthulhu 区当成缺失），再把边界向内腐蚀 15 像素去掉毛刺；缺失处在边界附近用附近已知像素的平均（归一化卷积，高斯半径为图宽的 2%），远处渐变到整幅图的平均亮度。所以这些区域显示为平淡的表面，不是真实的地貌。2020 年代太阳照着天王星卫星的北半球，所以现在看到的正好是填补的部分。

木星的大气特征会在经度上漂移，所以 `jupiter.jpg` 的经度对齐没有物理意义。

## 天空背景：`milky_way.jpg`

| 文件 | 来源 | 许可 |
|---|---|---|
| `milky_way.jpg` | NASA/Goddard Space Flight Center Scientific Visualization Studio，Deep Star Maps 2020，Milky Way background 图层 `milkyway_2020_4k.exr`（4096×2048，34.7 MB），<https://svs.gsfc.nasa.gov/4851>；数据来自 Gaia DR2（ESA/Gaia/DPAC） | SVS 内容为公有领域（“unless otherwise noted”，见 <https://svs.gsfc.nasa.gov/help/>）。按页面要求署名：“NASA/Goddard Space Flight Center Scientific Visualization Studio. Gaia DR2: ESA/Gaia/DPAC.” |

原始下载地址：<https://svs.gsfc.nasa.gov/vis/a000000/a004800/a004851/milkyway_2020_4k.exr>

- 内容：Gaia DR2 的暗星和弥漫光，**不含** Hipparcos/Tycho 亮星，所以和程序里的 BSC5 星点（到 6.5 等）不会重复。缺了 6.5–11.5 等的 Tycho 星，影响不大。
- 投影：ICRF（J2000）赤经赤纬的简单圆柱投影。赤经 0h 在图像中央、向左递增（从球内看），上边缘是 +90° 赤纬，即 u = 0.5 − α/360°，v = (90° − δ)/180°。用银心（α = 266.4°，图中最亮的核球）和大麦哲伦云（α ≈ 80°，δ ≈ −69°）的位置验证过。
- 处理：`tools/sky/convert_milkyway.py`（只依赖标准库和 Pillow；EXR 为 ZIP 压缩的 half 浮点）。原图是线性值，最大值恰好是 1.0（1% 分位 7.6e-4，99.99% 分位 0.72），所以直接按 sRGB 编码存为 JPEG（质量 92，不做色度抽样）。亮度由场景里的 `milky_way_brightness` 决定。
