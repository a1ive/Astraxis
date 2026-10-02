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

原始下载地址：

- <https://assets.science.nasa.gov/content/dam/science/psd/photojournal/pia/pia07/pia07782/PIA07782.jpg>
- <https://planetarymaps.usgs.gov/mosaic/Io_GalileoSSI-Voyager_Global_Mosaic_ClrMerge_1km.tif>（189 MB）
- <https://planetarymaps.usgs.gov/mosaic/Europa_Voyager_GalileoSSI_global_mosaic_500m.tif>（184 MB）
- <https://planetarymaps.usgs.gov/mosaic/Ganymede_Voyager_GalileoSSI_global_mosaic_1km.tif>（131 MB）
- <https://planetarymaps.usgs.gov/mosaic/Callisto_Voyager_GalileoSSI_global_mosaic_1km.tif>（110 MB）

## 经度约定

USGS 元数据标注的是“positive west”，但图像的实际排布是东经向右递增。左边缘位置如下：

- Io（经度域 −180..180）：左边缘为 180°W，即 −180°E
- Europa、Ganymede、Callisto（经度域 0..360）：左边缘为 360°W，即 0°E

这一点用已知地貌验证过：
- Callisto：Valhalla 盆地（14°N 56°W）
- Io：Pele（18.7°S 255°W）、Loki（13°N 309°W）
- Ganymede：Galileo Regio

- 地球（Blue Marble）和月球（CGI Moon Kit）都以 0° 为中心、东经递增，左边缘为 −180°E。月球的方向用危海（17°N 59°E）和第谷环形山（43°S 11°W）验证过。

场景文件统一写成 `texture_left_lon_deg`（东经）。

只有 USGS 拼接图会填补极区无数据像素。完整的地图不填补，否则地球极地海洋这类本来就暗的像素会被误当成缺失数据。

木星的大气特征会在经度上漂移，所以 `jupiter.jpg` 的经度对齐没有物理意义。

## 天空背景：`milky_way.jpg`

| 文件 | 来源 | 许可 |
|---|---|---|
| `milky_way.jpg` | NASA/Goddard Space Flight Center Scientific Visualization Studio，Deep Star Maps 2020，Milky Way background 图层 `milkyway_2020_4k.exr`（4096×2048，34.7 MB），<https://svs.gsfc.nasa.gov/4851>；数据来自 Gaia DR2（ESA/Gaia/DPAC） | SVS 内容为公有领域（“unless otherwise noted”，见 <https://svs.gsfc.nasa.gov/help/>）。按页面要求署名：“NASA/Goddard Space Flight Center Scientific Visualization Studio. Gaia DR2: ESA/Gaia/DPAC.” |

原始下载地址：<https://svs.gsfc.nasa.gov/vis/a000000/a004800/a004851/milkyway_2020_4k.exr>

- 内容：Gaia DR2 的暗星和弥漫光，**不含** Hipparcos/Tycho 亮星，所以和程序里的 BSC5 星点（到 6.5 等）不会重复。缺了 6.5–11.5 等的 Tycho 星，影响不大。
- 投影：ICRF（J2000）赤经赤纬的简单圆柱投影。赤经 0h 在图像中央、向左递增（从球内看），上边缘是 +90° 赤纬，即 u = 0.5 − α/360°，v = (90° − δ)/180°。用银心（α = 266.4°，图中最亮的核球）和大麦哲伦云（α ≈ 80°，δ ≈ −69°）的位置验证过。
- 处理：`tools/sky/convert_milkyway.py`（只依赖标准库和 Pillow；EXR 为 ZIP 压缩的 half 浮点）。原图是线性值，最大值恰好是 1.0（1% 分位 7.6e-4，99.99% 分位 0.72），所以直接按 sRGB 编码存为 JPEG（质量 92，不做色度抽样）。亮度由场景里的 `milky_way_brightness` 决定。
