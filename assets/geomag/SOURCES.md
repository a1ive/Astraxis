# 地磁与极光数据来源

`Planetary aurora` 的 `Event setting` 由活动事件的 `aurora_visible`、`aurora_ultraviolet` 控制。普通事件和无活动事件时关闭所有极光。`Visible light` 手动开启地球、Europa 和土星的可见光极光；`UV aurora overlay` 开启紫外叠加，同时保留可见光发射（土星改用紫外卵）。地球极光还需要开启 `Atmospheres`。

## Ganymede 自身的紫外极光

木星场景末尾的两个 Ganymede 导览分别朝向前导与拖尾半球，也可在面板的 `Planetary aurora` 中选择 `UV aurora overlay`。太阳系场景通过木卫三模板使用同一份数据。

`ganymede_aurora.txt` 来自 Marzok、Schlegel、Saur、Roth、Grodent、Strobel 和 Retherford 的 [Zenodo 6564687 数据集](https://doi.org/10.5281/zenodo.6564687)，原文件为 `Fig5_BrightnessMap.txt`，许可 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)。引用：Marzok et al. 2022, *Mapping the Brightness of Ganymede's Ultraviolet Aurora Using Hubble Space Telescope Observations*, JGR Planets 127, e2022JE007256，[论文](https://doi.org/10.1029/2022JE007256)。再分发依据是 Zenodo 数据许可；论文正文另采用 CC BY-NC，未随程序分发。

该图合并了 1998 至 2017 年的 46 次 HST/STIS 观测，记录 OI 135.6 nm 的垂直柱亮度（Rayleigh）。数组为 80×40，像元中心西经 2.25° 至 357.75°、纬度 −87.75° 至 87.75°，两轴步长均为 4.5°。烘焙保留网格和正亮度，将负背景残差截为零；原图中零值区域仍为零。运行时把东经取负后采样，纬度从南向北，经度环绕、纬度夹紧。上传线性 UNORM 纹理的范围为 0 至 200 R，最大量化误差小于 0.4 R，远景使用 mipmap。

发光分布采用论文 §2.3 校正视线柱长时使用的均匀 70 km 壳层，垂直积分还原图中的柱亮度。地图固定于木卫三本体；蓝色及每 100 R 的显示增益是示意。数据的低信噪比像元也参与双线性插值，部分暗区来自缺测，不能解释为没有极光。当前没有按木星电流片位置切换分组图，也没有模拟摇摆、风暴或细小斑块。

重新生成（Python 标准库）：

```powershell
New-Item -ItemType Directory -Force build/ganymede_aurora | Out-Null
curl.exe -L --fail https://zenodo.org/api/records/6564687/files/Fig5_BrightnessMap.txt/content -o build/ganymede_aurora/Fig5_BrightnessMap.txt
python tools/geomag/ganymede_aurora_bake.py build/ganymede_aurora/Fig5_BrightnessMap.txt assets/geomag/ganymede_aurora.txt
```

原文件 SHA-256：`f0dd34cdb5fd0f43353631eb26c488866c78e4081cf829879bd419028646a5d5`。工具检查此哈希及坐标网格，数据修订需要重新核对。`appearance_tests.cpp` 检查经纬度、接缝、量化、跨场景一致性和导览朝向，并用论文报告的前导低纬、拖尾高纬形态独立检查地图方向。

## 土星紫外主卵

`src/scene/saturn_aurora.cpp` 独立实现 [Carbary 2012, *The morphology of Saturn's ultraviolet aurora*](https://doi.org/10.1029/2012JA017670) 的式 (2)、表 1。输入是太阳地方时，纬度是行星中心纬度。论文统计 Cassini UVIS 2007–2009 年观测；模型在任意日期沿用这一时期的南北差异。

余纬（度）为 `A0 + Σ An cos[n π (LT − φn) / 12]`，LT 和相位单位为小时：

| 半球 | A0 | A1 | A2 | A3 | φ1 | φ2 | φ3 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 北 | 15.1 | 1.4 | 0.2 | 0.1 | 4.0 | 20.7 | 1.9 |
| 南 | 15.9 | 2.4 | 0.8 | 0.3 | 4.3 | 21.9 | 2.8 |

CPU 每半小时取样，GPU 周期线性插值。半高全宽用 `6 − 2 sin(π LT / 12)` 度；竖直柱亮度用 `10 + B exp(3[cos(π(LT−5)/12)−1])` kR，北 B=30、南 B=50。晨昏宽度和峰值取论文统计量，平滑函数及高斯横剖面由项目选择，并非逐地方时数据的拟合。不另画论文保留的日辉背景。

[Gérard et al. 2009](https://doi.org/10.1029/2008GL036554) 的 HST 夜侧边缘观测给出约 1100 km 的 H₂ 发光峰值。模型采用这一峰值，假设 Chapman 尺度高 200 km，在 500–3000 km 壳层内积分。高度按从一巴椭球沿地心射线增加的距离近似；纬度先撤销单位球坐标的极向拉伸。Carbary 的 1000 km 是投影高度，不当作发光峰值。

蓝色 `sRGB(0.3, 0.55, 1)` 和线性曝光增益 2（对应 40 kR 柱）是显示选择，昼夜都叠加在可见光云层和环上，不是肉眼颜色或绝对测光。使用既有 RSS 环剖面，以 `exp(−τ/μ)` 衰减相机到发光样本之间穿环的光线；这是灰色近似，RSS 射电光学厚度不等于紫外实测光学厚度。环为零厚度轴对称平面，不计算波长相关散射或极光对环的照明。

导览日期仅设置几何关系。没有太阳风、地球 Kp、风暴、SKR 振荡、螺旋、卫星足点、季节演化或随机闪烁。来源于 2026-11-03 核对：只使用公开发表的数值和公式，代码按项目 GPL-3.0-or-later 分发；未再分发论文图像、全文或 UVIS 原始观测。期刊的 Free Access 标记不当作 CC 或 Public Domain 许可。

## 土星可见光极光

面板 `Planetary aurora` 选择 `Visible light`，或在 `Event setting` 下打开土星场景的 `Saturn's visible aurora` 事件，可看夜侧粉红、紫色极光。该事件通过 `aurora_visible = true` 开启可见光发射。

依据 Dyudina, Ingersoll, Ewald & Wellington 2015, *Saturn's aurora observed by the Cassini camera at visible wavelengths*, [arXiv:1506.00664](https://arxiv.org/html/1506.00664v1)，[Icarus DOI](https://doi.org/10.1016/j.icarus.2015.05.022)，§1、§4 和图 8。2010-11-27 的观测显示距地平线几百公里为粉红色，1000 至 1500 km 为紫色，Hα 发射明显。可见光极光亮度约为日面的百万分之一至十万分之一。

渲染用两个 Chapman 层近似：峰值/尺度高分别为 450/140 km、1100/180 km，积分范围约 200 至 2000 km。低、高层的 sRGB 色值取 (1, 0.35, 0.55)、(0.65, 0.22, 1)，转成线性 RGB 后各乘 0.35。这些参数是显示假设，未拟合滤光片测光，也未把 UV kR 换算成可见光辐亮度。论文沿视线测到的色差包含可能的大气散射影响，模型将其近似为局部发光颜色。

横向复用 Carbary 统计卵及其相对亮度，固定于太阳地方时。表面太阳高度从 0° 到 −12° 平滑增强夜侧显示；此处理表达日光掩盖极光，未模拟昼夜发射率的变化。环透射采用与紫外模式相同的灰色近似。事件日期对应论文观测日，相机为示意视点；没有重建 Cassini 轨迹、螺旋亮弧、近共转结构或每小时增亮。仅引用公开结论，自行实现代码，未打包论文图像或观测数据。

## 木星紫外主卵

`jupiter_oval.txt` 来自 Linus Head 发布的 [Zenodo 12067692](https://doi.org/10.5281/zenodo.12067692)，许可为 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)。引用：Head et al. 2024, “Effect of magnetospheric conditions on the morphology of Jupiter's ultraviolet main auroral emission as observed by Juno-UVS”, A&A 688, A205，[论文](https://doi.org/10.1051/0004-6361/202450253)。数据作者为 Linus Head；Astraxis 将原始经纬度轮廓转换成等方位角采样的角半径表。

原始文件于 2026-11-03 下载，每个半球各 10,000 点：

| 文件 | SHA-256 |
|---|---|
| `UVSRefOval_N.csv` | `96023e6b935814acbb76c6dfa694a12e9b64731e338b2e18a22c7590e364e033` |
| `UVSRefOval_S.csv` | `5ae737ae83515b1e6d4563209424680af82259b501f9fed670ea519149bb95a4` |

复现（在仓库根目录运行，只需 Python 标准库）：

```powershell
New-Item -ItemType Directory -Force build/jovian_aurora
Invoke-WebRequest https://zenodo.org/records/12067692/files/UVSRefOval_N.csv -OutFile build/jovian_aurora/UVSRefOval_N.csv
Invoke-WebRequest https://zenodo.org/records/12067692/files/UVSRefOval_S.csv -OutFile build/jovian_aurora/UVSRefOval_S.csv
python -E tools/geomag/jupiter_oval_bake.py build/jovian_aurora assets/geomag/jupiter_oval.txt
```

CSV 的 `longitude` 是 System III 西经，`latitude` 是行星中心纬度。论文 4.1 节的参考点为北半球 (185°W, 74°N)、南半球 (32°W, 82°S)。烘焙工具先转成东经为正的本体单位向量，再绕各自参考点计算球面角距离与切平面方位角；切平面的 x 指东、y 指北。每圈保留 256 个等方位角样本，首尾周期插值。相对全部原始点的最大径向误差为北极 90.455 km、南极 82.410 km（按 71,492 km 换算），工具要求低于 100 km。

文本首行为 `AXJOVIAN1 256`，后面每行是北、南两个角半径（弧度）。第 k 行对应方位角 2πk/256。运行时将椭球单位坐标还原成行星中心方向后查表，南北轮廓各自固定在 System III 中。

论文 2 节使用距 1 bar 面 400 km 的投影高度，主发射区高度约 250 至 400 km；渲染的 Chapman 层峰值取 400 km，标高 80 km 为显示假设。高度按同一行星中心射线到基准椭球的距离计算，未做法线高度或折射修正。壳层覆盖约 150 至 1,000 km。横向高斯宽度 σ=200 km 为显示假设。论文报告的平均昏侧/晨侧亮度比约 2，用随地方时平滑变化的系数表示；蓝色、总显示亮度和噪声变化均为示意。主卵没有逐日期的活动序列；卫星足点另按下节计算。

## 木星卫星主足点

`jupiter_footprints.txt` 是离线计算的磁力线映射，不是观测坐标。内部场用 [Connerney et al. 2022 的 JRM33](https://doi.org/10.1029/2021JE007055)，取 13 阶；外部场用 [CON2020](https://doi.org/10.1029/2020JA028138)。参数依照 [Rabia et al. 2024 第 2 节](https://arxiv.org/html/2401.10115v1)：内外半径 7.8/51.4 RJ、半厚度 3.6 RJ、倾角 9.3°、倾向东经 155.8°，两项电流参数 139.6/16.7 nT，RJ=71,492 km。采用工具的 Edwards 解析近似，关闭平滑，方位分量为 `connerney`；未采用 L-MIC、磁层顶或随日期变化的电流参数。

计算工具为 Matthew Knight James 等人的 [JupiterMag 1.5.1](https://github.com/mattkjames7/JupiterMag) / libjupitermag，仅用于离线烘焙，程序不链接或分发它的 Python 包和 DLL。该版本包元数据写 GPL-3.0-or-later，随包 `LICENSE` 写 MIT、Copyright (c) 2021 Dr Matt James；这里保留两处原始声明，许可证副本见 [`tools/geomag/JupiterMag-LICENSE.txt`](../../tools/geomag/JupiterMag-LICENSE.txt)。这项声明涉及工具，不将数值输出标成第三方观测数据。

各卫星名义距离取场景的 421,800、671,100、1,070,400 km，径向采样为 97%、100%、103%，行星中心纬度为 −2°、0°、2°，东经每 1° 一个样本。每条线追踪南北两个方向，与半轴 (71,492+900, 71,492+900, 66,854+900) km 的椭球求交。表内存单位方向，顺序为卫星、半径、纬度、经度，每行六数是北 xyz、南 xyz。运行时做三线性插值并归一化；缺少对应的直属卫星或超出表域时隐藏足点。

领先角取 [Hue et al. 2023 表 1](https://arxiv.org/html/2304.14949v1)，DOI [10.1029/2023JA031363](https://doi.org/10.1029/2023JA031363)，该预印本标注 CC BY 4.0。作者为 V. Hue、G. R. Gladstone、C. K. Louis 等。输入是卫星的 System III **西经**；输出在赤道磁通管起点向**东**加，再查相应半球的表，不能直接加在足点经度上。Europa 北、Ganymede 南北的拟合频率 0.99/0.98 使原式首尾不闭合，0° 两侧各 15° 内以 smoothstep 混合两端延拓，最大改动分别为 0.047514°、0.252394°、0.146597°。其余经度保持原式，Io 的整数频率项不变。

复现（仓库根目录；依赖只装入构建目录）：

```powershell
python -m venv build/jovian_footprints/venv
build/jovian_footprints/venv/Scripts/python -m pip install JupiterMag==1.5.1
build/jovian_footprints/venv/Scripts/python tools/geomag/jupiter_footprint_bake.py
build/jovian_footprints/venv/Scripts/python tools/geomag/jupiter_footprint_check.py
```

生成文本按 LF 换行的 SHA-256 为 `83849f6e271485cc712d28f91da4278bc3cb10ca2a900f5116f00eb3709ac95c`。2026-11-03 本机烘焙检查：

| 卫星 | 未参与烘焙的 120 个随机位置，查表最大误差（北/南 km） | 同组位置若固定名义距离和零纬度的最大偏差（北/南 km） |
|---|---|---|
| Io | 5.618 / 7.274 | 647.241 / 638.088 |
| Europa | 5.688 / 5.179 | 445.319 / 406.575 |
| Ganymede | 6.310 / 5.137 | 379.744 / 371.000 |

比较基准为 4 倍细步长的直接追踪，随机种子固定，包含 0° 接缝两侧。另用 SciPy DOP853 和独立的椭球终止事件核对 12 个起点、共 24 个交点，和烘焙步长直接追踪的差异最大 1.073 km；这项检验共用磁场求值器，不是独立磁场模型。每颗卫星另抽查东经 0°、纬度 0.37° 的一个起点，解析近似与工具 hybrid 积分模式的差异（北/南）为 Io 0.292/0.218 km、Europa 16.790/19.504 km、Ganymede 42.146/30.666 km；这三个点不能给出全域误差上限。原生工具的足点插值曾返回空值或不收敛的坐标，烘焙脚本因此直接使用轨迹线段与椭球求交。

以上检验衡量数值误差，没有重新拟合 Juno 的观测点。磁场、电流片和领先角都是平均模型，原论文与观测的偏差是数百公里量级，不能把查表的公里级精度当成逐时预测精度。渲染以 900 km 为 Chapman 峰值，标高 366 km 是 Hue 引用的 Io 参考值，复用于另外两颗卫星；σ=250 km 的圆形斑点、6/2/4 的相对显示系数和蓝色均为示意。光柱沿径向，保留 150–3500 km，省略尾迹、反射翼和多重亮斑。

## 地球 Kp

| 文件 | 来源 | 许可 |
|---|---|---|
| `kp.bin` | GFZ Helmholtz Centre for Geosciences，行星际地磁指数 Kp（每 3 小时一个值，1932-01-01 起）。Matzka, J., Stolle, C., Yamazaki, Y., Bronkalla, O. and Morschhauser, A., 2021, The geomagnetic Kp index and derived indices of geomagnetic activity, Space Weather 19, e2020SW002641；数据集 Matzka et al. 2021, Geomagnetic Kp index, V. 1.0, GFZ Data Services, <https://doi.org/10.5880/Kp.0001> | CC BY 4.0（请注明 GFZ Helmholtz Centre for Geosciences） |

## 处理方式

`tools/geomag/kp_bake.py` 读 GFZ 的年度 WDC 格式文件（`kp1932.wdc` .. `kp2026.wdc`），只保留 Kp（以 1/3 为单位，0–27，每个三小时区间一个字节），从第一年的 1 月 1 日 00:00 UT 起连续排列，到最后一个文件的最后一天为止。WDC 文件里的太阳黑子数（SILSO，CC BY-NC 4.0）和 10.7 cm 射电流量不保留。

```
python -I tools/geomag/kp_bake.py <放 kpYYYY.wdc 的目录> assets/geomag/kp.bin
```

下载：GFZ 自己的服务器（`kp.gfz.de`、`ftp.gfz-potsdam.de`）在本机 TLS 握手失败，用 GWDG 的镜像 <https://ftp.gwdg.de/pub/geophys/kp-ap/wdc/yearly/>（格式说明在同目录的 `wdc_fmt.txt`）。当前文件到 2026-09-30 为止（2026-10-22 下载，没有缺值），更新时重新下载全部年份并重跑工具：GFZ 会修订最近的值。

## 程序里的其他地磁、极光常数

都写在 `src/scene/aurora.cpp`，注明了出处：

- IGRF-14 的偶极项 g10、g11、h11（1900–2025 每 5 年，加 2025–30 的长期变化）：IAGA Working Group V-MOD，<https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt>
- 极光卵边界：Starkov 1994 的系数，取自 Sigernes et al. 2011, Opt. Pura Apl. 44, 599, Table 1（开放获取，<https://folk.uib.no/st07350/papers/sigernes2011OPA599.pdf>）。表里 poleward 边界 A2 的 b0 印成“-4-44”，按 −4.44 用
- 各谱线的颜色（每 kR 的线性 RGB）：`tools/geomag/aurora_colors.py` 用 CIE 1931 2° 配色函数（CVRL，<http://cvrl.ucl.ac.uk/database/data/cmfs/ciexyz31_1.csv>，不入库）算出，以 5772 K 黑体代表的阳光为白色

## Europa 可见光氧极光

来源：de Kleer et al. (2023), *The Optical Aurorae of Europa, Ganymede, and Callisto*, PSJ 4, 37，doi:10.3847/PSJ/acb53c，[论文](https://arxiv.org/abs/2302.08487)，表 1、2；空间分布参考 de Kleer & Brown (2018), [Europa’s Optical Aurora](https://arxiv.org/abs/1809.04617)。只转录以下测量事实并自行实现模型，不分发论文、图像或地图，也不把论文素材声明为公有领域。

- 固定采用 2021-05-20 的盘平均氧谱线：630.0 nm 为 573±14 R，636.4 nm 为 190±7 R，557.7 nm 为 37±4 R。未把上限当作探测；未渲染近红外和 Hα。
- 沿用 `aurora_colors.py` 的 630.0、557.7 nm 每 kR 线性 RGB，636.4 nm 暂用 630.0 nm 色度并按 630/636.4 修正单光子能量；这是近邻红线合并近似，非完整光谱色度拟合。
- `satellite_aurora.hpp` 的 `aurora_disk_factor` 将包含临边外发射的可见球壳总通量除以实体盘面积，再将盘平均亮度换为竖直柱亮度。发光率随高度指数下降，H=50 km，上界 300 km，均为假设。测量孔径、PSF 和斑块未拟合。
- 本体方向 n 的横向权重为 `1 + 0.5*n.y`；东经 90°（西经 270°）的尾随半球最亮。3:1 对比是显示假设，球面平均为 1。归一化基准为朝木星半球；从其他方向看总亮度会变化，不强行维持盘平均值。
- 显示增益 `3e7` 独立于 Rayleigh 数值。夜侧可见，日照侧淡出；逐样本检查母行星的球形遮挡并平滑近似半影，使被木星遮日的日照半球也能显示。没有模拟木星扁率阴影、反照光、真实眼睛暗适应和极光照亮地面。
- 木星场景事件 18（从 0 计数）为 2021-05-20 13:40 TDB，测试确认整个发光层处于木星本影。观测表使用地球接收 UTC，而导览使用当地发射时刻的模拟 TDB；镜头不是望远镜视点，也不逐帧复现观测曝光。任意日期重复使用同组谱线强度，未实现电流片驱动或进出影后的大气变化。

模板开关为 `[aurora].europa_visible = true`，与其他极光字段互斥。可见光和紫外叠加视图中均保留这层可见光；没有把它改画成紫外伪彩。

## Io 与 Callisto 的可见光极光

`src/scene/satellite_aurora.hpp` 保存谱线的波长、Rayleigh 亮度、统计误差及每 kR 颜色。模板分别用 `[aurora].io_visible`、`[aurora].callisto_visible`，与其他极光字段互斥。普通事件默认不显示，导览的 `aurora_visible = true` 或手动可见光/UV 叠加模式启用；UV 模式仍保留这些可见光发射。

Io 来源为 Milby, de Kleer & Schmidt (2025), *Detection of New Auroral Emissions at Io and Implications for Its Interaction with the Plasma Torus*, doi:10.3847/PSJ/ade0b1，[论文](https://arxiv.org/abs/2508.04807)，表 2 的 2023-08-09 列。只采用这一天检测到的八条线，不拼接其他日期或填补未覆盖波段：

| 波长 nm | 发射 | 亮度 R | 统计误差 R |
|---|---|---:|---:|
| 406.9 | [S II] | 310 | 18 |
| 407.6 | [S II] | 87 | 11 |
| 458.9 | [S I] | 100 | 8 |
| 557.7 | [O I] | 461 | 12 |
| 589.0 | Na I | 12100 | 300 |
| 589.6 | Na I | 5070 | 140 |
| 630.0 | [O I] | 6090 | 110 |
| 636.4 | [O I] | 2200 | 20 |

这些是食中曝光的加权平均，误差不包含全部系统误差，尤其钠线。模型固定使用这组值，不随时间拟合钠消耗。Schmidt et al. (2023), *Io’s Optical Aurorae in Jupiter’s Shadow*, doi:10.3847/PSJ/ac85b0，[论文](https://arxiv.org/abs/2302.10849)，指出钠在入影后衰减、出影后恢复，时间尺度不同。当前没有实现这些化学过程，也没有用 SO₂ 凝结比例直接缩放氧发射。

Callisto 来源为 de Kleer et al. (2023), doi:10.3847/PSJ/acb53c，[论文](https://arxiv.org/abs/2302.08487)，表 2 的 2021-07-04 测量：630.0 nm 为 6.1±1.7 R，636.4 nm 为 3.3±1.2 R。其余谱线只有上限，不加入颜色。

空间模型均为自行实现的近似。Io 的氧/硫柱用 `1 + 0.75*(n.y² - n.z²)` 表示宽阔的赤道两侧增强，钠柱均匀；n 是本体单位方向。这不是磁场追踪，也不复现实际亮斑或火山发光。氧/硫指数标高 50 km、钠标高 100 km，上界均为 600 km。Callisto 使用均匀球壳，标高 50 km、上界 300 km，不添加极区形状。全部高度和横向对比为假设。

每层先按可见球壳（包含实体盘之外的临边发射）积分，换算为除以实体盘面积的孔径亮度。Io 横向权重在朝木星方向的投影积分中均值为 1，因 y²、z² 项抵消。没有拟合真实望远镜孔径和 PSF。独立显示增益为 Io `1e6`、Callisto `3e9`；不能按屏幕亮度比较它们的物理强弱。两者沿用 Europa 的逐样本球形母行星遮挡、近似半影和日照淡出。

导览事件 19、20（从 0 计数）分别为 Io 的 2023-08-09 13:40 TDB、Callisto 的 2021-07-04 14:00 TDB。测试确认整个发光层在木星本影。它们是相应观测日期的示意视点，不是地球望远镜位置，也不把接收 UTC 当作当地发射时刻。

### 配色复算与许可

`tools/geomag/aurora_colors.py` 增加上述钠、硫及 636.4 nm 波长，继续按 CIE 1931 2° 配色函数、5772 K 阳光白平衡和保持亮度的色域压缩换算。Io/Callisto 的红双线分别计算；Europa 原有的近邻红线近似保持不变。

CVRL 下载不可用时，本次从 Colour 的 `DATA_CMFS_STANDARD_OBSERVER["CIE 1931 2 Degree Standard Observer"]` 提取 1 nm 网格，写成 `wavelength,x,y,z` CSV 后输入原脚本。固定来源：[cmfs.py](https://github.com/colour-science/colour/blob/9072a98d04e048fe0ea2b1e3aa8429d6e95d37b6/colour/colorimetry/datasets/cmfs.py)，该文件 SHA-256 为 `a53076b9e6be4843443067986c925898d40a3c32c4eef9cb32ff08ca42e6493d`。Colour 为 BSD-3-Clause，声明见 `tools/geomag/Colour-LICENSE.txt`。不分发完整配色表、论文或图像，只保留计算后的常数和测量事实；不将论文素材声明为公有领域。


## 天王星与海王星的红外 H₃⁺

`uranus_h3p.txt`、`neptune_h3p.txt` 由 `tools/geomag/ice_giant_h3p_bake.py` 从下面的公开数据烘焙。太阳系场景的两个事件用 `aurora_infrared = true` 开启，也可以在 `Planetary aurora` 选 `UV/IR false color`。都是伪彩色，固定在本体上，不是实时模型。

天王星，CC BY 4.0：

- Melin et al. 2025, GRL 52, e2025GL118301（JWST GTO 1248，2023-01-08/09），[Zenodo 16085710](https://doi.org/10.5281/zenodo.16085710) 的 `uranus_superim_fits_2.5.fits.zip`：第 5 片是 H₃⁺ 总发射（W m⁻² sr⁻¹），2.5° 网格，已做视线修正，投影时去掉了黄昏亮区。列是 ULS 西经，行由南向北（IAU 意义，受照、被观测的极是 IAU 北、ULS 南）：图上约 60°、217°W 的暗斑和紧邻的亮区对应论文的“南极光和它极向旁边的暗斑”。约 5°N 以南的格子掉向边缘，噪声也大。
- Tiranti et al. 2026, GRL 53, e2025GL119304（JWST GO 5073，2025-01-19），[Zenodo 17036023](https://doi.org/10.5281/zenodo.17036023) 的 `results_lon_*.fits`：36 个经度、350 km 一格的 H₃⁺ 总发射（只在 ±25° 纬度的临边）。5025 km 以上论文说不可靠。文件里的数值和 Melin 的柱亮度差约 70 倍（见 CLAUDE.md），只用形状。
- Lamy et al. 2025, Nature Astronomy 9（doi:10.1038/s41550-025-02492-z）：自转周期 17.247864 ± 0.000010 h，`W = 203.81 − 500.9315936164617 d`（Zenodo 16085710 的 `uranus_lamy_IAU_new_lon.tpc`），写在场景里替换 IAU 的 −501.1600928。Melin 表 1 的 12 个中央经线在 1° 以内对上。

处理：先按距离平滑（σ ≈ 3°，约等于 1400 km 的空间分辨率，经度方向按 cos 纬度放宽，否则极区的投影台阶抹不掉）；2–18°N 过渡到 Tiranti 的临边经度曲线，按 15–35°N 的中位数定标（两者的经度曲线相关系数 0.73，低谷都在 195°W 附近）；35–25°S 淡到零，再往南没有观测，不发光。转到场景的本体系（`IAUnew 东经 = ULS 西经 − 50.6`）。高度表是 36 个经度的平均，上下各加一格零，归一成每 km 的柱份额。475 km 以下没有观测，线性降到 125 km 为零（均质层顶在 240–296 km）；顶部的零只是不让壳层有硬边。

海王星，论文和图的数值 CC BY 4.0（Melin et al. 2025, Nature Astronomy 9, 666），代码和 `figure_data` 在 [neptune-h3p](https://github.com/henrikmelin/neptune-h3p)（GPL-3.0，提交 `7c1c2b8`）：

- JWST GTO 1249，2023-06-22，两个视角（中央经线 131.0°W、303.2°W，时刻 06:20:54、14:01:13 UT）。`figure2c/d` 是每个像素的 H₃⁺ 辐射亮度，`figure2a/b` 里 3 μm 超过 20 μW m⁻² sr⁻¹ μm⁻¹ 的是云，`figure4_*` 是完全在盘内的像素的四角经纬度（跨 0° 经线的减过 3600°）。
- 单个像素信噪比约 3，不直接用。拟合每个视角一个平面的底光（晨昏渐变）加一个高斯亮区，得到中心 255°W、50°S，σ 12°（纬度）、34°（经度），峰值是底光的 1.64 倍；论文的范围是 60–30°S、200–280°W，约 2 倍，柱密度 1.7 倍。两个视角约束不住经向宽度（放开会到 60°），所以上限取论文的 80°。
- 图上只有底光和亮区：晨昏渐变随太阳固定，不放进本体图。经度差 −25.08° 是观测时刻（减光行时）场景的中央经线和 303.2°W 之差，`test_ice_giant_aurora` 打印；重烘时作为工具的第 2 个参数。

高度剖面没有观测，1500 km、标高 700 km 的 Chapman 层是假设（论文说模型预期在 1000 km 以上）。曝光、琥珀色和把行星反射光压到 3% 都只是显示选择。

下载和烘焙（放在 `build/`，不入库）：

```
mkdir -p build/h3p/lon build/h3p/nep
curl -L -o build/h3p/uranus_fits.zip https://zenodo.org/api/records/16085710/files/uranus_superim_fits_2.5.fits.zip/content
for l in $(seq 5 10 355); do curl -L -o build/h3p/lon/results_lon_$l.fits https://zenodo.org/api/records/17036023/files/results_lon_$l.fits/content; done
for f in figure2a_neptune_longitude1_clouds figure2b_neptune_longitude2_clouds figure2c_neptune_longitude1_h3p figure2d_neptune_longitude2_h3p figure4_coords_longitude1 figure4_coords_longitude2 figure4_polygons_longitude1 figure4_polygons_longitude2; do curl -L -o build/h3p/nep/$f.fits https://raw.githubusercontent.com/henrikmelin/neptune-h3p/7c1c2b80ee04e411d35a7ceaa7cd56bef5e5cb85/figure_data/$f.fits; done
python -E -P tools/geomag/ice_giant_h3p_bake.py uranus build/h3p/uranus_fits.zip build/h3p/lon assets/geomag/uranus_h3p.txt
python -E -P tools/geomag/ice_giant_h3p_bake.py neptune build/h3p/nep -25.08 assets/geomag/neptune_h3p.txt
```

## 行星内部磁场（磁力线）

`field_*.txt` 是施密特准归一化的内部场高斯系数（nT），每行 `g|h n m 值`。天王星、海王星、木星、土星的系数取自 [libinternalfield](https://github.com/mattkjames7/libinternalfield)（MIT，Copyright (c) 2022 Dr Matt James，提交 `9a6b6e5`），地球取自 IAGA 的 IGRF-14。程序只读这些数，不链接那个库。

| 文件 | 模型 | 参考半径 | 坐标 |
|---|---|---|---|
| `field_uranus_ah5.txt` | AH5，[Herbert 2009](https://doi.org/10.1029/2009JA014394)，4 阶 | 25,600 km | ULS（北是 IAU 的南） |
| `field_neptune_o8.txt` | O8，[Connerney et al. 1991](https://doi.org/10.1029/91JA01165)，取约束较好的 3 阶 | 24,765 km | NLS |
| `field_jupiter_jrm33.txt` | JRM33，[Connerney et al. 2022](https://doi.org/10.1029/2021JE007055)，取约束较好的 13 阶（文件到 30 阶） | 71,492 km | System III |
| `field_saturn_cassini11.txt` | Cassini 11，[Dougherty et al. 2018](https://doi.org/10.1126/science.aat5434)，11 阶，轴对称 | 60,268 km | 自转轴 |
| `field_mercury_anderson2012.txt` | [Anderson et al. 2012](https://doi.org/10.1029/2012JE004159)，4 阶，轴对称 | 2,440 km | 自转轴 |
| `field_ganymede_weber2022.txt` | [Weber et al. 2022](https://doi.org/10.1029/2022GL098633) 的偶极（Galileo 和 Juno），加四极项只略有改善，所以不用 | 2,631 km | 经度 0 朝木星，东向是后随面 |
| `field_earth_igrf14.txt` | IGRF-14 的 2025.0，13 阶 | 6,371.2 km | 地心 |

**许可**：libinternalfield 的 MIT 只覆盖它自己的代码和整理，不覆盖模型本身。这里入库的是论文发表的高斯系数（几十到不足两百个数），和场景里其他论文表格的数值一样按事实引用，注明出处。原文的条款：[JRM33](https://doi.org/10.1029/2021JE007055) 是 CC BY-NC 4.0，AH5、O8、Cassini 11、Anderson 2012、Weber 2022 是期刊的订阅条款。需要完全开放的木星模型时，可以换 Sharan et al. 2022（16 阶，[Zenodo 6564162](https://doi.org/10.5281/zenodo.6564162)，CC BY 4.0）。IGRF 是公开模型。

这些模型的经度按东经（右手）写：Connerney 系的模型西经 λ 对应 φ = −λ。依据是系数本身：O8 的 atan2(h11, g11) = −72°，论文写磁极倾向西经 72°。

核对：
- 从 libinternalfield 转来的数没有逐项对照论文的表（AH5 的原文没有拿到）。它们的总体性质对照过：O8 的偶极 0.142 G、倾角 46.9°，偏心 (0.05, 0.48, 0.00) R 对上 [A&A 2024 的 OTD-O8 表](https://www.aanda.org/articles/aa/full_html/2024/06/aa46699-23/T1.html)；Cassini 11 严格轴对称；IGRF-14 和 `scene/aurora.cpp` 的地磁极一致。
- AH5 和 Cassini 11 的参考半径没有从原文核实。天王星 25,600 与 25,559 km 只差 0.16%，画面上看不出来。

放到场景的本体系（`[bodies.magnetic_field]`）：
- **天王星**：`flip_pole` 把 ULS 的北翻到 IAU 的南，`longitude_offset_deg = -50.6` 用 Melin 2025 的 `IAU new 东经 = ULS 西经 − 50.6`。AH5 的南磁极落在 235.7°W、纬度 30.2°（IAU），Melin 的南极光在 232°W、37.5°，相距约 8°。
- **海王星**：`prime_meridian_deg = [253.18, 536.3128492]` 是 IAU 2009 以前的 W，周期 16.11 h，和 NLS 同一个系统；场景的 541.14°/天不是磁场的自转率。[PDS 的 Voyager 2 轨迹说明](https://pds.nasa.gov/data/vg2-n-pra-4-summ-browse-48sec-v1.0/vg_1001/catalog/pos/dataset.cat)规定 1989-08-25 03:56:00 飞船在西经 167.7°，按场景的轨迹和这个 W 算出 168.8°。现在的相位未知。
- **木星、土星、地球、水星**：模型的坐标系就是场景的本体系。木星只用内部场时，Io 的足点和 `jupiter_footprints.txt` 相差不到 3°。水星的模型算出偏心 479 km，[Anderson et al. 2011](https://doi.org/10.1126/science.1211001) 是 484 km。
- **Ganymede**：模型的坐标系也是本体系，这是推断出来的，不是论文明写的。[Kivelson et al. 2002](https://doi.org/10.1006/icar.2002.6834) 说南半球的磁极从朝木星的子午面向后随面转 24°，Weber 的系数给出 25.7°，只有 x 朝木星、y 朝后随面时才对得上（如果按 GPhiO 的 x 沿等离子体流方向读系数，磁极就差了约 65°）。场景的 x 轴朝木星、y 轴朝后随面，`test_magnetic_field` 检查。纯偶极的外层线落在纬度 48°，Marzok 的极光在后随面是 42–56°，在前导面只有 16–25°：极光在开闭磁力线的边界上（[Duling et al. 2022](https://doi.org/10.1029/2022GL101688)），而这条边界被木星的等离子体推得不对称，纯偶极表现不出来。

水星的线画到 L = 1.4，都在 [Winslow et al. 2013](https://doi.org/10.1002/jgra.50237) 的平均日下点磁层顶（1.45 个半径，太阳风动压变化时 1.35–1.55）以内。

重新生成（只用标准库）：

```bash
mkdir -p build/field
for f in uranus/ah5 neptune/gsfco8 jupiter/jrm33 saturn/cassini11 mercury/anderson2012 ganymede/weber2022dip; do curl -L -o build/field/$(basename $f).dat https://raw.githubusercontent.com/mattkjames7/libinternalfield/9a6b6e562b57a95146b7d608ee9cdf229b674a63/data/coeffs/$f.dat; done
curl -L -o build/field/igrf14coeffs.txt https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt
python -I tools/geomag/field_coeffs.py build/field assets/geomag
```

## 径向磁场图

`*.fieldmap` 是固定在本体上的磁场图，由 `tools/geomag/field_map_bake.py` 烘焙（需要 numpy）。径向图的文件头 `AXFIELDMAP1 宽 高 上限_nT` 之后是 宽×高 个有符号字节，值为 127·sign(Br)·√(|Br|/上限)；强度图的文件头 `AXFIELDMAG1 宽 高 下限_nT 上限_nT` 之后是 127·(2t − 1)，t = (|B| − 下限)/(上限 − 下限)。投影和天体贴图相同：行由北向南按网格的参数纬度，列从 180°W 向东，所以每格在椭球面上的真实位置求值。

| 文件 | 模型 | 位置 | 范围 |
|---|---|---|---|
| `mars_br.fieldmap` | [Morschhauser et al. 2014](https://doi.org/10.1002/2013JE004555)，110 阶，参考半径 3393.5 km | 椭球面上方 150 km | −709 .. 1278 nT |
| `jupiter_br.fieldmap` | JRM33，13 阶 | 1 bar 椭球面 | −13.9 .. 21.7 G |
| `jupiter_br_deep.fieldmap` | JRM33，13 阶 | 0.81 RJ 的球面，方向和云顶图相同（两张图在网格上对齐） | −79.4 .. 81.3 G |
| `moon_br.fieldmap` | [Tsunakawa et al. 2015](https://doi.org/10.1002/2014JE004785) 的球谐展开，Wieczorek 的 T2015_449 取到 240 阶，参考半径 1737.4 km | 月面上方 30 km | −25 .. 14 nT |
| `earth_b.fieldmap` | IGRF-14 的 2025.0，13 阶，总强度 |B| | 椭球面 | 22,073 .. 66,953 nT |

**月球**：系数取自 Wieczorek 在 Zenodo 的 [T2015_449.sh.gz](https://doi.org/10.5281/zenodo.3873648)（CC BY 4.0，单位是特斯拉）。取 30 km 而不是月面：直接画月面时，强异常周围有大致沿南北向的条纹（硬截断的振铃加上看来是轨道条带的噪声，在 120–240 阶之间做余弦衰减也去不掉），向上延拓到接近观测的高度后条纹消失。地表的峰值约 −279 nT，30 km 处只有 −25 nT，差得很多，展签按 30 km 写。

**地球**：画总强度而不是 Br。Br 基本只是偶极的南北两半，南大西洋异常看不出来。最弱处在 26.5°S、59.5°W。

**待做的地壳图**：[Thébault 的地球岩石圈模型](https://doi.org/10.5281/zenodo.5546528)（CC BY 4.0，约 40 km 分辨率）比主磁场弱得多，要单独一张图。

**火星模型的许可**：系数取自 Mark Wieczorek 在 Zenodo 的 [Morschhauser2014.txt.gz](https://doi.org/10.5281/zenodo.3876495)（CC BY 4.0，与论文附表 ts01 相同）。没有用分辨率更高的 Langlais 2019 和 Gao 2021：两篇论文都是 CC BY-NC（-ND），和项目的 GPL 不相容。Zenodo 上也有 Langlais 模型的 CC BY 版本，但它的许可和论文的 NC-ND 对不上，所以不用。

取 150 km 而不是地表：这里接近观测的高度（Mars Global Surveyor 的低轨和大气制动段），往下延拓到地表时小尺度放大，整张图布满高频细纹。

核对：
- 木星两张图和 `field_jupiter_jrm33.txt` 逐格一致（`test_magnetic_maps`，到一个字节的量化）。
- Great Blue Spot 在云顶约 −12.8 G、6.5°S、84.5°W；[Moore et al. 2018 的会议摘要](https://meetingorganizer.copernicus.org/EPSC2018/EPSC2018-556.pdf)写它在赤道附近、System III 约 90°。在扁椭球面上它不是全球最小值（南极附近 −13.9 G），在 1 RJ 的球面上才是。
- 火星的峰值在 48°S、175°E，属于 [Acuña et al. 1999](https://doi.org/10.1126/science.284.5415.790) 说的 Terra Cimmeria 和 Terra Sirenum；Hellas、Argyre 内的场不到峰值的 10%。
- 月球的峰值在 20.2°S、236.8°E，紧邻 Gerasimovich 环形山（22.9°S、122.6°W）；Reiner Gamma（7.5°N、59°W）处有明显的异常，Mare Crisium 处很弱；远月面的场比近月面强。
- 地球图和 `field_earth_igrf14.txt` 逐格一致（到一个字节的量化）。

重新生成（先按上一节取得 `field_jupiter_jrm33.txt`、`field_earth_igrf14.txt`）：

```bash
mkdir -p build/field_map
curl -L -o build/field_map/Morschhauser2014.txt.gz https://zenodo.org/api/records/3876495/files/Morschhauser2014.txt.gz/content
gunzip -k build/field_map/Morschhauser2014.txt.gz
curl -L -o build/field_map/T2015_449.sh.gz https://zenodo.org/api/records/3873648/files/T2015_449.sh.gz/content
python -E -P tools/geomag/field_map_bake.py build/field_map/Morschhauser2014.txt build/field_map/T2015_449.sh.gz assets/geomag
```

## Mars visible diffuse aurora

Enable with `[aurora] mars_visible = true` in the Mars body template. The Mars
scene includes **Mars visible aurora (2024-03-18)**; manual Visible light and
UV overlay also show this visible emission.

- Knutsen et al. (2025), *Detection of visible-wavelength aurora on Mars*,
  Science Advances 11, eads1563, https://doi.org/10.1126/sciadv.ads1563
  (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC12077521/).
- SuperCam measured OI 557.7 nm at 101 +/- 14 R, elevation 20 degrees,
  beginning 2024-03-18 06:46 UTC. Vertical dust optical depth was 1.1.
  The paper's radiative transfer model estimates **2.1 kR at the orbital limb**.
  That modeled limb estimate sets our reference intensity; 101 R is a
  dust-attenuated ground measurement and is not used as a zenith column.
- `src/scene/mars_aurora.hpp` assumes a horizontally uniform Gaussian layer,
  peak 60 km, sigma 15 km, truncated at 20 and 140 km. These heights approximate
  diffuse SEP emission; they are not a retrieved 2024 profile. A numerical
  tangent-ray integral normalizes its peak to 2.1 kR for radius 3396.19 km.
  The renderer approximates constant physical radial altitude on Mars' ellipsoid.
- 557.7 nm linear RGB per kR is the existing `aurora_colors.py` spectral conversion.
  Display gain is independently set to 3e7. Daylight masking fades between
  solar elevations -12 and 0 degrees; it represents visibility, not excitation.
- The event time is rounded to 06:47 TDB; its orbital camera is illustrative.
  The reference emission stays constant while enabled. There is no storm time
  series, terrestrial Kp driver, crustal-field structure, red-line upper-limit
  emission, dust transfer, multiple scattering or surface illumination.
- No paper images or datasets are redistributed; the implementation uses the
  cited scalar measurements and an independently written illustrative profile.
