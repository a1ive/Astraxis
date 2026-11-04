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
- `europa_aurora.hpp` 将包含临边外发射的可见球壳总通量除以实体盘面积，再将盘平均亮度换为竖直柱亮度。发光率随高度指数下降，H=50 km，上界 300 km，均为假设。测量孔径、PSF 和斑块未拟合。
- 本体方向 n 的横向权重为 `1 + 0.5*n.y`；东经 90°（西经 270°）的尾随半球最亮。3:1 对比是显示假设，球面平均为 1。归一化基准为朝木星半球；从其他方向看总亮度会变化，不强行维持盘平均值。
- 显示增益 `3e7` 独立于 Rayleigh 数值。夜侧可见，日照侧淡出；逐样本检查母行星的球形遮挡并平滑近似半影，使被木星遮日的日照半球也能显示。没有模拟木星扁率阴影、反照光、真实眼睛暗适应和极光照亮地面。
- 木星场景事件 18（从 0 计数）为 2021-05-20 13:40 TDB，测试确认整个发光层处于木星本影。观测表使用地球接收 UTC，而导览使用当地发射时刻的模拟 TDB；镜头不是望远镜视点，也不逐帧复现观测曝光。任意日期重复使用同组谱线强度，未实现电流片驱动或进出影后的大气变化。

模板开关为 `[aurora].europa_visible = true`，与其他极光字段互斥。可见光和紫外叠加视图中均保留这层可见光；没有把它改画成紫外伪彩。
