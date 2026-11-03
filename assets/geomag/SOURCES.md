# 地磁与极光数据来源

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
