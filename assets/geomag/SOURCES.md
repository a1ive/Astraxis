# 地磁与极光数据来源

`Planetary aurora` 的 `Event setting` 由活动事件的 `aurora_visible`、`aurora_ultraviolet` 控制。普通事件和无活动事件时关闭所有极光。`Visible light` 手动开启地球、Europa 和土星的可见光极光；`UV aurora overlay` 开启紫外叠加，同时保留可见光发射（土星改用紫外卵）。地球极光还需要开启 `Atmospheres`。

## Ganymede 自身的紫外极光

木星场景末尾的两个 Ganymede 导览分别朝向前导与拖尾半球（前导那个按 12 分钟/秒走，能看到电流片扫过），也可在面板的 `Planetary aurora` 中选择 `UV aurora overlay`。太阳系场景通过木卫三模板使用同一份数据。

`ganymede_aurora.txt` 来自 Marzok、Schlegel、Saur、Roth、Grodent、Strobel 和 Retherford 的 [Zenodo 6564687 数据集](https://doi.org/10.5281/zenodo.6564687)，原文件为 `Fig5_BrightnessMap.txt` 和 `Fig6_{Above,Inside,Below}CurrentSheet.txt`，许可 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)。引用：Marzok et al. 2022, *Mapping the Brightness of Ganymede's Ultraviolet Aurora Using Hubble Space Telescope Observations*, JGR Planets 127, e2022JE007256，[论文](https://doi.org/10.1029/2022JE007256)。再分发依据是 Zenodo 数据许可；论文正文另采用 CC BY-NC，未随程序分发。

该图合并了 1998 至 2017 年的 46 次 HST/STIS 观测，记录 OI 135.6 nm 的垂直柱亮度（Rayleigh）。数组为 80×40，像元中心西经 2.25° 至 357.75°、纬度 −87.75° 至 87.75°，两轴步长均为 4.5°。烘焙保留网格和正亮度，将负背景残差截为零；原图中零值区域仍为零。运行时把东经取负后采样，纬度从南向北，经度环绕、纬度夹紧。文件有四层（`AXAURORAMAP4`，平均图、电流片上方、内部、下方），上传线性 UNORM 纹理的四个通道，范围为 0 至 260 R，量化误差约 0.5 R，远景使用 mipmap。

**分组图**：论文第 3.2 节按木卫三在木星磁层里的磁纬分组，±6° 以内算电流片内部（19 次曝光），以上、以下各 14、13 次。磁纬是 VIP4 的偶极（倾斜 9.5°，朝东经 159.2°），论文表 1 的磁纬只由 System III 经度决定，这个偶极对得上（`test_ganymede_aurora`，8 行在 0.08° 内）。各组不能直接画：

- 覆盖不一致：表 1 的 46 次曝光里只有 6 次朝向亚木、反木一侧，180–210°W 几乎只有上方组看过（内部组在 180–225°W 只覆盖 11% 的格子），上方、内部、下方组的非零格分别缺 9.5%、10.8%、20% 的平均图格子；
- 噪声：每组只有约三分之一的光子，单格的差别大多不显著。论文也说下方组 0° 附近的孤立亮格几乎不显著。

所以每组 = 平均图 × 比值：比值是组和平均图在“这一组看过的格子”上按 σ = 8° 平滑后的比（`(Σ K G + 5 R·ΣK) / (Σ K M + 5 R·ΣK)`，负残差留在求和里，否则暗区会被抬高），没有数据处趋于 1；再逐格缩放，让三组按木卫三在各组停留的时间（0.282、0.437、0.282，由偶极和 ±6° 算出；HST 的曝光份额是 0.304、0.413、0.283）加权平均后严格回到平均图。这样缺测区没有接缝，组之间切换时也不会有噪声闪烁，代价是单格的结构全部来自平均图，比 σ 小的组间差别被抹平。比值范围 0.48–1.96。

**核对亮度衔接**：按论文式 (8) 的卵带（表 2 的拟合纬度上下各三格）算内部/外部的亮度比，未经处理的组图是下游（前导，40–150°W）1.21、上游（拖尾，240–340°W）0.85，烘焙后是 1.20、0.86，论文是 1.3 ± 0.31、0.78 ± 0.19。卵向赤道的移动（论文下游 −5.4° ± 3.2°）在北半球保留了约 3°（原图约 4°）。

运行时 `Simulation::moon_aurora` 按此刻的磁纬选层，在 5–7° 之间平滑过渡（`ganymede_sheet_shares`，真实的边界不是硬的，论文的 ±6° 只是分组）；倍速快到每秒超过 0.25 个电流片周期（约 5.25 小时）开始退回平均图，1 个周期以上完全是平均图（`plasma_sheet_steadiness`），免得闪烁。

发光分布采用论文 §2.3 校正视线柱长时使用的均匀 70 km 壳层，垂直积分还原图中的柱亮度。地图固定于木卫三本体；蓝色及每 100 R 的显示增益是示意。数据的低信噪比像元也参与双线性插值，部分暗区来自缺测，不能解释为没有极光。没有模拟摇摆、风暴或细小斑块。

重新生成（需要 numpy）：

```powershell
New-Item -ItemType Directory -Force build/ganymede_aurora | Out-Null
foreach ($f in "Fig5_BrightnessMap", "Fig6_AboveCurrentSheet", "Fig6_InsideCurrentSheet", "Fig6_BelowCurrentSheet") { curl.exe -L --fail "https://zenodo.org/api/records/6564687/files/$f.txt/content" -o "build/ganymede_aurora/$f.txt" }
python -E -P tools/geomag/ganymede_aurora_bake.py build/ganymede_aurora assets/geomag/ganymede_aurora.txt
```

工具检查四个原文件的 SHA-256（写在工具开头）及坐标网格，数据修订需要重新核对。`appearance_tests.cpp` 检查经纬度、接缝、量化、跨场景一致性和导览朝向，用论文报告的前导低纬、拖尾高纬形态独立检查地图方向，并检查上面的磁纬、时间加权平均、内外亮度比和缺测区的回退。

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

渲染用两个 Chapman 层近似：峰值/尺度高分别为 450/140 km、1100/180 km，积分范围约 200 至 2000 km，高度是假设（图 8 只有一次观测的三个采样框，没有剖面）。

颜色和两层的比例来自图 8，由 `tools/geomag/saturn_visible_colors.py`（`python -I`，只用标准库）算出，数值是按 300 dpi 渲染后定位误差棒读的，约 ±10%：

- 颜色按论文附录 A.2 的定义：RED、GRN、BL1 换算成 I/F 就当作真实颜色，和渲染器“线性 RGB 是 I/F”一致。低层取“brightest”框（RED 1.00、GRN 0.48、BL1 0.54），高层取“high”框（0.93、0.41、1.08），从 R/nm 换到 I/F 时除以 5772 K 黑体在滤光片中心的光子数，得到每单位 clear 亮度的线性 RGB (1.000, 0.491, 0.692) 和 (0.930, 0.420, 1.384)。不从滤光片重建光谱再积分 CIE：连续谱的插值方式会让高层的色相从蓝变到品红，而浅粉（绿约为红的一半）在各种做法里都成立。
- 论文图上采样框的颜色 (242, 0, 29)、(204, 0, 134)、(148, 0, 255) 绿色为 0，和同一张图 GRN 约为 RED 一半的测光对不上，看来是拉伸过的，所以这里比那些框淡得多。
- 比例按图上标的 clear 亮度：侧看时高层在 1100 km 处是低层峰值的 11/150，高层的柱量因此是低层的 0.075。“high”框横跨更暗的经度，另外两个框压在亮斑上，所以真实比值可能大几倍，这是高层的下限。
- 曝光 0.35 是显示选择，未把 UV kR 换算成可见光辐亮度。低处“deep”框更红，论文说可能是大气散射，这里不画散射，低层用“brightest”的颜色。

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

以上检验衡量数值误差，没有重新拟合 Juno 的观测点。磁场、电流片和领先角都是平均模型，原论文与观测的偏差是数百公里量级，不能把查表的公里级精度当成逐时预测精度。渲染以 900 km 为 Chapman 峰值，标高 366 km 是 Hue 引用的 Io 参考值，复用于另外两颗卫星；σ=250 km 的圆形斑点、6/2/4 的相对显示系数和蓝色均为示意。光柱沿径向，保留 150–3500 km，省略尾迹。

### 次级亮斑

Io 的足点最多有三个亮斑（[Bonfond et al. 2008](https://doi.org/10.1029/2007GL032418)，GRL 35, L05107；[2009](https://doi.org/10.1029/2009JA014312)，JGR 114, A07224；读的是 ORBi 仓库 `2268/3862`、`2268/28544` 的作者稿）。`JovianAurora::Spot` 按论文的解释各用一条规则，由上面的南北领先角推出，没有新的数据：

- **主斑（MAW）**：卫星朝这个半球发出的阿尔文翼到达处，就是上一节的足点。
- **电子束亮斑（TEB）**：对面半球的阿尔文翼把电子反向加速，电子沿磁力线几十秒内穿过等离子体片，落在这个半球。它们不受高密度的环拖慢，所以落点在“对面主斑所在的那条磁力线”上：领先角取对面半球的值。卫星离哪个半球远（那一侧的翼要穿过更多等离子体，领先角大），电子束亮斑就在那个半球的主斑**前面**（论文的 leading spot），在另一个半球的主斑后面（论文的第一个次级亮斑）。两个半球的间距因此在赤道上相等，论文第 10 段的观测正是如此。
- **反射亮斑（RAW，只画 Io 的）**：朝对面半球去的翼在环的远侧边界反射回来。延迟几乎都花在环里，所以它比主斑多走对面那一半两次：领先角取 `本半球 + 2 × 对面半球`。

和论文对照（`test_jovian_footprints`）：

| 量 | 这里 | 论文 |
|---|---|---|
| Io 主斑与电子束亮斑在赤道上的最大间距 | 4.93° | 约 4°（2008 第 4 段） |
| Io 反射亮斑落后主斑的最大值 | 12.6°（北）、14.4°（南） | 最多 12°（2008 第 16 段） |
| Io 北半球的电子束亮斑领先主斑的经度范围 | System III 西经 285–105° | 0–100°（2008 图 2） |
| 南半球的 | 105–285° | 130–300° |
| 2007-03-05 09:06 UT（Io 在西经 27.2°，JPL Horizons，扣除 44.26 分钟光行时） | 间距 4.9° | 3.1°、3.4° |
| Ganymede 两斑在行星上的最大间距 | 8,100 km | 所见最大 4,300 km（[Bonfond et al. 2013](https://doi.org/10.1002/grl.50989)，GRL 40, 4977；读的是 HAL `hal-00867977` 的作者稿） |

没有合上的地方：

- 2007-03-05 的例子里论文写的是“南半球的领先亮斑、北半球的次级亮斑”，而按 Horizons 那时 Io 在西经 27°，照同一篇的图 2（和这里的规则）领先亮斑应在北半球。这里只比较了间距，方向的出入没有查明。
- Ganymede 的间距是观测的近两倍：Hue 的南北两条拟合曲线相差最多 14°，Bonfond 2013 按线性理论预期约 10°，HST 在 2007 年覆盖的经度也有限，论文还指出间距随等离子体片的密度在几天内变化。两条各有几度不确定度的拟合相减，误差会放大。
- Europa 的次级亮斑没有观测报道，这里按同一规则画出（Bonfond 2013 认为 Io 的机制是普遍的）。

亮度：电子束亮斑取主斑显示系数的 0.3，反射亮斑 0.2，都是示意。Io 的领先亮斑实测 0.6–1.9 GW、典型 0.7 GW（2008 第 7 段），通常比主斑暗，但在南半球 0–100° 处可以和主斑一样亮（2009）；这些随经度的变化没有画，亮斑也一直画出，间距小时并入主斑。更远的多次反射和尾迹（Io 的 e 折长度约 21,000 km，[Bonfond et al. 2017](https://doi.org/10.1002/2017JA024370)）没有画。

## 卫星到足点的线

`jupiter_flux_tubes.bin` 由 `tools/geomag/jupiter_flux_tube_bake.py` 离线烘焙，磁场、工具和参数都和上一节的足点表相同（JRM33 13 阶 + CON2020，JupiterMag 1.5.1），用同一个 venv。每颗卫星 72 个东经（每 5°）、名义半径、赤道面上起一条线，南北各到 900 km 的椭球，重采样成 40 个点（`∫ds/r` 均匀，线之间逐点对应，近行星处更密）。文件头 `AXFLUXTUBE2 72 40`，之后按卫星依次是：小端 float32 的 x、y、z（单位 71,492 km）、w，按经度、半球（北、南）、点排列；再是该卫星名义位置上各经度的磁场 x、y、z（nT，同一模型，由 `jm.Internal.Field` 与 `jm.Con2020.Field` 相加），供 Ganymede 自身磁场叠加用（见“行星内部磁场”）。每条线绕 z 转了负的起点经度，所以相邻的线插值时不抄近路。线的脚和足点表逐条一致（`test_flux_tubes`，不到 0.1 km）。

**瞬时磁力线和亮斑的差别**：磁场随木星转，比卫星快，卫星扰动的是不断扫过它的磁通管，扰动再以阿尔文波沿磁通管传下去，要走几分钟到半小时（Hue 2023 表 2：Io 3.9–15.4、Europa 4.1–19.5、Ganymede 9.3–35.1 分钟）。到达电离层时，这条磁通管已经扫到卫星前面（东边）了，所以亮斑落在“卫星此刻所在的磁力线”的脚前面，赤道上相差的就是领先角。场景画两种线：

- **淡色**：此刻穿过卫星的磁力线，`flux_tube(..., false)`，落在 `instant_footprint`，那里没有亮斑。
- **橙色**：波此刻的位置，`flux_tube(..., true)`。走了时间份额 w 的那部分波，是卫星在这条线东边 δ·w 处（δ 是 Hue 的领先角）时发出的，所以第 k 个点取“经度 + δ·w_k 那条线上的第 k 个点”。它横穿磁力线，终点正好是 `footprint()` 的亮斑。

**w 的来源**：时间 ∝ ∫√ρ / B ds，B 取自追踪，ρ 取 Bagenal & Delamere 2011（JGR 116, A05209）的等离子体片：电子密度 `n0(ρ) exp(−(z/H(ρ))²)`，ρ 是到自转轴的距离，z 是离心赤道以上的高度（`centrifugal_equator_height_rj` 的切线近似，和 Europa、Callisto 极光用的是同一个）。

- `n0 = 1987 (ρ/6)^−8.2 + 14 (ρ/6)^−3.2 + 0.05 (ρ/6)^−0.65` cm⁻³，`H = 10^(−0.116 + 2.14 x − 2.05 x² + 0.491 x³ + 0.126 x⁴)` RJ，x = log₁₀(ρ/6)。2026-11-11 对过原文：密度是式 (1)，标高是式 (6)，系数是表 1 的 `1987, −8.2, 14, −3.2, 0.05, −0.65` 和 `−0.116, 2.14, −2.05, 0.491, 0.126`，垂直分布是式 (4) 的 `n0 exp(−(z/H)²)`。**式 (6) 下面印的是 `r = log10(R)`，漏了 /6**：照字面算 Io 处的标高是 3.9 RJ，而原文图 5 的曲线在 Io 处约 0.75 RJ、10 RJ 处约 1.8 RJ、50 RJ 以外趋于 3.9 RJ，只有 `log10(R/6)` 才对得上（0.77、1.83、3.9 RJ）；在三颗卫星处是 0.74、1.68、2.78 RJ，Callisto 处 3.66 RJ，和本页别处引用的 1.7、3.7 RJ 一致。原文取平均离子质量 20 amu（硫、氧离子的混合，式 (3) 下）。
- 拟合 Hue 2023 的领先角（t = P_syn δ / 360，P_syn 为 12.89、11.22、10.53 h）用两类参数：三颗卫星**共用**一个乘在 H 上的厚度系数，取使三者相对均方差之和最小的 1.55；每颗卫星一个“每个电子对应的质量”。

| 卫星 | 对 Hue 的均方根差 | 相关系数 | Hue 的范围 | 模型的范围 | 每个电子的质量 | H 不放大时 |
|---|---:|---:|---|---|---:|---|
| Io | 0.25° | 0.99 | 1.8–7.2° | 2.0–6.8° | 26.6 amu | 0.54°，50 amu |
| Europa | 1.10° | 0.91 | 2.2–10.4° | 3.6–9.9° | 18.4 amu | 1.13°，35 amu |
| Ganymede | 1.77° | 0.94 | 5.3–20.0° | 6.2–21.5° | 23.1 amu | 2.33°，40 amu |

和以前的做法比：以前 ρ 是沿磁力线**弧长**的高斯，每颗卫星各拟合一个宽度，Io 0.24°、Europa 1.13°、Ganymede 2.95°，Europa 和 Ganymede 的宽度约束不住（Ganymede 停在扫描上限）。电流片里的磁力线被拉得很长，弧长和离开赤道的高度不成比例，所以那个宽度没有物理含义。现在的高度是真正的 z，标高取自文献，一个共用的系数就让 Ganymede 的差降了四成，Io、Europa 持平；每个电子 18–27 amu 和以硫、氧离子为主、平均一到两价的等离子体相符（S⁺ 32、O⁺ 16、S²⁺ 16）。系数 1.55 说明要么实际的片比 Bagenal & Delamere 的厚，要么模型缺了什么（例如切线近似的离心赤道、卫星附近的局部减速），这里不区分。Europa 的相关系数仍只有 0.91：Hue 北半球的拟合频率是 0.99，曲线本身不闭合。离心赤道换成每条线离自转轴最远的点，结果几乎不变（差 0.01–0.04°）。

这只决定橙线的弯法，端点不受影响：起点是卫星，终点是 Hue 的领先角和足点表给的亮斑。

**运行时**：表只有名义半径和赤道面上的线。实际卫星的线从表里按经度插值后拉到卫星（偏移量乘 (r − r_foot)/(r_0 − r_foot)，偶极的线上偏移正比于 r），再拉到足点表给的脚（1 − 上面的比值，再取 8 次方），最后用 Catmull-Rom 曲线每段细分 3 份（外侧的点相距约 2 RJ）。工具用未参与烘焙的 60 个位置（随机的经度，半径和纬度在各卫星轨道的范围内：Io ±0.6%、±0.1°，Europa ±1.2%、±0.6°，Ganymede ±0.3%、±0.3°，取自根数加一点余量）和 2 倍细步长的直接追踪对照，最大差距（点到折线）：

| 卫星 | 拉过之后 | 不拉（只按经度插值） |
|---|---:|---:|
| Io | 2,045 km | 2,760 km |
| Europa | 9,546 km | 10,913 km |
| Ganymede | 18,286 km | 11,469 km |

不拉时差距主要在卫星一端，而且线穿不过卫星；拉过之后差距在线的中段。电流片两翼的磁场几乎是径向的（Ganymede 处），纬度差一点，整条线的 L 就差很多，这是非线性的，线性地拉消不掉：曾把测试范围放到足点表的 ±3%、±2° 时差到 8 万 km。卫星走不到那里，真要到那里（如另一个星历）时线会明显偏，但足点仍然准。

重新生成（在上一节的 venv 里，约 10 秒）：

```powershell
build/jovian_footprints/venv/Scripts/python tools/geomag/jupiter_flux_tube_bake.py
```

生成文件的 SHA-256：`e860c3f91e5885daabf2e578cd55309fe1bf442412e156e652ae598ab1783e0b`。

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
- `satellite_aurora.hpp` 的 `aurora_disk_factor` 将包含临边外发射的可见球壳总通量除以实体盘面积，再将盘平均亮度换为竖直柱亮度。发光率随高度指数下降，H = 150 km（Roth et al. 2016 第 4.2.4 节：临边以上 135.6 nm 的轮廓偏向更大的标高，150 km 最合，接近 O₂ 的标高），上界 600 km 是假设。测量孔径、PSF 和斑块未拟合。
- 横向分布见下面的“亮度图案”。图案按球面平均归一化，从任何方向看的总亮度都由等离子体决定，个别方向相差几个百分点。
- 显示增益 `3e7` 独立于 Rayleigh 数值。夜侧可见，日照侧淡出；逐样本检查母行星的球形遮挡并平滑近似半影，使被木星遮日的日照半球也能显示。没有模拟木星扁率阴影、反照光、真实眼睛暗适应和极光照亮地面。
- 木星场景的事件 “Europa's oxygen aurora in eclipse” 为 2021-05-20 13:40 TDB，测试确认整个发光层处于木星本影。观测表使用地球接收 UTC，而导览使用当地发射时刻的模拟 TDB；镜头不是望远镜视点，也不逐帧复现观测曝光。没有进出影后的大气变化：紫外亮度在日照和食中差不多（Roth et al. 2016），Milby 2026 也认为食中 O₂ 柱密度稳定。
- **等离子体片调制**：亮度 ∝ 离心赤道处的电子，`exp(−(z/H)²)`，H = 1.7 RJ（Bagenal & Delamere 2011 的模型值）。de Kleer 2023 单次日食拟合 1.67 RJ，Roth 2016 紫外 1.6 RJ；[Milby 2026 博士论文](https://thesis.caltech.edu/18581/)（Caltech，版权保留，只引用数值；也见 Milby et al. 2026, PSJ 7, doi:10.3847/PSJ/ae47f4）第 3 章十个夜晚中单夜的拟合为 1.24–1.70 RJ，图 3.10 把各夜归一后用 1.7 RJ 拟合得最好，散布 ±11%，他也把 1.7 RJ 当长期平均。z 用论文式 1.6、1.11（Phipps & Bagenal 2021 的离心赤道的切线近似，`centrifugal_equator_height_rj`），符号和他的表 3.1、Milby 2025 的表 1 一致，8 行对照在 0.025 RJ 内（`test_plasma_sheet`）。上面的 573 R 等是 2021-05-20 整段观测的平均，当时 z 从 0.003 到 −0.683 RJ，所以赤道上是它们除以这段的平均密度 0.95（`europa_reference_density`）。这样从赤道到最远处（约 1.15 RJ）只差到 0.63，和他图 3.5 的虚线相同。论文的单夜平均差到两倍，这里没有：夜间变化和每次 ±11% 的方位不均匀都没有画。de Kleer 指出 2021-05-20 这一天恰好看不出相关（6 月、7 月才有），这是把它当作归一化的已知代价。

**亮度图案**：来自 Roth et al. 2016, *Europa's far ultraviolet oxygen aurora from a comprehensive set of HST observations*, JGR Space Physics 121, doi:10.1002/2015JA022073（71 张 HST/STIS 图，20 次观测），只引用数值。可见光线和 135.6 nm 都来自 O₂ 的解离激发，Europa 的光学观测也说亮度偏向下午一侧（de Kleer & Brown 2018、de Kleer 2023），所以直接借用。

`europa_pattern`（`satellite_aurora.hpp`，`aurora.frag.hlsl` 有同样的一份）：`(1 + a n·d) exp(s n.z (α + β n.y) + κ ((n·b)² − n.z²))`，n 是本体方向（x 朝木星，y 朝拖尾一侧），d 是下午方向（z × 太阳），b 是穿过卫星的磁力线（`jovian_field_line_at_moon`，VIP4 偶极），s 是在等离子体片里的位置（磁纬 / 9.5°，离赤道最下方 +1，最上方 −1）。

- **下午一侧更亮**（`a`）：几乎所有图都是右半边亮，1.59 ± 0.33（图 8a），右边总是下午一侧。它跟着太阳，不跟着本体；在食中太阳在木星后面，下午一侧正好是拖尾一侧（以前的 `1 + 0.5 n.y` 只在食中才对）。论文说这个不对称的原因还不清楚。
- **面向等离子体片的半球更亮**（`α`、`β`）：离赤道上方时南边亮，下方时北边亮。拖尾（上游）一侧明显：图 7c 虚线是磁通管电子能量比的代理，10^(0.655 cos(λIII − 16°))，也就是 4.5 和 0.22，分箱约 3.6 和 0.26；前导（下游）一侧弱：图 7d 约 1.8 和 0.7。论文解释为南北磁通管的电子能量不等。
- **沿磁力线集中**（`κ`）：磁极（穿过卫星的磁力线和球面的交点）亮、磁赤道暗，象限比在投影磁场倾斜 20° 时约 1.6（图 7a）。这里只取随磁力线倾斜的部分 `(n·b)² − n.z²`，不是磁极本身的亮斑：薄壳的临边增亮让均匀分布的“极区/赤道”比就有约 1.4，已经超过观测的 0.9–1.2（图 8c），再加亮极区只会更差。论文也说球对称的壳既解释不了盘内亮度也解释不了临边轮廓。

系数由 `tools/geomag/europa_pattern_fit.py` 按图 7、8 拟合：前向模型是光学薄的壳加 σ = 0.08 RE 的 PSF，在论文的三个视向、论文的区域里量比值，比值都是从图上读的（`ln` 差的权重写在工具里，食中那张权重 0.25）。结果 a = 0.425、α = 0.575、β = 0.838、κ = 1.925：

| 量 | 视向 | 模型 | 论文（读图） |
|---|---|---|---|
| 南北比（下方、上方） | 拖尾 | 3.40、0.29 | 约 3.6、0.26 |
| 南北比（下方、上方） | 前导 | 1.56、0.64 | 约 1.8、0.7 |
| 象限比（\|β\| ≈ 20°） | 拖尾、前导 | 1.9、1.6 | 约 1.6 |
| 右/左 | 过赤道时 | 1.57 | 1.59 ± 0.33 |
| 南北比（下方、上方） | **亚木（食中）** | **1.81、0.55** | **约 0.33、3** |
| 极区/赤道 | 拖尾、前导 | 1.28、1.28 | 约 0.9、1.2 |

后两行没有合上。食中（图 7e，8 个短曝光点）的南北比是反的，论文猜是几何造成的：面向等离子体片的亮区转到了卫星背后。把磁力线的倾斜放大能往这个方向走（κ 加大时亚木视向的比值降到 1.28、0.78），但象限比就过头了；这里优先对上数据多的图 7a。事件正好是亚木视向，所以画面里的南北方向未必和观测一致。拖尾、前导的极区/赤道差异（图 8c、8d，1.3 倍）和盘内的斑块没有画。

**没有用上的**：[Roth 2021](https://doi.org/10.1029/2021GL094289)（GRL 48, e2021GL094289）说日照的拖尾半球中心 H₂O/O₂ = 12–22，H₂O 约 3×10¹⁵ cm⁻²、临边没有，而且多半来自升华。它改变的是紫外线比，而 135.6 nm 的亮度还是 O₂ 为主（表 2：中心 61 R 里 O₂ 占 55 R），这里只画夜侧和食中的可见光，Milby 2026 也没有在亚木半球看到 H₂O 的 Hα，所以图案里没有它。同一篇的 visit 17 用来核对等离子体片：从食中（146–161°W）到出影后（194–209°W），135.6 nm 从 40.2 ± 3.4 降到 36.3 ± 3.5 R（10 ± 11%），这里的模型降 23%，差 1.2σ（`test_europa_pattern`）。Roth 2016 拟合的等离子体片标高 1.8 ± 0.1 RJ（135.6 nm）、1.6 ± 0.1 RJ（130.4 nm）也和 1.7 RJ 相符，他们 1999 年的亮度高 50%，这里没有长期变化。

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

这些是食中曝光的加权平均，误差不包含全部系统误差，尤其钠线。

**钠的入影衰减和出影恢复**：Schmidt et al. (2023), *Io’s Optical Aurorae in Jupiter’s Shadow*, doi:10.3847/PSJ/ac85b0，[论文](https://arxiv.org/abs/2302.10849)：入影后 e 折叠约 10 分钟、约 30 分钟到稳定值（一小时以上的中值 1.7 kR，是前 10 分钟的 27%），出影后 e 折叠约 2 小时（他们认为是 NaCl 的光化学被打断）。`io_sodium_level` 就是这些：入影时 1，稳定值 0.19（由“27%”和 10 分钟的衰减反推），表里的 17.17 kR 对应 7 次曝光的平均水平 0.2073（Milby 2025 表 7 的文件名给出开始时刻，减去光行时；场景的入影在 12:47:37 TDB，曝光在其后 25–68 分钟，`test_io_eclipse_levels`）。Milby 的钠比 Schmidt 的亮一个量级（入影很久以后 Schmidt 是 1.7 kR，Milby 的 2022-11-24 出影前是 22.5 kR），而氧只差约 1.6 倍，Milby 也说钠有较大的系统误差；这里只取 Schmidt 的形状，绝对值跟着这张表，不调和两者。钠柱和亮度当作成正比，食初光学厚的影响（Schmidt 第 5 节）没有计入。食外（夜侧）用入影前的水平，约 83 kR：这只是外推，手动打开可见光时才看得到。

**SO₂ 的蓝色连续谱**：[Milby 2025](https://arxiv.org/abs/2508.04807) 表 3 把 Geissler et al. (2004) 的 Cassini 图换算成亮度，BL1 滤光片（390–500 nm）25 ± 10 kR，而这个波段的原子线只有 0.5 kR，Geissler 认为差额主要是 SO₂ 的电子碰撞准连续谱（可能还有 S₂）。这里把 24.5 kR 平摊在波段上（`aurora_colors.py` 的 `so2_bl1`），和氧、硫同一层和同一图案（Milby 表 3：BL1 只有赤道亮区）。它随 SO₂ 变：[de Pater et al. 2020](https://doi.org/10.3847/PSJ/abb93d)（ALMA）的盘积分流量入影时指数下降，食中比日照时低 2–3 倍，取 2.5，e 折叠约 10 分钟（Schmidt 第 5.1.1 节）；出影后约 10 分钟线性恢复（`io_so2_level`）。Cassini 的 BL1 图在入影后 38 分钟（2001-01-05 11:03 UTC 入影，11:41:15），所以 24.5 kR 对应这时的水平 0.413。**这是假设发射跟着大气主体走**：Milby 2025 由原子线反推的 SO₂ 在入影后并不马上下降，30–75 分钟处有些下降、2 小时前又回升，但这和离木星边缘的距离相关，他们怀疑是背景扣除的伪影，并提出发射可能来自更高的高度。两者还没有同时观测过。

**氧和硫不变**：Schmidt：氧对光照“无动于衷”，只跟着等离子体密度；但他们的相关系数只有 0.53–0.55，Milby 2025 拟合常数和 0.77 RJ 的高斯一样好（图 10），且论文说夜间变化到 2–3 倍。所以 Io 的氧不做等离子体片调制，和 Europa 不同，这是观测的差别，不是遗漏。

**亮斑的位置**：Milby 表 3 下面的讨论和 Retherford 2000 等：Io 的赤道亮区在亚木和反木一侧（等离子体流的两翼），不在前导和拖尾一侧（亚木亮斑可能再偏 20° 尾流方向，没有画）。图案是 `1 + 0.75*(n.x² - n.z²)`，x 轴朝木星（测试检查）；从前写成 n.y² 是错的。

**Callisto**：来源为 Milby 2026 博士论文第 5 章（上面 Europa 的同一篇）：2021-07-04 和 2025-11-29 的 HIRES，单次曝光的 630.0 nm 第一次够看出时间变化（表 5.2）。用 `exp(−(z/H)²)` 拟合 2025 年的曝光，赤道上 19 ± 2 R、H = 3.6 ± 0.7 RJ，几乎就是 Bagenal & Delamere 模型的 3.7 RJ；2021 年的数据也合得上。636.4 nm 取 630.0 nm 的 1/3.09（Wiese et al. 1996，表 5.3 平均 8.2、3.2 R 合得上）。其余谱线只有上限，不加入颜色。de Kleer et al. (2023) 表 2 的 6.1、3.3 R 是同一次 2021-07-04，背景扣除后来改进了（论文表 5.3 为 7.9、2.8 R）。论文指出最靠近赤道的三个点偏高，且都在食中较晚，也可能是电离层在食中塌陷使电子进入更多。2026-11 查过：没有空间分辨的观测（HST 17739 是无分辨的 COS 光谱，还没有发表）。

离心赤道的表 5.2 距离：场景算出 2021 年差 0.03 RJ，2025 年差约 0.3 RJ（相当于木星转 8 分钟），原因没有查明；亮度上影响很小。

空间模型均为自行实现的近似。Io 的氧/硫柱见上，钠柱均匀；n 是本体单位方向。这不是磁场追踪，也不复现实际亮斑或火山发光。氧/硫指数标高 50 km、钠标高 100 km，上界均为 600 km。Callisto 使用均匀球壳，标高 50 km、上界 300 km，不添加极区形状。全部高度和横向对比为假设。

每层先按可见球壳（包含实体盘之外的临边发射）积分，换算为除以实体盘面积的孔径亮度。Io 横向权重在半球上均值为 1，因 x²、z² 项抵消。没有拟合真实望远镜孔径和 PSF。独立显示增益为 Io `1e6`、Callisto `1.5e9`；不能按屏幕亮度比较它们的物理强弱。两者沿用 Europa 的逐样本球形母行星遮挡、近似半影和日照淡出。

运行时 `Simulation` 在时钟离开上一个周期时才重新搜索 Io 的入影和出影（Io 中心进入木星等效球的本影，5 分钟一步再二分）；超过 2 小时/秒不再搜索，钠和 SO₂ 取食外的水平，Europa 和 Callisto 的调制超过约 1 小时/秒开始退回长期平均。

木星场景的两个事件：Io 的 2023-08-09 12:57 TDB（入影后约 9 分钟，1 分钟/秒，钠的衰减正好看得见），Callisto 的 2021-07-04 14:00 TDB（2 分钟/秒，展签结束时到达赤道）。测试确认开场时整个发光层在木星本影。它们是相应观测日期的示意视点，不是地球望远镜位置，也不把接收 UTC 当作当地发射时刻。

### 配色复算与许可

`tools/geomag/aurora_colors.py` 增加上述钠、硫及 636.4 nm 波长，继续按 CIE 1931 2° 配色函数、5772 K 阳光白平衡和保持亮度的色域压缩换算。Io/Callisto 的红双线分别计算；Europa 原有的近邻红线近似保持不变。

CVRL 下载不可用时，本次从 Colour 的 `DATA_CMFS_STANDARD_OBSERVER["CIE 1931 2 Degree Standard Observer"]` 提取 1 nm 网格，写成 `wavelength,x,y,z` CSV 后输入原脚本。固定来源：[cmfs.py](https://github.com/colour-science/colour/blob/9072a98d04e048fe0ea2b1e3aa8429d6e95d37b6/colour/colorimetry/datasets/cmfs.py)，该文件 SHA-256 为 `a53076b9e6be4843443067986c925898d40a3c32c4eef9cb32ff08ca42e6493d`。Colour 为 BSD-3-Clause，声明见 `tools/geomag/Colour-LICENSE.txt`。不分发完整配色表、论文或图像，只保留计算后的常数和测量事实；不将论文素材声明为公有领域。


## 天王星与海王星的红外 H₃⁺

`uranus_h3p.txt`、`neptune_h3p.txt` 由 `tools/geomag/ice_giant_h3p_bake.py` 从下面的公开数据烘焙。天王星、海王星场景的事件用 `aurora_infrared = true` 开启（天王星的 1986 和 2007 年两个事件是同一张图换成当时的太阳，不是观测），也可以在 `Planetary aurora` 选 `UV/IR false color`。都是伪彩色，固定在本体上，不是实时模型。

天王星，CC BY 4.0：

- Melin et al. 2025, GRL 52, e2025GL118301（JWST GTO 1248，2023-01-08/09），[Zenodo 16085710](https://doi.org/10.5281/zenodo.16085710) 的 `uranus_superim_fits_2.5.fits.zip`：第 5 片是 H₃⁺ 总发射（W m⁻² sr⁻¹），2.5° 网格，已做视线修正，投影时去掉了黄昏亮区。列是 ULS 西经，行由南向北（IAU 意义，受照、被观测的极是 IAU 北、ULS 南）：图上约 60°、217°W 的暗斑和紧邻的亮区对应论文的“南极光和它极向旁边的暗斑”。约 5°N 以南的格子掉向边缘，噪声也大。
- Tiranti et al. 2026, GRL 53, e2025GL119304（JWST GO 5073，2025-01-19），[Zenodo 17036023](https://doi.org/10.5281/zenodo.17036023) 的 `results_lon_*.fits`：36 个经度、350 km 一格的 H₃⁺ 总发射（只在 ±25° 纬度的临边）。5025 km 以上论文说不可靠。文件里的数值和 Melin 的柱亮度差约 70 倍（见 CLAUDE.md），只用形状。
- Lamy et al. 2025, Nature Astronomy 9（doi:10.1038/s41550-025-02492-z）：自转周期 17.247864 ± 0.000010 h，`W = 203.81 − 500.9315936164617 d`（Zenodo 16085710 的 `uranus_lamy_IAU_new_lon.tpc`），写在场景里替换 IAU 的 −501.1600928。Melin 表 1 的 12 个中央经线在 1° 以内对上。
- Zenodo 16085710 的 `JupiterMag_Uranus_L_shells.txt`、`JupiterMag_Uranus_dip_equator.txt`：三个磁场模型（Q3、AH5、UH5，文件里写作 gsfcq3、ah5、umoh）的 L 壳等值线（3、5、10、20，南北各一圈）和偶极赤道。这里取 Q3，不是场景磁力线用的 AH5：论文第 4 节说暗 V 形的下半段和 Q3、UH5 的偶极赤道对得上，和 AH5 对得较差；图上也如此，偶极赤道上比旁边 10–15° 暗 3.9%（Q3）、2.9%（UH5）、0.9%（AH5）。表用的是 2 阶的 Q3（[libinternalfield](https://github.com/mattkjames7/libinternalfield) 提交 `9a6b6e5` 的 `gsfcq3.dat`，8 个系数写在 `test_ice_giant_aurora` 里；到 3 阶的 `gsfcq3full.dat` 对不上）。**文件里 L 的标签是反的**：`uranus_magnetic_field.ipynb` 按升序的等值线 `[3, 5, 10, 20]` 配了 `['20', '10', '5', '3']`，所以标 3 的是紧贴磁极的最内圈。用 Q3 系数从这些点追踪磁力线，标 3、5、10、20 的 Rmax 中值约 20、10、5、2.9（ULS 北），论文也说南极光“就在 L = 20 以内”。工具里按真实的 L 用，并断言这个顺序。

处理：先按距离平滑（σ ≈ 3°，约等于 1400 km 的空间分辨率，经度方向按 cos 纬度放宽，否则极区的投影台阶抹不掉）；2–18°N 过渡到 Tiranti 的临边经度曲线，按 15–35°N 的中位数定标（两者的经度曲线相关系数 0.73，低谷都在 195°W 附近）；35–25°S 过渡到下面的模型，得到“受照时”的一层；另有一层是没有阳光时剩下的，场景按此刻的太阳混合（见下面的“极夜和季节”）。转到场景的本体系（`IAUnew 东经 = ULS 西经 − 50.6`）。高度表是 36 个经度的平均，上下各加一格零，归一成每 km 的柱份额。475 km 以下没有观测，线性降到 125 km 为零（均质层顶在 240–296 km）；顶部的零只是不让壳层有硬边。

25°S 以南没有观测，图上是模型，不是数据：

- **磁坐标**：等值线和偶极赤道的点（相距至少 4°，共 397 个）取 ±√(1 − 1/L)，赤道取 0，ULS 北为正，用球面上的薄板样条插值成全球的光滑坐标，再取反正弦当不变纬度。对照逐格追踪 Q3：25°S 以北均方根差 0.02–0.06，以南 0.2，分歧集中在第二个偶极赤道回路外的一道新月形（越过它 L 会突变，四条等值线描述不了），对亮度的影响最多几个百分点。Q3 只来自一次飞掠，那里本来就约束很弱，所以没有为此另写追踪器。
- **亮度**：对 8°N 以北的观测格子拟合 `ln E = 常数 + a(|纬度|) + b(不变纬度)`（纬度 7.5° 一档，不变纬度 5° 一档，交替取平均，b 再做一次 [¼, ½, ¼] 平滑），得 a 为 0.76–1.23、b 为 0.92–1.17。a 按纬度的绝对值，也就是假设南半球和观测到的北半球对称；论文说观测到的极区偏暗，温度从 25°（ULS 南）往极区每度降 0.40 K，另一个极没有数据。ULS 北磁极的内层壳没有被观测到，b 取另一个磁极相同的壳。
- **依据**：在观测图上留出几个扇区（8–25°N、25–40°N、四个经度段）做预测（这次的 b 按 √(1 − 1/L) 分 20 档，形式相同），只用常数的对数均方根误差 0.14–0.34，加 a 降到 0.10–0.19（270–360°W 那一段没有改善），再加 Q3 的 b 总体不变：北极光一侧变好（0.132 → 0.115），270–360°W 变差（0.154 → 0.188）。所以受照时亮度的大尺度分布主要由纬度决定，b 只带来弱的结构：L = 2–6 的壳比 L = 1.1–1.6 亮约 9%，偶极赤道的暗带只有 2–3%。
- **对照论文**：b 是沿 L 壳的平均，不是真实的形态。论文（和 HST）看到的极光是有限范围的斑点，不是完整的卵或环；南极光比周围亮 64–94%，北极光平均高 0.12 µW m⁻² sr⁻¹，这些局部亮区按壳平均后被冲淡了。论文把 L < 3 算作非极光区，南半球最亮的格子也落在追踪得到的 L ≈ 3.7 上。
- **没有用共轭点映射**：论文第 3、4 节明确说南北极光在三个模型下都不共轭（也提醒两者不是同时观测的）。图上也如此：两端都有观测的磁力线（用 AH5 追踪），亮度残差的相关系数在足点相距 5–30° 时是 0.70（只是相邻的格子本来就相似），30–60° 时 0.15，60° 以上 −0.36；待补的格子里，另一端有观测的也只占 12%。
- **极夜和季节**：上面得到的是“受照时”的图 V。论文说非极光的 H₃⁺ 来自太阳光电离，寿命约 1 小时（Moore 2019），比 17.2 小时的自转短得多，所以它跟着阳光走；而 Tiranti 的曲线按 Melin 的受照亮度定标，也没有纬度信息。于是 `E = S·V + (1 − S)·D`，文件是两层（`AXAURORAMAP2`，纹理的红、绿通道），着色器逐样本按此刻的太阳混合（`AuroraPattern::BodyMapSeasons`，`aurora_map_pattern` 是同样的公式，测试用它）：
  - S 只看纬度，取一天里最亮的时候，也就是正午的太阳天顶角 χ = |日下点纬度 − 纬度|。按 χ 线性插值三个点（常数在 `aurora_model.cpp`）：χ = 50.83° 时 S = 1，这是观测可信的最南处 8°N 在 Melin 观测时的值（日下点 58.83°N，`test_ice_giant_aurora` 检查；论文的 55.1° 是 JWST 的视线纬度）；χ = 90° 时 S = 0.29，是 Tiranti 2026 接近晨昏线的柱密度 6.6×10¹⁴ m⁻² 比 Melin 2025 受照圆面的中值 2.27×10¹⁵ m⁻²；χ = 120.2° 时 S = 0，太阳落到 4000 km 高度的地平线以下，4000 km 是 Tiranti 引的 JWST 时代晨昏线附近光电子加热的主要高度。
  - D 是只有粒子电离时剩下的，只在 L ≥ 3 的壳上（L = 2–3 之间淡入）：有观测处（2–18°N 过渡）取图上超出背景（常数 + a + 非极光水平，L < 3 的 b 平均 1.004）的部分，所以论文的南北极光都留下来，最高 0.17 µW m⁻² sr⁻¹；没有观测处取 b 超出非极光水平的部分乘全图水平，这是沿壳的平均，最高约 0.05。
  - 混合后，2023 年的图和之前固定压暗的版本最多差 0.002。那时 31°S 是正午的晨昏线，61°S 以南除了壳上的极光就是黑的，偶极赤道的暗带也只留在还有阳光的纬度；正午晨昏线处的亮度是北半球同纬度的 0.31，L = 2–6 比 L = 1.1–1.6 亮 3.2 倍（磁场转 180° 后 0.43）。1986 年（场景的日下点 −82.1°）反过来，北半球极夜里只剩 Melin 的南极光（场景里约 0.1）。2007 年 12 月（昼夜分点）两个半球差不多亮，极区只被太阳掠过。

海王星，论文和图的数值 CC BY 4.0（Melin et al. 2025, Nature Astronomy 9, 666），代码和 `figure_data` 在 [neptune-h3p](https://github.com/henrikmelin/neptune-h3p)（GPL-3.0，提交 `7c1c2b8`）：

- JWST GTO 1249，2023-06-22，两个视角（中央经线 131.0°W、303.2°W，时刻 06:20:54、14:01:13 UT）。`figure2c/d` 是每个像素的 H₃⁺ 辐射亮度，`figure2a/b` 里 3 μm 超过 20 μW m⁻² sr⁻¹ μm⁻¹ 的是云，`figure4_*` 是完全在盘内的像素的四角经纬度（跨 0° 经线的减过 3600°）。
- 单个像素信噪比约 3，不直接用。拟合每个视角一个平面的底光（晨昏渐变）加一个高斯亮区，得到中心 255°W、50°S，σ 12°（纬度）、34°（经度），峰值是底光的 1.64 倍；论文的范围是 60–30°S、200–280°W，约 2 倍，柱密度 1.7 倍。两个视角约束不住经向宽度（放开会到 60°），所以上限取论文的 80°。
- 本体图上只有底光和亮区，经度差 −25.08° 是观测时刻（减光行时）场景的中央经线和 303.2°W 之差，`test_ice_giant_aurora` 打印；重烘时作为工具的第 2 个参数。
- 晨昏渐变随太阳固定，不放进本体图：论文说亮度从晨到昏（图上从左到右）大范围增加，可能是向阳面的累积，没有给数字。每个像素的晨昏坐标 x = cos(纬度) sin(中央经线 − 西经)（东边是下午，太阳和 JWST 的方向只差约 2°，日下点 20.8°S、视线 20.2°S），只由经纬度算，没有用“从左到右”这句话，算出来也是沿图的第二个下标增大。拟合每个视角一个底光、一个共同的梯度和上面的亮区，得 k = 0.228（`aurora_model.cpp` 的 `kNeptuneDuskGradient`，工具打印），避开亮区后两个视角单独拟合是 0.21、0.27，和 x 的相关系数 0.46、0.73。亮区在这个模型下是 1.62 倍（平面的拟合是 1.64 倍），所以图不用改。场景里是 `图 + k x`（`AuroraPattern::BodyMapDusk`）：渐变只加在底光上，亮区不跟着变。

高度剖面没有观测，1500 km、标高 700 km 的 Chapman 层是假设（论文说模型预期在 1000 km 以上）。曝光、琥珀色和把行星反射光压到 3% 都只是显示选择。

下载和烘焙（放在 `build/`，不入库）：

```
mkdir -p build/h3p/lon build/h3p/nep
curl -L -o build/h3p/uranus_fits.zip https://zenodo.org/api/records/16085710/files/uranus_superim_fits_2.5.fits.zip/content
curl -L -o build/h3p/L_shells.txt https://zenodo.org/api/records/16085710/files/JupiterMag_Uranus_L_shells.txt/content
curl -L -o build/h3p/dip_equator.txt https://zenodo.org/api/records/16085710/files/JupiterMag_Uranus_dip_equator.txt/content
for l in $(seq 5 10 355); do curl -L -o build/h3p/lon/results_lon_$l.fits https://zenodo.org/api/records/17036023/files/results_lon_$l.fits/content; done
for f in figure2a_neptune_longitude1_clouds figure2b_neptune_longitude2_clouds figure2c_neptune_longitude1_h3p figure2d_neptune_longitude2_h3p figure4_coords_longitude1 figure4_coords_longitude2 figure4_polygons_longitude1 figure4_polygons_longitude2; do curl -L -o build/h3p/nep/$f.fits https://raw.githubusercontent.com/henrikmelin/neptune-h3p/7c1c2b80ee04e411d35a7ceaa7cd56bef5e5cb85/figure_data/$f.fits; done
python -E -P tools/geomag/ice_giant_h3p_bake.py uranus build/h3p/uranus_fits.zip build/h3p/lon build/h3p/L_shells.txt build/h3p/dip_equator.txt assets/geomag/uranus_h3p.txt
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
| `field_earth_igrf14.txt` | IGRF-14 的全部历元：1900–2025 每 5 年（1945–2020 为 DGRF）加 2025–30 的长期变化率，13 阶（1995 年以前和变化率到 10、8 阶，其余为零） | 6,371.2 km | 地心 |

**许可**：libinternalfield 的 MIT 只覆盖它自己的代码和整理，不覆盖模型本身。这里入库的是论文发表的高斯系数（几十到不足两百个数），和场景里其他论文表格的数值一样按事实引用，注明出处。原文的条款：[JRM33](https://doi.org/10.1029/2021JE007055) 是 CC BY-NC 4.0，AH5、O8、Cassini 11、Anderson 2012、Weber 2022 是期刊的订阅条款。需要完全开放的木星模型时，可以换 Sharan et al. 2022（16 阶，[Zenodo 6564162](https://doi.org/10.5281/zenodo.6564162)，CC BY 4.0）。IGRF 是公开模型。

**地球随年份变化**：`field_earth_igrf14.txt` 的格式和其他文件不同，一行 `epochs` 列出历元，一行 `rate_until 2030.0`，之后每行 `g|h n m` 加各历元的值和最后的变化率（nT/年），模板里用 `coefficients_by_year` 读。按 IGRF 的规定，历元之间线性内插，2025 年以后按变化率外推到 2030 年；更早和更晚保持两端的场（`MagneticFieldSeries`）。`Scene::update` 按时钟最近的整年换模型、重新追踪磁力线并重画强度图（`follow_magnetic_field_year`），刚加载、还没更新时是 2025.0。一年内偶极轴移动不到 0.1°，整年的台阶看不出来。`test_magnetic_maps` 检查：1900 年的系数、各年的偶极和 `scene/aurora.cpp` 的地磁极一致（那里是同一批偶极项的另一份表）、2030 年等于 2025 年加五年的变化率、回到 2025 年时强度图和入库的文件逐字节相同。

这些模型的经度按东经（右手）写：Connerney 系的模型西经 λ 对应 φ = −λ。依据是系数本身：O8 的 atan2(h11, g11) = −72°，论文写磁极倾向西经 72°。

核对：
- 从 libinternalfield 转来的数没有逐项对照论文的表（AH5 的原文没有拿到）。它们的总体性质对照过：O8 的偶极 0.142 G、倾角 46.9°，偏心 (0.05, 0.48, 0.00) R 对上 [A&A 2024 的 OTD-O8 表](https://www.aanda.org/articles/aa/full_html/2024/06/aa46699-23/T1.html)；Cassini 11 严格轴对称；IGRF-14 和 `scene/aurora.cpp` 的地磁极一致。
- AH5 和 Cassini 11 的参考半径没有从原文核实。天王星 25,600 与 25,559 km 只差 0.16%，画面上看不出来。

放到场景的本体系（`[bodies.magnetic_field]`）：
- **天王星**：`flip_pole` 把 ULS 的北翻到 IAU 的南，`longitude_offset_deg = -50.6` 用 Melin 2025 的 `IAU new 东经 = ULS 西经 − 50.6`。AH5 的南磁极落在 235.7°W、纬度 30.2°（IAU），Melin 的南极光在 232°W、37.5°，相距约 8°。
- **海王星**：`prime_meridian_deg = [253.18, 536.3128492]` 是 IAU 2009 以前的 W，周期 16.11 h，和 NLS 同一个系统；场景的 541.14°/天不是磁场的自转率。[PDS 的 Voyager 2 轨迹说明](https://pds.nasa.gov/data/vg2-n-pra-4-summ-browse-48sec-v1.0/vg_1001/catalog/pos/dataset.cat)规定 1989-08-25 03:56:00 飞船在西经 167.7°，按场景的轨迹和这个 W 算出 168.8°。现在的相位未知。
- **木星、土星、地球、水星**：模型的坐标系就是场景的本体系。木星只用内部场时，Io 的足点和 `jupiter_footprints.txt` 相差不到 3°。水星的模型算出偏心 479 km，[Anderson et al. 2011](https://doi.org/10.1126/science.1211001) 是 484 km。
- **Ganymede**：模型的坐标系也是本体系，这是推断出来的，不是论文明写的。[Kivelson et al. 2002](https://doi.org/10.1006/icar.2002.6834) 说南半球的磁极从朝木星的子午面向后随面转 24°，Weber 的系数给出 25.7°，只有 x 朝木星、y 朝后随面时才对得上（如果按 GPhiO 的 x 沿等离子体流方向读系数，磁极就差了约 65°）。场景的 x 轴朝木星、y 轴朝后随面，`test_magnetic_field` 检查。极光在开闭磁力线的边界上（[Duling et al. 2022](https://doi.org/10.1029/2022GL101688)）：Marzok 的极光在后随面是 42–56°，在前导面只有 16–25°。
- **Ganymede 在木星磁场里**（模板的 `parent_field = true`）：偶极之外叠加木星在 Ganymede 处的磁场，当作均匀场。这个场来自 `jupiter_flux_tubes.bin`（JRM33 + CON2020，名义轨道、赤道面上每 5° 经度一个值，按 Ganymede 的经度插值）：72–117 nT，沿自转轴的分量 −70 至 −79 nT，径向分量随电流片扫过在 ±85 nT 间摆动；和 Galileo 的实测相符：[Kivelson, Khurana & Volwerk 2002](https://doi.org/10.1006/icar.2002.6834)（Icarus 157, 507）表 I 列出六次飞掠时 Ganymede 周围的背景场（由飞掠前后的数据拟合，GphiO 坐标：X 沿共转流，Y 指向木星，Z 沿自转轴），场景在这六个时刻（飞船 UT）算出的场与之相差 4、15、6、16、7、20 nT（G1、G2、G7、G8、G28、G29），均方根 13 nT，朝向或背离木星的分量符号全部一致（`test_flux_tubes`）：

  | 飞掠 | UT | Ganymede 的 System III 西经 | 实测 (X, Y, Z) nT | 模型 |
  |---|---|---:|---|---|
  | G1 | 1996-06-27 06:29:07 | 175.1° | (6, −79, −79) | (2, −79, −79) |
  | G2 | 1996-09-06 18:59:34 | 158.1° | (17, −73, −85) | (7, −65, −78) |
  | G7 | 1997-04-05 07:09:58 | 20.0° | (−3, 84, −76) | (2, 84, −77) |
  | G8 | 1997-05-07 15:56:09 | 286.5° | (−11, 11, −77) | (−13, −4, −70) |
  | G28 | 2000-05-20 10:10:10 | 351.6° | (−7, 78, −76) | (−4, 72, −74) |
  | G29 | 2000-12-28 08:58 | 239.2° | (−9, −83, −79) | (−13, −63, −76) |

  差得最多的是 G8（穿过电流片中心，径向分量变号处）和 G29。该文的图 4 就是偶极与均匀场的真空叠加，这里的画法与之相同；它采用的偶极是 719 ± 2 nT、倾角 176°（含四极项的拟合是 713 nT），和这里用的 Weber 2022 一致。木星的场在赤道上朝南，和 Ganymede 偶极在赤道上的场相反，所以两者在赤道面上 1.85–2.1 个半径处抵消（随场强）：赤道附近的线闭合（画 1.15、1.3、1.45 三圈，木星转一圈它们始终闭合），从磁纬 40°、55°、72° 出发的线离开，画到 4 个半径，末端和木星磁场的夹角在 10° 以内；木星穿过 Ganymede 的那条线（上面的磁通管，起始方向取同一个场）从它们中间接下去，事件“Ganymede's own magnetic field”把两者画在一起。开闭边界在纬度 31.8–35.7°（木星转一圈的平均，前导和后随子午线上；场纯沿偶极轴、强 113–75 nT 时的解析值是 30–35°），落在两处新月之间。
- 这是**真空叠加**：没有等离子体的流动，所以没有上下游的不对称（真实边界在后随面被推向两极、在前导面拖向赤道，正是两处新月的差别），开放的线也不向下游后掠成阿尔文翼；感应场没有加。场每变化 3% 重新追踪一次（约 1 ms），只在磁力线被画出时进行。

水星的线画到 L = 1.4，都在 [Winslow et al. 2013](https://doi.org/10.1002/jgra.50237) 的平均日下点磁层顶（1.45 个半径，太阳风动压变化时 1.35–1.55）以内。

重新生成（只用标准库）：

```bash
mkdir -p build/field
for f in uranus/ah5 neptune/gsfco8 jupiter/jrm33 saturn/cassini11 mercury/anderson2012 ganymede/weber2022dip; do curl -L -o build/field/$(basename $f).dat https://raw.githubusercontent.com/mattkjames7/libinternalfield/9a6b6e562b57a95146b7d608ee9cdf229b674a63/data/coeffs/$f.dat; done
curl -L -o build/field/igrf14coeffs.txt https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt
python -I tools/geomag/field_coeffs.py build/field assets/geomag
```

## 径向磁场图

`*.fieldmap` 是固定在本体上的磁场图，由 `tools/geomag/field_map_bake.py` 烘焙（需要 numpy）。径向图的文件头 `AXFIELDMAP1 宽 高 上限_nT` 之后是 宽×高 个有符号字节，值为 127·sign(Br)·√(min(|Br|/上限, 1))；强度图的文件头 `AXFIELDMAG1 宽 高 下限_nT 上限_nT` 之后是 127·(2t − 1)，t = (|B| − 下限)/(上限 − 下限)。行由北向南按行星中心纬度（`body.vert.hlsl` 就是这样采样的），列从 180°W 向东，在椭球面上该方向的真实位置求值。

| 文件 | 模型 | 位置 | 范围 |
|---|---|---|---|
| `mars_br.fieldmap` | [Morschhauser et al. 2014](https://doi.org/10.1002/2013JE004555)，110 阶，参考半径 3393.5 km | 椭球面上方 150 km | −710 .. 1268 nT |
| `jupiter_br.fieldmap` | JRM33，13 阶 | 1 bar 椭球面 | −13.9 .. 21.7 G |
| `jupiter_br_deep.fieldmap` | JRM33，13 阶 | 0.81 RJ 的球面，方向和云顶图相同（两张图在网格上对齐） | −79.3 .. 81.2 G |
| `moon_br.fieldmap` | [Tsunakawa et al. 2015](https://doi.org/10.1002/2014JE004785) 的球谐展开，Wieczorek 的 T2015_449 取到 240 阶，参考半径 1737.4 km | 月面上方 30 km | −25 .. 14 nT |
| `earth_b.fieldmap` | IGRF-14 的 2025.0，13 阶，总强度 |B| | 椭球面 | 22,073 .. 66,953 nT，字节范围 21,000 .. 70,000 nT |
| `earth_crust_br.fieldmap` | [Thébault et al. 2021](https://doi.org/10.5281/zenodo.5546528)，16–1050 阶，2048×1024 | 椭球面 | −6,830 .. 1,855 nT，600 nT 饱和 |
| `uranus_br.fieldmap` | AH5，4 阶，在场景的本体系里 | 1 bar 椭球面 | −14.9 .. 6.9 G |
| `neptune_br.fieldmap` | O8，3 阶，在 NLS 里 | 1 bar 椭球面 | −8.9 .. 5.9 G（单位见下） |

**月球**：系数取自 Wieczorek 在 Zenodo 的 [T2015_449.sh.gz](https://doi.org/10.5281/zenodo.3873648)（CC BY 4.0，单位是特斯拉）。取 30 km 而不是月面：直接画月面时，强异常周围有大致沿南北向的条纹（硬截断的振铃加上看来是轨道条带的噪声，在 120–240 阶之间做余弦衰减也去不掉），向上延拓到接近观测的高度后条纹消失。地表的峰值约 −279 nT，30 km 处只有 −25 nT，差得很多，展签按 30 km 写。

**地球**：画总强度而不是 Br。Br 基本只是偶极的南北两半，南大西洋异常看不出来。最弱处在 26.5°S、59.5°W。文件是 2025.0；其他年份由程序用同一份系数在同样的格子上重画（`field_strength_map`，和烘焙工具逐字节一致），所以字节范围取覆盖 1900–2030 年所有历元的 21,000–70,000 nT（一个字节 193 nT），各年的颜色可以直接比较。工具打印每个历元的范围：最弱处从 1900 年的 25,423 nT（21.5°S、36.5°W）降到 2025 年的 22,073 nT，按变化率到 2030 年是 21,909 nT（26.5°S、60.5°W），一百多年向西移了 23°；最强处从 69,435 降到 66,953 nT。地球–月球场景的两个事件分别是 2025 年和 1900 年的同一视角。

**地壳**：系数取自 Zenodo 5546528 的 `Spherical_HarmonicModel_GRL.zip`（Thébault, Hulot, Langlais & Vigneron 2021, *A Spherical Harmonic Model of Earth's Lithospheric Magnetic Field up to Degree 1050*, GRL 48, e2021GL095147，[doi:10.1029/2021GL095147](https://doi.org/10.1029/2021GL095147)；数据和论文都是 CC BY 4.0），由 CHAMP、Swarm 的磁场和几十年的近地标量异常合成。文件的 MD5 `2cf0d201b662ef5404aab58fd14df828`（Zenodo 给的），解压后的文本 SHA-256 `a9e0b40f01eefeaf54843539100e8299802748f1d96a02d72eaf90474481f8ad`，文本里夹着没有 `#` 的注释行，工具只取四列都是数的行（552,690 行）。按 Schmidt 准归一化的 Br 在 WGS84 椭球上求值，参考半径取 6371.2 km（记录里没写，是地磁模型的惯例）。面积的 50%、90%、99%、99.9% 分别在 45、147、327、589 nT 以下，峰值 −6,830 nT 在 50.5°N、37.0°E，正是 Kursk 磁异常。一个字节按平方根编码，再按峰值归一的话几乎全图都看不见，所以 600 nT 处饱和（`scale_nt = 600`，`value_nt` 在强异常处读到的是 600，不是真值）。地表直接画没有月球那样的条纹问题，海上的条带是海底扩张留下的。

**天王星、海王星**：内部场的系数和 `[magnetic_field]` 相同，所以图和线一定一致（`test_magnetic_maps` 逐格对照）。4 阶和 3 阶的模型本来就平滑，360×180 足够。天王星按 `flip_pole`、`longitude_offset_deg = -50.6` 直接烘在本体系里（模型的纬度是本体的相反数，东经是 −50.6° 减本体的东经；Br 是标量，翻转不变号）。海王星烘在 NLS 里：NLS 和场景的 IAU 自转率每天差 4.8°，不能烘死在本体上，渲染时按 `magnetic_map_turn_deg` 转动采样的 u。两张图的单位是 nT，表里写成高斯只是方便看（1 G = 10⁵ nT）。最强处都在偏心偶极靠近的一侧：天王星 −1.49 G、32.5°N、185.5°E（偏离偶极偏心方向 46°），对面的另一种符号最多 0.68 G；海王星 −0.89 G、34.5°S、NLS 东经 87.5°（偏离 32°），对面 0.59 G。各自的偶极单独只能给出 0.46、0.28 G。

**火星模型的许可**：系数取自 Mark Wieczorek 在 Zenodo 的 [Morschhauser2014.txt.gz](https://doi.org/10.5281/zenodo.3876495)（CC BY 4.0，与论文附表 ts01 相同）。没有用分辨率更高的 Langlais 2019 和 Gao 2021：两篇论文都是 CC BY-NC（-ND），和项目的 GPL 不相容。Zenodo 上也有 Langlais 模型的 CC BY 版本，但它的许可和论文的 NC-ND 对不上，所以不用。

取 150 km 而不是地表：这里接近观测的高度（Mars Global Surveyor 的低轨和大气制动段），往下延拓到地表时小尺度放大，整张图布满高频细纹。

核对：
- 木星两张图和 `field_jupiter_jrm33.txt` 逐格一致（`test_magnetic_maps`，到一个字节的量化）。
- Great Blue Spot 在云顶约 −12.8 G、6.5°S、84.5°W；[Moore et al. 2018 的会议摘要](https://meetingorganizer.copernicus.org/EPSC2018/EPSC2018-556.pdf)写它在赤道附近、System III 约 90°。在扁椭球面上它不是全球最小值（南极附近 −13.9 G），在 1 RJ 的球面上才是。
- 火星的峰值在 48.2°S、174.8°E，属于 [Acuña et al. 1999](https://doi.org/10.1126/science.284.5415.790) 说的 Terra Cimmeria 和 Terra Sirenum；Hellas、Argyre 内的场不到峰值的 10%。
- 月球的峰值在 20.2°S、236.8°E，紧邻 Gerasimovich 环形山（22.9°S、122.6°W）；Reiner Gamma（7.5°N、59°W）处有明显的异常，Mare Crisium 处很弱；远月面的场比近月面强。
- 地球图和 `field_earth_igrf14.txt` 逐格一致（到一个字节的量化）。
- 地壳图的 Kursk（51°N、37°E）周围 ±2° 的均方根 393 nT，Bangui（6°N、18°E）周围 ±6° 是 162 nT，南太平洋（30°S、160°W）是 49 nT。

重新生成（先按上一节取得 `field_*.txt`）：

```bash
mkdir -p build/field_map
curl -L -o build/field_map/Morschhauser2014.txt.gz https://zenodo.org/api/records/3876495/files/Morschhauser2014.txt.gz/content
gunzip -k build/field_map/Morschhauser2014.txt.gz
curl -L -o build/field_map/T2015_449.sh.gz https://zenodo.org/api/records/3873648/files/T2015_449.sh.gz/content
curl -L -o build/field_map/crust.zip https://zenodo.org/api/records/5546528/files/Spherical_HarmonicModel_GRL.zip/content
unzip -o build/field_map/crust.zip -d build/field_map
python -E -P tools/geomag/field_map_bake.py build/field_map/Morschhauser2014.txt build/field_map/T2015_449.sh.gz build/field_map/Spherical_HarmonicModel_GRL.txt assets/geomag
```

## 磁层顶与弓激波

`src/scene/magnetosphere.cpp` 用经验模型把磁层顶（蓝）和弓激波（橙）画成线框曲面，系数都是论文里的数值；数据文件只有地球的逐小时太阳风 `solar_wind.bin`（见下面的“逐小时太阳风”）和海王星的边界剖面表 `neptune_boundaries.txt`。行星有地球、木星、土星、天王星、水星和海王星。面板的 `Magnetosphere` 选 `Always` 时按模板里的典型太阳风画；`Event setting` 下由事件的 `magnetosphere` 列表给出各自的太阳风，`faint = true` 的一项是用来对比的淡色面（不画激波）。坐标轴：x 指向太阳（地球、水星另加光行差），z 是轨道面北（木星按 Joy 取自转轴），y 指向昏侧。

| 天体 | 磁层顶 | 弓激波 | 典型太阳风 |
|---|---|---|---|
| 地球 | Shue et al. 1998, JGR 103, 17691, doi:10.1029/98JA01103：`r0 = (10.22 + 1.29 tanh[0.184 (Bz + 8.14)]) Dp^(−1/6.6)`，`α = (0.58 − 0.007 Bz)(1 + 0.024 ln Dp)`，`r = r0 (2 / (1 + cos θ))^α`，半径 6371.2 km | Fairfield 1971（JGR 76, 6700）的平均二次曲线，取 SPEDAS `bshock_2`、`mpause_2` 的系数：鼻点 14.33、磁层顶 10.78 RE，晨昏线处 27.2、−22.9 RE。这里取对称的、焦点在地心的圆锥曲线：鼻点为磁层顶的 1.33 倍，偏心率 0.75（晨昏线半径 25.1 RE） | 2 nPa、Bz = 0，得 10.25 RE |
| 木星 | Joy et al. 2002, JGR 107, 1309, doi:10.1029/2001JA009146：`z² = A + B x + C x² + D y + E y² + F x y`，长度以 120 RJ 为单位，A、B、C 为 `c0 + c1 P^(−1/4)`，D、E、F 为 `c0 + c1 P`，z 沿自转轴 | 同一篇的另一组系数 | 0.071 nPa（Rutala et al. 2025 的中位数），得 83 RJ |
| 天王星 | Jasinski et al. 2024, Nature Astronomy 9, 66, doi:10.1038/s41550-024-02389-3（CC BY 4.0）的方法：`R = (B0² / (2 μ0 P))^(1/5.7)`，B0 = 23,000 nT，形状同 Shue、张开参数 0.6（Cao & Paty 的 MHD 模拟给 0.45–0.75），半径 25,559 km | 抛物线（偏心率 1，假设），鼻点为磁层顶的 1.275 倍（Farris & Russell 1994 的 `1 + 1.1 (γ − 1)/(γ + 1)`，高马赫数极限，γ = 5/3） | 平均 22.2 RU |

| 土星 | Kanani et al. 2010, JGR 115, A06207, doi:10.1029/2009JA014262：`r0 = 10.3 P^(−0.2)` RS，张开参数 `0.73 + 0.4 P`，形状同 Shue，半径 60,268 km | Went et al. 2011, JGR 116, A07202, doi:10.1029/2010JA016349：焦点在土星的圆锥曲线，偏心率 0.84，鼻点 `15 P^(−1/5.4)` RS | 0.02 nPa（取的整数值），得 22.5、31.0 RS |
| 水星 | Winslow et al. 2013, JGR Space Physics 118, 2213, doi:10.1002/jgra.50237：形状同 Shue，`Rss = 2.15 P^(−1/6.75)` RM，张开参数 0.5，半径 2440 km。中心在偶极中心（北移 0.196 RM），轴线带光行差（400 km/s） | 同一篇表 2（马赫数修正后的拟合）：焦点在轴上 0.5 RM，p = 2.90 RM，偏心率 1.02，鼻点 1.96 RM | 14.3 nPa（平均冲压） |
| 海王星 | Mejnertsen et al. 2016 的 MHD 模拟（图 4 的曲线，`neptune_boundaries.txt`），随偶极转动和季节变化。乘 1.007 以通过旅行者 2 号 1989-08-24 的穿越点（26.5 RN，Ness et al. 1989, Science 246, 1473，“tentatively”），半径 24,765 km | 同一模拟，乘 1.088 以通过 34.9 RN 的穿越点（同上） | 飞掠时的弱太阳风 |

**核对**（`test_magnetospheres`）：

- **Joy 的系数**是原文图 3 上印的 24 个数（式 2，长度以 RJ/120 为单位），2026-11-10 逐个对过；Rutala et al. 2025（arXiv:2502.09186）的公开代码 [JovianBoundaries](https://github.com/mjrutala/JovianBoundaries) 的 `JoyBoundaryCoords.py` 里是同一组数。原文的两个最可能的磁层顶是 63 RJ（0.306 nPa）和 92 RJ（0.039 nPa），这里得 63.1、92.0；激波平均 84 RJ（σ = 16），中位压强 0.258 nPa 的面在 78.2 RJ。Rutala 的表给出 Joy 模型在 0.02、0.071、0.22 nPa 的日下点距离 101、83、68 RJ（磁层顶）和 132、108、82 RJ（激波），这里得 101.7、83.4、67.6 和 131.8、107.8、82.1。原文说拟合只用到反日方向 250 RJ 以内（激波在行星后约 100 RJ 就出了模拟区），事件的磁尾画到 250 RJ，激波到 87.5 RJ。Kurth et al. 2002（Nature 415, 991）说 Cassini 在 2000-12-28 04:19 于 140.2 RJ 第一次遇到激波，紧挨 Joy 的 0.02 nPa 激波面：场景里那一刻 Cassini 距木星 140.2 RJ，是该面距离的 0.93。Kurth 写的 0.02 nPa 来自付印前的稿子，正式发表的表 2 里“10%”激波面是 0.063 nPa。坐标按原文：x 指向太阳，z 沿自转轴（取垂直于 x 的分量），y 指向昏侧；Kurth 图 1 的 z 垂直于轨道面，差 3°。
- **Cassini 与 Galileo（2001-01-10）**：Kurth 2002 给出两者离开磁层的时刻 20:35、20:52（探测器时间），x 相差 100 RJ 以上；Galileo 的穿越点在 Joy 的“25%”面上，Cassini 的在“10%”面之外，没有一个稳态面能同时通过两点，所以磁层正在收缩，压强锋面从 Galileo 到 Cassini 约 5 小时。场景里两点是 (40.7, 93.6, 1.1) 和 (−69.1, 192.4, −2.3) RJ；让 Joy 的面通过它们的压强是 0.0452、0.0097 nPa（日下点 89.9、112.2 RJ），事件就画这两个面，后者为淡色。Joy 表 2 的“25%”“10%”面分别是 0.048、0.030 nPa，和 Kurth 的说法一致：Galileo 的 0.0452 nPa 在 25% 面上（日下点差 1 RJ），Cassini 的在 10% 面之外。这是两个稳态面，论文图 1 中连接它们的弯折没有画，面也不随时间移动。
- **天王星**：论文的三组压强和日下点距离（0.0184 nPa → 17.3 RU，0.031 → 15.8，0.00043 → 33.5）都能复现。论文给出旅行者 2 号 10:07 UT 穿越磁层顶的位置 (15.784, 8.684, −0.632) RU，场景的轨迹在那一刻是 (15.818, 8.695, −0.633)，是 17.3 RU 模型面距离的 1.005。激波在 07:28（Treumann & Jaroschek 2008, arXiv:0808.1701 引 23.5 RU），场景里探测器距行星 23.67 RU，假设的抛物线在该方向是 23.15 RU。事件把 17.3 RU（飞掠时，出现的时间不到 5%）和 22.2 RU（平均）画在一起；论文估计一周前约 28 RU。
- **地球磁尾**：Slavin et al. 1985（JGR 90, 10875，摘要）用 ISEE-3 的 756 次磁层顶穿越得出 |X| = 130–225 RE 处磁尾直径 60 ± 5 RE，约 120 ± 10 RE 以外不再张开。Shue 的公式在那里没有约束，会一直张开，所以模板用 `flare_end_radii = 120`：更远处是半径等于该处 Shue 值（29.6 RE）的圆柱。轴线按地球公转速度和 430 km/s 的太阳风偏转约 4°（Fairfield 1971 的假设）。ISEE-3 在 1982-10 至 1983-12 位于 60 RE 以外尾向的时间里，约 55% 在这个平均面内，最远 236.2 RE 处离轴 16 RE；Hones et al. 1984（GRL 11, 5，摘要）说它在尾内时约三分之一在等离子体片里，磁尾摆动的幅度和自身宽度相当。
- **2024-05-10/11**：事件按下面的逐小时记录画，从 5 月 10 日 13:00 开始（每秒 15 分钟）。流压在 17 UT 那一小时从约 1 nPa 升到 17.75 nPa，20 UT 达 55.89 nPa，11 日 02–03 UT 为 48.76 nPa、Bz（GSM）−20.5 nT、速度 738 km/s，这一小时 Shue 给出 4.97 RE。事件开场时磁层顶在 11.3 RE，最小 4.86 RE（10 日 20:30 UT），此后 36 小时里有 19 小时在地球同步轨道（6.6 RE）以内。Fu et al. 2025（GRL 52, e2024GL114040）由地面磁强计、THEMIS-A 和 GOES-16/18 得出 5 月 10 日压缩到 5.12–5.20 RE，可作旁证。激波按磁层顶等比例缩小，没有计入风暴时较低的马赫数。

**逐小时太阳风**（`solar_wind.bin`，地球模板的 `solar_wind`）：NASA GSFC SPDF 的 OMNI2 小时值，[low_res_omni](https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/) 的 `omni2_1963.dat` 至 `omni2_2026.dat`（King & Papitashvili 2005, JGR 110, A02104；美国政府作品，SPDF 要求致谢 OMNIWeb 服务和 OMNI 数据），2026-11-11 下载。数据已由各探测器的位置平移到弓激波鼻点，每行一小时（UT）。`tools/geomag/solar_wind_bake.py` 只保留三个量，每小时 4 字节：流压（第 29 个字，0.01 nPa 一档）、GSM 的 Bz（第 17 个字，0.5 nT 一档，限在 ±63.5 nT，记录里没有超出的小时）、流速（第 25 个字，`150 + 8n` km/s，只用于约 4° 的光行差）。格式写在工具开头。

```
mkdir -p build/omni
for y in $(seq 1963 2026); do curl -s --fail -o build/omni/omni2_$y.dat https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_$y.dat; done
python -I tools/geomag/solar_wind_bake.py build/omni assets/geomag/solar_wind.bin
```

当前文件 559,418 小时，1963-01-01 至 2026-10-26 02:00 UT，SHA-256 `595097c2a3df7f6f3fed0fb57d9e521b06aaa8187f25b385c1b9e6ea5aba9cdb`；流压中位数 1.79 nPa、速度 414 km/s，最大的小时流压 76.06 nPa。有流压的小时占比：1960 年代 36%、1970 年代 66%、1980 年代 52%、1990 年代 68%、2000 年以后 96–99%（工具按年代打印）。更新时重新下载最近几年并重跑工具：最近的值可能修订。

运行时（`SolarWindRecord::mean`、`Scene::magnetosphere_wind`）：取以当前时刻为中心的时间窗内的平均，窗宽为屏幕上 1.5 秒对应的时间、至少 1 小时，所以正常速度下就是相邻两个小时均值之间的线性内插，高倍速时面缓慢涨缩而不闪烁。缺测的小时和记录以外按典型太阳风（2 nPa、Bz = 0、430 km/s）计入平均，三个量各自判断缺测。1979 年有几个小时的流压记为 0.00，平均后的流压不低于 0.02 nPa。面板的 `Magnetosphere` 选 `Always` 时地球按记录画；事件里写 `recorded = true`。其他行星没有这样的记录：从 1 au 外推到木星以远要用传播模型（如 Tao 2005、mSWiM），只在冲日前后可靠，这里没有做。`test_magnetospheres` 检查上面的小时值、内插、缺测（1983-10-27 整天没有数据）和记录两端。

**日球层**（太阳系场景的太阳，`model = "heliosphere"`，长度用 au）：同一套线框画日球层顶（蓝）和终端激波（橙，在日球层顶以内）。两者没有随太阳风变化的公式，大小直接给定。

- **方向**：鼻点朝星际来流。McComas et al. 2015（ApJS 220, 22，IBEX 六年的星际氦，“working values”）：流速 25.4 km/s，流向黄经 75.7°、黄纬 −5.1°，所以鼻点在 255.7°、+5.1°。z 轴取黄道北。
- **形状**：日球层顶用 Parker 1961（ApJ 134, 20）的不可压缩解，即点源加均匀来流的分界流面 `r = R / cos(θ/2)`，等于 Shue 形式取张开参数 0.5，远尾半径趋于 2R。终端激波取球面，这也是该解的极限情形。
- **大小**：两个事件各让日球层顶通过探测器的穿越点，激波取该探测器当年穿越的距离。场景里的旅行者号在四个日期的日心距与发表值都差不到 0.03 au（`test_magnetospheres`）：

| 探测器 | 终端激波 | 日球层顶 | 离鼻点 | Parker 的 R |
|---|---|---|---:|---:|
| Voyager 1 | 2004-12-16，94.0 au（Stone et al. 2005, Science 309, 2017） | 2012-08-25，121.6 au（Stone et al. 2013, Science 341, 150） | 29.9° | 117.5 au |
| Voyager 2 | 2007-08-30，83.7 au（Stone et al. 2008, Nature 454, 71） | 2018-11-05，119.0 au（Stone et al. 2019, Nature Astronomy 3, 1013） | 52.6° | 106.7 au |

- **这只是示意形状**：两次穿越放不进同一个 Parker 面（Voyager 2 的 R 是 Voyager 1 的 0.91 倍），所以每个事件只画自己的面。真实的日球层被星际磁场压得南北不对称，边界随太阳活动周进退几个 au，两次穿越又相隔六年；日球尾是长是短仍有争议，这里画到 250 au 为止。没有画弓激波或弓波（IBEX 的流速下不一定存在），也没有画日球层鞘里的结构。不选事件时（`Always`）用 120 au 的整数值。四个穿越日期和距离、Parker 的公式是凭已知文献写的，本次只下载核对了 McComas 2015 的原文；穿越距离由场景的星历独立复算。

**土星、水星**（`model = "power_law"`）**和海王星**（`"tabulated"`）。Kanani 2010、Went 2011 和 Winslow 2013 的原文在 2026-11-10 对过；海王星的几篇没有读到原文，数值来自摘要和转引：

- **土星**：Kanani 表 1 的新模型是 a1 = 10.3 ± 1.7、a2 = 0.20 ± 0.03、a3 = 0.73 ± 0.07、a4 = 0.4 ± 0.5（`R0 = a1 P^(−a2)`，`K = a3 + a4 P`），坐标是 KSM（x 指向太阳）；a4 的误差比数值大，张开程度随压强的变化并不显著。Went 的式 5 是 `r = (1 + ε) c1 P^(−1/c2) / (1 + ε cos θ)`，焦点固定在土星中心，表 1 的新模型为 ε = 0.84 ± 0.06、c1 = 15 ± 1 RS、c2 = 5.4 ± 0.5（203 个平均后的穿越，径向均方根偏差 4.6 RS），RS = 60,268 km。它在带光行差的 KSM 里拟合（原文说土星处不到 2°），Kanani 用的是不带光行差的 KSM；这里两个面都不加光行差。注意同年另有一篇 Went et al. 2011（JGR 116, A04224，木星与土星外磁层比较），不是这个模型的出处。Sulaiman et al. 2016（arXiv:1605.04844）说 Went 的模型 `R_SN ∝ P^(−1/5.4)`，转引里另有“约 0.04 nPa 时日下点 27 RS”，这里得 27.2。事件“Titan at the edge of Saturn's magnetosphere”画典型面（22.5 RS）和淡色的 0.05 nPa 面（18.8 RS）：那天（2013-12-01）土卫六离日下方向 21°、距土星 20.0 RS，在前者之内、后者之外、仍在激波之内。日期只用来把土卫六放到昼侧；当天实际的太阳风没有查。
- **水星**：原文给出 `Rss = (2.15 ± 0.10) P^[(−1/6.75) ± 0.024]`（RM、nPa），14.3 nPa 时 1.45 RM、张开参数 0.5，8.8–21.6 nPa 时 1.55–1.35 RM。激波是式 3 的圆锥曲线 `sqrt((X − X0)² + ρ²) = p ε / (1 + ε cos θ)`，表 2 有两组拟合：穿越中点的 (X0, p, ε) = (0.5, 2.75, 1.04)，鼻点 1.90 RM；按马赫数修正后的概率拟合 (0.5, 2.90, 1.02)，鼻点 1.96 RM。这里用后一组。坐标的中心在偶极上（Zd = 0.196 RM），x 逆着相对水星的太阳风，平均光行差约 7°（3.5–10.2°）；场景按 400 km/s 得 5.5–8.4°。原文说激波位置主要随阿尔文马赫数变（4.12–11.8 对应 2.29–1.89 RM），形状不变，这里激波只按磁层顶等比例移动。场景里 MESSENGER 在事件当天（2013-01-01，三圈）进出两个面各六次。
- **海王星**：Lepping et al. 1992 给出激波 14:38 UT、磁层顶 18:00 UT；场景里的旅行者 2 号那两刻距行星 34.89、26.45 RN，离日下方向 13.7°、13.4°，和 Ness 的 34.9、26.5 RN 一致。两个面各自定标到通过对应的穿越点（见下）。Cairns 1991 由远程探测得出激波鼻点在穿越前从 38.5 移到 34.2 RN，并认为位置受自转相位控制；磁轴倾斜 47°，磁层每 16 小时在“极对日”和“类地球”之间变化。
- **海王星的面随自转**（`model = "tabulated"`）：Mejnertsen et al. 2016（JGR Space Physics 121, 7497, doi:10.1002/2015JA022272，CC BY 4.0）用 Gorgon 做了飞掠时的全球 MHD 模拟，偶极按 O8 的强度绕自转轴进动，太阳风取旅行者 2 号入境前的实测（比平均低约 75%）。图 4 是偶极朝日（“极对日”）和半圈之后（“类地球”）两个位形下，激波和磁层顶在 y = 0、z = 0 两个平面里的截线，一共八条曲线。PDF 里是矢量图，`tools/geomag/neptune_boundaries_fig4.py`（需要 PyMuPDF，文章不入库）取路径坐标，线宽相当于 0.4 RN：

  | | 类地球 | 极对日 |
  |---|---:|---:|
  | 激波，轴上 | 31.3 RN | 32.8 RN |
  | 激波，z = ±20 RN 处的 x | 24.8、27.1 RN | 30.0、30.8 RN |
  | 激波，y = ±20 RN 处的 x | 26.5、28.0 RN | 30.2、30.5 RN |
  | 磁层顶，轴上 | 26.3 RN | 23.2–23.5 RN |
  | 磁层顶，最靠日的一点 | 26.5、26.8 RN（轴上） | 27.0 RN（北 9.9 RN）、25.0 RN（昏侧 9.7 RN） |

  激波鼻点在两个位形之间差 1.5 RN；正文第 3.2 节写的“类地球时约近 4 RN”在图上出现在偏轴处（上表第二、三行差 3–5 RN）。磁层顶在极对日时沿轴线内凹约 3 RN（日侧的极尖），最靠日的点移到偏轴约 10 RN 处。模拟的激波比实测近约 2 RN、磁层顶“在同一位置”（第 3.5 节）。

  **烘焙**：每条曲线在日向上切成两半，各自成为“到行星的距离随离日向夹角”的剖面，所以每个面、每个位形有四条：朝 +y（昏侧，逆着公转）、+z（轨道面北）、−y、−z，和场景的磁层轴相同（图的 NSO 坐标）。按 1° 分箱取描边两条边的中点，虚线的空隙沿曲线补上，做 σ = 2° 的平滑，每 2.5° 存一个值；四条剖面在日向上取同一个值（两个平面的平均，修正在 20° 内淡出）。图框只到横向 60 RN，所以激波存到 80°、磁层顶到 115°，更远处按 `r (1 + cos θ)^α = 常数` 外推，α 取末端 10° 的斜率（限在 0.3–1.2，实际 0.43–0.96）。文件 `AXBOUNDARY1`，格式写在工具开头。工具另外直接从轮廓打印各曲线在横向 ±10、±20 RN 处的 x，28 个值和表（经过分箱、平滑和下面的绕向相接）相差不到 0.6 RN（`test_magnetospheres`）。

  ```
  python -I tools/geomag/neptune_boundaries_fig4.py <文章的 PDF> assets/geomag/neptune_boundaries.txt
  ```

  **运行时**（`BoundaryTable`）用下列近似把论文的两个位形、两个平面扩展到各个方位和季节。

  - 四个方位之间取 `a + b cos φ + c sin φ + d cos 2φ`。图 4 的两个相位都取在偶极位于日向与自转轴所张的平面时。磁层顶以这个平面为零转角，按 O8 偶极投影到垂直日向平面上的方位转动，凹陷的非轴对称部分随之转动。弓激波的方位用 NSO 坐标，论文第 3.2 节指出其不对称也可能受行星际磁场影响。
  - 两个面的南北差由 `c sin φ` 给出。日下纬度低于 −5° 时用飞掠季节的剖面，高于 +5° 时用关于 `z = 0` 的镜像，中间用三次平滑插值。这个镜像是对另一季节的近似，论文第 4 节和结论提出还需模拟不同季节与行星际磁场。
  - `dipole_phase` 用偶极轴与日向夹角的余弦混合两个位形，偶极最接近日向时为 1（极对日），半圈之后为 0（类地球）。春分和秋分附近，两个磁极先后达到相近的日向夹角，混合权重逐渐收向 0.5。飞掠时的 O8 相位在激波穿越（14:38）为 0.35，在磁层顶穿越（18:00）为 0.92，在最近点（次日 03:56）为 0.01，和论文第 3.5 节的位形顺序一致。自转快到看不清时（每秒 0.25–1 圈）两个面渐变到中间相位，磁层顶同时渐变到各离日角上的方位平均。

  按穿越时的相位和方向取样，定标前的激波和磁层顶分别在 32.07、26.25 RN，实测为 34.89、26.45 RN。乘 `shock_scale = 1.0877` 和 `standoff_scale = 1.0074` 后，两个面通过各自的穿越点。飞掠季节一圈内沿日向的距离分别为 34.1–35.7 RN 和 23.3–26.5 RN。Cairns 由远程探测得出的激波鼻点范围为 34.2–38.5 RN。

线框：12 条经线加若干圈，磁尾画到 `tail_radii`，激波只画到它的 0.35 倍。这些都是平均位置：真实边界随太阳风起伏，磁层顶上还有波动和重联，激波的形状也随马赫数变。Joy 的面来自 MHD 模拟（IMF 为零的几次），鼻部权重大，磁尾远处会闭合，不能当真。

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
- Height: Schneider et al. (2018), *Global aurora on Mars during the September
  2017 space weather event*, GRL 45, 7391, https://doi.org/10.1029/2018GL077772,
  Fig. 2 (CO2+ UVD limb curve, orbit 5730; read from the plot). That curve is
  intensity edge-on, not emission at that height, so
  `tools/geomag/mars_aurora_profile.py` (`python -I`, standard library) fits the
  emission through the line-of-sight integral: `src/scene/mars_aurora.hpp` has a
  horizontally uniform layer with a Gaussian per side, peak 63 km, sigma 12 km
  below and 17 km above, cut at 20 and 140 km. Edge-on it peaks at 55 km (observed
  57.5) with an rms of 0.3 kR on the 20 kR curve. The layer below 23 km is not
  constrained, and the curve is not the green line: Knutsen shaped their model on
  the 297.2 nm line of the same event, which the paper does not plot; UVD is the
  paper's profile of the excitation, and O(1S) may be quenched lower down (their
  Cameron bands are). A numerical tangent-ray integral normalizes the peak to
  2.1 kR for radius 3396.19 km. The renderer approximates constant physical
  radial altitude on Mars' ellipsoid.
- September 2017 (`MarsAurora::limb_rayleigh_at`): Schneider Fig. 3, Cameron bands
  at the limb, limb scans and apoapse images read from the plot. Dates check
  against Fig. 1 (orbit 5731's marker falls at 13 Sep 05:34). The sites differ by
  up to 5 times within hours, so the timeline keeps the brightest of each 6 hours
  and raises what is under the 0.04 kR noise floor to it, joined in log brightness.
  Its peak is the limb scan of orbit 5730. Knutsen estimate 130 R of 297.2 nm for
  2024, "10 times fainter" than this event, so the green line at this peak is
  taken as 10 x 2.1 = 21 kR. This is inferred from the ultraviolet with the same
  line ratio as 2024; no visible observation exists, and the Cameron bands may not
  track the green line through the event. Over the day before and after the record
  it joins the reference level, so manual viewing has no jump.
- 557.7 nm linear RGB per kR is the existing `aurora_colors.py` spectral conversion.
  Display gain is independently set to 3e7. Daylight masking fades between
  solar elevations -12 and 0 degrees; it represents visibility, not excitation.
- The event time is rounded to 06:47 TDB; its orbital camera is illustrative.
  Outside September 2017 the reference emission stays constant while enabled.
  There is no storm time series for 2024, terrestrial Kp driver, crustal-field
  structure (including the discrete wisps of orbit 5738 in Schneider Fig. 4),
  red-line upper-limit emission, dust transfer, multiple scattering or surface
  illumination.
- No paper images or datasets are redistributed; the implementation uses the
  cited scalar measurements and an independently written illustrative profile.
