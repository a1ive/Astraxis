# White-light solar corona

The corona in `solar_eclipse.toml` and `parker.toml` is computed by this project's
own code. No solar photograph, magnetic map, paper figure or third-party shader
is included. The streamers and plumes are illustrative geometry.

## Scattering and density

- Electron density: Leblanc, Dulk & Bougeret (1998),
  `n_e(r) = 3.3e5/r^2 + 4.1e6/r^4 + 8.0e7/r^6` cm⁻³, with `r` in solar radii.
  The equation is reproduced in Section 3.1, page 8 of
  [Reginald, Newmark & Rastaetter (2023), Solar Physics 298:73](https://ntrs.nasa.gov/api/citations/20230015141/downloads/paper1.pdf?attachment=true),
  DOI 10.1007/s11207-023-02160-3. The cited comparison covers 1.5–3.5 solar radii;
  extending this profile down to 1 and out to 32 is a modeling choice here.
- Thomson scattering follows the finite solar disk treatment described in
  [Inhester (2016), Thomson Scattering in the Solar Corona](https://arxiv.org/abs/1512.00651).
  `src/scene/corona.cpp` numerically integrates the incident disk into coefficients
  `A(r) + B(r) cos²(chi)`, in mean solar brightness per solar radius of path.
  It uses 256 angular samples, a linear limb-darkening coefficient of 0.6 and
  the [CODATA 2022](https://physics.nist.gov/cuu/pdf/wall_2022.pdf) Thomson cross
  section, 6.6524587051e-29 m². The solar radius is 695,700 km.
- The GPU interpolates 128 radial coefficient pairs and integrates each ray
  with 96 samples in `z = b sinh(t)`. The geometric contrast multiplies the
  radial density; it is not normalized to preserve its spherical average.
  The volume fades from 24 to 32 solar radii. Half-resolution integration is
  recomputed at full resolution near opaque silhouettes.
- Solar axes use the IAU pole (286.13°, 63.87°) and prime meridian
  `84.176° + 14.1844° d`, as in
  [NAIF pck00011.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc).

These are independently implemented scientific formulas and published constants;
the code is GPL-3.0-or-later. The papers retain their authors' and publishers'
rights. Their availability through NASA or arXiv does not make their text or
figures public domain.

## Cycle and appearance

The geometric activity parameter is smoothly interpolated through three dates:
April 2014 maximum, December 2019 minimum, and the October 15, 2024 announcement
that the Sun had reached its maximum phase. Sources:
[NASA, Solar Cycle 25 (2020)](https://www.nasa.gov/missions/sdo/what-will-solar-cycle-25-look-like/)
and [NASA/NOAA maximum-phase announcement (2024)](https://science.nasa.gov/science-research/heliophysics/nasa-noaa-sun-reaches-maximum-phase-in-11-year-solar-cycle/).
The 2024 announcement is not an exact peak date. Values outside the anchor span
are held at their endpoints; this is not an eleven-year cycle forecast.
No SILSO sunspot series or GFZ sunspot column is bundled or read.

Low activity emphasizes a wrinkled equatorial sheet and polar plumes. Higher
activity adds a tilted sheet and stronger warping. The angles, widths, density
contrasts and fine striations are project design choices, not reconstructions
of the 2017 or 2024 corona. No PFSS extrapolation or CME is simulated.

The displayed white-light brightness is enhanced: default `display_gain=8e5`,
`display_gamma=0.65`. Ordinary eclipse views reduce the gain according to the
uncovered, limb-darkened solar disk; the Parker event explicitly enables full
enhancement. This contrast compression is separate from the scattering model
and is not a calibrated exposure or naked-eye brightness prediction. The tint
is the project's 5772 K solar continuum color. F-corona, EUV emission,
polarization and instrument response are omitted. The Parker view is not a
WISPR image reconstruction.

## Eclipse geometry

`tools/corona/eclipse_bake.py` generates the Earth and Moon tables and prints
the two observer offsets used by the scene. See
[ephemeris sources](../ephem/SOURCES.md#solar-eclipses-2017-and-2024) for sampling,
verification errors and data attribution. These viewpoints are 1000 km above
a spherical Earth, fixed in geocentric ICRF at each event, not ground sites or
propagated satellite orbits. The events use UTC; the ephemerides use TDB.
The scene opens the 2024 event; the event menu also offers 2017 and the
2022 lunar viewpoint described below.

## 地球边缘的大气折射像

运行 `astraxis --scene solar_eclipse --event 3`，从月球近地侧上方 200 km 看 2022-11-08 的地球遮日。相机在 10:59 UTC 的位置由 `tools/corona/eclipse_bake.py` 计算，随后保持月心 ICRF 偏移；它是设计的观测点。星历取样和独立核对见 [星历来源](../ephem/SOURCES.md#solar-eclipses-2017-and-2024)。

物理参数来自地球模板 `assets/bodies/earth.toml`：气体标高 8 km，550 nm 的垂直 Rayleigh 光学厚度 0.1085，按波长的负四次方换算到项目的 RGB 代表波长。气溶胶消光厚度 0.005328、标高 1.2 km，与既有大气散射采用同一组参数。直射光沿弯曲路径的衰减采用完整消光厚度，不使用单次散射显示中的前向峰折减。

底部折射率采用 Ciddor 的标准干空气（15 °C、101325 Pa、450 ppm CO2）公式：

```text
sigma = 1 / 0.55                 # vacuum wavenumber, um^-1
n0 - 1 = 1e-8 * (5792105 / (238.0185 - sigma^2)
                 + 167917 / (57.362 - sigma^2))
       = 0.000277837635429306
```

系数核对 [NIST Appendix A, A22](https://emtoolbox.nist.gov/Wavelength/EquationA22.asp)，计算式见 [A32](https://emtoolbox.nist.gov/Wavelength/EquationA32.asp)。密度和折射率都采用指数剖面：`n(r)=1+(n0-1) exp(-(r-R)/H)`。这是本项目选择的晴空简化，未使用当日温压观测。

`src/scene/refraction.cpp` 从球对称介质的 Snell 不变量 `b=n(r) r sin(psi)` 计算整条光线，其中 `psi` 是光线与径向的夹角。在切点 `r_t`，`b=n(r_t) r_t`。由此得到：

```text
bending = -2 b integral[r_t..infinity] ((d ln n / dr) / sqrt((n r)^2 - b^2)) dr
optical_depth = 2 integral[r_t..infinity] (extinction(r) n r / sqrt((n r)^2 - b^2)) dr
transmission = exp(-optical_depth)
```

积分用 `r=r_t+x^2` 去掉切点奇性，以 512 个中点积分到切点以上 32 个最大标高。加载器拒绝任意高度的 `n(r)r` 非递增、表顶低于最低出射瞄准距离或非有限参数。128 个表项均匀分布在渐近瞄准距离 `b` 上，分别存 RGB 透射率和折射角；天体半径、气体标高与折射率决定表的实际范围。没有逐帧重算积分。

着色器把相机射线朝天体中心弯折，再检查它是否落在有限日面上；日面采用线性临边昏暗系数 0.6。对观察方向积分辐亮度和沿曲线路径使用 Beer-Lambert 衰减的物理依据，见 [García Muñoz et al. (2012), Section II](https://arxiv.org/html/1206.4344)。这里独立实现了球对称直接光的积分，没有移植论文的代码、图像或观测光谱。论文的多次散射和完整大气光谱模型不在本实现内。

环常常窄于一个像素，GPU 按像素的径向覆盖宽度取 24 个样本并平均，以保留环的覆盖率。太阳色取现有 5772 K 连续谱，基础辐亮度使用恒星表面的 6 个线性 HDR 单位。`refraction_gain=64` 是日食场景的显示取值，其他场景默认 1；它不改变积分出的折射角或透射率。关闭大气显示也会关闭红环。

几何仍有近似：球对称表应用在现有椭球单位坐标中，折射使用 550 nm 的同一个角度，未计色散、臭氧、云层、地形或天气。只补充直射日面外的折射像，直射日面和既有单次散射仍沿直线绘制；不计算月面收到的折射照明或星空畸变。径向像素平均未做完整二维像素积分。图像适合展示地球边缘的红环，不能用来预测实际曝光或当日环上的细节。

新增代码由本项目独立编写，按 GPL-3.0-or-later 发布。引用的是公式和数值事实；没有新增第三方纹理或代码依赖。测试用另一条数值路径（笛卡尔坐标中的 eikonal 方程，按弧长 RK4 积分）核对折射、守恒的瞄准距离和三个通道的消光，并检查月球视点的本影位置。
