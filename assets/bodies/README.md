# 复用天体属性

场景的 `[[bodies]]` 可以用 `template` 引用一个天体模板。模板文件顶层直接写属性，字段名与场景内的天体相同。

例如，在 `assets/scenes/example.toml` 中创建一颗使用地球属性的根天体，并调整大气显示增益：

```toml
name = "Earth example"

[[bodies]]
name = "Earth"
template = "bodies/earth.toml"

[bodies.atmosphere]
gain = 0.5
```

`template` 和模板中的贴图、形状、自转表等资源路径都相对于 `assets/`。模板路径使用 `/`，绝对路径、盘符和 `..` 会报错。模板可以引用普通资源文件；模板内再写 `template` 或 `unset` 会报错。每个天体最多引用一个模板。

场景也可以直接写完整的天体属性。`name` 由场景提供，父天体仍需在子天体之前定义。

## 覆盖与清除

加载顺序是读取模板、执行场景中的 `unset`、应用场景属性，再按合并后的内容构造天体并校验。未指定的属性采用加载器默认值。

- 标量由场景值覆盖，包括 `0` 和空字符串。
- 数组整体替换。例如 `radii_km = [100.0]` 会替换模板的整组半径。
- 子表逐字段合并。空表沿用模板中的所有字段。
- `unset` 删除模板中的顶层字段，也可以删除整个子表。随后在场景中重写该字段，就得到一次完整替换。

关闭地球的大气：

```toml
name = "Earth without atmosphere"

[[bodies]]
name = "Earth"
template = "bodies/earth.toml"
unset = ["atmosphere"]
```

`unset` 只接受模板中已有的顶层字段名。重复、拼错、点分子路径和非字符串条目都会报错。未引用模板的天体使用 `unset` 也会报错。关闭子表用 `unset`；`atmosphere = false` 不是合法的大气配置。

自转模型存在互斥关系。例如让月球跟随轨道同步转动时，需要清除模板中的 IAU 自转字段：

```toml
name = "Synchronous Moon example"

[[bodies]]
name = "Earth"
template = "bodies/earth.toml"
gm_km3_s2 = 403503.235625 # Earth + Moon, DE440; see scenes/earth_moon.toml.

[[bodies]]
name = "Moon"
parent = "Earth"
template = "bodies/moon.toml"
unset = ["pole_ra_dec_deg", "prime_meridian_deg", "pole_rate_deg_per_century"]
synchronous = true

[bodies.orbit]
type = "ephemeris"
file = "ephem/moon_geo.eph"
```

## 模板允许的属性

模板按以下清单校验。未知字段、场景字段和错误的字段类型会报错；子表中的字段也按对应模型的清单校验。场景覆盖模板里已有的属性时（例如在继承来的 `[bodies.atmosphere]` 里改一项），合并后的结果按同一份清单再查一次，报错信息会注明错在场景的覆盖里。清单与加载器的解析代码是分开写的：加载器给这些属性或子表新增字段时，要同时加进 `scene_loader.cpp` 的 `validate_body_template` 和下表。

| 属性 | 字段 |
|---|---|
| 类型与大小 | `kind`、`radii_km` |
| 自身引力 | `body_gm_km3_s2`、`gravity_radius_km`、`zonal_harmonics` |
| 恒星与黑洞自身参数 | `temperature_k`、`luminosity_solar`、`spin` |
| IAU 自转 | `pole_ra_dec_deg`、`prime_meridian_deg`、`pole_rate_deg_per_century`、`nut_prec_angles`、`nut_prec_ra`、`nut_prec_dec`、`nut_prec_pm` |
| 其他独立自转模型 | `spin_table`、`free_precession` |
| 表面与形状 | `color`、`style`、`shape`、`texture`、`texture_left_lon_deg`、`texture_lon_direction` |
| 外观子表 | `atmosphere`、`night_lights`、`relief`（起伏图：`texture`，和 `texture` 同一投影的坡度和水面，见 `textures/SOURCES.md`；可选 `exaggeration` 放大坡度，默认 1）、`clouds`（云层：`texture` 是和 `texture` 同一投影的灰度云量图，`albedo` 是云的反照率，可选的 `tropical_drift_deg_per_day`、`midlatitude_drift_deg_per_day` 让云图沿纬线慢慢流动（赤道和 45° 处每天漂过的经度，东为正，示意）；随大气层的开关显示，没有大气层的天体不画）、`sea`（起伏图里水面的太阳反光：`wind_m_s` 定 Cox–Munk 的波面坡度，`refractive_index` 定反射率；需要 `relief`）、`aurora`、`magnetic_field`、`magnetic_map`、`rotor`、`shape_after`（从某个时刻起换用的形状模型：`shape` 和 `time` 或 `time_tdb`，与 `shape` 共用纹理图集）、`rings`（含 `[[rings.bands]]`）、`[[plumes]]`、`[[dark_spots]]`（大气里的暗斑：按十进制年份列出的纬度、宽度和对比度，经度锚点和漂移，淡入淡出，见 `DarkSpot`）、`[[impacts]]`（撞击点：`name`、`class`、`time_tdb`、`lat_deg`、`west_lon_deg`，见 `scene/impact.hpp`） |

`atmosphere.refractivity` 是大气底部的 `n-1`，默认 0；设为 0 可关闭折射像。地球模板采用 550 nm 标准干空气的值，随 Rayleigh 气体的指数密度变化。`atmosphere.refraction_gain` 默认 1，独立调整折射像的显示亮度，日食场景采用 64。加载时计算弯折和透射率表；不支持 `n(r)r` 非单调的波导大气。公式与适用范围见 [日冕与折射来源](../corona/SOURCES.md#地球边缘的大气折射像)。

`[bodies.rotor]` 的 `period_s` 是绕本体 +z 旋转一周的仿真秒数。省略子表里的 `shape` 时，整份本体网格旋转，快进时每秒最多转半圈。指定子表里的 `shape` 时，该网格作为独立转动的部件，每秒最多转四分之一圈。

表数组（`[[rings.bands]]`、`[[plumes]]`）和其他数组一样整体替换：场景写了 `[[bodies.rings.bands]]`，就替换模板里的全部环带；只改 `[bodies.rings]` 里的 `gain` 等字段时，环带沿用模板。

`nut_prec_ra/dec/pm` 可以使用祖先天体的 `nut_prec_angles`，和 PCK 一样：土星的卫星（`mimas.toml`、`tethys.toml`、`rhea.toml`）只写周期项的系数，角度取自 `saturn.toml`。加载器检查本体或最近的祖先给出的角度够不够用，不够时报错。

`name`、`parent`、`orbit`、`orbit_color`、`trail*`、标签和事件关联写在场景中。`synchronous`、`roche`、`accretion_disk`、`pulsar` 等涉及父天体或其他天体的模型也留在场景。当前允许清单以外的模型（包括彗尾和脉动）继续内联定义。

`[bodies.accretion_disk]` 写在受体星上，`companion` 指向为它供质的洛希面恒星。盘绕受体、与轨道共面，`outer_radius_km` 和 `half_height_km` 给出不透明表面的大小，半厚至少为盘内恒星的半径（`radii_km`）。`rim_temperature_k`、`temperature_power` 和 `inner_temperature_k` 控制盘面温度，最后一项是盘心温度上限。可用 `rim_heating_k` 和 `rim_heating_radius_km` 设定朝供体一侧的温升和范围。`[light_curve].stars` 列出这对恒星时，加载场景会按远方观察者的方向计算一圈双食光变。

`illustrative_torus = true` 会在画面中积分厚盘气体的发光和吸收，并画出受体星、气流和喷流。光变曲线按不透明 slab 计算。

`body_gm_km3_s2` 是天体自身 GM，可以放进模板。`gm_km3_s2` 用于轨道计算，可能包含卫星质量，由场景提供。例如地月场景中 Earth 的轨道计算 GM 是地球与月球之和。加载器仍会在缺少 `body_gm_km3_s2` 时采用 `gm_km3_s2` 等已有回退规则；模板给出了自身 GM 时，场景只覆盖 `gm_km3_s2` 不会同时修改自身 GM。

## 已有模板与来源

- [earth.toml](earth.toml)：地球的半径、自身 GM、自转、表面、大气、夜光和极光。物理参数、散射模型和显示取值的来源在文件头。
- [moon.toml](moon.toml)：月球的半径、自身 GM、简化 IAU 自转和贴图。地月与空间望远镜场景通过 `unset` 省略长期极轴漂移。
- [sun.toml](sun.toml)：太阳的半径、自身 GM、温度和光度。各场景的自转和显示颜色单独定义。
- 行星：[mercury.toml](mercury.toml)、[venus.toml](venus.toml)（含大气）、[mars.toml](mars.toml)（含大气和 J2–J4）、[jupiter.toml](jupiter.toml)（含 J2、J4、环和紫外统计极光卵）、[saturn.toml](saturn.toml)（含 J2、J4、环和卫星周期项的角度）、[uranus.toml](uranus.toml)（含环、红外 H₃⁺、磁场和卫星周期项的角度）、[neptune.toml](neptune.toml)（含环弧、红外 H₃⁺、磁场和卫星周期项的角度）。
- 火卫：[phobos.toml](phobos.toml)、[deimos.toml](deimos.toml)（形状模型）。
- 伽利略卫星：[io.toml](io.toml)（含羽流）、[europa.toml](europa.toml)、[ganymede.toml](ganymede.toml)、[callisto.toml](callisto.toml)。
- 土星的大卫星：[mimas.toml](mimas.toml)、[enceladus.toml](enceladus.toml)（含羽流）、[tethys.toml](tethys.toml)、[dione.toml](dione.toml)、[rhea.toml](rhea.toml)、[titan.toml](titan.toml)（含大气）、[iapetus.toml](iapetus.toml)。
- 天王星的卫星：[miranda.toml](miranda.toml)、[ariel.toml](ariel.toml)、[umbriel.toml](umbriel.toml)、[titania.toml](titania.toml)、[oberon.toml](oberon.toml)。
- 海王星的卫星：[proteus.toml](proteus.toml)（形状模型）、[triton.toml](triton.toml)（含大气和羽流）。

只在一个场景里出现的天体（Janus 等）仍内联定义，第二个场景用到时再移进模板。

属性移入模板时，将它的来源与简化说明一起放入模板。贴图、形状和数据文件的许可证与处理步骤仍在各资源目录的 `SOURCES.md` 中。场景保留轨道、事件与覆盖值的依据。

木星模板的 `[aurora]` 用 `jovian_oval` 指定统计主卵、`footprints` 指定卫星磁场映射表。后者需要木星的 1 bar 半径，只计算场景中名为 Io、Europa、Ganymede 的直属卫星，缺少某颗时跳过它。两项数据用于紫外叠加视图，与地球的 `kp` 模型互斥；格式、领先角和显示假设见 [地磁与极光来源](../geomag/SOURCES.md)。`flux_tubes` 是从这三颗卫星到足点的线（同一磁场，需要 `footprints`），由事件的 `flux_tubes = true` 或手动开关（`B`）打开。

天王星、海王星用 `[aurora].uranus_map`、`neptune_map` 读取 JWST 的 H₃⁺ 图，要求各自的半径，与其他极光字段互斥，只在红外视图显示，见 [极光来源](../geomag/SOURCES.md#天王星与海王星的红外-h₃)。

Ganymede 模板的 `[aurora].ganymede_map` 读取 HST 的柱亮度图（平均图和电流片上方、内部、下方三组），要求木卫三的球形半径，与其他极光字段互斥。它随本体旋转，按木卫三此刻在木星磁层里的磁纬选组，通过紫外叠加开关显示；数据和 70 km 壳层的说明见 [极光来源](../geomag/SOURCES.md#ganymede-自身的紫外极光)。

`[magnetic_field]` 是内部磁场的磁力线（水星、地球、木星、土星、天王星、海王星、Ganymede 在模板里）：`coefficients` 和 `reference_radius_km` 必填，随年份变化的模型改写 `coefficients_by_year`（地球的 IGRF，模型、磁力线和 `[magnetic_map]` 的强度图按整年跟随时钟）；`l_shells`、`lines_per_shell`、`max_radius` 决定画哪些线，`parent_field = true` 把母行星在该处的磁场当作均匀场叠加进来（Ganymede，需要母行星的 `flux_tubes` 表），同时用 `open_latitudes_deg` 给出从表面哪些磁纬出发画离开的线，`l_shells` 里只保留闭合的线；`color` 是颜色；模型坐标系和本体系不同时用 `flip_pole` 加 `longitude_offset_deg`，或用 `prime_meridian_deg` 给模型自己的 W（两种写法互斥）。默认不显示，由事件的 `magnetic_lines = true` 或手动开关（`B`）打开，见 [地磁与极光来源](../geomag/SOURCES.md#行星内部磁场磁力线)。

`[magnetosphere]` 是磁层顶和弓激波的线框（水星、地球、木星、土星、天王星、海王星在模板里）：`model` 必填，`"shue1998"`（地球，用 `pressure_npa` 和 `bz_nt`）、`"joy2002"`（木星，只用 `pressure_npa`）或 `"pressure_balance"`（需要 `surface_field_nt`，可用 `compressibility`、`flaring`，典型状态可写 `standoff_radii` 代替压强）。`radius_km` 是模型的长度单位（默认赤道半径），`tail_radii` 是磁尾画多远，`flare_end_radii` 以外磁尾不再张开，`aberration = true` 时轴线按 `wind_speed_km_s` 偏向行星运动方向，`shock_standoff_ratio`、`shock_eccentricity` 是前两种模型的激波，`color`、`shock_color` 是颜色。`"power_law"` 是由探测器穿越点拟合的幂律（土星、水星）或固定大小：`standoff_at_1npa` 加 `compressibility`，或只给 `standoff_radii`；`flaring`、`flaring_slope` 是张开参数及其随压强的斜率；激波可用自己的 `shock_at_1npa` 加 `shock_compressibility`，`shock_focus_radii` 把激波的焦点沿轴前移，`center_offset_radii` 把两个面的中心从天体中心移开（水星的偶极在北）。`"tabulated"` 是查表的面（海王星，需要写在前面的 `[magnetic_field]`）：`boundaries` 指定剖面表，两个面各有两个自转相位、四个方位的剖面，按偶极此刻的相位混合，所以不轴对称、随自转变化、不随太阳风变；`standoff_scale`、`shock_scale` 是乘在表上的系数。`solar_wind` 是逐小时的太阳风记录（地球的 `geomag/solar_wind.bin`），面板选 `Always` 时或事件项写 `recorded = true` 时按当时的记录画，缺测和记录以外用典型值。`"heliosphere"` 是太阳的日球层顶和终端激波（写在太阳系场景的太阳上，不在 `sun.toml` 里，因为别的恒星也用这个模板）：需要 `nose_ecliptic_deg = [黄经, 黄纬]` 和 `standoff_radii`，`radius_km` 取 1 au；事件用 `standoff_radii`、`shock_radii` 给出两个面的大小。默认不显示，由面板的 `Magnetosphere` 或事件的 `magnetosphere = [{ body = "...", pressure_npa = ..., bz_nt = ..., wind_speed_km_s = ..., standoff_radii = ..., tail_radii = ..., faint = true }]` 打开：每项画一个状态，没写的量取模板的典型值，`recorded = true` 改用天体的太阳风记录（不能再写自己的压强、Bz 和速度），`faint` 的一项是对比用的淡色面，见 [磁层顶与弓激波](../geomag/SOURCES.md#磁层顶与弓激波)。

`[magnetic_map]` 是固定在本体上的磁场伪彩色图（火星、月球、木星、天王星、海王星是径向场，地球的 `surface` 是总强度，由文件决定）：`surface` 必填，`deep` 是更深处的第二张图（木星），`crust` 是只有地壳的图（地球）。行按行星中心纬度，列从 180°W 开始，和贴图的 `texture_left_lon_deg` 无关；`[magnetic_field]` 用 `prime_meridian_deg` 时图在模型的经度里，渲染时跟着转（海王星）。默认不显示，由事件的 `magnetic_map = "surface"`、`"deep"` 或 `"crust"`、手动开关（`G`）打开，见 [径向磁场图](../geomag/SOURCES.md#径向磁场图)。

## 加载接口与离线工具

`load_scene_file` 按 `<assets>/scenes/<file>.toml` 推导资源根目录。使用 `load_scene_string` 加载带模板的内容时，需要传入 `asset_root`。错误信息包含场景来源、天体名和模板路径。

读取天体属性的 Python 工具使用 `tools/scene_templates.py` 的 `load_scene(path)` 展开模板；`philae_bake.py` 和 `fit_moon_phase.py` 已使用此接口。这个辅助函数负责引用、清除和合并，完整的属性与物理模型校验由 C++ 加载器执行。只读取场景级 `[sky]` 或独立数据文件的工具可以继续使用 `tomllib`。

测试命令：

```powershell
cmake --build --preset win-msvc-release --target astraxis_tests
ctest --test-dir build/win-msvc -C Release --output-on-failure
python -m unittest discover -s tools -p test_scene_templates.py
```

Europa 模板的 `[aurora].europa_visible = true` 在引用它的场景中启用可见光氧极光，要求球形半径 1560.8 km，与其他极光字段互斥。亮度的分布跟着太阳、等离子体片和磁力线，见 `Simulation::moon_aurora`。木星遮日按当前父天体位置计算；没有父天体时只用当地昼夜可见度。数据和显示假设见 `../geomag/SOURCES.md`。

Io 和 Callisto 模板分别使用 `[aurora].io_visible = true`、`[aurora].callisto_visible = true`，要求各自的球形参考半径（1821.49、2410.3 km）。与其他极光字段互斥。导览或手动极光视图启用显示，谱线强度和显示增益分别保存，亮度由 `Simulation::moon_aurora` 按此刻的情况调整（Io 的钠和 SO₂ 跟着日食，Europa、Callisto 跟着等离子体片）；模型与来源见 `../geomag/SOURCES.md`。
