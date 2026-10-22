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
| 外观子表 | `atmosphere`、`night_lights`、`aurora`、`rotor`、`rings`（含 `[[rings.bands]]`）、`[[plumes]]` |

表数组（`[[rings.bands]]`、`[[plumes]]`）和其他数组一样整体替换：场景写了 `[[bodies.rings.bands]]`，就替换模板里的全部环带；只改 `[bodies.rings]` 里的 `gain` 等字段时，环带沿用模板。

`nut_prec_ra/dec/pm` 可以使用祖先天体的 `nut_prec_angles`，和 PCK 一样：土星的卫星（`mimas.toml`、`tethys.toml`、`rhea.toml`）只写周期项的系数，角度取自 `saturn.toml`。加载器检查本体或最近的祖先给出的角度够不够用，不够时报错。

`name`、`parent`、`orbit`、`orbit_color`、`trail*`、标签和事件关联写在场景中。`synchronous`、`roche`、`pulsar` 等涉及父天体或其他天体的模型也留在场景。当前允许清单以外的模型（包括彗尾和脉动）继续内联定义。

`body_gm_km3_s2` 是天体自身 GM，可以放进模板。`gm_km3_s2` 用于轨道计算，可能包含卫星质量，由场景提供。例如地月场景中 Earth 的轨道计算 GM 是地球与月球之和。加载器仍会在缺少 `body_gm_km3_s2` 时采用 `gm_km3_s2` 等已有回退规则；模板给出了自身 GM 时，场景只覆盖 `gm_km3_s2` 不会同时修改自身 GM。

## 已有模板与来源

- [earth.toml](earth.toml)：地球的半径、自身 GM、自转、表面、大气、夜光和极光。物理参数、散射模型和显示取值的来源在文件头。
- [moon.toml](moon.toml)：月球的半径、自身 GM、简化 IAU 自转和贴图。地月与空间望远镜场景通过 `unset` 省略长期极轴漂移。
- [sun.toml](sun.toml)：太阳的半径、自身 GM、温度和光度。各场景的自转和显示颜色单独定义。
- 行星：[mercury.toml](mercury.toml)、[venus.toml](venus.toml)（含大气）、[mars.toml](mars.toml)（含大气和 J2–J4）、[jupiter.toml](jupiter.toml)（含 J2、J4 和环）、[saturn.toml](saturn.toml)（含 J2、J4、环和卫星周期项的角度）。
- 火卫：[phobos.toml](phobos.toml)、[deimos.toml](deimos.toml)（形状模型）。
- 伽利略卫星：[io.toml](io.toml)（含羽流）、[europa.toml](europa.toml)、[ganymede.toml](ganymede.toml)、[callisto.toml](callisto.toml)。
- 土星的大卫星：[mimas.toml](mimas.toml)、[enceladus.toml](enceladus.toml)（含羽流）、[tethys.toml](tethys.toml)、[dione.toml](dione.toml)、[rhea.toml](rhea.toml)、[titan.toml](titan.toml)（含大气）、[iapetus.toml](iapetus.toml)。

只在一个场景里出现的天体（天王星、海王星及其卫星、Janus 等）仍内联定义，第二个场景用到时再移进模板。

属性移入模板时，将它的来源与简化说明一起放入模板。贴图、形状和数据文件的许可证与处理步骤仍在各资源目录的 `SOURCES.md` 中。场景保留轨道、事件与覆盖值的依据。

## 加载接口与离线工具

`load_scene_file` 按 `<assets>/scenes/<file>.toml` 推导资源根目录。使用 `load_scene_string` 加载带模板的内容时，需要传入 `asset_root`。错误信息包含场景来源、天体名和模板路径。

读取天体属性的 Python 工具使用 `tools/scene_templates.py` 的 `load_scene(path)` 展开模板；`philae_bake.py` 和 `fit_moon_phase.py` 已使用此接口。这个辅助函数负责引用、清除和合并，完整的属性与物理模型校验由 C++ 加载器执行。只读取场景级 `[sky]` 或独立数据文件的工具可以继续使用 `tomllib`。

测试命令：

```powershell
cmake --build --preset win-msvc-release --target astraxis_tests
ctest --test-dir build/win-msvc -C Release --output-on-failure
python -m unittest discover -s tools -p test_scene_templates.py
```
