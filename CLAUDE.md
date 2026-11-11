# Astraxis

Astraxis 是一个天文动态背景程序。目标是安静、低功耗，并让画面尽量符合已知的物理和观测。交互服务于观看。C++20、CMake、SDL3 和 SDL_GPU，支持 Windows 和 Linux。场景、安装和操作见 [README](README.md)，模型假设和可见偏差见 [已知局限](limitations.md)。

## 构建与验证

Windows 使用 VS 2026 自带的 CMake。以下是本机 Enterprise 的路径，其他安装可按 [CI 的 vswhere 用法](.github/workflows/build.yml) 查找。

```powershell
$vsCMakeBin = 'C:\Program Files\Microsoft Visual Studio\18\Enterprise\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin'
& "$vsCMakeBin\cmake.exe" --preset win-msvc
& "$vsCMakeBin\cmake.exe" --build --preset win-msvc-release
& "$vsCMakeBin\ctest.exe" --test-dir build/win-msvc -C Release --output-on-failure --no-tests=error
```

Linux 依赖和安装命令见 [README](README.md)。WSL 中把源码同步到 Linux 文件系统后再构建，排除 `build/` 和 `.git/`。

```bash
cmake --preset linux
cmake --build --preset linux-release
ctest --test-dir build/linux -C Release --output-on-failure --no-tests=error
```

- 调试时把构建预设的 `-release` 改为 `-debug`，CTest 配置改为 `Debug`。N 体积分让 Debug 测试慢得多。GCC 和 Clang 均应零警告。WSL 用不同编译器时给每个构建目录指定独立的 `FETCHCONTENT_BASE_DIR`，避免 SDL 的预编译头互相覆盖。
- 修改 `core/`、`ephem/`、`scene/` 后运行测试。测试从源码树读 `assets/`，入口是 [test_main.cpp](tests/test_main.cpp)。按领域把回归检查放到 `tests/*_tests.cpp`，共用 `test_util`。
- 渲染改动要实际运行程序看画面。Windows 可用 `SDL_GPU_DRIVER=direct3d12` 或 `SDL_GPU_DRIVER=vulkan` 检查两种着色器输出。WSL 的软件 Vulkan 可检查功能和画面，性能以硬件实测为准。
- 用 `--scene <名称> --event <n>` 打开指定事件，编号从 1 开始。运行时读取 exe 旁的资源副本，改动 `assets/` 后重新构建。场景改名或删除时清理构建目录中对应的旧副本。
- [截图脚本](tools/capture/capture_window.ps1) 的 `-Config` 可传入一个空配置文件以复现默认设置。镜头未就位时调大 `-WaitSeconds`，留出空闲 60 秒自动导览的余量。检查颜色时留意展签遮挡和 ACES 色调映射。
- Release 和系统启动的屏保用 `ASTRAXIS_LOG=<文件>` 收集日志。测试 `.scr` 参数时用 `ProcessStartInfo` 并设置 `UseShellExecute = false`，否则文件关联会插入 `/S`。`astraxis_layer_probe.exe` 可用 GDI 单独测试 Windows 壁纸挂载。

## 架构与约定

| 路径 | 职责 |
|---|---|
| [src/app](src/app)、[src/platform](src/platform) | 窗口、输入、设置、屏保、壁纸与平台集成 |
| [src/view](src/view)、[src/render](src/render) | `SceneRenderer` 构造绘制项并编排 pass，渲染层使用 `CameraView` 和绘制项 |
| [src/scene](src/scene) | 场景、天体、参考系、相机、导览与 `Simulation` |
| [src/ephem](src/ephem)、[src/core](src/core) | 运动来源、星历、积分、时间、单位和双精度数学 |
| [shaders](shaders)、[assets](assets)、[tools](tools) | HLSL、运行时资源、离线资源生成 |

- 依赖按表中顺序向下。`astraxis_sim` 包含 `core`、`ephem`、`scene`，`astraxis_settings` 负责设置，两者均独立于 SDL 和 GPU。图形后端统一使用 SDL_GPU，保持现有分层和小型程序的范围。
- `MotionSource::eval(t_tdb)` 返回相对父天体的 ICRF 状态。`Scene::update` 沿父子链得到 `icrf_position`，再经 `FrameTransform` 得到供相机、导演、标签和渲染使用的 `world_position`。测试中先更新场景再读位置或光照。`Scene::icrf_state_at` 返回相对场景根的状态，求相对轨道时先减去中心天体。
- 宿主每帧推进一次 `Simulation`，再为每个输出计算视图并渲染。GPU 设备共享，交换链和 `SceneTargets` 属于各输出，各自等待交换链。屏保 `.scr` 由主程序复制，`/p` 的跨进程预览使用 GDI。壁纸和屏保设置的优先级为代码默认值 → `[view]` → 模式小节 → 命令行。
- 仿真时间为自 J2000 起的 TDB 秒数，使用 `double`。长度用 km、时间用 s、GM 用 km³/s²。积分型运动生成节点表或使用滚动窗口，以支持跳时和倒放。
- CPU 使用右手系和双精度位置。先以双精度减相机位置，再转为 GPU 的 float。深度用 reversed-Z。渲染在线性 HDR 空间计算，配置和代码中的颜色常量用 sRGB，颜色纹理用 sRGB 格式。
- 类型用 `PascalCase`，函数和变量用 `snake_case`，成员加 `m_`，常量加 `k`。代码和注释使用英文，项目内开发说明使用中文。

## 着色器

- [Shaders.cmake](cmake/Shaders.cmake) 下载并校验固定版本的 DXC，生成嵌入程序的 DXIL 和 SPIR-V。需要覆盖编译器时用 `ASTRAXIS_DXC` 指向支持 SPIR-V 的 DXC。
- 把 `<name>.<vert|frag|comp>.hlsl` 加到 [CMakeLists.txt](CMakeLists.txt) 的 `astraxis_add_shaders`。包含生成的 `<shaders/<name>.<stage>.h>`，把 DXIL 和 SPIR-V 一起交给 `create_shader`。
- 按 SDL_GPU 的 `register`/`space` 约定绑定资源。片元纹理用 `FRAGMENT_TEXTURE` 配对纹理和采样器，输入和输出按 `TEXCOORD0`、`TEXCOORD1` 等顺序声明，系统语义照常使用。读存储缓冲区用 `ByteAddressBuffer`，与 D3D12 的 RAW SRV 对齐。SDL 的 Vulkan 后端每个 uniform 绑定限 4 KiB。用 `static_assert` 检查大小，大表用纹理或存储缓冲区。
- D3D12 的 ImGui 后端需要 DXBC，创建设备时声明 `DXIL | DXBC`。各宿主都调用 `ImGui_ImplSDLGPU3_NewFrame()` 以初始化后端资源。全局快捷键用 `WantTextInput` 判断是否正在输入文本。
- 随时间变的磁力线（Ganymede 叠加木星磁场）只在被画出时更新：绘制方设 `BodyMagneticField::lines_in_view`，`Scene::update` 下一帧才追踪，测试里要先设这个标志。
- 磁层顶、弓激波等线框曲面走轨道线的 pass（`OrbitLinePass::add_line`），每帧由 `boundary_wireframe` 生成。画面用的太阳风经 `Scene::magnetosphere_wind` 取得，不直接读事件里的值。
- 表面叠加图按 `body.vert.hlsl` 的行星中心纬度采样。光线和小天体的垂距用叉积求，避免两个大平方数相减。极薄的壳层先裁剪射线再积分，大面积加性光斑先画到半分辨率单采样目标再合成。
- `SceneRenderer::render` 在检查黑洞 pass 的 `ready()` 之前更新星表资源。透镜雅可比用有限差分求，轨道线使用几何位置。修改黑洞合成时同时检查屏幕追踪、星表和立方体天空回退。

## 场景与数据

- 场景在 [assets/scenes](assets/scenes)，共用天体字段见 [Body Templates](assets/bodies/README.md)，解析入口是 [scene_loader.cpp](src/scene/scene_loader.cpp)。顶层 TOML 键放在首个表头之前。切换互斥字段时用 `unset` 清除模板值。`body_gm_km3_s2` 为天体自身 GM，`gm_km3_s2` 可含卫星质量。
- 没有太阳天体的场景可用 `from_body = "Sun"` 引用太阳方向。没有地球时用 `observer` 和 `observer_offset_km` 锚定地球方向的设计视点。`flyby` 需要最近点两侧都有数据，撞击事件用视点字段构图。`elevation_deg` 覆盖俯仰角，`phase_deg` 调整方位角。
- 视点超过星表完整距离 `max_viewer_pc` 时，用 [模型星空](assets/sky/SOURCES.md)。近地航天器若用平根数画一圈轨道，用 `trail = "orbit"`。换参考系或时间倍率后检查轨迹密度和事件节奏。
- 用原始观测或论文表格标注物理量来源。场景 TOML 文件头列出引用标签，正文引用这些标签。能用独立数据交叉验证的位置、方向、距离或时刻，写进相关测试。核对时间尺度、光行时、参考中心、经纬度正方向以及距离是到中心还是到表面。
- 星历工具在 [tools/bake](tools/bake)。压缩容差约束位置，点火和数据接缝要另留节点。Δv 按原表节点计算，靠近扁行星时扣除带状谐项。分段 SPK 按实际记录覆盖取值，详见 [星历资料](assets/ephem/SOURCES.md)。
- 航天器转换用 [make_spacecraft.py](tools/shapes/make_spacecraft.py)。导入时应用节点变换，用已知部件标定尺寸，核对 OBJ 导出轴和关键部件法线。模型用硬边标志，删除游离顶点。`point_at` 指定 +z 方向，`fill_light` 是示意补光。烘焙自阴影时要核对模型的 +z 与场景中的入射阳光方向。
- 依赖由 CMake FetchContent 固定到版本或提交。外部数据选用许可明确的来源，优先公有领域。原始大文件按资料页下载。新增依赖、移植代码或外部资源时更新 [第三方声明](THIRD_PARTY_NOTICES.md) 和相关资料页。新增或改名场景时更新 README 的场景表。

| 资料 | 复现入口 |
|---|---|
| [星历](assets/ephem/SOURCES.md)、[星表与光变](assets/stars/SOURCES.md)、[模型星空](assets/sky/SOURCES.md) | `tools/bake/`、任务目录、`tools/stars/`、`tools/sky/` |
| [纹理与地形](assets/textures/SOURCES.md)、[形状与航天器](assets/shapes/SOURCES.md) | `tools/textures/`、`tools/shapes/` |
| [地磁与极光](assets/geomag/SOURCES.md)、[日冕](assets/corona/SOURCES.md)、[行星环](assets/rings/SOURCES.md) | 数据处理命令和模型依据见各资料页 |

### 数据核对备注

- 提取 PDF 表格时用 `pdftotext -raw`，并与原页核对列和负号。Wiley、Nature、ADS 拒绝脚本抓取时，Nature 的开放获取文章可在浏览器里读正文，作者主页的 PDF 可下载后转文本，DSpace 机构仓库（如 Imperial Spiral）可从 `/server/api/core/items/<uuid>/bundles` 的 TEXT 包取全文文本（没有图）。系数印在图里时 `pdftotext` 取不到，用 PyMuPDF 把那一页渲染成图再读。曲线图先用 `page.get_drawings()` 看是不是矢量：是的话直接取路径坐标（描边常被转成填充的轮廓，颜色在 `fill` 里），例见 `tools/geomag/neptune_boundaries_fig4.py`。老文章的 ADS PDF 地址若只返回第一页，在浏览器打开 `ui.adsabs.harvard.edu/scan/manifest/<bibcode>` 查看扫描件。需要单页时，可用页面的 `/v1/accounts/bootstrap` 返回的访问令牌请求 `/v1/scan/image/iiif/2/bitmaps-~seri-~<刊>-~<卷>-~600-~<页>.000/full/1600,/0/default.jpg`。
- Tiranti 2026 的 FITS 中 `INTENSITY_PROFILES` 约为 `3e-12`，按同温度和密度换算的总发射约 `4e-14 W m⁻³ sr⁻¹`。与 Melin 2025 的柱亮度差约 70 倍，因此仅用剖面形状，定标见 [冰巨星红外资料](assets/geomag/SOURCES.md)。

## 待办

- 定期刷新地球的逐小时太阳风记录（OMNI2，现到 2026-10），评估用传播模型给外行星接入太阳风。按 Mejnertsen 2016 图 4 的两个截面给海王星磁层顶加极对日时的极尖凹陷，激波补偏轴处的变化（都需要非轴对称的面）。卫星足点补尾迹，核对 Bagenal & Delamere 2011 的密度剖面原文，给 Ganymede 的开闭边界加上下游不对称和阿尔文翼的后掠。
- 定期刷新 Roman、TESS、Euclid、JWST、Chandra 的预报段。随 ESA 新核更新 Hera 抵近弧和两颗小行星的位置，补 Milani、Juventas 分离后的模型和轨迹，以及撞击后的 Dimorphos 形状。
- 按 Dyudina 2015 的旋转周期试做土星可见光亮弧，先核对 Leicester 相位表许可。等待 Callisto 空间分辨观测、Io 的 SO₂ 食中时序、Europa 更多食中图像、火星新绿线观测与高度剖面。冰巨星待补夜侧、地方时和垂直剖面的约束。
- 补海王星暗斑伴云、GDS 形变及其他暗斑，核对 NDS-2018 的后期消退。待正式资料更新 Sincérité 的长度、光学深度和 Égalité 的演化。
- 补 PDS 70 的光学厚 ¹²CO 和 Mie 波长颜色。HD 206893 采用更高分辨率缺口拟合，并跟进 B I 候选体。旋臂优先参考已有图案速度测量的 MWC 758。
- 核对嫦娥五号 T1 和天问一号尺寸，补 Orion 再入姿态与防热罩曲面、ISEE-3 整体自旋和筒身带、DART 电池片图案。评估运行时自阴影和镜面高光。
- 为 Castor 补 TESS 第 44、45、47、71 扇区，核对 2021-12-08 耀斑归属，绘制轨道法线，等 C 的观测弧延长后更新轨道。
- 支持 macOS 的 Metal 着色器和 CI。J1141−6545 待几何解进一步约束后加入自转轴进动。

提交前更新已完成的待办，把可复用的排障经验放到对应章节或代码旁，把新物理或数据近似写进 [已知局限](limitations.md)。开发过程记在提交信息里。
