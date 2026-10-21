# Astraxis

用来“发呆”的天文动态背景程序：太阳系与各探测器轨迹、JWST L2 晕轨道、金星、火星、木星和土星系统、彗星与哈雷舰队、Alpha Centauri、系外行星共振链、多星系统、食变星 Algol 和造父变星 δ Cephei（带从地球看的光变曲线）、脉冲星双星、Sgr A* 黑洞外观、《星际穿越》的虚构黑洞 Gargantua（按论文参数，电影版与带多普勒频移的真实样子对比）、致敬《坎巴拉太空计划》的 Kerbol 系统（游戏里读出的轨道、Kerbin 历法、程序生成的表面）等，共 25 个场景。
C++20 + CMake + SDL3 + SDL_GPU。Windows（D3D12）和 Linux（Vulkan）都能构建运行；Windows 上还有全屏、屏保和动态壁纸。将来做 macOS。

- **目标**：安静、好看、物理上可信，主要用来发呆，交互只是辅助；功耗低，可以长时间挂在后台当壁纸；单个 exe，启动快，依赖少
- **不做**：游戏性、飞船操控、通用引擎、编辑器
- **工作结束前**：新踩的坑写进下面的“已知的坑”；做完的待办从“待办”里删掉，新的下一步加进去；物理或数据上的新简化写进 [limitations.md](limitations.md)。不另写进度日志，过程写在提交信息里

## 构建（Windows）

工具链来自 VS 2026 Enterprise（cmake 4.3 和 ninja 是 VS 自带的，不在 PATH 上）：

- CMake: `C:\Program Files\Microsoft Visual Studio\18\Enterprise\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe`
- dxc：不用 Windows SDK 里的（它不能输出 SPIR-V）。`cmake/Shaders.cmake` 在配置时下载固定版本的微软官方 DXC 发行版（v1.9.2609，校验 SHA-256，放在 `build/_deps/dxc-src`），同一个 dxc 输出 DXIL 和 SPIR-V

```
cmake --preset win-msvc
cmake --build --preset win-msvc-debug
build/win-msvc/Debug/astraxis_tests.exe     # 单元测试
```

建议用 Bash 工具运行构建命令：在 PowerShell 里，msbuild 的 `-v:minimal` 参数会被误解析（报 MSB1016）。改动 `core/`、`ephem/`、`scene/` 之后要跑测试；涉及渲染的改动要实际运行程序看画面。依赖源码在 `build/_deps/`，可以直接查 SDL 和 ImGui 的源码与示例。本机没有 Vulkan SDK，也没有 glslc。

## 构建（Linux，本机用 WSL Ubuntu）

```
cmake --preset linux                    # Ninja Multi-Config，build/linux
cmake --build --preset linux-debug
build/linux/Debug/astraxis_tests
```

- GCC 和 Clang 都要求 0 警告
- 在 WSL 里构建时，先把仓库同步到 WSL 自己的文件系统（如 `~/astraxis`，不带 `build/`、`.git/`）再构建。直接在 `/mnt/c` 上构建很慢。
- WSL 只有 llvmpipe（CPU 软件 Vulkan），只能检查功能和画面，不能看性能。WSLg 窗口会出现在 Windows 桌面上（属于 `msrdc.exe`，标题是“Astraxis (Ubuntu)”）。要在真 GPU 上检查 SPIR-V 着色器，在 Windows 上用 `SDL_GPU_DRIVER=vulkan` 运行即可。

## 运行与调试

- 命令行 `--scene <名称> --event <n>` 可以直接打开某个场景的第 n 个事件（从 1 开始），方便验证画面。还有 `--mode window|fullscreen|screensaver|wallpaper`、`--display <n>`、`--fps <n>`、`--config <path>`；命令行只覆盖本次运行，不写进 `config.toml`。
- Windows 上构建会在 exe 旁边复制一份 `astraxis.scr`（屏保就是同一个程序，识别 `/s`、`/p <hwnd>`、`/c`）。测试 .scr 不要用 `Start-Process`：ShellExecute 对 .scr 的关联是 `"%1" /S`，会在你的参数前面插入 `/S`。要用 `ProcessStartInfo` 并设 `UseShellExecute = false`。
- `astraxis_layer_probe.exe`（仅 Windows）：不用 GPU 检查壁纸挂载（用 GDI 窗口代替），可以在虚拟机或旧版 Windows 上测 `WallpaperLayer`，参数是保持的秒数。
- Release 版和 Windows 启动的屏保都没有控制台；设置环境变量 `ASTRAXIS_LOG=<文件>` 会把日志追加到这个文件，`SDL_LOGGING=app=verbose` 还会每 10 秒记录一次屏保各输出的帧率。
- 截图看画面：`tools/capture/capture_window.ps1` 启动 astraxis（可带 `--scene`/`--event`），按键（如 `OLH` 隐藏轨道、标签和界面）后用 PrintWindow 截取 1280×720 客户区存成 PNG。D3D12、Vulkan 和 WSLg 窗口都能抓。

依赖全部用 CMake `FetchContent` 拉取并**固定 tag**（stb 上游没有 tag，固定到 commit）：SDL3、glm、Dear ImGui（使用 `imgui_impl_sdl3` + `imgui_impl_sdlgpu3` 后端）、toml++、stb_image。不使用 vcpkg。Windows 版默认用 VC-LTL（Debug 不用，见 `cmake/VCLTL.cmake`）。

## 架构

依赖方向严格自上而下，下层不得 include 上层：

```
src/app/       宿主：main、窗口（App）、屏保（Screensaver）、壁纸（Wallpaper，仅 Windows）；后两者共用 SceneHost
src/platform/  SDL 窗口、路径；Win32 的壁纸挂载（WallpaperLayer）和托盘图标
src/view/      SceneRenderer：把 Simulation 转成 draw item，编排各 pass 和后处理
src/render/    SDL_GPU 封装与各渲染 pass（星空、天体、大气层、彗星、行星环、黑洞、轨道线等）
src/scene/     天体、参考系、相机、场景定义；Simulation（场景 + 时钟 + 相机 + 导览）
src/ephem/     运动来源：解析开普勒 / 星历插值 / N 体（Yoshida）/ GR（Kerr 测地线，IMR；赤道圆轨道用解析解 kerr_circular，任意倾角的轨道用 kerr_apoapsis 在 Kerr-Schild 坐标里 3D 积分）
src/core/      时间、单位、double 数学 —— 不依赖 SDL 和 GPU
shaders/       HLSL 源码，构建时编译并嵌入 exe（共享代码放 *.hlsli）
assets/        运行时资源（场景 TOML、星表、纹理、形状模型、字体），构建后复制到 exe 旁
tools/         离线工具（星表转换、纹理缩放、银河图转换、形状模型转换、Horizons 星历烘焙），不进入运行时
```

CMake 目标：
- `astraxis_sim`：静态库，包含 core、ephem、scene，不依赖 SDL 和 GPU
- `astraxis_settings`：静态库，读写 exe 旁边的 `config.toml`（`app/settings.*`），不依赖 SDL，将来的启动器也会用它
- `astraxis`：可执行程序（Windows 上另复制一份 `astraxis.scr`）
- `astraxis_wallpaper`：仅 Windows，壁纸的设置对话框和启动器，和 astraxis.exe 放在同一目录
- `astraxis_tests`：单元测试

核心抽象：所有天体运动都实现 `MotionSource::eval(t_tdb) -> State`，返回相对父天体的 ICRF 状态，场景层不关心来源是解析、星历还是积分。`Scene::update` 先按父子链算出 ICRF 位置（`icrf_position`），再经当前显示参考系（`FrameTransform`：惯性、主–次天体连线旋转、随天体自转三种）变换成 `world_position`。渲染、相机、标签、导演都只用 `world_position`。

渲染层的 pass 只接收 `CameraView` 和简单的 draw item，由 `src/view/` 的 `SceneRenderer` 把 scene 转换成这些数据。宿主无关的部分：
- `Simulation`（`scene/simulation.hpp`，不依赖 SDL）持有 Scene、时钟、相机和导览；`compute_view` 给每个输出算一份 `OutputView`（相机视图、卫星淡出、content scale）。导览的普通镜头由 `CameraDirector` 选；每隔一两个镜头插一个事件镜头，因为要动时钟，由 `Simulation` 执行：按时间顺序取下一个事件，淡出黑场后把时钟跳到事件前（快到的事件不跳），倍速按剩余时间指数降到事件自己的倍速（飞掠降到最近点前 20 s），显示展签，飞掠镜头持续到最近点之后。飞掠体由事件的 `flyby` 指定（自动找最近的探测器在多探测器场景里常挑错，比如“MRO arrives”找到 Mars Express），`--event` 和导览用同一个视角（`Simulation::event_shot`）。事件可以带幽灵（`[events.ghost]`，格式见 `mars.toml` 的天问一号入轨）：加载时合成一个 `Ghost` 天体（`Body::event`），从点火前的状态外推（`CoastMotion`，`ephem/coast.hpp`）：中心天体有 `zonal_harmonics` 时计入扁率（木星、土星），`tides` 列出的天体计入潮汐（月球场景的地球），否则是二体；有摄动时在加载时积分成节点表。只在该事件激活时显示，轨迹整段画出；探测器是太阳的子天体、只是星历相对行星时（信使号、破晓号），用 `center` 指定绕哪个天体外推。点火的 Δv（`unexplained_delta_v`，按星历节点到节点算）由 `Scene::burn_delta_v` 给出，画在探测器标签下。事件还可以用 `elevation_deg` 指定相机在 up 轴法平面上方的仰角，用 `disk_shifts = false` 在事件激活时关掉吸积盘的多普勒和引力频移（Gargantua 的电影版）。黑场用 `Simulation::fade()`，在后处理合成和标签里乘上
- `LabelLayout`（`scene/label_layout.hpp`）算标签和标记放在哪里（优先级、避让、渐隐、拾取），不负责绘制；`App` 用 ImGui 画
- 光变曲线（场景的 `[light_curve]`，`scene/light_curve.hpp`）：从某个天体方向看，各恒星按黑体圆盘、线性临边昏暗、前面的星挡住后面的星算总光度。`Simulation` 每帧按固定时间网格更新一次采样（已算的保留），各宿主和标签一起画在角落（`app/light_curve_overlay.*`），倍速快到几秒就换掉整条曲线时淡出。星等从“每颗星最亮时、不遮挡”的固定参考光度量起。脉动星（`[bodies.pulsation]`，`Pulsation`）按相位表查半径、温度、光度，`Scene::update` 每帧写回天体；光变曲线按任意时刻查表（`star_disk_at`）。充满洛希瓣的星（`[bodies.roche]`，`scene/roche.hpp`，Algol B）加载时生成临界等势面的网格，逐顶点温度带引力昏暗和伴星的照射；它随轨道转（`Body::synchronous_with`），按形状模型画（星的着色乘上顶点亮度），光变曲线对它按面元求和，它挡在别的星前面时用网格投影的凸包
- `SceneRenderer` 持有 pass 和每个场景的 GPU 资源，显示选项通过 `ViewOptions` 传入
- GPU 分两层：`GpuDevice` 只有一个；`RenderOutput` 每个窗口一个，持有交换链和本输出的 `SceneTargets`。`Frame` 的 `width/height` 是场景尺寸，`output_width/output_height` 是交换链尺寸，两者在有渲染缩放（`render_scale`，各模式都有）时不同。MSAA 采样数（`msaa` 设置：1/2/4，`GpuDevice::set_msaa` 按硬件支持往下取）存在 `GpuDevice` 的 `SceneTargetFormat` 里，各 pass 的管线和 `SceneTargets` 都依赖它：运行时切换要调 `SceneRenderer::refresh_pipelines()` 重建管线，`SceneTargets::ensure` 发现采样数变了会自己重建。新 pass 要实现 `set_format`，并加进 `refresh_pipelines`。所有和输出尺寸相关的纹理（MSAA/HDR/深度、泛光链、黑洞追踪目标）都放在 `SceneTargets` 里，pass 内部不要再缓存按尺寸分配的纹理，否则多个不同尺寸的输出轮流渲染时会反复重建。各输出都用 SDR 交换链，格式相同，共用一个 `SceneRenderer`

宿主（`App`、屏保、壁纸）只负责窗口、输入、UI 和闲置策略，不要把仿真或渲染逻辑放回 `src/app/`。宿主每帧只 `sim.update` 一次，再对每个输出调用 `SceneRenderer::render`；与相机无关的工作（例如尘埃彗尾的传播）不要每个输出重做一遍。

渲染流程：
1. 黑洞光追（有黑洞时）：在屏幕外以半分辨率追踪，记下每条光线的天空方向和盘的光（`shaders/black_hole.hlsli` 是两边共用的追踪代码）；第 2 步的黑洞合成按全分辨率逐像素插值方向，相邻光线不一致处（阴影边缘、光子环附近）在合成着色器里重新追踪；天空逐像素查：星按局部线性化的透镜映射画成点（`StarfieldPass` 把星按方向排进网格，供合成着色器查找），银河直接采样它自己的贴图，天空被压缩到星点超出网格查找范围处（光子环附近）改用带 mip 的立方体贴图
2. 场景 pass，渲染到 HDR（RGBA16F + MSAA）：银河 → 星空 → 天体 → 大气层 → 羽流 → 彗发/离子彗尾 → 尘埃彗尾 → 行星环 → 小行星带/柯伊伯带 → 黑洞合成 → 恒星 → 脉冲星光束 → 轨道线
3. `PostProcess`：泛光 → ACES → sRGB → 交换链
4. ImGui 叠加 pass

场景用 `assets/scenes/*.toml` 描述（格式见 `jupiter.toml` 的注释），新增场景不用改代码。

### 屏保、壁纸与设置

- 所有模式都由 astraxis.exe 在本进程内渲染。`astraxis.scr` **就是 astraxis.exe 本身**，不是启动器：`/s` 下 Windows 把屏保进程退出当作屏保结束（锁屏也跟着它）；`/p <HWND>` 如果跨进程建子窗口，两个线程的输入队列会关联起来，容易卡住或留下孤儿进程
- `astraxis_wallpaper.exe` 是 Win32 设置对话框，选好场景、显示器、帧率等以后拉起 `astraxis.exe --mode wallpaper`；屏保的 `/c` 也用这个对话框。开机自启就是一行 `astraxis.exe --mode wallpaper`
- 设置放在 exe 旁边的 `config.toml`（便携），目录不可写时不保存。`[view]` 是各模式共用的显示选项，`[window]`、`[wallpaper]`、`[screensaver]` 放各模式自己的项，也可以覆盖 `[view]`。优先级：代码默认值 < `[view]` < 模式那一节 < 命令行
- 壁纸要同时支持两种桌面层级：Win11 24H2 起 WorkerW 是 Progman 的子窗口；旧层级要先发 `0x052C` 让 Explorer 拆出 WorkerW。Explorer 重启（`TaskbarCreated`）后重新挂载，退出时用 `SPI_SETDESKWALLPAPER` 刷出原壁纸
- 功耗：壁纸被最大化或全屏窗口盖住、锁屏、显示器关闭时暂停；用电池时按 `battery` 设置限帧或暂停；`render_scale` 降低内部分辨率（黑洞追踪跟着一起降），`msaa` 降低采样数（窗口模式的面板里也有这两项，可以运行时切换）

## 约定

- **时间**：仿真时间 = 自 J2000 起的 TDB 秒数，`double`。真实帧间隔 × 时间倍率推进；单帧真实 dt 上限 0.1 s，防止拖动窗口后跳变。
- **单位**：长度 km，时间 s，GM 单位 km³/s²。GR 模块内部可以用 SI 制或几何单位制，但接口边界统一换算成 km/s。
- **坐标**：右手系。CPU 端位置一律用 `double`（`glm::dvec3`）。
- **精度**：相对相机渲染——每帧在 double 下先减去相机位置，再转成 float 上传 GPU。深度用 reversed-Z（近 1 远 0，`GREATER` 比较）。不引入浮动原点或网格系统。
- **着色器**：HLSL，构建时用 dxc 同时编译为 DXIL（D3D12）和 SPIR-V（Vulkan），生成头文件嵌入 exe。
  - SDL_GPU 对 HLSL 的 `register`/`space` 有固定约定（例如顶点 uniform 用 `space1`，片元资源用 `space2`，片元 uniform 用 `space3`），写着色器前先查 `SDL_CreateGPUShader` 文档。
  - 新增着色器只需加到 `CMakeLists.txt` 的 `astraxis_add_shaders(...)` 里，文件命名为 `<name>.<vert|frag|comp>.hlsl`。然后 `#include <shaders/<name>.<stage>.h>`，就能得到 `k<Name><Stage>Dxil` 和 `k<Name><Stage>Spirv`，两个都传给 `create_shader`。
  - 片元纹理一律用 `common.hlsli` 里的 `FRAGMENT_TEXTURE(type, tex, sampler, n)` 声明。Vulkan 后端要求纹理和采样器是同一个 binding 的 combined image sampler，所以第 n 个纹理只能配第 n 个采样器。
  - 阶段输入/输出（含顶点输入）的语义一律写 `TEXCOORDn`，并在结构体里按 n = 0, 1, 2… 的顺序声明。SPIR-V 的 location 按声明顺序分配，这样 location 才等于 n，顶点属性和前后两个阶段才能对上（D3D12 后端也要求顶点输入用 `TEXCOORD{location}`）。
  - 着色器里读存储缓冲区用 `ByteAddressBuffer`，不用 `StructuredBuffer`：SDL 的 D3D12 后端建的是 RAW SRV（`StructureByteStride = 0`），声明对不上时在部分显卡上满屏花纹（GTX 750 Ti 上的轨道线）。
- **命名**：类型 `PascalCase`，函数/变量 `snake_case`，成员变量 `m_` 前缀，常量 `k` 前缀。代码与注释全部用英文（文档 .md 用中文）。
- **物理常数/轨道根数**：必须注明来源（如 JPL SSD），不凭记忆写数字。论文数值要从原文表格中摘取（标注表号），能用独立数据交叉验证的就写成测试（例如 α Cen 根数对照 ALMA 实测位置）。场景 TOML 文件头用 `[标签]` 列出处，正文引用标签。
- **积分型运动**（N 体、GR）：在场景加载时预先积分成节点表（或用滚动窗口），不逐帧积分，以支持任意跳转和倒放。
- **颜色**：渲染全程在线性 HDR 空间。TOML 和代码里的颜色常量按 sRGB 书写，着色器里用 `srgb_to_linear` 转换；纹理用 `*_UNORM_SRGB` 格式。
- **外部资源**：只用许可证清楚的数据（公有领域优先），并在对应目录的 `SOURCES.md` 里写明出处、许可证和处理方式。原始大文件（例如约 600 MB 的 USGS TIFF）不入库，需要时按 `SOURCES.md` 重新下载。项目本身是 GPL-3.0-or-later；新增第三方库、移植代码或数据时同步更新 `THIRD_PARTY_NOTICES.md`。
- **README**：新增或改名场景时同步更新 README 的场景表。

## 测试

不依赖测试框架，覆盖 `core/`、`ephem/`、`scene/` 和设置文件。CTest 也注册了（`ctest --test-dir build/win-msvc -C Debug`）。测试直接读源码树里的 `assets/`（通过 `ASTRAXIS_ASSET_DIR`）。Debug 下 N 体积分慢，整套测试要一分多钟；Release 只要几秒。

按领域分文件，每个文件导出一个 `run_<领域>_tests()`，由 `tests/test_main.cpp` 依次调用；`check`、`load_scene_or_die` 等公用函数在 `tests/test_util.*`。新增测试放进对应领域的文件，只有这个文件用到的辅助函数放在它的匿名命名空间里：

| 文件 | 内容 |
|---|---|
| `core_tests.cpp` | 日历（含 Kerbin 历法）/时间尺度、星表（含重投影）与黑体颜色、银道坐标系 |
| `ephem_tests.cpp` | 开普勒方程与传播、幽灵的外推（带状谐项、潮汐）与 Δv、星历表、目视/凌星轨道约定、N 体 |
| `planet_scene_tests.cpp` | 场景里的行星和卫星：太阳系、木星、地月、土星；Kerbol 系统对照游戏里的位置和自转角、6 小时太阳日、潮汐锁定、Laythe:Vall:Tylo ≈ 1:2:4、示意的 Mun 转移（停泊轨道、两次点火的 Δv、捕获轨道） |
| `mission_tests.cpp` | 探测器：先驱者/新视野、Parker、伽利略/朱诺、Artemis/CAPSTONE/嫦娥五号 T1 及其火箭末级、卡西尼/惠更斯、信使/贝皮科伦坡、火星的八个轨道器（含天问一号 2021 年的重建）、赛丁泉彗星与 3I/ATLAS、太阳系场景的 1I/3I、Europa Clipper/Hera 飞掠、金星快车/破晓号/Zoozve 与金星借力飞掠、ISEE-3、哈雷舰队、Rosetta/Philae、JWST |
| `shape_tests.cpp` | 小天体形状模型 |
| `appearance_tests.cpp` | 彗尾、尘埃彗尾、大气层、羽流（含各场景参数一致） |
| `stellar_system_tests.cpp` | α Cen、TRAPPIST-1、Kepler 系统、TIC 168789840、Algol（CHARA 位置与 B 的拉长、TESS 极小、GCVS 食深）与光变曲线模型、洛希面（Kim 1989 的半径、Eggleton 公式）、δ Cephei（SPIPS 半径、Engle 极大、视向速度轨道）、PSR B1620-26 |
| `relativity_tests.cpp` | Kerr 测地线（类时/类光）、后开普勒双脉冲星、Sgr A*、Gargantua（[JA15] 的盘和相机、最内稳定圆轨道上的示例行星、米勒星的时间膨胀、Mann 星的轨道对照 [TH14]、假想的 Edmunds 星看到的阴影大小和用光线追踪积分的盘的辐射） |
| `scene_tests.cpp` | 场景加载、Simulation、随天体自转的参考系、标签布局、跨场景一致性（同一天体在各场景的轨道参数）、相机导演、事件说明 |
| `settings_tests.cpp` | `config.toml` 读写与场景列表 |

## 资源与工具

- 星历烘焙：`python tools/bake/horizons_bake.py tools/bake/ephem.toml [name ...]`，联网访问 Horizons，输出到 `assets/ephem/`，说明见 `assets/ephem/SOURCES.md`。卫星相位拟合：`tools/bake/fit_moon_phase.py`
- `assets/stars/hyg.csv`：`tools/stars/convert_hyg.py`（HYG v4.4，带距离，供近邻系统重投影）；M4 星团天空：`tools/stars/make_m4.py`（Gaia DR3）
- `assets/sky/`：远处场景和银心的模型星空，`tools/sky/make_galaxy_sky.py`（需要 numpy；先 `calibrate` 对照 `milky_way.jpg` 标定，结果在 `tools/sky/galaxy_calibration.json`，再 `render <场景.toml>`，一个场景约 10 分钟，银心约 20 分钟）。改了模型或标定要重新生成全部场景。工具用 multiprocessing：Windows 上子进程重新导入模块，在主进程里改的模块全局变量传不过去，参数要经 `_init` 传
- `assets/belts/*.bin`：`tools/belts/make_belts.py`（JPL SBDB，主带 H < 15、全部 TNO）
- `assets/textures/*.jpg`：`tools/textures/prepare_textures.py` 从原图缩放；`uranus.jpg`、`neptune.jpg` 由 `tools/textures/make_ice_giant_textures.py` 从哈勃 OPAL 的 FITS 图合成（纬向剖面加一张真实图，需要 numpy）；Kerbol 场景 16 个天体的纹理（`kerbin.jpg` 等）由 `tools/textures/make_kerbol_textures.py` 程序生成（球面三维噪声、自带哈希，需要 numpy 和 Pillow，全部约 10 分钟）
- `tools/kerbol/`：Kerbol 场景的数据。`kerbol_dump.ks` 是 kOS 脚本（复制到 KSP 的 `Ships/Script/`，在带 kOS 处理器的飞船终端里 `SWITCH TO 0. RUNPATH("kerbol_dump").`），导出 `kerbol_bodies.txt`（已入库）；`kerbol_scene.py`（只用标准库）核对约定，打印场景的天体块、测试用的位置和自转角、大气参数；`kerbol_transfer.py`（只用标准库，借用 `horizons_bake.py` 的 Kepler 传播、抽稀和写表）按拼接圆锥曲线设计示意的 Kerbin→Mun 转移，烘成 `kerbol_mun_transfer.eph`。交会的最近距离由远拱点半径定，不是由点火相位定（飞船沿 Mun 的轨道方向到达、被 Mun 从后面追上，只改相位时都近乎迎头相撞）
- `assets/shapes/*.mesh`（AXMESH2）：`tools/shapes/` 下按来源格式分：`make_arrokoth.py`、`make_grid_table_shape.py`（Thomas / Stooke 经纬网格表）、`make_dtm_shape.py`（DLR 全球 DTM）、`make_plate_shape.py`（板块模型、OBJ、VRML、Gaskell 的 ICQ）
- `tools/bake/rejoin_burn.py`：把表里的入轨点火接到其后接缝处的下一段弧上（希望号，见 `assets/ephem/SOURCES.md`），烘焙之后运行
- `tools/isee3/`、`tools/armada/`、`tools/ce5t1/`：ISEE-3、哈雷舰队、嫦娥五号 T1 服务舱的轨道重建（需要 numpy，下载缓存不入库），做法见 `assets/ephem/SOURCES.md`
- `tools/tianwen1/`：天问一号环绕器。`momag_fetch.py` 从 USTC 的 MOMAG 磁强计数据下载位置（约 15 GB，每天一个文件，逐个下载后只留每分钟一个样本，`cache/` 不入库）；`tianwen1_reconstruct.py <jupyter_notebooks/Tianwen 目录>`（需要 numpy）把 D. Estévez 的 GMAT 脚本里的遥测向量外推、在机动处接起来，得到 2021-02-08 .. 11-16 的路径（`cache/tw1_2021.txt`；纯 Python 的 10 阶重力场积分半年的 8 小时轨道要二十多分钟，弧段缓存在 `cache/arcs/`）；`tianwen1_bake.py` 把它和 MOMAG 烘成 `tianwen1_mars.eph`，只用标准库
- `tools/philae/`：67P 的自转和 Philae（只用标准库，`daf.py`、`spk.py`、`ck.py` 读 SPICE 的 DAF、SPK 5/18 型段和 CK 5 型段）。`c67p_spin.py <目录>` 把 ESA 的三个彗核姿态 CK 拟合成 `c67p_spin.toml`（场景里 67P 的 `spin_table`）；`philae_bake.py <目录>` 读两个着陆器 SPK，用这张表换到 67P 本体系烘成 `philae_67p.eph`，所以改了自转表要重新烘 Philae。下载地址和做法见 `assets/ephem/SOURCES.md`
- `tools/pdftext.py`：只用标准库从 arXiv PDF 提取文字（ar5iv 经常超时）
- `tools/algol/algol_ephemeris.py <TESS *_lc.fits>`（需要 numpy，自带最小的 FITS 读取）：从 TESS 光变曲线量主极小、扣掉外轨道的光行时效应，拟合 `algol.toml` 内轨道用的历表；下载地址见 `assets/ephem/SOURCES.md`
- `tools/cepheid/spips_table.py <delcep2.dat> <pc> <out.csv>`（只用标准库）：把 SPIPS 模型（VizieR J/A+A/584/A80）转成脉动表，见 `assets/stars/SOURCES.md`
- `res/astraxis.ico`：`tools/icon/make_icon.py`（需要 Inkscape）从 `res/astraxis.svg`（≥ 48 px）和 `res/astraxis-small.svg`（≤ 40 px）生成；.ico 入库，CI 不依赖 Inkscape

## 已知的坑

各任务数据集自己的坑（SSCWeb、ICE 导航文件、卡西尼号最后一段等）写在 `assets/*/SOURCES.md` 里；实现上的取舍写在代码注释里。

### 构建、着色器与 SDL / ImGui

- VS 生成器会对加进 target 的 `.hlsl` 自动执行 FXC（报 X3721），必须设 `VS_TOOL_OVERRIDE "None"`（`HEADER_FILE_ONLY` 不起作用）
- CMakeCache 里的旧变量会盖过新逻辑：dxc 路径因此改名为 `ASTRAXIS_DXC_EXECUTABLE`
- GCC 的 `-Wmissing-field-initializers`：SDL 结构体用指定初始化器时漏写 `props` 也会报警告（MSVC 不报），要写 `.props = 0`
- HLSL 2021 不支持向量三元运算，用 `select()`；`m[i]` 是第 i 行，取列用 `mul(m, float4(1, 0, 0, 0))`；`fwidth` 不能放在非一致的分支里
- dxc 默认不是 IEEE 严格模式，`isnan`/`isinf` 可能被优化掉，排查 NaN 时要用 `asuint` 看指数位
- 体渲染里很薄的壳层（几 km 厚，而步长十几 km）会显出抖动噪声的规则斜纹：先把视线解析地裁剪到壳层内再采样
- 大的时钟值进着色器前要先 `frac` 再放大，否则 float 精度不够会闪烁
- 在片元着色器里按像素足迹解析地滤波细线（如天王星的窄环）：足迹用 `ddx_fine`/`ddy_fine`，`fwidth` 是 2×2 粗导数，低角度看的细线会变成 2 像素一节的珠串；滤波核只有 1 像素宽的盒子时边缘呈阶梯，用 2 像素宽的三角核（与纹理 mip 三线性过滤相当）
- ImGui 的 sdlgpu3 后端在 D3D12 上用 DXBC，建设备时要同时声明 `DXIL | DXBC`
- 不要用 `SDL_GenerateMipmapsForGPUTexture` 给星空立方体贴图生成 mip：各面是按 Direct3D 面布局、用镜像投影画的，SDL 生成的 mip 和本层对不上，天空在面的交界处出现亮度台阶。`StarfieldPass::render_cubemap` 直接画每一级 mip
- 不用平台后端时也要调用 `ImGui_ImplSDLGPU3_NewFrame()`（第一次调用时创建采样器）；漏了 Debug 能跑，Release 在 `SDL_BindGPUFragmentSamplers` 里崩
- ImGui 开了 `NavEnableKeyboard` 后，有窗口获得焦点时 `WantCaptureKeyboard` 就是 true，全局快捷键只看 `WantTextInput`
- 多输出时每个输出都要等自己的交换链（`SDL_WaitAndAcquireGPUSwapchainTexture`）
- 不要让 SDL 包装别的进程窗口下的子窗口（屏保 `/p` 预览）：对话框在第二块屏上时 SDL 会按顶层窗口的逻辑移动并放大它，预览区只剩黑底。预览窗口现在用 GDI 画
- TOML 的表头会“吞掉”后面的键：顶层键（如 `name`）必须写在第一个 `[table]` 之前
- 以行星为根、没有太阳天体的场景（木星、土星、火星）里，事件的 `from_body = "Sun"` 指太阳方向（`SceneEvent::from_sun`，用 `Scene::sun_position()`），不是天体；这类场景里不要再加名叫 Sun 的天体
- 撞击不能当飞掠事件（`flyby`）：`find_flyby` 要找前后都有数据的距离极小值，星历在撞击点结束就找不到（场景测试报“the named flyby is found”）。撞击事件用 `from_body` + `phase_deg` 定视角（撞在背面时 `from_body = "Earth"`、`phase_deg = 180`）
- 彗星可以是行星根天体的子天体（火星场景的赛丁泉），但尘埃彗尾要按太阳的 GM 传播（`Scene::sun_gm`），不能用父天体的 GM；彗星的星历要相对那颗行星烘焙，日心表减去近似根数的行星位置会差几千 km
- 改了 `assets/` 里的场景要重新构建才生效：程序读的是构建时复制到 exe 旁的副本

### 运行、截图与 shell

- 截图不要用 GDI+ 的 `Graphics.GetHdc()`（`CopyFromScreen` 内部也用它）：GDI+ 会把恰好等于 RGB(13, 11, 12) 的像素变成透明黑，在暗的光晕里看起来像一圈白色坏点。要把 PrintWindow 画进 `CreateCompatibleDC` + `CreateCompatibleBitmap`，`capture_window.ps1` 就是这么做的
- 截图脚本偶尔丢最后一个按键：需要切换状态的键放在前面，截完核对。等待别超过 60 秒，空闲 60 秒后自动导览会接管相机；`--event` 超出范围时只打日志、停在默认视图
- Release 版是 GUI 程序，后台进程无法把它切到前台；computer-use 的 `open_application` 会再开一个实例，要用 `SetForegroundWindow`/`ShowWindow`。computer-use 和 `SendKeys` 发的 Esc 进不了 SDL，用 `PostMessage(hwnd, WM_KEYDOWN, VK_ESCAPE)`
- PowerShell 的变量名不区分大小写：`$S` 和循环变量 `$s` 是同一个变量
- PowerShell 里的 `/tmp` 是 `C:	mp`，不是 Git Bash 的 `/tmp`（在用户的 Temp 目录下）：两边交换文件用 Windows 绝对路径，`-ErrorAction SilentlyContinue` 会把拷贝失败藏起来
- 截图要特定视角时，可以在 exe 旁的场景副本末尾临时追加 `[[events]]`（`focus` + `from_body` + `phase_deg` + `distance_km`），再用 `--event` 打开；用完从 `assets/` 拷回原文件
- 在后台跑烘焙工具并把输出重定向到文件时，Python 会缓冲 stdout，几十分钟看不到进度：用 `python -u`
- Claude Code 的沙箱里访问 Horizons 会反复报 `SSL: UNEXPECTED_EOF_WHILE_READING`（不是服务器的问题），烘焙要在沙箱外运行
- 从 Bash 用 `&` 启动的 exe 会随 shell 一起退出，长时间运行用 PowerShell 的 `Start-Process`；后台命令里 `cd dir && a & b &` 只有 a 在 dir 里执行，下载用 `curl -o 绝对路径`
- Bash 工具的 heredoc 即使写成 `<<'EOF'` 也会吃掉反斜杠（`\\n` 变成 `\n`）：含反斜杠的编辑用 Edit 工具，或先用 Write 写成脚本文件。heredoc 里的中文传给 `python -` 或 `cat >> 文件` 也会变成乱码（写成 GBK，不是 UTF-8），改中文文档同样用 Edit；几百行、夹着 C++ 代码的 heredoc 还可能整段报 `unexpected EOF while looking for matching`，批量改代码的脚本一律先 Write 到 scratchpad 再 `python <文件>`
- WSL 里用 Clang 再建一个构建目录时，`linux` 预设把 `FETCHCONTENT_BASE_DIR` 固定在 `build/linux/_deps`，两个编译器会共用依赖的构建目录（SDL 的预编译头互相覆盖，报 `cmake_pch.h has been modified`）。Clang 构建要另给 `-DFETCHCONTENT_BASE_DIR`，再用 `-DFETCHCONTENT_SOURCE_DIR_<SDL3|GLM|IMGUI|STB|TOMLPLUSPLUS|DXC>` 指向已下载的源码
- 用 `python -I` 跑读下载文件的工具时，脚本自己的目录不在 `sys.path` 上：同目录的模块（如 `tools/philae/spk.py`）要在脚本里显式 `sys.path.insert`。本机的 numpy 装在用户 site-packages 里，`-I`（含 `-s`）找不到它：要 numpy 的工具用 `python -E -P`
- Windows 上 Python 的控制台编码是 GBK，打印 PDF 提取的文字会 `UnicodeEncodeError`：加 `PYTHONIOENCODING=utf-8`（grep 这类文件要 `grep -a`）。`python -I` 会忽略 `PYTHON*` 环境变量，这时在脚本里 `sys.stdout.reconfigure(encoding='utf-8')`
- 标签、信息面板用 ImGui 的默认字体，只有 Latin-1：天体名里的 ʻ（U+02BB）等字符会显示成“?”，名字用 ASCII（`1I/'Oumuamua`，和 Horizons 一样）；展签字体 Jost 另有希腊字母
- 从 Git Bash 调 `wsl.exe`：`/mnt/c/...` 参数要加 `MSYS_NO_PATHCONV=1`；`bash -lc '...'` 里的 `$变量` 会被提前展开，脚本最好写成文件再执行

### 数据源

- **Horizons**：
  - 航天器轨道不全是重建值：看说明里的轨道段名，`*_PREDICT` 是预报
  - 用 curl 查向量要加 `REF_PLANE='FRAME'`，默认是黄道坐标
  - 覆盖范围以实际查询为准（说明里写的结束时间可能偏晚）；恰好结束在 2100-01-01 时，结束时间写成 2099-12-31
  - 日心坐标和以行星为中心的坐标可能不自洽（旅行者号、先驱者号早期的 rough 轨道在飞掠时差上万公里），比较时用日心查询相减
  - 小天体可能有两套轨道：`486958;`（地面定轨）和 `2486958`（任务组轨道），飞船轨道只和后者自洽；以彗星为中心的数据（如 `@1000012`）相对 ESA 自己的彗星解，场景里的彗星要用同一个解
  - 飞掠可能落在每日样本之间，留一法看不出来：用 `refine_near` 按距离强制加密
  - 日心坐标里有太阳的反射运动，容差太严会让节点处处很密（行星 50 km、航天器 20 km）；行星本体中心绕质心抖动，场景里的行星用系统质心
  - 短时点火期间的状态可能是假的（THEMIS-B 2008-10-18 有 3 分钟速度到 120 km/s，位置偏离平滑路径上千公里，只有到地心的距离看着平滑），会让插值溢出；烘焙工具按邻近速度的中位数剔除这种样本，用前后的节点跨过去。轨道文件接缝处的位置跳变（LCROSS 2009-08-25 达 656 km）不剔除，要让 1 分钟样本把它夹在两个节点之间
  - 留一法只在样本点上估误差，而且假设运动光滑：基础步长太粗时会漏掉快速的近心点（ARTEMIS-P2 20 分钟步长差 15 km），也看不见两个样本之间的机动或轨道文件接缝（ARTEMIS-P1 2012-07-05 一分钟内跳 44 km）。短的表直接用 1 分钟步长；长的表设 `verify_step_minutes`，用错开的样本复查超差处再加密
  - 万有变量 Kepler 方程的纯 Newton 迭代在偏心轨道上会在两个值之间来回跳（e≈0.83、约 0.2 周期时偏 2 万 km，只在个别时刻出现），`propagate_kepler` 和烘焙工具都用带区间的 Newton（f 单调，根在 0 和 √μ·dt/q 之间）
  - 小天体页头的根数是解的历元（如 3I/ATLAS 的 2026-02-19）的吻切值，带非引力力时 QR、Tp 和近日点当时的值不同（3I 的 q 差 5,000 km、Tp 差 0.013 天）：对照近日点要用 `EPHEM_TYPE=ELEMENTS` 在近日点附近另查一次
  - 观测表（`EPHEM_TYPE=OBSERVER`）的距离 `delta`（量 20）带光行时改正：天体取光离开时的位置，地球取观测时刻，和场景里的几何距离差得不少（2I 离地球最近时 1.93673 vs 1.93690 au，差 2.5 万 km）。对照几何距离用 `VECTORS` 加 `VEC_CORR='NONE'`
  - Kepler 相对插值的参考 GM 要准，GM 偏差会让节点暴增（卡戎用 DE440 的系统 GM 要 10390 个节点，拟合的等效值只要 277 个）；节点也可能隔好几圈，`history_times` 和烘焙工具都按平均运动细分
  - 火星低轨道每圈都要几个节点（扁率摄动），10 km 容差下 MRO 每月 43 KB，全程烘焙前先拿一个月试算大小；`encounters` 在 Python 里逐点求值，几十年的表按分钟步长扫要几十分钟
  - 绕火星的探测器吻切倾角会周期变化几度（MAVEN 在 73°–75.7° 之间），对照名义倾角时容差要放宽
- **按容差抽稀的表**只保证位置，速度和加速度在节点之间可能很差：天问一号 15 分钟的入轨点火原来整段落在两个节点之间，三次插值把推力曲线扭弯，按表求出的 Δv 少 2%，从点火开始时刻的状态外推的幽灵 84 分钟后偏 110 km。要从表里读机动的就让点火期间每个样本都成为节点（`ephem.toml` 的 `keep_all`；已有的表用 `horizons_bake.py --patch` 只改点火那一段，其他节点不动；天问一号是 `tianwen1_bake.py` 的 `KEEP_ALL`）。读数窗口要在点火前后各留一点
- 从表里读 Δv 要按节点到节点算（`unexplained_delta_v` 用 `history_times` 取节点）：节点是精确样本，三次插值在点火突然开始、结束的那一分钟里过冲，按 4 秒一步累加会把来回摆的两边都算上（天问一号多 1%）。ESA 的轨道还常在点火结束处换文件、位置对不上（火星快车 0.9 km、金星快车 3.8 km、TGO 99 km 和 50 m/s），按 4 秒累加时金星快车多 13%；窗口终点落进台阶那一分钟时连按节点算也不行（TGO 按 ESA 新闻稿的 13:05 开始会多 1.1 km/s），这时用更早的出处时间或在台阶前截止
- Horizons 的“tag-up fit”轨道（希望号、MOM）是一段段拟合的弧，带着入轨点火的那一段可能是计划中的点火，下一段才是实际到达的轨道：希望号点火 8 小时后两段相差 2,397 km、126 m/s，场景里轨迹跳一下。用 `tools/bake/rejoin_burn.py` 以后一段为准，重新拟合点火（先核对哪一段与公布的捕获轨道一致：MOM 点火那一段的 Δv 与 ISRO 公布的一致，没有改）。往回积分要带上前一段里二体解释不了的加速度（火星扁率），只按二体积分 8 小时差 2.5 km、2.3 m/s
- 在巨行星近心点附近从星历读推力，要先扣掉扁率：朱诺号入轨时木星 J2、J4 的加速度（0.3–0.45 m/s²）比发动机（0.26 m/s²）还大，只按点质量算会把它当成推力，Δv 多出几百 m/s。点质量的 GM 也要用行星本体的（系统 GM 把卫星算在里面，近心点在所有大卫星的轨道以内）
- Horizons 的 ARTEMIS（THEMIS-B/C）轨道在机动处是文件接缝：位置在一分钟里跳 200–380 km（P1 2011-06-27 16:32、P2 2011-07-18 00:26 TDB），整点附近还有一来一回的假样本，读不出推力，入月球轨道也就没做幽灵
- 新闻稿里点火的开始时刻口径不一：JPL 新闻资料（奥德赛号、MRO、朱诺号、卡西尼号）和 ESA 的新闻资料时间表（TGO）给地面接收时间，ESA 新闻稿、ISRO、APL、MBRSC、JAXA 给航天器上的时间（多数只到分钟）。先用 1 分钟的 Horizons 样本找出星历里推力开始的时刻（二体引力解释不了的加速度），再对照出处判断是哪一种；幽灵的起点要在星历的推力开始之前。主发动机前可能还有一段弱推力（火星快车点火前 2.5 分钟起 11.5 mm/s²，看来是沉底推进剂），1 分钟样本对照二体时淹没在扁率的残差里（近火点前每分钟几十 mm/s），要按 10 秒样本对照带扁率的引力找
- 幽灵在事件开始时和真实路径的偏差（测试里打印的 `info`）不全是外推的误差：真实路径来自按 5–10 km 容差抽稀的表，事件开始处节点稀，插值本身就偏 1–1.6 km（奥德赛号、MAVEN、TGO）。检验引力模型要对照 Horizons 的原始样本
- **天问一号**：
  - GMAT 脚本里的 `MarsInertial` 是 `BodyInertial` 轴（z 为 J2000 火星自转轴，x 为火星赤道与地球 J2000 平赤道的升交点），不是 ICRF 轴；火星极要用 GMAT 的 IAU 2000 值（317.68143°、52.88650°），用 pck00011 的极换算会差约 80 km。当成 ICRF 用时倾角会算成 74° 而不是 87°
  - MOMAG 2021 年的文件位置单位是火星半径（3389.5 km，三位小数，3.4 km 一级），之后是 km：单位要按文件头的 `Columns` 行判断；用相邻 1 Hz 样本差分求速度只会得到量化噪声，要用前后几分钟的最小二乘拟合
  - MOMAG 的位置是 MSO（x 指向太阳、z 沿火星公转角动量），换到惯性系时速度要加坐标系转动项 ω × r（ω = h / r²）
  - USTC 的下载接口不支持 Range、不压缩，偶尔返回截断的文件：按文件头检查，失败就重试
- **JPL 卫星平根数表**：
  - 各表的 P 不是同一种周期（JUP365 是平近点角周期，SAT441 是恒星周期），每张表都要对照 PCK 自转速率核实；进动周期只给绝对值，符号要自己判断（Ω+ω+M 的变化率要和同步自转速率一致）
  - 周期位数不够时（火卫、木星内侧卫星），几十年就失相：改由 PCK 的自转速率反推。拱线、交点周期也只有两三位（MAR099 火卫一的 1.1 年、2.3 年，平根数因此差 2°）：PCK 里摆动卫星极的角就是交点（`BODY4_NUT_PREC_ANGLES` 的 M1、M6），火卫一的 M5 是平近点角，拱线速率 = W′ − M5′ − 交点速率
  - 火卫一的受迫天平动（W 里的 −1.143° sin M5）和光学天平动（2e = 1.72°）反相，长轴偏离火星方向最大约 2.9°，不是两者相减
  - 天王星卫星的“equatorial”平面以天王星的正极为极，与 PCK 的北极相反；海卫一的交点周期约是实际值的一半；SAT441 的平根数相位对不上 SAT441 本身，土星卫星的相位用 `fit_moon_phase.py` 对 Horizons 拟合
  - 远处的不规则卫星（Phoebe）受太阳摄动，椭圆在一个行星年内被拉着转，平根数的进动椭圆跟不上（试过的一种进动方向组合对 Horizons 拟合后均方根差 11°，在表的两端交接会跳几百万 km）：烘长时间段的表，表外用 `extrapolate = "kepler"`。`fit_moon_phase.py` 每 2 天取一个样本，对一百多年的表一种符号组合要几十分钟
- 以行星本体中心（如 `500@699`）为中心的卫星表包含行星绕系统质心的晃动（土星被 Titan 带着每 16 天晃约 290 km）：外侧卫星的表因此每几天就要一个节点（Phoebe 128 年 1 万个），这是真实的相对位置，不要为省节点改用质心，否则和同样以本体中心为中心的探测器表对不上
- **ESA SPICE 的 Rosetta CK**：`ck/` 下的最新 CATT（彗核姿态）只覆盖 2015-12 以后，早的在 `former_versions/`，该用哪几个、谁优先看 `mk/ROS_OPS.TM`（后加载的优先）；同一文件里相邻段共用端点但值不同（导航解交接处跳到 0.5°），要按“后一段优先”去重，否则插值时除以零的时间差。CATT 的时钟是人造 SCLK（`ROS_V38.TF`），计数单位是 1/65536 s
- **IAU 本体坐标系**：W 从本体赤道与 ICRF 赤道的升交点（RA = α₀+90°）起量；极轴在 Dec=90° 时叉积退化，用 `iau_pole_frame(α₀, δ₀)`。pck00011 的火星极有一项周期约 7 万年的“章动”（`BODY499_NUT_PREC_RA/DEC` 的最后一项：RA 加 0.419057° sin(79.398797° + 0.5042615° T)，Dec 加 1.591274° cos(166.325722° + 0.5042615° T)），J2000 的极是 (317.681°, 52.886°)，只取常数项 (317.269°, 54.433°) 会偏 1.5°；W 也有对应的 0.584542° 项，漏掉时本初子午线偏 1–1.4°。场景里用 `nut_prec_*` 写了这三项
- **USGS 地图**：元数据写的“positive west”不可信，图像实际是东经向右递增。拼接图服务器跳转到 S3，用 `?list-type=2&prefix=mosaic/<名称>` 列文件
- **PDS Rings Node** 拒绝 Python urllib 的默认 User-Agent（403）
- **Stooke 形状模型 V2.0**：小行星和彗星的经度表改成了东经（360° 减 V1.0），但标签里的说明文字照抄 V1.0，其中的经度数值（如 Halley 的 270°）仍是西经。Vega 电视图像的时间有 UT（`TIM--OBS`、`TIME-ENC`）和开机后的星上计数（`TONBHRDW`，4.766485 计数/秒）两种，Stooke 标签里的“2:00:30”看来是后者
- **论文**：PDF 表格用 `pdftotext -layout` 提取时列会错位，用 `-raw`；`tools/pdftext.py` 会丢掉表格里的负号（Iess 2018 表 1 的 J4 提取成 586.609，实际是 −586.609），符号要核对；有的书的 PDF（《The Science of Interstellar》）用 `tools/pdftext.py` 只提取出空格，Git 自带的 `pdftotext`（`/mingw64/bin`）可以；本机没有 `pdftoppm`，Read 工具不能渲染 PDF 页面，看插图用 Inkscape：`inkscape --pdf-poppler --pages=<页> --export-type=png <pdf>`；凌星拟合的 ω 约定各代码不同（NbodyGradient 凌星时 ω+f = 270°，Phodymm 是 90° 且状态取反），场景里用 `transit_u_deg` 写明；共振链对积分相位误差很敏感，用 6 阶 Yoshida
- **轨道重建**：平动点轨道在刀刃上，按稳定流形二分出的机动方向只对同一条积分路径成立（换一种分段方式调用积分器，1e-9 的差别一个月后就是上万 km），求解和生成星历要走同一条路径。瞄准月球低空飞掠或入轨时，数值近月点的时刻在掠过月面附近会跳变，数“第 n 个近月点”也会整圈跳：先用十几小时前的密切二体轨道瞄准，再用各时刻的密切平近点角细调
- **共振轨道**：星下点每几圈重复的火星轨道（天问一号的 2 sol 停泊轨道、1/3 sol 中继轨道）受田谐项共振影响，几十天后的相位对初值极其敏感（02-23 的 Δv 变 1.5 m/s，80 天平均周期在 48.8–49.6 h 之间乱跳），不能把一次机动拟合到几十天后的向量上；用前后两个向量的弧段在中间交接
- **食双星的历元**：论文给的 T0 不一定是它说的那个时刻。Baron 2012 表 6 注明 Algol 的 T0 是“主极小”，实际是内轨道过升交点（主极小在其后 3/4 周期）。圆轨道的目视轨道 T0 常是过节点，要用观测的极小时刻（GCVS、TESS）和他们自己的位置表核对；几十年的外推还要考虑 O−C 漂移，用最近的观测拟合（`tools/algol/`）
- **论文里的数自相矛盾**：摘要和表格可能不一致（Anderson 2015 摘要的 e = 0.647，表 3 是 0.674）。用表里别的量反算核对（a₁ sin i 只和 0.674 相符），以表为准
- **半径的口径**：充满洛希瓣的星，论文给的“半径”可能是某个方向的（Kolbas 2015 的 Algol B 3.43 R☉ 是侧向半径，按他们自己的 q 和 a 算），不是等体积半径；由 v sin i 加同步自转假设推出的半径是上限（Algol A 的 2.87 R☉）。拿来和别的来源的质量、倾角混用前，先用洛希面换算到场景自己的 a。光变解（Wilson-Devinney）的倾角和半径是一起拟合的，单独换掉一个（例如用干涉测量的倾角）会让食深偏离很多
- **ADS 的老文章扫描**：`articles.adsabs.harvard.edu/pdf/<bibcode>` 用 curl 只能拿到第一页，其余页要令牌。在内置浏览器打开 `ui.adsabs.harvard.edu/scan/manifest/<bibcode>`，在页面里用 `/v1/accounts/bootstrap` 返回的 `access_token` 取 `/v1/scan/image/iiif/2/bitmaps-~seri-~<刊>-~<卷>-~600-~<页>.000/full/1600,/0/default.jpg`（如 `ApJ__`、`0342`、`0001065`），放进页面后截图读
- **TESS 上没有亮星的光变曲线**：TIC 里亮星的 Tmag 可能是错的（δ Cep 的 TIC 415741355 写成 14.0 等），没有任何管线给它做光变曲线；MAST 按位置搜到的同名附近星（HD 213307）不是它，先按 HIP 查 TIC 号
- **星表距离和场景距离不一致**：场景采用的视差（例如轨道视差）和 Hipparcos 不同时，星表里的“自己”会离观察者超过 `kOwnSystemPc`，在天空里成为一颗极亮的星（Algol 差 1.25 pc，−4.6 等）：用 `[sky] hide_hip` 去掉
- **KSP / kOS**：kOS 脚本的变量和参数名不能遮住内置函数（`v`、`r`、`q` 等，报 CLOBBERBUILTINS）。kOS 的向量在 Unity 的左手系里，y 轴朝系统北极，原点在飞船，低空时坐标轴还会转：换成场景坐标要 y、z 互换，x 轴取同一物理帧里采样的 `SOLARPRIMEVECTOR`。天体的 `VELOCITY:ORBIT` 是相对飞船所在 SOI 的速度，不是相对母天体，要和母天体的相减。游戏里的根数多是 float32（e = 0.2 读出来是 0.200000002980232，历元平近点角的配置值是 3.14 rad 而不是 π），照抄读出的值。历法常数在 `Assembly-CSharp.dll` 的 KerbinDay/KerbinYear getter 里（`ldc.i4 21600; ret` 和 `ldc.i4 9201600; ret`）
- **时间和距离**：MJD 转 JD 要加 2400000.5；新闻稿里的时间是地面接收时间（差光行时），对照测试只比较距离；NASA 的“离地球多远”和飞掠高度从地表算；地面成像给出的环弧、卫星经度是光离开天体时的位置，对应的时刻要从观测 UT 里减去光行时（海王星约 4 h，Adams 环弧段转 137°），论文按 UT 算的时间间隔还会少掉其间的闰秒
- **借力飞掠的“速度改变”**：新闻稿的口径不一。卡西尼号 1998 年的“about 7 km/s”是相对金星速度的转折 |v_out − v_in|，伽利略号、信使号的是日心轨道能量的变化折合成金星距离上的速度（信使号 2007 年转折 6.95 km/s，能量口径 6.64 km/s，即“nearly 15,000 mph”）；对照前两种都算一下

### 测试

- 场景测试里查光照或位置前要先 `scene.update(t)`，否则所有天体都在原点
- GR 的近心点方向要对 r(φ) 做抛物线拟合（节点角间距和每圈进动同一量级）；IMR 是辛积分器，能量误差有界但不精确守恒
- 黑洞光线是从相机往回追踪的（时间反演），测试里自旋不对称性的符号要反过来看
- Kerr-Schild 坐标里协变动量 p_i 不是光线方向（dx^i = p_i − f (l·p) l_i，150 M 处相差约 2/r），屏幕方向也不是角度，要经观者标架（`static_observer_frame`）。a = 0 时它的空间坐标就是史瓦西的面积半径坐标，弱场“从 R 到无穷远偏折 (2/b)(1 − cos)”是各向同性坐标里的分法，拿来接积分段会在 b = R 处错开约 0.013 rad；用 `weak_leg_angle`
- 光变曲线在食外也不为 0（洛希形状和照射带来的起伏）：判断是否在食中，要比较有遮挡和无遮挡时的光；食深从最亮处量起，不要从 0 量
- 导演测试用多个随机种子长时间跑，抓相机穿过天体、掠过黑洞这类问题；带滚动窗口积分的场景要用低倍速，否则 Debug 下很慢

## 待办

按容易程度排：

1. Dear ImGui 等 1.93.x 补丁版出来再从 v1.92.9b 升级（1.93.0 重写了 ImDrawList，刚发布）。SDL3/SDLGPU3 后端和用到的 API 都兼容，要截图核对的是：`label_overlay.cpp` 的空心圆点（`AddCircle` 半径不再减 0.5、默认内侧描边）和标记菱形（`AddQuad`）、`captions.cpp`/`app_ui.cpp` 里非整数坐标的 `AddRectFilled`（现在有抗锯齿边缘）、`StyleColorsDark()` 改过的面板配色
2. B1913+16 自转轴的测地线进动：用 Weisberg & Taylor 2002（astro-ph/0205280）表 1 的解 1（i = 47.2°，与场景一致；δ = 21.1°、α = 157.5°、T0 = 1981.0）；Kramer 1998 表 1 的解 I 是同一几何（δ = 22°，α = 27° 即自转反向）。速率由质量按 Barker & O'Connell 算（1.215°/年，对照 WT02 的 1.213°/年写测试），自转轴绕轨道法线走半角 δ 的锥；画面上只有 ζ(t)（自转轴与视线的夹角）可见，它只依赖 cos Φ，不受进动方向和假设的天空节点影响，测试对照 WT02 的公式。`pulsar_beams` 要多传仿真时间。再加一个从地球方向（视线）看的事件，事件视角目前只有 `from_orbit_normal`、`from_body`，要新增一种：现有的“光束扫过相机时闪一下”会一直闪到约 2025 年（Kramer 的预测）后停下；停下的时刻取决于光束宽度（现在是示意的 10°，WT02 的光束是轴比 1.59 的沙漏形），要调或注明。2019 年以后没找到观测文献，展签写成“预测”。[WH16] 4.4 节说各家几何解“不一致”，是指精度不够扣除光行差，不影响画面
3. 变星：Polaris（复用 `[bodies.pulsation]`：振幅很小，周期每年变长约 4.5 s；Aa–Ab 约 29.6 年的轨道，Evans 等的 HST/CHARA 解；远处的 B）。要找有相位表的脉动模型（像 SPIPS 那样），或从视向速度积分半径（p 因子要注明）
4. Kerbol 场景：高度图和法线（晨昏线上的立体感，也可用于月球、火星的 USGS DEM）与海面高光（Kerbin、Eve、Laythe、地球）；还可以用 `kerbol_transfer.py` 的办法加别的示意航线（Minmus、Duna 的转移窗口等，标明是设计的）
5. 黑洞附近的天体做透镜：Gargantua 的示例行星和 Mann 星现在按真实位置画，转到黑洞背后时被挡住（真实情况下光会弯到阴影上方）。可以对每个天体解相机到天体的类光测地线，把球和标记画在主像（和次像）的位置
6. 远期：macOS（Metal，需要 SPIR-V → MSL 和 macOS CI）；J1141-6545 的自转轴进动（Venkatraman Krishnan 2019 表 1：α ≈ 89°，δ 各模型从 35° 到 155° 不等，倾向的模型 3 是 91° ± 60°；Manchester 2010 约 110°），等几何解定下来再做

## 不要做

- 不要手写多套图形后端（GL/DX/VK）—— 只用 SDL_GPU
- 不要为了“通用”引入 ECS 或复杂的引擎层；保持小而直接
- 不要让 `core/`、`ephem/` 依赖 SDL 或 GPU
