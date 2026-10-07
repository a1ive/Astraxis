# Astraxis

用来“发呆”的天文动态背景程序：太阳系与各探测器轨迹、JWST L2 晕轨道、木星和土星系统、彗星与哈雷舰队、Alpha Centauri、系外行星共振链、多星系统、脉冲星双星、Sgr A* 黑洞外观等，共 19 个场景。
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
src/render/    SDL_GPU 封装与各渲染 pass（星空、天体、轨道线、标记）
src/scene/     天体、参考系、相机、场景定义；Simulation（场景 + 时钟 + 相机 + 导览）
src/ephem/     运动来源：解析开普勒 / 星历插值 / N 体（Yoshida）/ GR（Kerr 测地线，IMR）
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

核心抽象：所有天体运动都实现 `MotionSource::eval(t_tdb) -> State`，返回相对父天体的 ICRF 状态，场景层不关心来源是解析、星历还是积分。`Scene::update` 先按父子链算出 ICRF 位置（`icrf_position`），再经当前显示参考系（`FrameTransform`）变换成 `world_position`。渲染、相机、标签、导演都只用 `world_position`。

渲染层的 pass 只接收 `CameraView` 和简单的 draw item，由 `src/view/` 的 `SceneRenderer` 把 scene 转换成这些数据。宿主无关的部分：
- `Simulation`（`scene/simulation.hpp`，不依赖 SDL）持有 Scene、时钟、相机和导览；`compute_view` 给每个输出算一份 `OutputView`（相机视图、卫星淡出、content scale）
- `LabelLayout`（`scene/label_layout.hpp`）算标签和标记放在哪里（优先级、避让、渐隐、拾取），不负责绘制；`App` 用 ImGui 画
- `SceneRenderer` 持有 pass 和每个场景的 GPU 资源，显示选项通过 `ViewOptions` 传入
- GPU 分两层：`GpuDevice` 只有一个；`RenderOutput` 每个窗口一个，持有交换链和本输出的 `SceneTargets`。`Frame` 的 `width/height` 是场景尺寸，`output_width/output_height` 是交换链尺寸，两者在有渲染缩放（壁纸、屏保的 `render_scale`）时不同。所有和输出尺寸相关的纹理（MSAA/HDR/深度、泛光链、黑洞追踪目标）都放在 `SceneTargets` 里，pass 内部不要再缓存按尺寸分配的纹理，否则多个不同尺寸的输出轮流渲染时会反复重建。各输出都用 SDR 交换链，格式相同，共用一个 `SceneRenderer`

宿主（`App`、屏保、壁纸）只负责窗口、输入、UI 和闲置策略，不要把仿真或渲染逻辑放回 `src/app/`。宿主每帧只 `sim.update` 一次，再对每个输出调用 `SceneRenderer::render`；与相机无关的工作（例如尘埃彗尾的传播）不要每个输出重做一遍。

渲染流程：
1. 黑洞光追（有黑洞时）：在屏幕外以半分辨率追踪
2. 场景 pass，渲染到 HDR（RGBA16F + MSAA）：银河 → 星空 → 天体 → 大气层 → 羽流 → 彗发/离子彗尾 → 尘埃彗尾 → 行星环 → 黑洞合成 → 太阳 → 脉冲星光束 → 轨道线
3. `PostProcess`：泛光 → ACES → sRGB → 交换链
4. ImGui 叠加 pass

场景用 `assets/scenes/*.toml` 描述（格式见 `jupiter.toml` 的注释），新增场景不用改代码。

### 屏保、壁纸与设置

- 所有模式都由 astraxis.exe 在本进程内渲染。`astraxis.scr` **就是 astraxis.exe 本身**，不是启动器：`/s` 下 Windows 把屏保进程退出当作屏保结束（锁屏也跟着它）；`/p <HWND>` 如果跨进程建子窗口，两个线程的输入队列会关联起来，容易卡住或留下孤儿进程
- `astraxis_wallpaper.exe` 是 Win32 设置对话框，选好场景、显示器、帧率等以后拉起 `astraxis.exe --mode wallpaper`；屏保的 `/c` 也用这个对话框。开机自启就是一行 `astraxis.exe --mode wallpaper`
- 设置放在 exe 旁边的 `config.toml`（便携），目录不可写时不保存。`[view]` 是各模式共用的显示选项，`[window]`、`[wallpaper]`、`[screensaver]` 放各模式自己的项，也可以覆盖 `[view]`。优先级：代码默认值 < `[view]` < 模式那一节 < 命令行
- 壁纸要同时支持两种桌面层级：Win11 24H2 起 WorkerW 是 Progman 的子窗口；旧层级要先发 `0x052C` 让 Explorer 拆出 WorkerW。Explorer 重启（`TaskbarCreated`）后重新挂载，退出时用 `SPI_SETDESKWALLPAPER` 刷出原壁纸
- 功耗：壁纸被最大化或全屏窗口盖住、锁屏、显示器关闭时暂停；用电池时按 `battery` 设置限帧或暂停；`render_scale` 降低内部分辨率（黑洞追踪跟着一起降）

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
| `core_tests.cpp` | 日历/时间尺度、星表与黑体颜色、银道坐标系 |
| `ephem_tests.cpp` | 开普勒方程与传播、星历表、目视/凌星轨道约定、N 体 |
| `planet_scene_tests.cpp` | 场景里的行星和卫星：太阳系、木星、地月、土星 |
| `mission_tests.cpp` | 探测器：先驱者/新视野、Parker、伽利略/朱诺、Artemis/CAPSTONE、卡西尼/惠更斯、信使/贝皮科伦坡、ISEE-3、哈雷舰队、Rosetta、JWST |
| `shape_tests.cpp` | 小天体形状模型 |
| `appearance_tests.cpp` | 彗尾、尘埃彗尾、大气层、羽流（含各场景参数一致） |
| `stellar_system_tests.cpp` | α Cen、TRAPPIST-1、Kepler 系统、TIC 168789840、PSR B1620-26 |
| `relativity_tests.cpp` | Kerr 测地线（类时/类光）、后开普勒双脉冲星、Sgr A* |
| `scene_tests.cpp` | 场景加载、Simulation、标签布局、跨场景一致性（同一天体在各场景的轨道参数）、相机导演、事件说明 |
| `settings_tests.cpp` | `config.toml` 读写与场景列表 |

## 资源与工具

- 星历烘焙：`python tools/bake/horizons_bake.py tools/bake/ephem.toml [name ...]`，联网访问 Horizons，输出到 `assets/ephem/`，说明见 `assets/ephem/SOURCES.md`。卫星相位拟合：`tools/bake/fit_moon_phase.py`
- `assets/stars/bsc5.csv`：`tools/stars/convert_bsc5.py`（CDS V/50）；M4 星团天空：`tools/stars/make_m4.py`（Gaia DR3）
- `assets/belts/*.bin`：`tools/belts/make_belts.py`（JPL SBDB，主带 H < 15、全部 TNO）
- `assets/textures/*.jpg`：`tools/textures/prepare_textures.py` 从原图缩放
- `assets/shapes/*.mesh`（AXMESH2）：`tools/shapes/` 下按来源格式分：`make_arrokoth.py`、`make_grid_table_shape.py`（Thomas / Stooke 经纬网格表）、`make_dtm_shape.py`（DLR 全球 DTM）、`make_plate_shape.py`（板块模型、OBJ、VRML）
- `tools/isee3/`、`tools/armada/`：ISEE-3 和哈雷舰队的轨道重建（需要 numpy，下载缓存不入库），做法见 `assets/ephem/SOURCES.md`
- `tools/pdftext.py`：只用标准库从 arXiv PDF 提取文字（ar5iv 经常超时）
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
- ImGui 的 sdlgpu3 后端在 D3D12 上用 DXBC，建设备时要同时声明 `DXIL | DXBC`
- 不用平台后端时也要调用 `ImGui_ImplSDLGPU3_NewFrame()`（第一次调用时创建采样器）；漏了 Debug 能跑，Release 在 `SDL_BindGPUFragmentSamplers` 里崩
- ImGui 开了 `NavEnableKeyboard` 后，有窗口获得焦点时 `WantCaptureKeyboard` 就是 true，全局快捷键只看 `WantTextInput`
- 多输出时不要每个输出都阻塞等一次 VSYNC，用不阻塞的 `SDL_AcquireGPUSwapchainTexture`
- 不要让 SDL 包装别的进程窗口下的子窗口（屏保 `/p` 预览）：对话框在第二块屏上时 SDL 会按顶层窗口的逻辑移动并放大它，预览区只剩黑底。预览窗口现在用 GDI 画
- TOML 的表头会“吞掉”后面的键：顶层键（如 `name`）必须写在第一个 `[table]` 之前
- 改了 `assets/` 里的场景要重新构建才生效：程序读的是构建时复制到 exe 旁的副本

### 运行、截图与 shell

- 截图不要用 GDI+ 的 `Graphics.GetHdc()`（`CopyFromScreen` 内部也用它）：GDI+ 会把恰好等于 RGB(13, 11, 12) 的像素变成透明黑，在暗的光晕里看起来像一圈白色坏点。要把 PrintWindow 画进 `CreateCompatibleDC` + `CreateCompatibleBitmap`，`capture_window.ps1` 就是这么做的
- 截图脚本偶尔丢最后一个按键：需要切换状态的键放在前面，截完核对。等待别超过 60 秒，空闲 60 秒后自动导览会接管相机；`--event` 超出范围时只打日志、停在默认视图
- Release 版是 GUI 程序，后台进程无法把它切到前台；computer-use 的 `open_application` 会再开一个实例，要用 `SetForegroundWindow`/`ShowWindow`。computer-use 和 `SendKeys` 发的 Esc 进不了 SDL，用 `PostMessage(hwnd, WM_KEYDOWN, VK_ESCAPE)`
- 从 Bash 用 `&` 启动的 exe 会随 shell 一起退出，长时间运行用 PowerShell 的 `Start-Process`；后台命令里 `cd dir && a & b &` 只有 a 在 dir 里执行，下载用 `curl -o 绝对路径`
- Bash 工具的 heredoc 即使写成 `<<'EOF'` 也会吃掉反斜杠（`\\n` 变成 `\n`）：含反斜杠的编辑用 Edit 工具，或先用 Write 写成脚本文件
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
  - Kepler 相对插值的参考 GM 要准，GM 偏差会让节点暴增（卡戎用 DE440 的系统 GM 要 10390 个节点，拟合的等效值只要 277 个）；节点也可能隔好几圈，`history_times` 和烘焙工具都按平均运动细分
- **JPL 卫星平根数表**：
  - 各表的 P 不是同一种周期（JUP365 是平近点角周期，SAT441 是恒星周期），每张表都要对照 PCK 自转速率核实；进动周期只给绝对值，符号要自己判断（Ω+ω+M 的变化率要和同步自转速率一致）
  - 周期位数不够时（火卫、木星内侧卫星），几十年就失相：改由 PCK 的自转速率反推
  - 天王星卫星的“equatorial”平面以天王星的正极为极，与 PCK 的北极相反；海卫一的交点周期约是实际值的一半；SAT441 的平根数相位对不上 SAT441 本身，土星卫星的相位用 `fit_moon_phase.py` 对 Horizons 拟合
- **IAU 本体坐标系**：W 从本体赤道与 ICRF 赤道的升交点（RA = α₀+90°）起量；极轴在 Dec=90° 时叉积退化，用 `iau_pole_frame(α₀, δ₀)`
- **USGS 地图**：元数据写的“positive west”不可信，图像实际是东经向右递增。拼接图服务器跳转到 S3，用 `?list-type=2&prefix=mosaic/<名称>` 列文件
- **PDS Rings Node** 拒绝 Python urllib 的默认 User-Agent（403）
- **论文**：PDF 表格用 `pdftotext -layout` 提取时列会错位，用 `-raw`；凌星拟合的 ω 约定各代码不同（NbodyGradient 凌星时 ω+f = 270°，Phodymm 是 90° 且状态取反），场景里用 `transit_u_deg` 写明；共振链对积分相位误差很敏感，用 6 阶 Yoshida
- **时间和距离**：MJD 转 JD 要加 2400000.5；新闻稿里的时间是地面接收时间（差光行时），对照测试只比较距离；NASA 的“离地球多远”和飞掠高度从地表算

### 测试

- 场景测试里查光照或位置前要先 `scene.update(t)`，否则所有天体都在原点
- GR 的近心点方向要对 r(φ) 做抛物线拟合（节点角间距和每圈进动同一量级）；IMR 是辛积分器，能量误差有界但不精确守恒
- 黑洞光线是从相机往回追踪的（时间反演），测试里自旋不对称性的符号要反过来看
- 导演测试用多个随机种子长时间跑，抓相机穿过天体、掠过黑洞这类问题；带滚动窗口积分的场景要用低倍速，否则 Debug 下很慢

## 待办

按容易程度排：

1. 按事件设轨迹长度（事件覆盖天体的轨迹设置）：帕克场景金星视角里此前所有飞掠的轨迹都穿过金星
2. 导览按事件走：飞掠前自动降倍速、飞掠特写镜头、显示展签；让屏保和壁纸里也出现展签（现在只有窗口版有）
3. 标题字体加希腊字母后备（Jost 没有希腊字母，`α` 显示成 `?`）；面板打开时长标题两端会被挡住
4. 随天体自转的参考系；Rosetta 场景加 Philae（ESA SPICE 里有着陆轨迹，Horizons 没有）
5. 远期：macOS（Metal，需要 SPIR-V → MSL 和 macOS CI）；Hulse–Taylor 并合（3 亿年后 double 秒数只有约 2 s 分辨率，要做成以并合为零点的单独场景）；脉冲星自转轴的测地线进动（B1913+16 的几何解各论文不一致）

## 不要做

- 不要手写多套图形后端（GL/DX/VK）—— 只用 SDL_GPU
- 不要为了“通用”引入 ECS 或复杂的引擎层；保持小而直接
- 不要让 `core/`、`ephem/` 依赖 SDL 或 GPU
