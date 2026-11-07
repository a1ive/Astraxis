# Astraxis

用来“发呆”的天文动态背景程序：太阳系与各探测器轨迹、地球夜侧的城市灯光和极光（按 1932 年以来的 Kp）、空间望远镜（近地轨道上的 Hubble，高椭圆轨道上的 Chandra，L2 上的 WMAP、Herschel、Planck、JWST、Euclid、Roman、Gaia，与月球 2:1 共振的 TESS，落在地球后面的 Spitzer 和 Kepler）、金星、火星、木星和土星系统、彗星与哈雷舰队、Shoemaker-Levy 9（木星场景里的解体、“一串珍珠”和 21 个碎片的撞击，按 JPL 的最终解，有尘埃尾、闪光、弹道羽流、环和扩散成一条带的疤痕）、DART 撞击 Dimorphos（幽灵是没被撞的轨道）、Alpha Centauri、β Pictoris（第一个被拍到的碎屑盘，侧看：倾斜的副盘、行星 b 和 c、毫米波的星子带、随时间转的 CO 团、JWST 的“猫尾”和星际介质吹弯的尾尖、TESS 凌星日期上的示意系外彗星）、PDS 70（不透光的原行星盘，Mie 的颗粒和多次散射，缺口里的行星 b、c 和 c 的环行星盘，CO 视图）、HD 206893（ALMA 看到的带缺口的碎屑带，里面的褐矮星 B、接近氘燃烧极限的 c 和暂定的卫星）、系外行星共振链、多星系统、食变星 Algol 和造父变星 δ Cephei、北极星（带从地球看的光变曲线）、WR 104 的尘埃风车（几何模型，倾斜 42° 的轨道从地球看仍像正对的螺线）和 WR 140 每次近星点产生的尘埃壳（JWST 看到的 17 层）、η Carinae 的 Homunculus（ESO 的三维模型，1847 年起均匀膨胀，有红外视图的两层尘埃；里面是 1890 年的 Little Homunculus、5.54 年的大质量双星和每次近星点留下的化石星风壳，星和星云的亮度跟着 1820 年以来的光变曲线）、脉冲星双星、Sgr A* 黑洞外观、《星际穿越》的虚构黑洞 Gargantua（按论文参数，电影版与带多普勒频移的真实样子对比；还有致敬的 Endurance 号，环在转）、致敬《坎巴拉太空计划》的 Kerbol 系统（游戏里读出的轨道、Kerbin 历法、程序生成的表面）等，共 34 个场景（含 2017 / 2024 太空视点日食和 2022 月球附近视点的地球遮日）。
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
src/ephem/     运动来源：解析开普勒 / 星历插值 / 根数表（进动椭圆的根数随时间插值，Janus、Epimetheus）/ N 体（Yoshida）/ GR（Kerr 测地线，IMR；赤道圆轨道用解析解 kerr_circular，任意倾角的轨道用 kerr_apoapsis 在 Kerr-Schild 坐标里 3D 积分）
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
- `Simulation`（`scene/simulation.hpp`，不依赖 SDL）持有 Scene、时钟、相机和导览；`compute_view` 给每个输出算一份 `OutputView`（相机视图、卫星淡出、content scale）。导览的普通镜头由 `CameraDirector` 选；每隔一两个镜头插一个事件镜头，因为要动时钟，由 `Simulation` 执行：按时间顺序取下一个事件，淡出黑场后把时钟跳到事件前（快到的事件不跳），倍速按剩余时间指数降到事件自己的倍速（飞掠降到最近点前 20 s），显示展签，飞掠镜头持续到最近点之后。飞掠体由事件的 `flyby` 指定（自动找最近的探测器在多探测器场景里常挑错，比如“MRO arrives”找到 Mars Express），`--event` 和导览用同一个视角（`Simulation::event_shot`）。远处场景的事件可设 `from_earth = true`，从 `sky.viewer_pc` 指向太阳系的方向近看目标（需要惯性参考系，导览时视线不漂移）；要像天文图像那样北在上，让事件切到一个带 `up = "north"` 的惯性参考系（`[[frames]]` 的 `up`，同 `[view]` 的写法，η Car 场景）。事件可以带幽灵（`[events.ghost]`，格式见 `mars.toml` 的天问一号入轨）：加载时合成一个 `Ghost` 天体（`Body::event`），从点火前的状态外推（`CoastMotion`，`ephem/coast.hpp`）：中心天体有 `zonal_harmonics` 时计入扁率（木星、土星），`tides` 列出的天体计入潮汐（月球场景的地球），否则是二体；有摄动时在加载时积分成节点表。只在该事件激活时显示，轨迹整段画出（普通天体也可以只在某些事件里显示：事件的 `show` 列出天体名，如 α Cen 场景里已被否定的 Bb；几个事件都可以 `show` 同一个天体，如 Hubble）；探测器是太阳的子天体、只是星历相对行星时（信使号、破晓号），用 `center` 指定绕哪个天体外推。幽灵也可以是卫星（Didymos 场景里没被 DART 撞的 Dimorphos）：外推的 GM 加上它自己的 `body_gm_km3_s2`，画成和它一样大的空心圆。点火的 Δv（`unexplained_delta_v`，按星历节点到节点算）由 `Scene::burn_delta_v` 给出，画在探测器标签下，总量不到 1 m/s 时用 mm/s。事件还可以用 `elevation_deg` 指定相机在 up 轴法平面上方的仰角，用 `disk_shifts = false` 在事件激活时关掉吸积盘的多普勒和引力频移（Gargantua 的电影版），用 `dust_visible = true` 把风车星云的尘埃从红外伪彩色换成肉眼看到的颜色，用 `nebula_infrared = true` 把膨胀星云换成热红外的伪彩色（只画有红外层的星云）。黑场用 `Simulation::fade()`，在后处理合成和标签里乘上。旋转参考系里星空跟着转：按真实时间每秒转过 45° 以上开始淡出、180° 以上隐去（`Simulation::sky_visibility()`，随 Janus 转的参考系在 2 天/秒时每秒转 1000° 以上）。天体自转可以写 `synchronous = true`（长轴始终指向母天体，Janus、Epimetheus）。探测器有 `shape` 时按网格画（否则只是标记），导演的跟拍镜头也拉近到能看清它（距离是 `radii_km` 的 4–7 倍，长吊杆不要算进半径）。NASA 的模型（`tools/shapes/make_spacecraft.py`，Voyager、New Horizons、Pioneer 10、Cassini、Galileo、Juno、Hubble、JWST、Parker）带纹理图集（网格第 2 位标志，`ShapeModel::map_v`，纯色材质在图集的色板格里）和烘焙的阳光（顶点“反照率”通道，太阳在 +z 时是否被自身挡住），场景写 `shape`、`texture`、`point_at`（+z 指向某个天体，+y 靠近黄道北，示意的姿态，“Sun”在没有太阳天体的场景里也行）和 `fill_light`（背光面的示意补光，大致从相机方向来）；`[bodies.rotor]` 是同一位置再绕本体 z 轴转的第二个网格（Endurance 的环），转角由 `Simulation::rotor_turns` 按倍速累加，画面上最快 0.25 圈/秒（和吸积盘动画限速同一个思路）
- `LabelLayout`（`scene/label_layout.hpp`）算标签和标记放在哪里（优先级、避让、渐隐、拾取），不负责绘制；`App` 用 ImGui 画
- 光变曲线（场景的 `[light_curve]`，`scene/light_curve.hpp`）：从某个天体方向看，各恒星按黑体圆盘、线性临边昏暗、前面的星挡住后面的星算总光度。`Simulation` 在曲线可见时每帧按固定时间网格更新一次采样（已算的保留，高倍速完全隐藏时停止更新，恢复可见时补齐），各宿主和标签一起画在角落（`app/light_curve_overlay.*`），倍速快到几秒就换掉整条曲线时淡出。星等从“每颗星最亮时、不遮挡”的固定参考光度量起。脉动星（`[bodies.pulsation]`，`Pulsation`）按相位表查半径、温度、光度，`Scene::update` 每帧写回天体；光变曲线按任意时刻查表（`star_disk_at`）。充满洛希瓣的星（`[bodies.roche]`，`scene/roche.hpp`，Algol B）加载时生成临界等势面的网格，逐顶点温度带引力昏暗和伴星的照射；它随轨道转（`Body::synchronous_with`），按形状模型画（星的着色乘上顶点亮度），光变曲线对它按面元求和，它挡在别的星前面时用网格投影的凸包
- 观测的光变曲线（场景的 `[observed_light_curve]`，`ObservedLightCurve`，η Car 从 1820 年到现在的 V 星等）：不是模型，是归并好的观测表（`year,v,...`），按 `from`/`to` 画出整段历史，标出时钟所在的年份和星等（`observed_magnitude` 插值，相距超过 `gap_years` 的点不连、不插值）。有它时角落的光变曲线图换成它（`draw_light_curve` 里分派），不随倍速淡出。可以用 `stars` 和 `reference_time` 让这几颗星、以及它们照亮的星云和风车星云按观测的亮度画（`observed_magnitude_offset`：观测星等减参考时刻的星等，空缺连起来、两端保持，恒星换算为按 0.56 次方压缩的限幅后显示增益，星云和风车星云的光乘上对应倍数再按各自的 `display_gamma` 压缩）
- 风车星云（场景的 `[pinwheel]`，`scene/pinwheel.hpp`，WR 104、WR 140）：White 2025 和 Han 2022 的几何模型，尘埃在锥面上成环，从双星质心飞出（可以按 Han 2022 式 1 被辐射压加速，`pinwheel_radius`）。偏心轨道可以只在一段真近点角里产尘（环按真近点角均匀分布，每环的权重含产尘率和 dt/dν），产尘率在轨道上和锥的方位上可以变化。尘埃是加载时生成的一组固定颗粒（方向、发射相位、圈号），轨道每个周期重复，所以位置只是时间的函数（`pinwheel_grain_offset`）。`PinwheelPass` 把颗粒作为静态缓冲区上传一次，在 GPU 上按同一公式算位置，每颗画成贴在膨胀球壳上的扁平三维高斯（EWA 投影，侧看锥壁时有临边增亮）。颜色和亮度系数按 ln r 查表（红外伪彩色一张、可见光一张，`pinwheel_visible_light` 算热辐射和散射光的份额）；倍速快到螺线每秒转 1 圈以上时淡出（`Simulation::pinwheel_visibility`）；导演的系统范围算上了内三圈。半径由 GPU 查一张以 √age 为横坐标的表。显示亮度以场景的 `display_reference_au` 为基准、按 `display_gamma` 压缩；最新的尘埃远亮于基准时整体调暗（`pinwheel_exposure`，像眼睛按最亮处适应，从外面看最多 20 倍，相机在最内一圈以内时放开），由 `Simulation::pinwheel_exposure` 按真实时间缓动，时钟跳变时直接取目标值。同一套几何也画碰撞星风留下的压缩星风壳（η Car 的化石星风壳）：`emission_color` 让它在两种视图里都按一种谱线的颜色发光（不需要尘埃温度和反照率），`ring_weights = "angle"` 让窗口里每环的量相同（不乘 dt/dν），`from` 之前发出的颗粒不画
- 膨胀星云（场景的 `[[nebulae]]`，最多 4 个，`scene/nebula.hpp`，η Car 的 Homunculus 和 Little Homunculus）：一次喷发抛出、此后匀速飞离的壳，形状按 (t − t_ejection) / (t_reference − t_ejection) 等比放大（均匀膨胀）。壳是观测做的三维网格（`shape`，模型自己的单位；`model_center` 是中心星在模型里的位置，`unit_au` 是 `reference_time` 时一个单位的长度，`polar_scale` 再沿模型的 z 轴缩放），朝向由 `axis_pa_deg`、`axis_tilt_deg`、`roll_deg` 定（`nebula_model_to_icrf`）。加载时按面积在三角形上撒固定的面片（R2 序列，法线取网格的平滑法线），所有星云的面片放在一个缓冲区里，`NebulaPass` 每个星云（红外视图里每一层）一次绘制，在 GPU 上缩放，按风车星云的办法画成 EWA 光斑。`light = "scattered"`（默认）：亮度 = 星光 1/d² × 柱密度 1/scale² × Henyey-Greenstein 相函数；`"emission"`：谱线发光，1/d²、各向同性。都按 `display_gamma` 压缩，年轻的小星云最多亮 2 倍（`nebula_exposure`）。`[[nebulae.infrared]]` 是热红外的尘埃层（大小份额、壁厚、相对柱密度、平均距离处的温度，温度按 d^−`temperature_index` 变），事件 `nebula_infrared` 时逐层画成伪彩色：光 = 柱密度（1/scale²）× 黑体在 8–25 μm 的光（T⁴ × 份额，`nebula_infrared_emission`），颜色查一张按 ln T 的表（温度乘 26.7 的黑体色）。风车星云和它共用 `SplatPass`（点缓冲区、半分辨率单采样的离屏目标 `SplatTargets`、合成；一次 render 可以有多次绘制，各自推 uniform，面片的起点由 uniform 给，不用 first_vertex），新的光斑类效果也用它
- 盘（场景的 `[disk]`，`scene/disk.hpp`，碎屑盘 β Pic 和原行星盘 PDS 70）：轴对称的尘埃，碎屑盘光学薄，透过去看，没有消光和阴影，最多 4 个组分（`[[disk.components]]`：各自的盘面 `node_deg`、`i_deg`，同目视轨道的约定；垂直是两个高斯，`flare_radius_au` 以内厚度不变、以外正比于 r；柱密度是 `radial_au` + `column` 的表或 `ring_center_au` + `ring_sigma_au`）。四种视图，由事件切换：默认散射光（柱密度 × 垂直形状 × 三个 Henyey-Greenstein 的相函数 / d²），`nebula_infrared` 是中红外（× 黑体在 8–25 μm 的光，伪彩色，和星云共用这个键和波段），`disk_millimeter` 是毫米波（× T），`disk_gas` 是 CO 线（柱亮度，颜色是开普勒速度沿视线的分量，`disk_line_of_sight_km_s`，蓝近红远，`gas_velocity_km_s` 是满色的速度，`gas_display_gain` 只改这个视图的亮度）；组分的 `scattered`、`infrared`、`millimeter`、`gas` 是它在各视图的增益，0 就不画。环可以有第二个高斯（`ring2_center_au`、`ring2_sigma_au`、`ring2_share`）。谱线组分可以写 `line_temperature_k` 和 `line_upper_k`，表乘上 exp(−E/T)/T。`[[disk.clumps]]`（最多 2 个）是气体里的团：柱密度乘 1 + 增益 × 径向高斯 × 方位形状（前面高斯，后面按尾巴的长度指数衰减），按自己的周期转（`disk_clump_factor`，着色器的 `disk_clumps`，两边一起改）。`DiskPass` 在半分辨率、单采样的离屏纹理里沿每条光线积分（`shaders/disk.hlsli`，每步取到盘面距离的一半和到恒星距离的 0.3 中小的那个，最多 48 步，薄盘正看和侧看都不用很多样本，CPU 的 `disk_radiance` 走同样的步子，测试用它），合成时在挡光体的轮廓附近按全分辨率重算（恒星和够大的行星，结构同 `CoronaPass`）。亮度相对于“过 `display_reference_au` 沿盘面的光线”按 `display_gamma` 压缩，相机离恒星不到 `display_adapt_au` 时按距离调暗。`thick = true` 的盘在散射光里光学厚（毫米、红外和 CO 仍光学薄）：有 `tau_per_column` 的组分在同一平面，表存垂向光学深度，每步散射 `消光 × (增益 × 相函数 × exp(−τ★) + D) / d²`，并按累积的透射率挡住后面。写了 `albedo`、`asymmetry` 时是 TMS 近似（`disk_apply_similarity`）：光学深度乘 1 − ωg，增益 ω/(1 − ωg)，D 是后面各次散射的各向同性源（`disk_build_diffuse`，加载时算，约 0.35 s，存在和 τ★ 同一张 RG16F 纹理的 g 通道）；不写就是单次散射。不透光度和相函数用 `tools/pds70/mie_opacity.py` 算。τ★ 是加载时烘的二维表（`disk_build_star_depth`，ln r 和仰角的平方根，256×128，存 ln(1+τ)，R16F 纹理），合成改成 `L + T × 场景`（薄盘 T = 1），恒星按相机到星的光学深度变暗（`disk_optical_depth`，`draw_stars`）。参考亮度取“垂直穿过盘”。步长是两段式的（`march_thick` 和 `disk_march_thick` 要一起改）：有尘埃处（一个 σ 的光学深度 > 1e-3）按高度、光学深度和到星距离取小步，其余跳到离中面 5σ 处，最多 96 步。`[[disk.sources]]` 是盘里发光的小东西（吸积的 Hα、环行星盘），不被盘挡，太小时按固定的像素尺寸画。从地球看的特写不要把相机放进盘里（透视和前向散射都不对，见下面的坑），用 `observer` + `observer_offset_km` + `fov_deg`。一次碰撞的尘埃流（`[nebulae.stream]`，`DustStream`）是星云的另一种来源：加载时按引力、辐射压和星际介质的拖曳（`ism_h_per_cm3` 等，`ism_beta_index` 是拖曳随 β 的指数）积分，`infrared_only` 的星云只在红外视图显示（猫尾）。彗星可以是任何恒星的子天体（β Pic 的示意系外彗星，恒星要有 `gm_km3_s2`），不一定是根天体的
- 极光的结构（`scene/aurora_model.hpp`）：天体最多有一个模型（`Body::aurora.model`，有数据的模型另带 `kp`、`jovian`、`ganymede`），`[aurora]` 表只写一个模型（木星的 `jovian_oval` 可以再带 `footprints`），`aurora_shown` 决定哪个视图显示它（可见光、紫外、红外三个开关，事件键 `aurora_visible`、`aurora_ultraviolet`、`aurora_infrared`，手动的伪彩色视图同时开紫外和红外）。`aurora_profile()` 给出一个发光壳的全部参数（壳层范围、两层颜色、层参数、卵表、高度表），物理和显示常数都在 `scene/`；`SceneRenderer::add_aurora_items` 只做坐标变换和打包。着色器的行为是三个独立的选择：`AuroraPattern`（横向）、`AuroraVertical`（垂直）、`AuroraLighting`（昼夜和遮日），数值和 `aurora.frag.hlsl` 的 `PATTERN_`、`VERTICAL_`、`LIGHTING_` 常量对应。新模型：加枚举值和 `aurora_profile` 的分支，能用已有的图案和剖面就不加新的，在 `appearance_tests.cpp` 的 `test_aurora_profiles` 里列上。固定在本体上的亮度图都用 `AuroraMap`（`scene/aurora_map.hpp`，文件头 `AXAURORAMAPn 宽 高 上限`，n 层依次排列，放在纹理的 r、g、b、a 通道，南到北的行、西经的列、都在场景的本体系里），着色器走 `PATTERN_BODY_MAP_SHARES`（四层按份额混合），要按太阳变化的走 `PATTERN_BODY_MAP_SEASONS`、`PATTERN_BODY_MAP_DUSK`（参数都是 `AuroraProfile::map_sun`，`aurora_map_pattern` 是 CPU 上同样的公式，测试用它，改的时候两边一起改）
- 冰巨星的红外 H₃⁺（`[aurora].uranus_map`、`neptune_map`，`tools/geomag/ice_giant_h3p_bake.py`）：天王星横向取 Melin 2025 的全球图，低纬接 Tiranti 2026 的临边经度曲线，25°S 以南没有数据，是外推（对观测部分拟合纬度剖面和 L 壳剖面，磁坐标由 Melin 的 Q3 L 壳和偶极赤道表插值），这是“受照时”的一层，另一层是没有阳光时剩下的极光（观测到的斑点和沿壳的平均），着色器按此刻日下点算出的正午太阳天顶角逐样本混合，所以有季节，垂直取 Tiranti 的总发射剖面（`VERTICAL_TABLE`）；海王星是 Melin 2025 的底光加一块拟合的南半球亮区，再加一条随太阳固定的晨昏渐变（昏缘亮、晨缘暗，±23%，只加在底光上），高度剖面是假设。都是伪彩色，图的结构是固定的快照，随本体转。3–4 μm 里行星本身是暗的，亮度差只有约 2 倍的 H₃⁺ 叠在可见光行星上会看不见，所以红外视图把它们的反射光压到 3%（`aurora_surface_gain`）。天王星的自转率用 Lamy 2025（不是 IAU 的），`test_ice_giant_aurora` 对照 Melin 表 1 的 12 个中央经线
- 地球的夜间灯光和极光：天体的 `[bodies.night_lights]`（灰度灯光图，和 `texture` 同一投影，在天体着色器里按地面的太阳高度角亮起）和 `[bodies.aurora]`（Kp 记录，`scene/aurora.hpp`）。`Simulation::aurora_kp()` 按倍速取 Kp 的滑动平均（画面上 3 秒，最少 3 小时），`aurora_emission(kp)` 由 Starkov 卵、半球功率和谱线比值算出卵的边界（按磁地方时 48 个采样）和一个竖直柱的颜色；`AuroraPass` 在 80–450 km 的壳层里光线步进两层 Chapman 发光剖面，磁极取 IGRF 偶极，红线按局部 427.8 nm 亮度逐点算。卵里的弧和射线查一张平铺噪声纹理（投影到垂直于偶极轴的平面，没有接缝）。灯光和极光都按“长曝光”放大（`aurora_model.cpp` 的 `kEarthExposure`、场景里的 `gain`），白天淡出
- 磁力线（`[bodies.magnetic_field]`，`scene/magnetic_field.hpp`）：`MagneticField` 是施密特球谐的内部场，坐标用模型自己的半径和右手系。加载时在偏心偶极赤道的 `l_shells` 环上取种子，RK4 追踪到行星椭球面或 `max_radius`，抽稀成折线（w 是淡出量）存在 `BodyMagneticField::lines`。线固定在模型坐标系里，`magnetic_model_to_body` 给出它在本体系的方向（`flip_pole` 翻极、`longitude_offset_deg` 或模型自己的 `prime_meridian_deg`），`SceneRenderer::build_magnetic_lines` 每帧变换后交给 `OrbitLinePass`（`tail_opacity = 0`），不另开 pass。默认不显示：事件的 `magnetic_lines = true`（只画焦点天体的线）、设置键 `magnetic_view`（0 跟随事件、1 总是，所有天体、2 从不）或 `B` 键打开。转得太快（`field_line_visibility`）或行星太小（`field_line_size_fade`）时淡出。事件的 `focus_orbits_only = true` 只画焦点天体的天然卫星和 `trails` 列出的天体的轨道，近景不被远处天体的轨道穿过；`hide_orbits = true` 全都不画
- 磁场图（`[bodies.magnetic_map]`，`FieldMap`，`surface`，可选的 `deep` 和 `crust`）：离线烘焙的图，行按行星中心纬度（和 `body.vert.hlsl` 的 v 一样），列从 180°W 起，有两种，由文件头决定。径向图（火星、月球、木星、天王星、海王星、地球的地壳）存 Br 的平方根，可以在某个值饱和（地壳的 600 nT），强度图（`magnitude`，地球的 `surface`）线性地存 |B|。`SceneRenderer` 把它作为线性纹理交给 `BodyPass`，`body.frag.hlsl` 在光照之后混合：径向图红向外、蓝向内，按 |value|^1.5 混合；强度图用 viridis，固定不透明度；整个天体统一压暗，夜面也显示。采样的 u 由 `field_map_u` 给，和贴图的 `texture_left_lon_deg` 无关；图默认固定在本体系，模型有自己的本初子午线时（海王星的 NLS）图烘在模型的经度里，`magnetic_map_turn_deg` 每帧转动。事件的 `magnetic_map = "surface"` / `"deep"` / `"crust"`、设置键 `magnetic_map_view`（0 跟随事件、1 表面、2 深层、3 从不、4 地壳）或 `G` 键打开；天体没有要的那张图时退回 `surface`
- 卫星到足点的线（模板 `[aurora].flux_tubes`，`JovianAurora::flux_tube`，`tools/geomag/jupiter_flux_tube_bake.py`）：Io、Europa、Ganymede 各两种线、南北各一条。淡色是此刻穿过卫星的磁力线（JRM33 + CON2020），落在没有亮斑的 `instant_footprint`；橙色是阿尔文波此刻的位置，走过时间份额 w 的点取“卫星东边 δ·w 那条线上的同一个点”（δ 是 Hue 2023 的领先角），横穿磁力线，终点就是 `footprint()` 的亮斑。表是名义半径、赤道面上每 5° 一条线，每条 40 个点，绕 z 转到起点经度为 0；运行时按经度插值，拉到卫星和足点表给的脚，再用 Catmull-Rom 细分。`Scene::update` 把卫星的本体系位置写进 `Body::jovian_moons`，`SceneRenderer::build_magnetic_lines` 交给 `OrbitLinePass`。事件的 `flux_tubes = true`（不画 JRM33 的线）或 `B` 键打开，w 和拟合的标高见 `assets/geomag/SOURCES.md`
- 彗星撞击（`[[bodies.impacts]]`，`scene/impact.hpp`，木星模板里有 SL9 的 21 个）：撞击点是一张表（名字、Hammel 表 2 的等级、到达木星的时刻、行星中心纬度、System III 西经，`tools/sl9/impact_views.py` 从碎片表的终点打印），其余都只是时间的函数。羽流用 Harrington & Deming 2001 的弹道模型：此刻在 p 处的物质以 v = (p + g t²/2 z) / t 出发，密度 f(v)/t³，什么都不用存，`ImpactPass` 在羽流的盒子里沿光线积分，地影以上被照亮，影子里有逐渐消失的辉光。羽流的亮度定在实测上（`kPlumeExtinctionKm2`、`kPlumeAlbedo`：Hammel 的 I/F 和 Chapman 的上限）。闪光是闭合形式积分的高斯，光变曲线是火流星加火球（`impact_flash`，比例按等级），峰值有测量的写在 `flash_share`（Martin、Chapman），颜色按黑体随温度变（`impact_flash_temperature_k`）；闪光不做深度测试，是否被挡只看火球中心。疤痕在 `impact_scar.frag.hlsl` 里逐像素算，用乘法混合画在天体之后：光学深度 = 每个点的图案（核、条纹、新月形、环，`impact_scar_shape`，加上线性风切）+ 一条带，乘到云层上的因子是 T + ρ(1 − T)，T = exp(−a τ (1/μ + 1/μ₀))（`impact_debris_factor`，薄处透出云层，所以越靠近边缘越暗，厚处只显出残骸自己的反照率，常数来自 `tools/sl9/debris_optics.py`）。`impact.cpp` 和两个着色器是同一套公式，要一起改。事件可以用 `aim_offset_km` 让相机看向焦点以外的地方（配 `observer`，observer 可以是焦点自己），`observer_body_fixed = true` 时两个偏移都在焦点的本体系里，相机跟着转（从正上方看撞击点）
- `SceneRenderer` 持有 pass 和每个场景的 GPU 资源，显示选项通过 `ViewOptions` 传入
- 木星紫外主卵：模板 `[aurora].jovian_oval` 读取 Head 2024 的南北统计轮廓（`JovianAurora`，每极 256 点），`AuroraPass` 在木星本体坐标中画成发光壳。`ViewOptions::aurora_view` / 设置键 `aurora_view`：0 跟随事件、1 可见光、2 紫外极光叠加；事件用 `aurora_visible = true` 开启可见光，`aurora_ultraviolet = true` 开启紫外叠加；普通事件和无活动事件时关闭所有极光，包括地球、Europa 和土星的可见光。手动模式覆盖事件，紫外叠加保留地球和 Europa 的可见光；地球极光还需要开启大气。叠加保留可见光云层，蓝色和明暗细节是示意，白天也画。木星不需要配置大气散射；加载与地球的 `kp` 互斥。烘焙命令及许可见 `assets/geomag/SOURCES.md`。
- 土星紫外主卵：模板 `[aurora].saturn_oval = true` 启用 Carbary 2012 的南北地方时统计轮廓，与 `kp` / `jovian_oval` / `footprints` 互斥。`SaturnAurora::profile` 按半小时采样形状、宽度、晨侧亮度；`noon` 从本体坐标的太阳方向投影得到，不能随 System III 转。仍用 `aurora_view` 和 `aurora_ultraviolet` 切换。先画环，再画土星极光；每个样本光路查 RSS 环剖面，只有环在样本前面时衰减。环纹理导数要在发散的壳层裁剪前算。参数、假设和许可边界见 `assets/geomag/SOURCES.md`。
- 土星可见光极光：`aurora_view` 的可见光模式和 `aurora_visible` 事件使用夜侧粉红/紫色双 Chapman 层（Dyudina 2015），紫外事件使用蓝色统计叠加。两种模式都在环后合成并逐样本处理透射。可见光的高度、颜色和增益是观测趋势约束下的显示假设，横向仍取 Carbary 地方时统计卵。土星末尾的 2010-11-27 导览给出夜侧示意视点；参数及局限见 `assets/geomag/SOURCES.md`。
- Ganymede 自身极光：模板 `[aurora].ganymede_map` 读取 Marzok 2022 的 80×40 HST 柱亮度图，西经、南至北纬度，和其他极光字段互斥。文件四层：平均图和木卫三在电流片上方、内部、下方（±6° 磁纬）的分组，分组 = 平均图 × 平滑的比值，缺测区趋于 1，三组按停留时间平均回到平均图（`tools/geomag/ganymede_aurora_bake.py`）。`Simulation::moon_aurora` 由 `Body::in_parent_axes`（卫星在母天体轴里的位置，木星的卫星就是 System III）算磁纬，给出四层的份额（`scene/jovian_plasma_sheet.hpp`，倍速高时退回平均图），`aurora_profile` 的最后一个参数 `MoonAuroraState` 把它带进着色器。`AuroraPass` 在 0 至 70 km 均匀壳层积分，使用独立线性纹理和 mipmap，紫外叠加开关控制显示。木星场景末尾两个导览看前导、拖尾半球。缺测零值和负残差截零不能解释为无极光，颜色、曝光和高度近似见 `assets/geomag/SOURCES.md`。
- 木星卫星足点：同表的 `footprints` 读取 JRM33（13 阶）+ CON2020 的离线三维表，覆盖 Io、Europa、Ganymede 的本体经度、纬度和距离。`Scene::update` 按卫星相对木星的 ICRF 位置转回 System III，先加 Hue 2023 的半球独立赤道领先角，再查南北交点；多个输出共用结果，位置不走限速动画时钟。`AuroraPass` 单独积分足点的 900 km 发光层，主卵仍用 400 km 层；光斑大小和显示强度为示意。
- GPU 分两层：`GpuDevice` 只有一个；`RenderOutput` 每个窗口一个，持有交换链和本输出的 `SceneTargets`。`Frame` 的 `width/height` 是场景尺寸，`output_width/output_height` 是交换链尺寸，两者在有渲染缩放（`render_scale`，各模式都有）时不同。MSAA 采样数（`msaa` 设置：1/2/4，`GpuDevice::set_msaa` 按硬件支持往下取）存在 `GpuDevice` 的 `SceneTargetFormat` 里，各 pass 的管线和 `SceneTargets` 都依赖它：运行时切换要调 `SceneRenderer::refresh_pipelines()` 重建管线，`SceneTargets::ensure` 发现采样数变了会自己重建。新 pass 要实现 `set_format`，并加进 `refresh_pipelines`。所有和输出尺寸相关的纹理（MSAA/HDR/深度、泛光链、黑洞追踪目标）都放在 `SceneTargets` 里，pass 内部不要再缓存按尺寸分配的纹理，否则多个不同尺寸的输出轮流渲染时会反复重建。各输出都用 SDR 交换链，格式相同，共用一个 `SceneRenderer`

宿主（`App`、屏保、壁纸）只负责窗口、输入、UI 和闲置策略，不要把仿真或渲染逻辑放回 `src/app/`。宿主每帧只 `sim.update` 一次，再对每个输出调用 `SceneRenderer::render`；与相机无关的工作（例如尘埃彗尾的传播）不要每个输出重做一遍。

渲染流程：
1. 黑洞光追（有黑洞时）：在屏幕外以半分辨率追踪，记下每条光线的天空方向和盘的光（`shaders/black_hole.hlsli` 是两边共用的追踪代码）；第 2 步的黑洞合成按全分辨率逐像素插值方向，相邻光线不一致处（阴影边缘、光子环附近）在合成着色器里重新追踪；天空逐像素查：星按局部线性化的透镜映射画成点（`StarfieldPass` 把星按方向排进网格，供合成着色器查找），银河直接采样它自己的贴图，天空被压缩到星点超出网格查找范围处（光子环附近）改用带 mip 的立方体贴图。有风车星云或膨胀星云时，尘埃也在这一步画进屏幕外半分辨率、单采样的纹理
2. 场景 pass，渲染到 HDR（RGBA16F + MSAA）：银河 → 星空 → 天体 → 日冕合成 → 大气层 → 极光 → 羽流 → 彗发/离子彗尾 → 尘埃彗尾 → 行星环 → 土星极光 → 小行星带/柯伊伯带 → 风车星云合成 → 膨胀星云合成 → 碎屑盘合成 → 黑洞合成 → 恒星 → 脉冲星光束 → 轨道线
3. `PostProcess`：泛光 → ACES → sRGB → 交换链
4. ImGui 叠加 pass

场景用 `assets/scenes/*.toml` 描述（格式见 `jupiter.toml` 的注释），新增场景不用改代码。

天体固有属性可以放在 `assets/bodies/`，由场景的 `template` 单层引用；覆盖、`unset` 和字段清单见 [Body Template](assets/bodies/README.md)。模板与资源路径相对 `assets/`。太阳、地球、月球、水星到土星的行星和多个场景共用的卫星，共享参数及来源都在模板中（含环、羽流、J2/J4），轨道计算的 `gm_km3_s2` 留在场景；只在一个场景出现的天体仍内联。读取天体属性的 Python 工具用 `tools/scene_templates.py` 展开模板。

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

不依赖测试框架，覆盖 `core/`、`ephem/`、`scene/` 和设置文件。CTest 也注册了（`ctest --test-dir build/win-msvc -C Debug`）。测试直接读源码树里的 `assets/`（通过 `ASTRAXIS_ASSET_DIR`）。Debug 下 N 体积分慢，整套测试要六分多钟；Release 约 20 秒。

按领域分文件，每个文件导出一个 `run_<领域>_tests()`，由 `tests/test_main.cpp` 依次调用；`check`、`load_scene_or_die` 等公用函数在 `tests/test_util.*`。新增测试放进对应领域的文件，只有这个文件用到的辅助函数放在它的匿名命名空间里：

| 文件 | 内容 |
|---|---|
| `core_tests.cpp` | 日历（含 Kerbin 历法）/时间尺度、星表（含重投影）与黑体颜色和波段份额、银道坐标系 |
| `ephem_tests.cpp` | 开普勒方程与传播、幽灵的外推（带状谐项、潮汐）与 Δv、星历表、目视/凌星轨道约定、N 体 |
| `planet_scene_tests.cpp` | 场景里的行星和卫星：太阳系（含 Proteus 对照 NEP098 的独立时刻和同步自转）、木星、地月、土星；Janus 和 Epimetheus 的根数表对照 Horizons、19 次交换（a 的跳变之比对照 GM 比、2006-01-21、最近距离）；Kerbol 系统对照游戏里的位置和自转角、6 小时太阳日、潮汐锁定、Laythe:Vall:Tylo ≈ 1:2:4、示意的 Mun 转移（停泊轨道、两次点火的 Δv、捕获轨道） |
| `mission_tests.cpp` | 探测器：先驱者/新视野、Parker、伽利略/朱诺、Artemis/CAPSTONE/嫦娥五号 T1 及其火箭末级、卡西尼/惠更斯、信使/贝皮科伦坡、火星的八个轨道器（含天问一号 2021 年的重建）、赛丁泉彗星与 3I/ATLAS、太阳系场景的 1I/3I、Europa Clipper/Hera 飞掠、金星快车/破晓号/Zoozve 与金星借力飞掠、ISEE-3、哈雷舰队、Rosetta/Philae、空间望远镜（WMAP 的月球借力与 L2 轨道、Chandra 的近远地点与周期、JWST、Euclid、Roman、Gaia、Herschel、Planck 在 L2 附近和弃置轨道、TESS 与月球的 2:1 共振、Spitzer 和 Kepler 的漂离、Hubble 示意轨道的周期与 J2 交点进动及事件显隐）、DART 与 Dimorphos（解 547 的撞前周期和撞后稳定周期、同步自转对照 Horizons、幽灵三天后仍在旧轨道上、2.65 mm/s 的读数、撞击速度和来向、LICIACube 的最近距离）、Shoemaker-Levy 9（`test_shoemaker_levy_9`：近木点和解体对照 Sekanina，21 次撞击的时刻和位置对照 Chodas & Yeomans 表 5，G 的速度和入射角对照表 6，发现时和 1994-05-17 的碎片列长度，Q1、Q2 的分裂，都从母体出发） |
| `shape_tests.cpp` | 小天体形状模型（含 Didymos、Dimorphos 的体积，Proteus 的体积和 PDS 西经约定）；Endurance 致敬模型（硬边、尺寸、环的转速对应 1 g）；探测器模型（`test_spacecraft_models`：图集坐标、`point_at` 的朝向、Voyager、New Horizons、Pioneer 10、Cassini 的天线直径、Voyager 的吊杆、JWST 遮阳罩的长度、Juno 太阳能板的跨度、Galileo 的总高、Hubble 的长和宽、烘焙的阴影：阳光照到 Hubble 太阳能板的前面、天线挡住机身而 New Horizons、Pioneer 10 的发电机在外面、遮阳罩挡住望远镜和 Galileo 的机身、隔热罩挡住 Parker 的机身而太阳能板在外面、阳光照到 Juno 的太阳能板） |
| `appearance_tests.cpp` | 碎片的尾巴（`test_bound_dust_tail`：15,000 km 处的 β、放出早晚的差别、1994 年 1 月和 7 月的位置角、碎片列的位置角和 30° 夹角、长度和宽度，都对照 Sekanina）；彗星撞击（`test_impact_sites`：撞击点对照碎片表的终点和 Hammel 表 2，羽流对照 HD01 和 Hammel（22% 的外壳、3000 km、6–9 分钟到顶、15–20 分钟落回、密度积分等于总质量、朝东南），羽流顶越过边缘的时刻对照 Chodas 表 1，羽流亮度对照 Hammel 的 I/F 和 Chapman 的上限，闪光对照 Martin、Chapman、Carlson（2 秒上升、平台、K 的半高宽 27 秒、N 的微光、七个峰值、K 在 Q1 和 H 之间、色温、火球 40 km 和上升速度），环的半径和消失，G 的新月形外缘 13,000 km，小撞击点和未探测到的，光学深度对照 Sanchez-Lavega 表 IV，核心和新月形的明暗对照 West 图 3，带的位置和消失）；卫星到足点的线（`test_flux_tubes`：表的脚对照足点表、波的时间份额单调、两种线的两端、波在磁力线东边不超过领先角、亮斑离线的脚的距离）；彗尾、尘埃彗尾、大气层、羽流（含各场景参数一致）；地球的夜间灯光和极光（各场景一致、大磁暴的 Kp、Starkov 卵的性质、Longyearbyen 在磁正午处于 Kp 3 的卵下、IGRF 偶极极、半球功率、谱线产额和比值、Kp 随倍速的平均）；各行星和卫星极光的模型（`test_aurora_profiles`：每个模型的图案、剖面和光照选择，哪个视图显示，壳层和颜色）；等离子体片（`test_plasma_sheet`：离心赤道的高度对照 Milby 的表和场景里的观测时刻，Callisto 的拟合对照曝光，Europa、Callisto 的调制幅度，高倍速退回长期平均）；Europa 的图案（`test_europa_pattern`：Roth 2016 的最上方和最下方、磁力线的倾斜方向和大小、平均归一化、下午一侧 1.59、面向等离子体片的半球和上下游的差别、象限比、visit 17 的亮度变化，本体 x 轴朝木星，高倍速平均掉）；Io 的日食（`test_io_eclipse_levels`：钠的 27%、10 分钟和 2 小时，SO₂ 的 2.5 倍和 10 分钟恢复，钠的参考水平对照 Milby 的 7 次曝光，入影的搜索，食长，x 轴朝木星）；Ganymede 的极光（`test_ganymede_aurora`：图的方向和接缝，磁纬对照 Marzok 表 1，分组按停留时间平均回到平均图，内外亮度比对照论文，缺测区的回退，按磁纬选层，高倍速退回平均图）；冰巨星的 H₃⁺（`test_ice_giant_aurora`：天王星中央经线对照 Melin 表 1、南极光和旁边暗斑的对比、临边的亮度低谷、外推的南半球（2023 年的日下点、往南变暗、晨昏线处的份额、极夜的极点全黑，沿 Q3 追踪 L 壳查极光的位置，磁场转 180° 作对照，接缝；季节：1986 年日下点在南、北半球的极夜里只剩 Melin 的南极光，2007 年昼夜分点两个半球差不多亮）、高度表归一；海王星亮区的位置和亮度、北半球没有亮区、晨昏渐变（昏缘亮、跟着太阳走、亮区不变，太阳在本体上向西走））；磁力线（`test_magnetic_field`：偶极的解析形状、无散无旋、O8 的偶极和偏心对照论文、Voyager 2 的 NLS 经度、AH5 的南磁极对照 Melin 的南极光、IGRF 偶极对照极光的地磁极、Io 的线对照足点表、水星的偏心和南北不对称及线都在磁层顶以内、Ganymede 的磁极方向和本体坐标系及线的落点对照极光、转动时线的摆动、线的数量和范围、事件开关和错误配置）；磁场图（`test_magnetic_maps`：火星的强场在 Terra Cimmeria 和 Sirenum、南北差异、Hellas 和 Argyre 几乎没有；木星两张图逐格对照 JRM33、Great Blue Spot 的位置和强度、深层更突出；月球的峰值在 Gerasimovich、Reiner Gamma 和 Mare Crisium、远近月面的差别；地球图逐格对照 IGRF、南大西洋异常的位置和强度；地壳图的 Kursk、Bangui 和开阔的太平洋；天王星、海王星逐格对照 AH5、O8（海王星在两个时刻，随 NLS 转动），最强处在偏心偶极一侧；事件开场时特征转到合适的位置） |
| `stellar_system_tests.cpp` | HD 206893（B、c 对照 GRAVITY 的位置，周期和相互倾角对照 Kral 2025，盘的朝向、近侧、两个峰和缺口对照 Marino 2020，候选卫星的周期和晃动）、α Cen（含 Ab 对照 S1/C1 的位置和相互倾角、只在事件里显示的 Bb、阿凡达虚构事件的显隐与轨道稳定范围粗筛）、β Pic（行星对照 GRAVITY 的位置、周期、2017 年的希尔球凌星；盘的脊线幂律和宽度对照 Golimowski、副盘的倾斜和亮度、近侧在 NW、毫米带的两个高斯、积分步子对照细分求和、猫尾的 45° 爬升、尾尖 21° 和 265 au；CO 团的 30% 流量、位置和 3–4 km/s、SW 更亮、转动和尾巴的方向；6 颗彗星在 TESS 的时刻穿过恒星、时长和核半径）、PDS 70（b、c 对照 GRAVITY，候选 d 对照 Hammond，周期比；标高对照 Portilla-Revelo，τ★ 表对照直接求和，散射亮环约 54 au 和毫米环 73.7 au，缺口暗、西侧亮，内盘占的流量，侧看挡住恒星，穿过盘轴的光线和步子对照细分求和；相似关系和相函数的归一，后面各次散射的表对照小蒙特卡洛；第二峰的肩对照 Keppler 图 2；气体的浅缺口和转向）、TRAPPIST-1、Kepler 系统、TIC 168789840、Algol（CHARA 位置与 B 的拉长、TESS 极小、GCVS 食深）与光变曲线模型、洛希面（Kim 1989 的半径、Eggleton 公式）、δ Cephei（SPIPS 半径、Engle 极大、视向速度轨道）、北极星（Anderson 的脉动视向速度、Evans 的极大时刻、光极大在视向速度极小前 0.21 d、V 振幅、Ab 对照 HST/CHARA/斑点干涉的位置、K1）、PSR B1620-26、WR 104（Hill 的合相和 K、Richardson 的质量、旋臂间距、从地球看的旋向和 Tuthill 的 θ0、宽锥让投影的旋臂接近圆、Soulain 的 r⁻³ 径向剖面、伴星 B 的位置、可见光里成核处以热辐射为主而一圈外只剩散射光）、WR 140（Thomas 的视向速度、Monnier 和 Thomas 的干涉位置、质量、产尘窗口和时长、开始时的星距、远星点时最新的尘埃、近星点产尘率的低谷和后缘的密度比、Han 式 1 的三段运动、Lau 的壳间距、壳的朝向对照 Han 的代码、最新的壳最亮）、η Carinae（轨道朝向对照 Madura 2012、伴星在近星点时位于主星后面、2025.6 的近星点、Homunculus 的轴和极半径对照 Smith 2006、与轨道轴对齐、两个凸起在远星点两侧 ±55°、2012 年的长度、红外两层的温度对照 Smith 2009、腰部更暖、暖层在 8–13 μm 更亮而冷层在更长波长、膨胀时变冷、Little Homunculus 的长度和帽盘半径对照 Smith 2005/2006 并始终在 Homunculus 里面、化石星风壳的间距对照 Gull 2016、每环等量、在近星点一侧和星的后面、1896 年以前没有、AAVSO 光变曲线的 1844 年峰值和 1916–1943 年的缺口、由它定的亮度） |
| `relativity_tests.cpp` | Kerr 测地线（类时/类光）、后开普勒双脉冲星、B1913+16 测地线进动（WT02 视线夹角、速率、示意光束退出）、Sgr A*、Gargantua（[JA15] 的盘和相机、最内稳定圆轨道上的示例行星、米勒星的时间膨胀、Mann 星的轨道对照 [TH14]、假想的 Edmunds 星看到的阴影大小和用光线追踪积分的盘的辐射、Endurance 假设的 30 M 轨道的周期和时间膨胀） |
| `scene_tests.cpp` | 场景加载、Simulation、随天体自转的参考系、旋转参考系转得快时星空淡出、风车星云在高倍速时淡出、WR 140 尘埃的曝光（调暗上限、缓动、跳变时直接到位）、标签布局、跨场景一致性（同一天体在各场景的轨道参数）、相机导演（含 WR 104）、事件说明 |
| `settings_tests.cpp` | `config.toml` 读写与场景列表 |

## 资源与工具

- 星历烘焙：`python tools/bake/horizons_bake.py tools/bake/ephem.toml [name ...]`，联网访问 Horizons，输出到 `assets/ephem/`，说明见 `assets/ephem/SOURCES.md`。卫星相位拟合：`tools/bake/fit_moon_phase.py`。离扁行星很近的卫星用根数表：`python -E -P tools/bake/element_table_bake.py [name ...]`（支持 Janus、Epimetheus、环旁小卫星及特洛伊卫星，名称小写；需要 numpy，下载每小时的样本缓存在 `tools/bake/cache/`），输出 `assets/ephem/*.elm`。独立时刻复查：`python -E -P tools/bake/verify_saturn_small_moons.py`
- Proteus 星历：`python -E tools/bake/proteus_bake.py <nep098 摘取文件.bsp>`（numpy、jplephem），直接读取 JPL NEP098，下载命令见 `assets/ephem/SOURCES.md`。SPK 的 808 和 899 都相对系统质心 8，需相减后才能作为海王星中心的卫星轨道；每四小时偏移 37 秒的独立时刻检查拟合误差。
- `assets/stars/hyg.csv`：`tools/stars/convert_hyg.py`（HYG v4.4，带距离，供近邻系统重投影）；M4 星团天空：`tools/stars/make_m4.py`（Gaia DR3）
- `assets/sky/`：远处场景和银心的模型星空，`tools/sky/make_galaxy_sky.py`（需要 numpy；先 `calibrate` 对照 `milky_way.jpg` 标定，结果在 `tools/sky/galaxy_calibration.json`，再 `render <场景.toml>`，一个场景约 10 分钟，银心约 20 分钟）。改了模型或标定要重新生成全部场景。工具用 multiprocessing：Windows 上子进程重新导入模块，在主进程里改的模块全局变量传不过去，参数要经 `_init` 传
- `assets/belts/*.bin`：`tools/belts/make_belts.py`（JPL SBDB，主带 H < 15、全部 TNO）
- `assets/geomag/kp.bin`：`python -I tools/geomag/kp_bake.py <kpYYYY.wdc 目录> assets/geomag/kp.bin`（GFZ 的年度 WDC 文件，下载地址见 `assets/geomag/SOURCES.md`），更新时全部重下；`tools/geomag/aurora_colors.py <ciexyz31_1.csv>` 打印 `scene/aurora.*` 里各谱线每 kR 的颜色；`tools/geomag/ganymede_aurora_bake.py <Zenodo 目录> <输出>`（需要 numpy）把 Marzok 的平均图和三张分组图烘成 `ganymede_aurora.txt`；`tools/geomag/europa_pattern_fit.py`（需要 numpy，用 `python -E -P`，约 10 秒）按 Roth 2016 拟合 Europa 图案的系数，打印模型和论文的比值；`tools/geomag/ice_giant_h3p_bake.py`（需要 numpy，用 `python -E -P`）烘焙 `uranus_h3p.txt`、`neptune_h3p.txt`，打印高度表（写进 `aurora_model.cpp`），下载命令见 `assets/geomag/SOURCES.md`；`python -I tools/geomag/field_coeffs.py <下载目录> assets/geomag` 从 libinternalfield 和 IGRF-14 写磁力线用的 `field_*.txt`（只用标准库）；`python -E -P tools/geomag/field_map_bake.py <Morschhauser2014.txt> <T2015_449.sh.gz> <Spherical_HarmonicModel_GRL.txt> assets/geomag`（需要 numpy）烘焙磁场图 `*.fieldmap`（含 1050 阶的地壳，约半分钟）；`tools/geomag/jupiter_flux_tube_bake.py`（JupiterMag，和足点表同一个 venv）烘焙 `jupiter_flux_tubes.bin`；`python -I tools/geomag/saturn_visible_colors.py` 由 Dyudina 2015 图 8 的读数打印土星可见光极光的颜色和两层比例（`saturn_aurora.hpp`），改了层的高度要重跑；`python -I tools/geomag/mars_aurora_profile.py` 由 Schneider 2018 图 2、图 3 的读数拟合火星极光的高度并打印 2017 年 9 月的时间线（约半分钟）
- `assets/textures/*.jpg`：`tools/textures/prepare_textures.py` 从原图缩放（地球夜间灯光 `earth_night.jpg` 在线性空间里平均缩放，需要 numpy）；`uranus.jpg`、`neptune.jpg` 由 `tools/textures/make_ice_giant_textures.py` 从哈勃 OPAL 的 FITS 图合成（纬向剖面加一张真实图，需要 numpy）；Kerbol 场景 16 个天体的纹理（`kerbin.jpg` 等）由 `tools/textures/make_kerbol_textures.py` 程序生成（球面三维噪声、自带哈希，需要 numpy 和 Pillow，全部约 10 分钟）
- `tools/kerbol/`：Kerbol 场景的数据。`kerbol_dump.ks` 是 kOS 脚本（复制到 KSP 的 `Ships/Script/`，在带 kOS 处理器的飞船终端里 `SWITCH TO 0. RUNPATH("kerbol_dump").`），导出 `kerbol_bodies.txt`（已入库）；`kerbol_scene.py`（只用标准库）核对约定，打印场景的天体块、测试用的位置和自转角、大气参数；`kerbol_transfer.py`（只用标准库，借用 `horizons_bake.py` 的 Kepler 传播、抽稀和写表）按拼接圆锥曲线设计示意的 Kerbin→Mun 转移，烘成 `kerbol_mun_transfer.eph`。交会的最近距离由远拱点半径定，不是由点火相位定（飞船沿 Mun 的轨道方向到达、被 Mun 从后面追上，只改相位时都近乎迎头相撞）
- `assets/shapes/*.mesh`（AXMESH2）：`tools/shapes/` 下按来源格式分：`make_arrokoth.py`、`make_grid_table_shape.py`（Thomas / Stooke 经纬网格表）、`make_dtm_shape.py`（DLR 全球 DTM）、`make_plate_shape.py`（板块模型、OBJ、VRML、Gaskell 的 ICQ）
- `tools/shapes/make_spacecraft.py ref [名称 ...]`（需要 numpy、Pillow 和 DracoPy，用带 DracoPy 的 venv 和 `-E -P`，约半分钟）：把 NASA 的 GLB（`ref/NASA-3D-Resources` 的克隆和 `ref/probes/` 里单独下载的）转成网格和图集，做法和署名见 `assets/shapes/SOURCES.md`。新探测器在文件开头的 `CRAFT` 表里加一项（哪根轴是 +z、原点、缩放、要跳过的材质、要涂掉的徽标、要换的颜色、要细分的长度），先把 GLB 里的图和页面上的署名看一遍
- `tools/bake/rejoin_burn.py`：把表里的入轨点火接到其后接缝处的下一段弧上（希望号，见 `assets/ephem/SOURCES.md`），烘焙之后运行
- `tools/isee3/`、`tools/armada/`、`tools/ce5t1/`：ISEE-3、哈雷舰队、嫦娥五号 T1 服务舱的轨道重建（需要 numpy，下载缓存不入库），做法见 `assets/ephem/SOURCES.md`
- `tools/tianwen1/`：天问一号环绕器。`momag_fetch.py` 从 USTC 的 MOMAG 磁强计数据下载位置（约 15 GB，每天一个文件，逐个下载后只留每分钟一个样本，`cache/` 不入库）；`tianwen1_reconstruct.py <jupyter_notebooks/Tianwen 目录>`（需要 numpy）把 D. Estévez 的 GMAT 脚本里的遥测向量外推、在机动处接起来，得到 2021-02-08 .. 11-16 的路径（`cache/tw1_2021.txt`；纯 Python 的 10 阶重力场积分半年的 8 小时轨道要二十多分钟，弧段缓存在 `cache/arcs/`）；`tianwen1_bake.py` 把它和 MOMAG 烘成 `tianwen1_mars.eph`，只用标准库
- `tools/philae/`：67P 的自转和 Philae（只用标准库，`daf.py`、`spk.py`、`ck.py` 读 SPICE 的 DAF、SPK 5/18 型段和 CK 5 型段）。`c67p_spin.py <目录>` 把 ESA 的三个彗核姿态 CK 拟合成 `c67p_spin.toml`（场景里 67P 的 `spin_table`）；`philae_bake.py <目录>` 读两个着陆器 SPK，用这张表换到 67P 本体系烘成 `philae_67p.eph`，所以改了自转表要重新烘 Philae。下载地址和做法见 `assets/ephem/SOURCES.md`
- `tools/dart/`：DART 场景（只用标准库）。`s547_bake.py <dimorphos_s547.bsp>` 读解 547 的 SPK（figshare，Horizons 现在不给这两个天体 2026 年以前的数据），烘成 `didymos_bary.eph`、`dimorphos.eph`，撞击前后各留一个节点；之后 `approach_bake.py` 联网取 DART 最后三天和 LICIACube（相对系统质心），把 DART 平移到撞击时刻正好到 Dimorphos 中心。系统质心和 DART 的巡航在 `ephem.toml` 里
- `tools/sl9/debris_optics.py`（需要 numpy，用 `python -E -P`，几秒，借用 `tools/pds70/mie_opacity.py`）：残骸的 Mie 和两流计算，对照 West 1996 图 3，打印 `impact.hpp` 里 `kDebris*` 的值
- `tools/sl9/impact_views.py`（只用标准库，用 `python -I`，联网）：从碎片表打印木星模板的 `[[impacts]]`，以及 G 从地球、Galileo 方向和正上方看的相机偏移（事件的 `observer_offset_km`、`aim_offset_km`）；重新烘焙碎片后要重跑
- `tools/sl9/sl9_bake.py`（只用标准库，用 `python -I`，联网）：Shoemaker-Levy 9。下载 21 个碎片和母彗星（K 的路径）相对木星的向量，缓存在 `tools/sl9/cache/`（不入库），在分裂后按时间从母体混合到 JPL 解，表在 100 mbar 面结束，烘成 `assets/ephem/sl9_*.eph`；打印对照 Chodas & Yeomans 表 5 的结果和场景的天体块。同期的卫星在 `tools/bake/ephem.toml`（`*_sl9`）
- `tools/eta_car/aavso_light_curve.js`：η Car 的光变曲线 `assets/stars/eta_car_aavso.csv`。AAVSO 的检索页有浏览器验证（curl 拿到的是验证页），要在浏览器里打开 `apps.aavso.org/v2/data/search/photometry/` 后在控制台运行（不用登录），`per_page=1000` 一页取一千条，在页面里去重、分箱后打印结果；后台标签页的定时器会被节流，窗口隐藏时要几分钟
- `tools/pds70/mie_opacity.py <.lnk 目录>`（需要 numpy，用 `python -E -P`，几秒）：按 Portilla-Revelo 的粒径分布和 DIANA 的材料用 Mie 理论算 1.25 μm 的不透光度、反照率、g，并拟合三个 Henyey-Greenstein，打印 `pds70.toml` 用的数；光学常数从 optool 的仓库下载，写在文件头
- `tools/hd206893/scene_numbers.py`（只用标准库）：打印 `hd206893.toml` 的盘表（Marino 2020 表 2 的参数模型）、雅可比坐标的尺度、近星点时刻和朝向地球的相机偏移
- `tools/betapic/exocomets.py`（只用标准库）：把 Lecavelier 2022 表里的几次凌星变成 `beta_pictoris.toml` 里的彗星天体块（近星点、半径来自观测，其余是假设，见文件头）
- `tools/betapic/fit_disk.py`（需要 numpy，用 `python -E -P`）：按 Golimowski 2006 的脊线幂律反推 β Pic 主盘和副盘的柱密度（对视线积分反复做乘性修正），打印 `beta_pictoris.toml` 的 `radial_au`、`column`；改了相函数、盘面或垂直形状要重跑
- `tools/pdftext.py`：只用标准库从 arXiv PDF 提取文字（ar5iv 经常超时）
- `tools/algol/algol_ephemeris.py <TESS *_lc.fits>`（需要 numpy，自带最小的 FITS 读取）：从 TESS 光变曲线量主极小、扣掉外轨道的光行时效应，拟合 `algol.toml` 内轨道用的历表；下载地址见 `assets/ephem/SOURCES.md`
- `tools/cepheid/spips_table.py <delcep2.dat> <pc> <out.csv>`（只用标准库）：把 SPIPS 模型（VizieR J/A+A/584/A80）转成脉动表，见 `assets/stars/SOURCES.md`；北极星没有这样的模型，`tools/cepheid/polaris_table.py <out.csv>`（只用标准库，数据写在脚本里）由论文里的视向速度、平均半径和温度拼出脉动表，并拟合场景用的极大历表
- `res/astraxis.ico`：`tools/icon/make_icon.py`（需要 Inkscape）从 `res/astraxis.svg`（≥ 48 px）和 `res/astraxis-small.svg`（≤ 40 px）生成；.ico 入库，CI 不依赖 Inkscape

## 已知的坑

- Daphnis 的 Horizons 覆盖只有 1990-06-15 至 2018-01-17，且 2010、2016 年元旦有轨道文件接缝；根数拟合窗口要分段，不能跨接缝平滑，处理见 `assets/ephem/SOURCES.md`。

各任务数据集自己的坑（SSCWeb、ICE 导航文件、卡西尼号最后一段等）写在 `assets/*/SOURCES.md` 里；实现上的取舍写在代码注释里。

### 构建、着色器与 SDL / ImGui

- 用解析的椭球画在天体网格上（`SV_Depth`，如撞击疤痕）：网格的顶点就在椭球上，那里两个深度只差舍入，比较时落败的像素会留下一个个亮点，排在顶点上。深度要往相机方向偏一点（乘 1.0005，反向 Z）
- 贴在天体上的图，行要和 `body.vert.hlsl` 一致，按行星中心纬度（`asin(z/|r|)`），不是网格的参数纬度：磁场图曾按后者烘焙，木星上 45° 处错开 3.7°，贴近赤道和两极时看不出来，所以 Great Blue Spot 的测试也没发现
- 带符号的伪彩色叠加：不要让压暗量随强度变化，过零处底图不被压暗，会画出一圈看起来像真实结构的亮线。在线性空间里按压缩后的值直接混合，弱场也会显得很浓（火星曾经全球都是同样饱和的斑块），透明度要再取幂
- **观测光变的显示增益**：η Car 在事件距离上已超过 SunPass 的亮度上限；观测变化通过 `observed_star_display_gain` 在限幅后乘到星点、光晕及恒星球面，不能只加到星等上。年轻发射星云的曝光补偿用 scale^(2γ)，散射光用 scale^(4γ)。观测点精确命中时先返回采样值，再检查相邻缺口。

- VS 生成器会对加进 target 的 `.hlsl` 自动执行 FXC（报 X3721），必须设 `VS_TOOL_OVERRIDE "None"`（`HEADER_FILE_ONLY` 不起作用）
- CMakeCache 里的旧变量会盖过新逻辑：dxc 路径因此改名为 `ASTRAXIS_DXC_EXECUTABLE`
- GCC 的 `-Wmissing-field-initializers`：SDL 结构体用指定初始化器时漏写 `props` 也会报警告（MSVC 不报），要写 `.props = 0`
- HLSL 2021 不支持向量三元运算，用 `select()`；`m[i]` 是第 i 行，取列用 `mul(m, float4(1, 0, 0, 0))`；`fwidth` 不能放在非一致的分支里
- dxc 默认不是 IEEE 严格模式，`isnan`/`isinf` 可能被优化掉，排查 NaN 时要用 `asuint` 看指数位
- 体渲染里很薄的壳层（几 km 厚，而步长十几 km）会显出抖动噪声的规则斜纹：先把视线解析地裁剪到壳层内再采样
- 大的时钟值进着色器前要先 `frac` 再放大，否则 float 精度不够会闪烁
- 在着色器里把线段裁剪到相机平面，不能 `lerp` 完就用：端点的 w 上百万时，算出的 w 舍入误差约 1，比目标值（`kMinW`，1e-4）大一万倍，符号也随机，端点飞到屏幕另一侧，线段每帧消失或翻转 180°。`orbit_line.vert.hlsl` 从 |w| 小的一端插值，并直接把 w 设成 `kMinW`。线越长、节点越稀越明显（Voyager 的节点相隔几百万 km），中等距离上也有，只是拉近到几十米时最显眼。判断这类闪烁，按同一公式用 float 和 double 各算一遍比较屏幕方向
- 从很远处判断两个天体是否在一条视线上，不要用 |c|² − (c·d)² 求垂距：两个很大的数相减全是舍入（4 au 外看相距 1.2 km 的 Didymos 和 Dimorphos，Didymos 的标签被当成挡住了），用叉积 |c × d|
- 在片元着色器里按像素足迹解析地滤波细线（如天王星的窄环）：足迹用 `ddx_fine`/`ddy_fine`，`fwidth` 是 2×2 粗导数，低角度看的细线会变成 2 像素一节的珠串；滤波核只有 1 像素宽的盒子时边缘呈阶梯，用 2 像素宽的三角核（与纹理 mip 三线性过滤相当）
- ImGui 的 sdlgpu3 后端在 D3D12 上用 DXBC，建设备时要同时声明 `DXIL | DXBC`
- 不要用 `SDL_GenerateMipmapsForGPUTexture` 给星空立方体贴图生成 mip：各面是按 Direct3D 面布局、用镜像投影画的，SDL 生成的 mip 和本层对不上，天空在面的交界处出现亮度台阶。`StarfieldPass::render_cubemap` 直接画每一级 mip
- 不用平台后端时也要调用 `ImGui_ImplSDLGPU3_NewFrame()`（第一次调用时创建采样器）；漏了 Debug 能跑，Release 在 `SDL_BindGPUFragmentSamplers` 里崩
- ImGui 开了 `NavEnableKeyboard` 后，有窗口获得焦点时 `WantCaptureKeyboard` 就是 true，全局快捷键只看 `WantTextInput`
- EWA 光斑的投影是在面片中心线性化的：相机在壳层里面时，中心在视野外很远（WR 140 约 50° 离轴）的大面片会被算成伸进画面中央的小椭圆，排成一串等距的红点。中心超出 1.5 倍 NDC 的光斑要剔除（`kMaxOffView`）
- 大面积重叠的加性光斑不要直接画进 MSAA 的 HDR 目标：每个样本都要混合，填充率按采样数翻倍（WR 104 的尘埃在 Arc A380、1280×720、4x MSAA 下只有 8 FPS）。光斑总数不影响这个开销（数量多了每个就小），起作用的是总覆盖面积、重叠系数和分辨率：画进半分辨率、单采样的离屏纹理再合成（`PinwheelPass`，60 FPS）
- 体渲染的每个样本里不要算哈希值噪声：三维值噪声一次 8 个哈希，极光每像素 48 个样本、每个样本 5 次，近看铺满屏幕时只有 45 FPS；改查一张平铺的噪声纹理（`SampleLevel`，循环里不能用导数）后到 60 FPS
- SDL 的 Vulkan 后端给每个 uniform 缓冲区的绑定范围固定是 4 KiB（`MAX_UBO_SECTION_SIZE`），D3D12 没有这个限制：超出的部分是未定义行为，有的驱动碰巧读得到，没有报错也看不出来（极光的 uniform 曾有 5 KB）。大的 uniform 结构要写 `static_assert(sizeof(...) <= 4096)`，放不下的表用纹理或存储缓冲区
- 多输出时每个输出都要等自己的交换链（`SDL_WaitAndAcquireGPUSwapchainTexture`）
- 不要让 SDL 包装别的进程窗口下的子窗口（屏保 `/p` 预览）：对话框在第二块屏上时 SDL 会按顶层窗口的逻辑移动并放大它，预览区只剩黑底。预览窗口现在用 GDI 画
- 视点超过星表的完整距离（`max_viewer_pc`，约 30 pc）时场景加载失败，要用模型星空（`milky_way`、`stars`，`make_galaxy_sky.py render`），即使只有 41 pc：HD 206893
- TOML 的表头会“吞掉”后面的键：顶层键（如 `name`）必须写在第一个 `[table]` 之前
- 以行星为根、没有太阳天体的场景（木星、土星、火星）里，事件的 `from_body = "Sun"` 指太阳方向（`SceneEvent::from_sun`，用 `Scene::sun_position()`），不是天体；这类场景里不要再加名叫 Sun 的天体
- 撞击不能当飞掠事件（`flyby`）：`find_flyby` 要找前后都有数据的距离极小值，星历在撞击点结束就找不到（场景测试报“the named flyby is found”）。撞击事件用 `from_body` + `phase_deg` 定视角（撞在背面时 `from_body = "Earth"`、`phase_deg = 180`）
- 有形状模型的彗核，彗发和离子尾的光线要按网格求交，并使用与表面相同的 `Body::orientation`；用等体积球遮挡会在不规则彗核周围显出圆形暗斑。网格 BVH 在加载场景时上传，自转时把光线转到本体系。
- 彗星可以是行星根天体的子天体（火星场景的赛丁泉），但尘埃彗尾要按太阳的 GM 传播（`Scene::sun_gm`），不能用父天体的 GM；彗星的星历要相对那颗行星烘焙，日心表减去近似根数的行星位置会差几千 km。长期绕行星转的彗星（Shoemaker-Levy 9）不适用：尘埃尾只绕太阳运动，木星的引力几天就把颗粒和彗星拉开，这时写 `dust_only = true`，只画尘埃彗发，大小由 `dust_speed_km_s` 定（气体彗发的尺度在 5 au 处有 200 万 km，不能用）；要尾巴就加 `bound_tail_release = [起, 止]`（`BoundDustTail`，`scene/comet.hpp`）：颗粒只放出一次，对彗星路径线性化后 d = P(t) v + β b(t)，每个放出时刻一张表，第一次用到时积分（一个碎片约 2 ms，运行时查表）
- 改了 `assets/` 里的场景要重新构建才生效：程序读的是构建时复制到 exe 旁的副本。复制用的是 `copy_directory_if_different`，不删文件：场景改名或删除后，旧副本还留在 `build/*/<配置>/assets/scenes/`，场景列表里会多出一项，要手动删掉
- 航天器默认画历史轨迹。给 Hubble 这样的近地航天器用 `mean_elements` 时，要写 `trail = "orbit"`；否则默认历史窗口会取几百年，密集的低轨路径叠成一圈亮带。
- 人造物体的网格（探测器）要写硬边标志（`write_mesh(..., hard_edges=True)`）：加载器默认把位置相同的顶点合并成一个法线（为了小行星的接缝和极点），按面复制顶点的方块会被抹成边缘带渐变的“枕头”
- 烘焙后留在网格里、没有被三角形用到的顶点，加载器给的法线是 (0, 0, 1)：不会被画出来，但按顶点统计时（测试）会被算成朝太阳的。合并顶点或删三角形之后要清掉
- 相机的碰撞检测让它离任何天体中心至少 1.1 倍半径（`kBodyClearance`）。Hubble 的 569 km 轨道就在这以内，相机曾被推到 69 km 外，模型退化成一个标记。现在目标本身在某个天体的这个范围内时，只对这个天体放宽到 1.01 倍（`kLowOrbitClearance`）。排查“镜头拉不近”时看信息面板的 From camera 距离
- 事件的 `elevation_deg` 是直接替换俯仰角，`phase_deg` 只加在方位角上。要让相机在探测器背对某个天体的一侧（背景里看到那个天体），要的是方位角差约 180°、俯仰角取那个方向的相反数，否则那个天体不在画面里。天线朝太阳的探测器，太阳几乎在相机背后（相位很小）时，天线正对镜头，把机身全挡住（Cassini、Pioneer 10）：挑相位 50–60° 的时刻，天线才是侧面
- NASA 的模型：GitHub 仓库里的 GLB 是 Draco 压缩，网站上单独下载的不一定（New Horizons 没压缩）；`_root` 节点是一个 1 m 的占位方块，要排除；金属材质的颜色（金箔近乎纯黄）在这里的漫反射着色下太亮，要换；页面上的署名可能是公司（Hubble (A) 是 DigitalSpace），不要用；徽标可能是单独的材质，纹理里也可能有一整条横幅（Hubble 的 NASA 和 ESA），把模型用到的纹素画出来再决定怎么涂，并看一眼生成的图集；太阳能板分清哪面是电池片（Hubble 深色的一面才朝太阳，装反了测试也发现不了，要看画面）；单位不一定是米（Galileo）；颜色也可能是假的（Galileo 叫 black 的材质是亮紫色）；完全透明的材质（Juno）里的面和可见的面重合，工具会跳过；低面数的大扇形（Galileo 的遮阳板）要先细分，否则烘焙的阳光从被挡住的中心一路插值到边上，整块发黑；多数材质是 `doubleSided`，薄片不封闭，这里的管线又剔除背面，要复制反向的面；同一架探测器的 (A)、(B) 差别很大（JWST (A) 有 49 万个三角形且配色发灰，(B) 只有 10 万个），先看 `.png` 预览；纹理是照片时可能带 NASA 徽标
- 探测器向太阳的一面亮、背面几乎全黑，从背面看就是剪影。没有 `fill_light` 时 JWST 的主镜看不见（它从来不见阳光），加了以后要在展签里说明是示意的
- 只被吸积盘照亮的天体（Gargantua 场景），从盘的对面看是全黑的剪影，侧着看又只有掠射光：事件视角要让亮面朝向相机。Endurance 的环轴在惯性空间里固定，事件按“环轴和黑洞方向夹 67.5°”的时刻选，再用 `phase_deg` 选相机在哪一侧

### 运行、截图与 shell

- exe 旁的 `config.toml` 会影响截图（`aurora_view = 2` 时天王星、海王星是橙色的）：要可重复的画面，用 `capture_window.ps1 -Config <空文件>`
- 截图不要用 GDI+ 的 `Graphics.GetHdc()`（`CopyFromScreen` 内部也用它）：GDI+ 会把恰好等于 RGB(13, 11, 12) 的像素变成透明黑，在暗的光晕里看起来像一圈白色坏点。要把 PrintWindow 画进 `CreateCompatibleDC` + `CreateCompatibleBitmap`，`capture_window.ps1` 就是这么做的
- 带 `--event` 的截图偶尔落在相机还没飞到位的时候，整幅画面都不对（土星只有一点、看不到目标）：对比差异很大时先加长 `WaitSeconds` 重截一次，再怀疑渲染
- 从截图量亮度比要算上色调映射：ACES 在暗处的斜率约 1.6–2（d ln f / d ln x），海王星 ±23% 的晨昏渐变在截图里量出约 ±47%。方向可以直接量，幅度对照 CPU 的公式（`aurora_map_pattern`）
- 截图脚本偶尔丢最后一个按键：需要切换状态的键放在前面，截完核对。等待别超过 60 秒，空闲 60 秒后自动导览会接管相机；`--event` 超出范围时只打日志、停在默认视图
- Release 版是 GUI 程序，后台进程无法把它切到前台；computer-use 的 `open_application` 会再开一个实例，要用 `SetForegroundWindow`/`ShowWindow`。computer-use 和 `SendKeys` 发的 Esc 进不了 SDL，用 `PostMessage(hwnd, WM_KEYDOWN, VK_ESCAPE)`
- PowerShell 的变量名不区分大小写：`$S` 和循环变量 `$s` 是同一个变量
- PowerShell 里 `[IO.File]::ReadAllText` 等 .NET 方法的相对路径按进程的当前目录解析，不是 `Set-Location` 的位置：先 `Resolve-Path`。`Get-Content -Raw` 读进来的场景文件是 CRLF，`.Replace` 里写 `` `n `` 的多行串匹配不上
- PowerShell 里的 `/tmp` 是 `C:	mp`，不是 Git Bash 的 `/tmp`（在用户的 Temp 目录下）：两边交换文件用 Windows 绝对路径，`-ErrorAction SilentlyContinue` 会把拷贝失败藏起来
- 事件镜头只停展签的阅读时间（25–45 秒，按字数算），倍速要让要看的那一刻落在这段时间里：Janus 和 Epimetheus 的交换按 1 天/秒还没到最近点镜头就结束了。看节奏要用 Release：Debug 只有几帧每秒，单帧真实 dt 又被限在 0.1 s，仿真时间走得很慢
- 从地球看的盘（β Pic）不能把相机放到盘里面再拉近：相机在 150 au 时，近处的尘埃前向散射铺满画面，离得近的尾巴还被透视压扁（猫尾 480 au 深的尖端看起来朝恒星折回去）。放到 6000 au 外用 `observer` + 小 `fov_deg`，视线才接近平行。这样放大后银河贴图是一团团模糊的斑，事件里不显眼，但别用更小的视场。横贯画面的亮带经泛光会抬亮整幅画面（尤其全屏的毫米带），亮度要比单个小天体低得多
- 光学厚的盘，步长有两个反着的坑：按“密度随高度指数上升”解析放大步长、或者空旷处走大步，都会跳过盘面薄薄的发光层，整片盘出现同心的莫尔纹（发光层是 τ★ 和柱密度的乘积，厚度只有约 σ²/z）；步子太碎又会在到盘面前耗尽步数，耗尽的光线不报错，只是没有光（盘轴附近的 σ 趋于 0，远侧短轴方向出现一个黑楔形）。改步长后要俯视和斜视都看一遍，并看恒星周围
- 光学薄的盘，发光的温度要按到恒星的距离算，不能查按柱半径 r 的表：极轴方向 r 很小但离恒星很远，会在恒星上下画出两条白亮的竖条
- 发光的星云（`light = "emission"`，各向同性，不像散射光那样背向时变暗）把相机包在里面时，四周的壁层层叠加、铺满画面，盖住里面的东西：η Car 化石星风壳的镜头放在 Little Homunculus 里（4000 au）时只看得见一片紫色，退到它外面（6700 au）才看清。新的近景事件先算一下相机是不是在哪个星云的壳里
- 光学厚的盘多了多次散射后，很薄的高空尘埃也在发光（D 不像星光那样被挡），光线在“空旷”分支里走大步时会画出一层层台阶：判断“有尘埃”的门槛要够低（`kThickDustDepth` 1e-4），和着色器里的 `1e-4` 一起改
- 从密度变化很快的表里积分光线（后面各次散射的表）：源函数要按对数插值（进入阴影时按 e 指数降），线性插值会把阴影里的尘埃点亮；步子也要按高度取够细（0.5 倍高度时和蒙特卡洛差 3 倍）。靠近恒星的亮区光线找不到，要当成中心的一个点光源，朝向按平面的兰伯特定律
- 展签的半透明底会把画在它后面的东西压成灰色，`H` 也不隐藏事件展签：判断颜色前先确认目标不在展签底下（土星可见光极光曾被它盖着），落在画面下方的目标要调相机仰角
- 截图要特定视角时，可以在 exe 旁的场景副本末尾临时追加 `[[events]]`（`focus` + `from_body` + `phase_deg` + `distance_km`），再用 `--event` 打开；用完从 `assets/` 拷回原文件
- 在后台跑烘焙工具并把输出重定向到文件时，Python 会缓冲 stdout，几十分钟看不到进度：用 `python -u`
- Claude Code 的沙箱里访问 Horizons 会反复报 `SSL: UNEXPECTED_EOF_WHILE_READING`（不是服务器的问题），烘焙要在沙箱外运行
- 从 Bash 用 `&` 启动的 exe 会随 shell 一起退出，长时间运行用 PowerShell 的 `Start-Process`；后台命令里 `cd dir && a & b &` 只有 a 在 dir 里执行，下载用 `curl -o 绝对路径`
- Bash 工具的 heredoc 即使写成 `<<'EOF'` 也会吃掉反斜杠（`\\n` 变成 `\n`）：含反斜杠的编辑用 Edit 工具，或先用 Write 写成脚本文件。heredoc 里的中文传给 `python -` 或 `cat >> 文件` 也会变成乱码（写成 GBK，不是 UTF-8），改中文文档同样用 Edit；几百行、夹着 C++ 代码的 heredoc 还可能整段报 `unexpected EOF while looking for matching`，批量改代码的脚本一律先 Write 到 scratchpad 再 `python <文件>`
- WSL 里用 Clang 再建一个构建目录时，`linux` 预设把 `FETCHCONTENT_BASE_DIR` 固定在 `build/linux/_deps`，两个编译器会共用依赖的构建目录（SDL 的预编译头互相覆盖，报 `cmake_pch.h has been modified`）。Clang 构建要另给 `-DFETCHCONTENT_BASE_DIR`，再用 `-DFETCHCONTENT_SOURCE_DIR_<SDL3|GLM|IMGUI|STB|TOMLPLUSPLUS|DXC>` 指向已下载的源码
- 用 `python -I` 跑读下载文件的工具时，脚本自己的目录不在 `sys.path` 上：同目录的模块（如 `tools/philae/spk.py`）要在脚本里显式 `sys.path.insert`。本机的 numpy 装在用户 site-packages 里，`-I`（含 `-s`）找不到它：要 numpy 的工具用 `python -E -P`
- Windows 上 Python 的控制台编码是 GBK，打印 PDF 提取的文字会 `UnicodeEncodeError`：加 `PYTHONIOENCODING=utf-8`（grep 这类文件要 `grep -a`）。`python -I` 会忽略 `PYTHON*` 环境变量，这时在脚本里 `sys.stdout.reconfigure(encoding='utf-8')`
- 标签、信息面板用 ImGui 的默认字体，只有 Latin-1：天体名里的 ʻ（U+02BB）等字符会显示成“?”，名字用 ASCII（`1I/'Oumuamua`，和 Horizons 一样）；展签字体 Jost 另有希腊字母
- WSL 里只有 `linux-debug` 这类预设有的配置才有构建预设，写错名字时 `cmake --build` 直接报错退出；用 grep 过滤输出会把这个错误一起滤掉，旧的测试程序还能运行，看起来像编译通过了。看测试程序的时间戳。很久没用的 `build/linux` 可能让 ninja 断言失败（`RefreshDyndepDependents`），另开一个目录，用 `-DFETCHCONTENT_SOURCE_DIR_<...>` 指向已下载的源码
- 从 Git Bash 调 `wsl.exe`：`/mnt/c/...` 参数要加 `MSYS_NO_PATHCONV=1`；`bash -lc '...'` 里的 `$变量` 会被提前展开，脚本最好写成文件再执行

### 数据源

- **磁场模型贴着地表时不可信**：火星、月球的模型往下延拓到地表，小尺度被放大（火星的高频细纹，月球强异常周围大致南北向的条纹，后者在 120–240 阶做余弦衰减也去不掉）。贴图要放在接近观测的高度（火星 150 km、月球 30 km），向上延拓会自然压掉噪声，代价是场强小很多，展签要写明高度
- **Ganymede 磁场系数的坐标系**：论文的模拟用 GPhiO（y 朝木星），但系数读作 x 朝木星、y 朝后随面才和 Kivelson 写的磁极方向对得上。模型是纯偶极，开闭边界和极光却被木星的等离子体推得不对称，不要拿线的落点去对全球平均的极光纬度
- **火星磁场模型的许可**：Langlais 2019 和 Gao 2021 是 CC BY-NC（-ND），不能进 GPL 项目；Zenodo 上 Langlais 模型标的 CC BY 和论文对不上，也不用。Morschhauser 2014 有 Wieczorek 的 CC BY 版本（Zenodo 3876495）。libinternalfield 带着这些系数不代表许可没问题，它的 MIT 只管它自己的代码。上万个系数的模型文件按数据集对待；几十到两百个系数的论文表格按事实引用（JRM33 的论文虽是 CC BY-NC，见 `assets/geomag/SOURCES.md`）
- **磁场模型的经度和相位**：Connerney 系模型的系数按东经写（φ = −西经），论文和 Voyager 的经度却是西经，先用 atan2(h11, g11) 对照论文写的磁极经度。海王星的 NLS 不是场景的 IAU 经度（536.31 与 541.14°/天），它的本初子午线由“Voyager 2 在 1989-08-25 03:56:00 位于西经 167.7°”定义（PDS 说明），03:56 是真实的最近点；场景的重建轨道晚 6 分钟、远 4,000 km，用同一时刻比较会差 28°，要取轨道最近点的方向配上真实时刻的自转
- **事件的太阳视点**：`from_body = "Sun"` 在加载后可表示为 `from_sun = true`、`from_body = -1`。检查视点方向时用 `Scene::sun_position()`，直接用 `from_body` 索引天体会越界。
- **天王星的两套坐标**：ULS 的北是 1986 年受照的半球，也就是 IAU 的南（Melin 2025）；Melin 的“南极光”在 IAU 的北半球，Tiranti 说的 67.7°N 太阳直射点却是 IAU 的。ULS 西经和 IAUnew 东经同向，只差 50.6°（`IAUnew 东经 = ULS 西经 − 50.6`）。IAU 的自转率 501.16°/天比 Lamy 2025 的 500.93°/天快，累积到 2025 年差几千度，对照观测的经度前要先换自转率。论文的中央经线是光到达时的值，算场景里的视角要减光行时（天王星 2.7 小时，约 56°）
- **Melin 2025 的 L 壳表**：Zenodo 16085710 的 `JupiterMag_Uranus_L_shells.txt` 里 L 的标签是反的（notebook 把升序的等值线配给了 `['20', '10', '5', '3']`），标 3 的是紧贴磁极的最内圈，追踪 AH5 得 Rmax 约 19。用论文附带的等值线前先拿模型核对一下。场景的磁力线是 AH5，但 Melin 说暗带和 Q3、UH5 更吻合，图上也如此，所以这张图用 Q3；AH5 的 L 壳和表都不能混用。外推时也别凭图的观感：L 壳看着和亮斑很贴，留出法却显示亮度的大尺度分布主要由纬度决定。亮度和磁场的相关要做对照（把磁场转过去再算一次）：Q3 的第二个偶极赤道回路里 Br 反号，只看 |sin dip| 时转不转都差 ±5%。`uranus_h3p.txt` 的西经列 = 50.6 − ULS 西经，方向是反的，写反了是绕 230.6°W 的镜像，而南极光恰好在对称点上，看着还是对的；对照北极光的亮弧（70°W）。Melin 的 55.1° 是 JWST 的视线纬度，不是日下点（场景在观测时刻算出 58.83°，看的方向和太阳的方向不一样）
- **海王星的自转相位**：16.108 ± 0.006 h 每年累积约 73° 的不确定，Melin 2025 明说现在的相位未知，场景的 IAU 自转率（541.14°/天）又和他们的经度系统对不上。亮区只能在观测时刻按中央经线对齐到场景（差 −25.08°，`test_ice_giant_aurora` 打印），离观测时刻越远越只是示意
- **论文 FITS 的单位**：Tiranti 2026 的 `INTENSITY_PROFILES` 约 3e-12，按同一温度和密度算总发射只有约 4e-14 W m⁻³ sr⁻¹，和 Melin 2025 的柱亮度对不上约 70 倍，这里只用它的形状
- **极光亮度地图与 mipmap**：Marzok 2022 的 Ganymede 图用西经、南至北纬度，已按 70 km 壳层校正成垂直柱亮度；体积分直接恢复视线柱长，勿再做一次临边校正。创建带 mipmap 的采样器时要设置 `max_lod`，零初始化会把采样限制在第 0 层。
- **离心赤道的符号**：Milby 2026 式 1.6 的 θ_ceq 按字面是离心赤道的纬度，但他的表（东经 88.9° 处 Io 在 +0.232 RJ）说明它是加在卫星纬度上的，也就是卫星相对赤道的高度 ≈ r tan(θ + θ_ceq)；式 1.7 的写法和文字都不对。这里按表，磁纬的方向也如此（偶极朝东经 159.2°）
- **Io 的钠**：Milby 2025 表 2 的钠（17–29 kR）比 Schmidt 2023（前 10 分钟约 6 kR，之后 1.7 kR）亮一个量级，出影前也一样亮，氧却只差 1.6 倍。Milby 说钠有明显的系统误差。只能取一边的绝对值，不要拿两篇的钠直接比。Milby 的 25 kR BL1 里只有 0.5 kR 是谱线，剩下的 SO₂ 连续谱在蓝端比钠还亮，Io 在食中的颜色因此是紫红色，不是单纯的橙黄
- **Roth 2016 的左右和南北**：图 8a 的“右半边亮”固定在图像上，因为 HST 的相角总小于 11°，右边总是下午一侧，这里按太阳方向画，不是固定在本体上；食中它碰巧是拖尾一侧，只看光学观测（都在食中）会以为是拖尾。图 8c 的“极区/赤道”不能直接拟合：均匀的薄壳因为临边增亮就有约 1.4，观测才 0.9–1.2，论文也说球壳解释不了盘内亮度。“磁极亮”用 `(n·b)² − n.z²` 只取随磁力线倾斜的部分，它同时给出图 7a 的象限比；图 7e 食中的南北反转只能部分靠它
- **Marzok 2022 的上下游**：论文说的 downstream 是前导半球（90°W），upstream 是拖尾半球（270°W），等离子体从拖尾一侧追上来。第 3.3 节的 1.3 和 0.78 也按这个分；别按字面把 downstream 当成拖尾
- **木星统计极光卵**：Head 2024 的 CSV 经度是 System III 西经，纬度是行星中心纬度。着色器里的椭球单位坐标要先乘回极半径/赤道半径之比再归一化，否则极区轮廓会偏。主卵的 400 km 投影高度也不同于卫星足点常用的 900 km。

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
- **离扁行星很近的卫星**（Janus、Epimetheus 在 2.5 个土星半径处）：拱线每天进动约 2°，AXEPH2 的二体弧几小时就偏开，放宽容差也没用（50 km 时一年仍要 1000 个节点，调参考 GM 只省一半），所以用根数表（`element_table_bake.py`）。拟合进动椭圆时 a 和 n 要各自独立（扁率让 n 比开普勒第三定律给的快），进动速率要先拟合再固定（不带进动的两天窗口残差 35–50 km，带上后 3–4 km）
- **在旋转参考系里看快的相对运动**：参考系每转一圈，星空和土星的晨昏线也跟着转一圈（Janus 系 16.8 小时一圈，2 天/秒时每秒约 3 圈）。星空已自动淡出；镜头要从高处往下看（环的相角几乎不变），并让行星本体不入画。俯视时 `from_orbit_normal` 的屏幕朝向是任意的，要控制朝向用 `from_body` + `phase_deg` + 接近 90° 的 `elevation_deg`
- **NASA 图片说明的数**不一定对得上星历：PIA08170 说 Janus 和 Epimetheus“从不近于约 15,000 km”，Horizons 的中心距最近 10,235 km（各次 10,150–10,540 km）。展签按星历写，出处里注明
- **按容差抽稀的表**只保证位置，速度和加速度在节点之间可能很差：天问一号 15 分钟的入轨点火原来整段落在两个节点之间，三次插值把推力曲线扭弯，按表求出的 Δv 少 2%，从点火开始时刻的状态外推的幽灵 84 分钟后偏 110 km。要从表里读机动的就让点火期间每个样本都成为节点（`ephem.toml` 的 `keep_all`；已有的表用 `horizons_bake.py --patch` 只改点火那一段，其他节点不动；天问一号是 `tianwen1_bake.py` 的 `KEEP_ALL`）。读数窗口要在点火前后各留一点
- 从表里读 Δv 要按节点到节点算（`unexplained_delta_v` 用 `history_times` 取节点）：节点是精确样本，三次插值在点火突然开始、结束的那一分钟里过冲，按 4 秒一步累加会把来回摆的两边都算上（天问一号多 1%）。ESA 的轨道还常在点火结束处换文件、位置对不上（火星快车 0.9 km、金星快车 3.8 km、TGO 99 km 和 50 m/s），按 4 秒累加时金星快车多 13%；窗口终点落进台阶那一分钟时连按节点算也不行（TGO 按 ESA 新闻稿的 13:05 开始会多 1.1 km/s），这时用更早的出处时间或在台阶前截止
- Horizons 的“tag-up fit”轨道（希望号、MOM）是一段段拟合的弧，带着入轨点火的那一段可能是计划中的点火，下一段才是实际到达的轨道：希望号点火 8 小时后两段相差 2,397 km、126 m/s，场景里轨迹跳一下。用 `tools/bake/rejoin_burn.py` 以后一段为准，重新拟合点火（先核对哪一段与公布的捕获轨道一致：MOM 点火那一段的 Δv 与 ISRO 公布的一致，没有改）。往回积分要带上前一段里二体解释不了的加速度（火星扁率），只按二体积分 8 小时差 2.5 km、2.3 m/s
- 在巨行星近心点附近从星历读推力，要先扣掉扁率：朱诺号入轨时木星 J2、J4 的加速度（0.3–0.45 m/s²）比发动机（0.26 m/s²）还大，只按点质量算会把它当成推力，Δv 多出几百 m/s。点质量的 GM 也要用行星本体的（系统 GM 把卫星算在里面，近心点在所有大卫星的轨道以内）
- Horizons 的 Spitzer 在一些月初有位置接缝，均匀的两小时复查会漏过。烘焙后运行 `tools/bake/spitzer_seams.py`，保留每个月初前后各一秒的原始向量，避免插值把跳变摊到几小时；中间两秒是源数据的连接，不是机动。
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
- **GFZ 的 Kp**：`kp.gfz.de` 在本机 TLS 握手失败（沙箱内外都一样），用 GWDG 的镜像 `ftp.gwdg.de/pub/geophys/kp-ap/wdc/yearly/`。WDC 文件里的 Kp 是“十分之一”码（0、3、7、10、13…，即 0o、0+、1-、1o、1+），换成三分之一要四舍五入；同一文件里的太阳黑子数是 CC BY-NC，不能一起带走
- **临边曲线不是发射剖面**：临边扫描的亮度是沿视线积分的，切点低于发光层的视线仍然穿过它。直接拿曲线和发射剖面比，会觉得层的下侧很宽（火星的 UVD 在 23 km 还有峰值的一半以上），沿视线积分后才知道发射的下侧比上侧陡（σ 12 与 17 km），临边的峰也比发射的峰低 8 km。要比就都换成临边
- **滤光片测光换颜色**：只有几个宽带滤光片时，不要重建光谱再积分 CIE，滤光片之间的插值假设会让色相差很多（土星极光的高层从蓝变到品红）；直接取 RED、GRN、BL1 的 I/F。论文图上的“真实颜色”也可能是拉伸过的（Dyudina 2015 图 8 的框绿色为 0，测光却是红的一半），要和图里的测光核对
- **猫尾的星际介质拖曳**：Rebollido 2024 给了 ISM 的密度、风速和颗粒密度，按 Q = 1 算拖曳只把尾尖推几 au（≈ 2% 的辐射压）；拖曳和辐射压都按面积/质量时，只把尾巴整体压平，弯不出尾尖。论文的文字说尾尖是最小的颗粒：Q 在小颗粒时变小，β 趋平而拖曳继续涨，这里用 β⁸ 模拟
- **从论文的图取数**：Portilla-Revelo 2023 的图是位图（PDF 转 SVG 后取 base64 的 PNG，1280×960），Keppler 2019 的图 2 是矢量（SVG 的 `stroke` 路径，属性写在标签上不在 `style` 里），两者都比目测好；Keppler 图 2 的横轴第一点是 0.05″，不是 0
- **论文先看 `ref/`**（不入库）：用户下载过的 PDF 放在这里（包括 Wiley 和 HAL 打不开的），上网找之前先列一下；用 Git 自带的 `pdftotext -raw` 提取。扫描件没有文字层（`pdftotext` 只输出分页符），用 PyMuPDF 把页面渲染成 PNG 再看（`python -E -P`，`fitz.open(...)[i].get_pixmap(dpi=110).save(...)`，先拼一张几十页的缩略图定位）；文件名里的弯引号在 Git Bash 里打不开，先用 PowerShell 拷到 scratchpad。
- **期刊网站的验证**：SWSC（swsc-journal.org）对 curl 和 WebFetch 返回 captcha 或 403，内置浏览器能打开 full_html，表格在单独的页面（`.../T2.html`）；AGU/Wiley（agupubs.onlinelibrary.wiley.com）连内置浏览器也过不了 Cloudflare 验证，换开放的转述或同作者的开放版本（Sigernes 2011 的 Starkov 系数用 uib 上的 Opt. Pura Apl. 版本，表 1 有一处印成“-4-44”）
- **电流片两翼的磁力线**几乎是径向的（Ganymede 处，电流片倾斜后离开卫星几个 RJ）：卫星的纬度差零点几度，整条线的 L 就差很多，足点却几乎不动。线的表不能当线性的，验证要按卫星实际的轨道范围（足点表的 ±3%、±2° 太宽，会差到 8 万 km），比较两条线也要用点到折线的距离，按下标比较会把沿线的滑动算成误差
- **Jena 的宇宙尘埃光学常数库**（OCDB，`www2.astro.uni-jena.de`）：curl 报 TLS 握手失败（沙箱内外都一样），内置浏览器也打不开。WR 104 的颗粒反照率改从 Harries 2004 图 3 读（他们用的就是 Zubko 1996 的无定形碳常数）
- **Stooke 形状模型 V2.0**：小行星和彗星的经度表改成了东经（360° 减 V1.0），但标签里的说明文字照抄 V1.0，其中的经度数值（如 Halley 的 270°）仍是西经。Vega 电视图像的时间有 UT（`TIM--OBS`、`TIME-ENC`）和开机后的星上计数（`TONBHRDW`，4.766485 计数/秒）两种，Stooke 标签里的“2:00:30”看来是后者
- **Chandra 的 Horizons 简介**：近远地点高度把英里的数值标成了 km（5,999 / 86,487）；1999 年入轨公告给的是约 9,700 / 139,200 km。用状态向量和任务中心原始公告核对，勿抄简介的单位。
- **带颜色的残骸层别用绝对 I/F 叠到贴图上**：贴图的亮度没有校准，木星南纬 44° 的棕色带线性值只有 0.2–0.4，和残骸层的 I/F（约 0.24）差不多，加上去疤痕就消失了（West 的 0.23 对的是 955 nm 里 0.36 的周围）。用相对云层的比值。比值也要按可见光算：残骸在 540 nm 的光学深度是 955 nm 的近 3 倍，955 nm 的比值（核 0.64、新月形 0.84）比可见光（0.50、0.74）淡。West 图 3 的两个连续谱波段是拟合 k 的依据，他的图 4 的 k 在 955 nm 要乘 2.7
- **位置角的东向**：天球上东 = 北 × 视线（从地球指向天体），写反了位置角变成 360° − PA，尾巴和碎片列的朝向都会差很多，但两者的夹角还是对的，不容易发现。对照 Sekanina 图 5 才看出来
- **SL9 的闪光和等级无关**：直觉上大撞击闪光也亮，实测不是：K（1 级）比 H（2a）暗，N（3 级）也有 L 的四分之一，七次只差 4.6 倍（Martin、Chapman）。别按撞击点等级算闪光。论文之间也有出入：Hammel 1995 引的 G 是 3%（预印本），Martin 的 L 是 6%；Carlson 的 t = 0 比 Martin 的闪光晚约 5 秒，要加上。羽流亮度不要拍脑袋加显示增益，Hammel 1995 给了 I/F（G 0.09，“相当暗”），Chapman 给了上限
- **从地球方向看 SL9 的撞击**：相机要放得很远、视场很窄（100 万 km、2.4°，`aim_offset_km` 对准撞击点）。放近了透视会把地平线往里拉，撞击点比 6.4° 藏得更深。场景里没有地球和 1994 年的 Galileo，方向是 `tools/sl9/impact_views.py` 算好写死的。论文的引文别想当然：常见的“比全世界核武器还多的能量”这类说法，这里的论文里没有，不写进展签
- **Shoemaker-Levy 9 的记录**：Horizons 不会让碎片停在木星上，向量穿过去继续走，页头 “Valid … to” 也不是撞击时刻（预报的、地球接收时间，比实际晚约 42 分钟）：要自己找入射穿过 100 mbar 面的时刻。按 `DES=D/1993 F2; CAP` 列出的只有 18 条，Q1、P2、Q2 要直接查，“P” 其实是 P1。Chodas & Yeomans 表 5 的时刻也是地球接收时间，表 6 的“垂直”是椭球法线（和径向差 3.2°），速度在随木星转动的参考系里。解往回推到 1992 年的近木点时碎片相距 9 万 km，不是真的解体时的样子。从表 4 的吻切根数只带太阳和 J2、J4 外推，71 天后差约 25,000 km（约 20 分钟），但碎片之间的时刻差只差几十秒
- **Horizons 不给的数据**：Didymos 本体（920065803）和 Dimorphos（120065803）的说明写着覆盖 2000–2030，可 2026-10 起查 2026-07 以前什么都不返回（没有数据也没有报错，大概是质心换成 ESA 的 Hera 文件后断了链）。解 547 的 SPK 在 figshare（`tools/dart/s547_bake.py`），是 ECLIPJ2000 坐标的 5 型段，要转到 ICRF
- **DART 的最后 48 秒**：Horizons 在导航重建后接了一段到撞击的 `spc` 段，起点偏 6.8 km、方向偏 1.3°。重建本身在撞击时刻离 Dimorphos 中心不到 2 km（早 0.3 s），场景用重建并整体平移（`tools/dart/approach_bake.py`）
- **卫星的幽灵**：要让外推的周期等于撞前的周期，按论文的 J2 和系统 GM 还差 0.1%（Dimorphos 每圈慢 49 s，三天后偏 213 m，是撞击效应的一成），因为撞前的轨道是拟合的圆，GM 本身有 4% 的误差。系统 GM 取 n0² a³ / (1 + 3/2 J2 (R/a)²)，在 GM 的误差内，三天只偏 1 m。Didymos 的 J2 让 n 比开普勒第三定律给的快 1.7%：画轨道用的 `gm_km3_s2` 要取 n0² a³，否则吻切周期和偏心率是错的（12.3 h、0.015）
- **PDS 70**：Trevascus 2025 的 i > 90° 是顺时针，和目视轨道的约定一致，近侧要算出来核对（Müller 2018：西侧近且亮）；后验中值拼不成一条轨道，同一列也可能差十几 mas，逐列对照 GRAVITY 再挑。星等 2404.11641 表 4 的三区模型不适合画（大颗粒区 r^−4.88，亮度在 R_in 而不在 ALMA 的 75 au；他们自己也说小颗粒质量多了上千倍），Portilla-Revelo 2023 只有图没有表，柱密度从图 A.2 读（用 Inkscape 渲染 PDF 页）。A&A 的网页（aanda.org）对 WebFetch 返回 403，用 arXiv 的 PDF
- **β Pic 的副盘**：Golimowski 2006 表 2 的 +5.9° 和 −4.7° 是两翼各自相对主盘的角，一正一负是因为两翼方向相反，整体是一条约 5° 倾斜的直线（Apai 2015 叫 CCW tilt），不是 V 形。Rebollido 2024 第 5.2 节的“位置角差 2.5°、离侧向 6°”是偏心盘模型的参数，换成圆盘时脊线只斜约 2.5°，要按量到的 5.3° 取；近侧反过来（i = 96°）才能让盘和 CO 气体同向转、猫尾从远侧出发。Matrà 2019 的 a 在 Eq. 2（窄高斯的质量份额）和表 1 的注释（峰值比）里定义不一致，按 Eq. 2
- **论文**：PDF 表格用 `pdftotext -layout` 提取时列会错位，用 `-raw`；`tools/pdftext.py` 会丢掉表格里的负号（Iess 2018 表 1 的 J4 提取成 586.609，实际是 −586.609），符号要核对；有的书的 PDF（《The Science of Interstellar》）用 `tools/pdftext.py` 只提取出空格，Git 自带的 `pdftotext`（`/mingw64/bin`）可以；本机没有 `pdftoppm`，Read 工具不能渲染 PDF 页面，看插图用 Inkscape：`inkscape --pdf-poppler --pages=<页> --export-type=png <pdf>`；凌星拟合的 ω 约定各代码不同（NbodyGradient 凌星时 ω+f = 270°，Phodymm 是 90° 且状态取反），场景里用 `transit_u_deg` 写明；共振链对积分相位误差很敏感，用 6 阶 Yoshida
- **轨道重建**：平动点轨道在刀刃上，按稳定流形二分出的机动方向只对同一条积分路径成立（换一种分段方式调用积分器，1e-9 的差别一个月后就是上万 km），求解和生成星历要走同一条路径。瞄准月球低空飞掠或入轨时，数值近月点的时刻在掠过月面附近会跳变，数“第 n 个近月点”也会整圈跳：先用十几小时前的密切二体轨道瞄准，再用各时刻的密切平近点角细调
- **共振轨道**：星下点每几圈重复的火星轨道（天问一号的 2 sol 停泊轨道、1/3 sol 中继轨道）受田谐项共振影响，几十天后的相位对初值极其敏感（02-23 的 Δv 变 1.5 m/s，80 天平均周期在 48.8–49.6 h 之间乱跳），不能把一次机动拟合到几十天后的向量上；用前后两个向量的弧段在中间交接
- **食双星的历元**：论文给的 T0 不一定是它说的那个时刻。Baron 2012 表 6 注明 Algol 的 T0 是“主极小”，实际是内轨道过升交点（主极小在其后 3/4 周期）。圆轨道的目视轨道 T0 常是过节点，要用观测的极小时刻（GCVS、TESS）和他们自己的位置表核对；几十年的外推还要考虑 O−C 漂移，用最近的观测拟合（`tools/algol/`）
- **直接成像的候选行星轨道**：论文只给后验的均值（Beichman、Sanghi 2025 的 Table 4 只有 a、e 和倾角），各量的均值拼不成一条轨道；完整的轨道在作者放到 Zenodo 的后验样本里（orbitize! 的根数，tau 从参考历元按周期的分数量起）。拿样本先算族的均值和 Table 4 核对，再挑一条：离均值最近的那条与观测差 2σ（OFTI 的样本按似然接受，不是最佳拟合），要把与观测的 χ² 也算进挑选标准
- **轨道表的 τ 历元和坏行**：Kral 2025 表 2 没写 τ 从哪天量起，用表 1 的位置试几个常见历元（MJD 58849、59000）立刻能分出来（差 25 mas 和 1.7 mas）。他们表 1 里 c 的第一行（59454.125）离任何轨道都差 27 mas，Hinkley 2023 同一时刻的值才对
- **论文里的数自相矛盾**：摘要和表格可能不一致（Anderson 2015 摘要的 e = 0.647，表 3 是 0.674）。用表里别的量反算核对（a₁ sin i 只和 0.674 相符），以表为准
- **半径的口径**：充满洛希瓣的星，论文给的“半径”可能是某个方向的（Kolbas 2015 的 Algol B 3.43 R☉ 是侧向半径，按他们自己的 q 和 a 算），不是等体积半径；由 v sin i 加同步自转假设推出的半径是上限（Algol A 的 2.87 R☉）。拿来和别的来源的质量、倾角混用前，先用洛希面换算到场景自己的 a。光变解（Wilson-Devinney）的倾角和半径是一起拟合的，单独换掉一个（例如用干涉测量的倾角）会让食深偏离很多
- **ADS 的老文章扫描**：`articles.adsabs.harvard.edu/pdf/<bibcode>` 用 curl 只能拿到第一页，其余页要令牌。在内置浏览器打开 `ui.adsabs.harvard.edu/scan/manifest/<bibcode>`，在页面里用 `/v1/accounts/bootstrap` 返回的 `access_token` 取 `/v1/scan/image/iiif/2/bitmaps-~seri-~<刊>-~<卷>-~600-~<页>.000/full/1600,/0/default.jpg`（如 `ApJ__`、`0342`、`0001065`），放进页面后截图读
- **TESS 上没有亮星的光变曲线**：TIC 里亮星的 Tmag 可能是错的（δ Cep 的 TIC 415741355 写成 14.0 等），没有任何管线给它做光变曲线；MAST 按位置搜到的同名附近星（HD 213307）不是它，先按 HIP 查 TIC 号
- **被尘埃包着的亮星，Gaia 视差不可信**：WR 104 的 DR3 视差 1.56 mas（0.64 kpc），旋臂间距和风速给的是 2.58 kpc（DR2 是 0.24 mas）。这类星用论文里的几何距离，Gaia 的值只在 TOML 注释里提一下
- **WR 140 的论文**：Han 2022 扩展数据表 2 印的 Ω = 53.6° 是排印错误，他们的代码（github.com/yinuohan/WR140 的 `Run_spiral_visualisation.m`）用的是 Monnier 的 353.6°；代码的轨道约定（X 北、Y 东、Z 远离观察者，R_z(Ω) R_x(i) R_z(ω)，锥轴是 O 相对 WR）和场景的目视轨道一致，可以直接拿来对照。Han 的速度、加速度和半径是按 1.67 kpc 从角运动换算的，换距离时一起缩放。Lau 2022 图 3 的壳在图上朝西南张开，但按柱密度加权的方向是 165°（Han 的代码同样），对照方向时不要凭图的观感
- **η Carinae 的轨道朝向**：Madura 2012 用 (i, φ, PAz) 定朝向（φ 是视线在轨道面上的投影与远星点一侧的夹角，PAz 是 +z 轨道轴在天上的位置角），换成场景的目视轨道根数是 i 不变、node = PAz − 90°、伴星的 ω = 90° − φ（他们写的 ω = 270° − φ 是主星的）。各家的 ω、e、近星点时刻互相差得不少（Strawn 2023 表 1），场景里写明各取自哪篇
- **观测做的三维星云模型**（ESO 的 Homunculus OBJ）：原点不一定在中心星（这个差了约半个瓣宽），单位也不一定是物理单位。中心按网格的点对称找，尺度按论文里的尺寸标定，在天空中的朝向对照论文的图拟合，做法写在 `assets/shapes/SOURCES.md`
- 检验分段的运动方程（WR 140 尘埃的加速度在 r_nuc 和 r_t 处跳变）不要用跨过跳变的 RK4 积分：跨跳变的那一步速度差 a·dt，之后位置误差线性增长（0.1 天步长 60 年后差 0.3 au）。在每一段内用中心差分检验加速度
- **WR 104 的论文约定**：Tuthill 2008 的 θ0 是阿基米德螺线原点处的位置角（拟合时把投影的旋臂当平面螺线）；Richardson 2026 图 6 的交点 78° 是 xenomorph 自己的约定，和场景的目视双星约定（16.12°）对不上，不能直接搬
- **星表距离和场景距离不一致**：场景采用的视差（例如轨道视差）和 Hipparcos 不同时，星表里的“自己”会离观察者超过 `kOwnSystemPc`，在天空里成为一颗极亮的星（Algol 差 1.25 pc，−4.6 等）：用 `[sky] hide_hip` 去掉
- **KSP / kOS**：kOS 脚本的变量和参数名不能遮住内置函数（`v`、`r`、`q` 等，报 CLOBBERBUILTINS）。kOS 的向量在 Unity 的左手系里，y 轴朝系统北极，原点在飞船，低空时坐标轴还会转：换成场景坐标要 y、z 互换，x 轴取同一物理帧里采样的 `SOLARPRIMEVECTOR`。天体的 `VELOCITY:ORBIT` 是相对飞船所在 SOI 的速度，不是相对母天体，要和母天体的相减。游戏里的根数多是 float32（e = 0.2 读出来是 0.200000002980232，历元平近点角的配置值是 3.14 rad 而不是 π），照抄读出的值。历法常数在 `Assembly-CSharp.dll` 的 KerbinDay/KerbinYear getter 里（`ldc.i4 21600; ret` 和 `ldc.i4 9201600; ret`）
- **SPK 摘取的分段范围**：jplephem 的 `excerpt` 可能把每一段的摘要时间都写成整个请求区间。读取跨段文件时，还要按 Chebyshev 记录的起点、步长和记录数限制有效范围，否则会在段切换处报 `OutOfRangeError`（Proteus / NEP098）。
- **时间和距离**：MJD 转 JD 要加 2400000.5；新闻稿里的时间是地面接收时间（差光行时），对照测试只比较距离；NASA 的“离地球多远”和飞掠高度从地表算；地面成像给出的环弧、卫星经度是光离开天体时的位置，对应的时刻要从观测 UT 里减去光行时（海王星约 4 h，Adams 环弧段转 137°），论文按 UT 算的时间间隔还会少掉其间的闰秒
- **借力飞掠的“速度改变”**：新闻稿的口径不一。卡西尼号 1998 年的“about 7 km/s”是相对金星速度的转折 |v_out − v_in|，伽利略号、信使号的是日心轨道能量的变化折合成金星距离上的速度（信使号 2007 年转折 6.95 km/s，能量口径 6.64 km/s，即“nearly 15,000 mph”）；对照前两种都算一下

- **脉冲星进动的倾角约定**：`visual_orbit` 的 north/east/depth 基底给出 `h·K = -cos(i)`（K 从地球指向脉冲星）；Weisberg & Taylor 2002 的 i 定义为 h 与 K 的夹角。B1913+16 的 WT02 小倾角解要在 TOML 写补角 132.8746°，并用 `cos ζ = -cos Φ sin δ sin i - cos δ cos i` 检查符号。有限光束的闪光需要在边缘降到零，高斯尾巴会让离开视线后的脉冲继续闪。近看的中子星已触及 `SunPass` 的恒星亮度上限，射电脉冲的示意增亮要在限幅后乘上，否则改星等也看不出闪烁

### 测试

- 木星足点烘焙：JupiterMag 1.5.1 的 `ionosphere` 插值在部分完整轨迹上返回 NaN 或不收敛的坐标，直接从 `x/y/z` 轨迹求与目标椭球的交点，并用独立积分器核对。Hue 2023 领先角拟合有 0.98/0.99 倍频项，不能简单取模后直接拼接首尾；在西经 0° 附近平滑衔接，其他经度保留原式。

- Body Template 在 TOML 层合并后才构造天体。切换自转模型时用 `unset` 清除互斥字段；`body_gm_km3_s2` 是自身质量，`gm_km3_s2` 可能含卫星质量，迁移时分别核对。跨场景存在的自转简化要在场景里保留（如地月场景省略月球的长期极轴漂移）。
- 场景测试里查光照或位置前要先 `scene.update(t)`，否则所有天体都在原点
- `Scene::icrf_state_at` 是 ICRF 里的绝对位置（相对场景的根，通常是太阳），卫星绕主星的角度和轨道面要先减去主星：用日心状态定轨道面时，退行的 Dimorphos 在其中一直倒转，“扫过 n 圈”的循环永远停不下来
- GR 的近心点方向要对 r(φ) 做抛物线拟合（节点角间距和每圈进动同一量级）；IMR 是辛积分器，能量误差有界但不精确守恒
- 黑洞光线是从相机往回追踪的（时间反演），测试里自旋不对称性的符号要反过来看
- Kerr-Schild 坐标里协变动量 p_i 不是光线方向（dx^i = p_i − f (l·p) l_i，150 M 处相差约 2/r），屏幕方向也不是角度，要经观者标架（`static_observer_frame`）。a = 0 时它的空间坐标就是史瓦西的面积半径坐标，弱场“从 R 到无穷远偏折 (2/b)(1 − cos)”是各向同性坐标里的分法，拿来接积分段会在 b = R 处错开约 0.013 rad；用 `weak_leg_angle`
- 光变曲线的倍速淡出和采样开关共用 `light_curve_visibility`；采样放在导览修改时钟之后，暂停或减速恢复显示的当帧就要补齐，不能显示旧时间窗口
- 光变曲线在食外也不为 0（洛希形状和照射带来的起伏）：判断是否在食中，要比较有遮挡和无遮挡时的光；食深从最亮处量起，不要从 0 量
- 导演测试用多个随机种子长时间跑，抓相机穿过天体、掠过黑洞这类问题；带滚动窗口积分的场景要用低倍速，否则 Debug 下很慢

- **大气折射红环**：地球模板的 `refractivity` 控制折射，`refraction_gain` 单独控制显示增益。加载时沿弯曲光线积分成 128 点表，`AtmospherePass` 查表并按像素径向宽度积分；只取像素中心会让不足一像素宽的红环断裂。月球旁的 `observer_offset_km` 也要满足镜头离表面至少 0.1 个半径的限制，2022 月食事件采用 200 km 高度。散射、夜光和极光继续使用已有大气模型
- **日冕范围**：K 日冕的几何结构、2017 / 2024 太空视点日食、Parker 增强白光和 2022 月球附近的折射红环已实现；不做 CME。若引入黑子数序列，先核实 SILSO 的许可证

- **日冕遮挡**：月球在太阳半径单位里很小，射线与球的判别式用 `r² - |cross(center, ray)|²`，不要用两个大数相减的 `b² - c`，否则月缘会因 float 精度破碎。半分辨率积分也要在遮挡轮廓附近全分辨率重算，避免日冕漏进黑色月面

`[corona]` 启用白光 K 日冕，`display_gain` / `display_gamma` 是显示增强参数。`CoronaPass` 先在本输出 `SceneTargets` 的半分辨率单采样纹理里积分，再叠加到主 HDR pass；MSAA 切换要同步重建合成管线。事件的 `observer` + `observer_offset_km` 把相机放在指定天体旁的固定 ICRF 偏移，`fov_deg` 控制视场；这类镜头滚轮改视场，拖动解除锚定，切换普通焦点恢复 45°。`corona_enhanced` 用于 Parker 的增强视图。`[view].event` 按名称选择初始事件。散射、遮挡、日食本影位置及相机行为在 `appearance_tests.cpp` 检查。

- **Europa 氧极光**：`[aurora].europa_visible` 启用 2021-05-20 光学氧线参考值，`satellite_aurora.hpp` 的 `aurora_disk_factor` 用可见球壳体积积分将孔径盘平均 Rayleigh 换为柱亮度。不能直接把盘平均数值当作天顶柱。着色器用指数层（150 km）、逐点母行星遮日（`LIGHTING_NIGHT_OR_ECLIPSED`）和 `PATTERN_EUROPA`（`europa_pattern`，CPU 上同样的公式，改的时候两边一起改）：下午一侧亮（跟着太阳）、面向等离子体片的半球亮（`MoonAuroraState::sheet`，拖尾一侧更明显）、沿穿过卫星的磁力线集中（`MoonAuroraState::field`，放在 `AuroraDrawItem::pole`），系数由 `tools/geomag/europa_pattern_fit.py` 拟合 Roth 2016 的图 7、8。亮度乘 `MoonAuroraState::plasma`：离心赤道（`centrifugal_equator_height_rj`）上的高斯，标高 1.7 RJ，以当晚的平均为 1；图案除以 `europa_pattern_mean`，只搬动光，不改变总量。可见光事件、手动可见光和 UV 叠加视图均显示其可见光发射，曝光与物理谱线分离。来源和假设见 `assets/geomag/SOURCES.md`。

- **Io / Callisto 可见极光**：模板 `io_visible` / `callisto_visible` 与其他极光字段互斥，`satellite_aurora.hpp` 将观测 Rayleigh、谱线颜色、假设层高和显示增益分开保存。三颗卫星共用遮日与指数层，但各有横向分布，Io 钠层是不受图案调制的第二层。Io 跟着日食变：`Simulation` 找出此刻所在的一次入影和出影（`io_eclipse`，时钟离开这个周期才重算，超过 2 小时/秒不算），钠（`io_sodium_level`，第二层）和 SO₂ 的蓝色连续谱（`io_so2_level`，第一层的一部分）乘相对参考观测的水平，氧和硫不变。Callisto 跟着等离子体片变（标高 3.6 RJ，参考是赤道上的亮度）。导览用 `aurora_visible`，手动 UV 叠加仍保留可见光；默认普通事件不显示。事件日期与发光层本影范围、模板冲突、各层孔径归一化在 `appearance_tests.cpp` 验证。

- **火星可见光极光**：模板 `[aurora].mars_visible`，与其他极光模型互斥。`mars_aurora.hpp` 用 Knutsen 2025 的轨道临边估计 2.1 kR 归一化两侧不等宽的高斯层（Schneider 2018 图 2 沿视线拟合）；101 R 是地面含尘埃消光的斜视观测值。2017 年 9 月的亮度按 Schneider 图 3 随日期变（`MarsAurora::limb_rayleigh_at`，经 `MoonAuroraState::mars_limb_rayleigh` 进 `aurora_profile`），峰值 21 kR 由紫外推算，其他日期是参考值。事件和手动可见光／UV 叠加均可启用；来源、层高和曝光假设见 `assets/geomag/SOURCES.md`。

## 待办

按容易程度排：

1. Kerbol 场景：高度图和法线（晨昏线上的立体感，也可用于月球、火星的 USGS DEM）与海面高光（Kerbin、Eve、Laythe、地球）；还可以用 `kerbol_transfer.py` 的办法加别的示意航线（Minmus、Duna 的转移窗口等，标明是设计的）
2. 黑洞附近的天体做透镜：Gargantua 的示例行星和 Mann 星现在按真实位置画，转到黑洞背后时被挡住（真实情况下光会弯到阴影上方）。可以对每个天体解相机到天体的类光测地线，把球和标记画在主像（和次像）的位置
3. 磁场可视化后续：① 磁层边界（只放进少数事件：天王星压缩与常态对比、Cassini 和 Galileo 测木星磁层收缩、ISEE-3 穿过磁尾、2024 年 Kp 9，水星最需要它，地球可用 Shue 1998 和 OMNI）；② 地球的 IGRF 随年份变化（南大西洋异常的漂移和分裂）；③ 卫星到足点的线：加二级亮斑（反射波），用有依据的等离子体片模型替换 Europa 和 Ganymede 拟合不出的标高，让 Ganymede 的线和它自己的磁层接上
4. 空间望远镜星历更新：Roman、TESS、Euclid、JWST 的表大半是预报，Chandra 的预报到 2027-03-30，过几个月重新烘焙
5. DART 场景加 ESA 的 Hera（−91）：2026 年 10 月底到达 Didymos，Horizons 现在只有到 10-17 的跟踪数据和到 11-04 的预报（离 Didymos 约 30 km）。等 ESA 发布近距离阶段的重建轨道再烘焙，同时加两颗立方星 Milani、Juventas，以及 Hera 测出的撞击后 Dimorphos 的形状和撞击坑（Naidu 2024 推断长轴比约 1.3；形状模型发布后替换或在撞击时刻切换）。届时 Didymos 和 Dimorphos 的表要延长到 2027 年：解 547 的 SPK 覆盖到 2024 年底，2026 年起要用 ESA 的新解（Horizons 的 2026-07 以后已经换成 ESA 的文件）；DART 撞击也改变了系统的日心轨道（量级很小），新解里应该有
6. 土星可见光极光后续：加随本体转动的亮弧（Dyudina 2015 的 10.65 h、10.8 h，1 万 km 尺度，持续 100 小时以上；Bader 2018 的统计调制只有 4–6 kR，真实结构的幅度论文没给，要标成示意），事件倍速要提到几百倍才看得出转动，紫外卵用同一套调制。相位表（Leicester 档案，2004–2017）的许可要先查
7. 伽利略卫星极光后续（都在等观测，2026-11 查过）：Callisto 有空间分辨的观测后换掉均匀壳层（HST 17739 是无分辨的 COS 光谱，发表后可以核对紫外是否在食中消失）；Io 的 SO₂ 连续谱有穿过日食的时间序列后，核对它是否跟着 ALMA 的大气主体走；Europa 食中亚木视向的南北比（Roth 2016 图 7e 是反的，图案复现不了）有更多食中图像后再定
8. 火星可见极光后续：Perseverance 的第二次探测等论文发表（EPSC 2025 摘要没有日期和亮度），有 297.2 nm 或绿线自己的高度剖面后换掉 UVD。尘埃消光不用做：从轨道看只影响临边约 40 km 以下，2.1 kR 也已经扣掉了尘埃
9. 冰巨星红外极光后续（都在等观测，2026-11 查过还没有）：海王星有高度剖面后换假设的 Chapman 层，有地方时的观测后把线性的晨昏渐变换成真实的曲线；天王星 25°S 以南是外推，有夜侧的观测或模型后核对压暗的三个定点
10. 盘的后续：PDS 70 的 ¹²CO（光学厚，像温度图）和 Mie 按波长的颜色；HD 206893 换成 Imaz Blanco 2023 更高分辨率的缺口，等 GRAVITY 确认或否定 B I。HD 100546 没有确认的行星和旋臂的图案速度，要旋臂先看 MWC 758（Ren 2020 测到 0.22°/年）
11. 探测器模型的下一批：NASA 3D Resources 里还有 Chandra、Roman、Spitzer、Kepler、TESS、WMAP、Rosetta、Dawn、MESSENGER、MRO、MAVEN 等；Huygens 有单独的模型；整个探测器自旋（Juno、Pioneer 10、Galileo）还没有做。烘焙的阳光只适用于一直用同一面朝太阳的，不是这样的要加不烘焙的选项，或者运行时算阴影。`ref/probes/` 里已经下载了 Europa Clipper、TGO、IBEX（`src.txt` 有出处，IBEX 那行的链接指向 Rosetta 的页面，要核对；Huygens 是 JPL SPACE 的 DXF，只有几何）。ESA 的望远镜、天问一号等没有找到模型，可以照 `make_endurance.py` 自己生成简化造型。观感上可以再加高光（金色主镜和隔热膜）
12. 远期：macOS（Metal，需要 SPIR-V → MSL 和 macOS CI）；J1141-6545 的自转轴进动（Venkatraman Krishnan 2019 表 1：α ≈ 89°，δ 各模型从 35° 到 155° 不等，倾向的模型 3 是 91° ± 60°；Manchester 2010 约 110°），等几何解定下来再做

## 不要做

- 不要手写多套图形后端（GL/DX/VK）—— 只用 SDL_GPU
- 不要为了“通用”引入 ECS 或复杂的引擎层；保持小而直接
- 不要让 `core/`、`ephem/` 依赖 SDL 或 GPU
