# Astraxis

A quiet, physically grounded astronomy visualizer meant to be left running in the background — something to stare at.

## Scenes

| Scene | What you see | Motion source |
|---|---|---|
| Jupiter | Jupiter, its faint rings, the Galilean and four inner moons, and the orbits of Galileo (1995–2003) and Juno (2016–) | JPL Horizons ephemerides, JPL mean orbital elements |
| Voyager | The Sun, planets and both Voyager probes from 1977 on | JPL Horizons ephemerides |
| JWST | JWST's halo orbit around Sun–Earth L2 | JPL Horizons ephemerides |
| Earth–Moon | Artemis II's free-return flight, Artemis I's distant retrograde orbit and CAPSTONE's near-rectilinear halo orbit, in the Earth–Moon rotating frame | JPL Horizons ephemerides |
| Saturn | Saturn's rings (Cassini radio-occultation optical depths), the seven major moons, Cassini's tour and Huygens' descent to Titan | JPL Horizons ephemerides, JPL mean orbital elements (phases fitted to Horizons) |
| Parker Solar Probe | Seven Venus gravity assists step the perihelion down to 9.86 solar radii; petals in the Sun–Venus rotating frame | JPL Horizons ephemerides |
| Alpha Centauri | A, B and Proxima | Newtonian N-body integration |
| Sgr A* | The Galactic Centre black hole and the star S2 | Kerr geodesics + GPU ray-traced black hole |
| TRAPPIST-1 | Seven Earth-sized planets in a resonant chain around an ultracool dwarf | Newtonian N-body integration from transit-timing fits |
| Kepler-223 | Four sub-Neptunes in a 3:4:6:8 resonant chain | Newtonian N-body integration from transit-timing fits |

Scenes are described in `assets/scenes/*.toml`; adding one needs no code changes. Physical constants and orbital elements cite their sources in the scene files, and data licenses are listed in each `assets/*/SOURCES.md`.

## Building

All dependencies (SDL3, glm, Dear ImGui, toml++, stb_image, and Microsoft's DirectX Shader Compiler for HLSL → DXIL / SPIR-V) are fetched by CMake at pinned versions.

**Windows** (Direct3D 12): Visual Studio 2026 with its bundled CMake.

```bash
cmake --preset win-msvc
cmake --build --preset win-msvc-debug
```

Run `build/win-msvc/Debug/astraxis.exe`. Unit tests: `build/win-msvc/Debug/astraxis_tests.exe`.

Windows MSVC non-Debug builds enable [VC-LTL5 v5.3.1](https://github.com/Chuyu-Team/VC-LTL5/tree/v5.3.1) by default, using the Windows 10 system `ucrtbase.dll` to reduce the runtime footprint. CMake downloads the pinned binary package and verifies its SHA-256. To build the CI Release configuration locally:

```bash
cmake --preset win-msvc -B build/win-vcltl -DASTRAXIS_USE_VC_LTL=ON -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded
cmake --build build/win-vcltl --config Release
ctest --test-dir build/win-vcltl -C Release --output-on-failure
```

`ASTRAXIS_USE_VC_LTL` defaults to `ON` for Windows MSVC and `OFF` elsewhere. Debug configurations always use the standard CRT, including its Debug heap support, even when this option is ON. Set `-DASTRAXIS_USE_VC_LTL=OFF` to also disable VC-LTL for non-Debug configurations. Existing CMake caches retain their previous setting; pass `-DASTRAXIS_USE_VC_LTL=ON` to enable it in an existing build directory. Linux builds do not use VC-LTL. This runtime choice does not extend the application's graphics support to older Windows versions.

**Linux x64** (Vulkan): GCC or Clang with C++20, CMake ≥ 3.28, Ninja, and the X11/Wayland development headers SDL3 needs, e.g. on Ubuntu:

```bash
sudo apt install ninja-build libx11-dev libxext-dev libxrandr-dev libxcursor-dev libxi-dev libxss-dev   libxtst-dev libxfixes-dev libxkbcommon-dev libwayland-dev libdecor-0-dev libdrm-dev libgbm-dev libegl-dev
cmake --preset linux
cmake --build --preset linux-release
```

Run `build/linux/Release/astraxis` (needs a Vulkan driver). Unit tests: `build/linux/Release/astraxis_tests`.

The Windows build also embeds the SPIR-V shaders: set `SDL_GPU_DRIVER=vulkan` to run it on Vulkan.

## Command line

```bash
astraxis --scene sgr_a --event 1
```

`--scene` takes a file name from `assets/scenes/` without `.toml` (default `jupiter`); `--event` jumps to the n-th entry of the scene's Events list.

## Controls

| Input | Action |
|---|---|
| Left/right drag | Rotate camera |
| Mouse wheel | Zoom |
| `1`–`9` | Focus on a body |
| `Space` | Pause / resume |
| `R` | Reverse time |
| `[` / `]` | Halve / double time warp |
| `N` | Jump to the current time |
| `A` | Toggle auto tour (starts by itself after 60 s idle) |
| `O` / `L` | Toggle orbits / labels |
| `H` / `F1` | Toggle UI |
| `Esc` | Quit |

## License

Astraxis is licensed under the [GNU General Public License v3.0 or later](LICENSE) (GPL-3.0-or-later).

Third-party libraries and data keep their own licenses; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
