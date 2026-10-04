<br />
<div align="center">
  <img src="res/astraxis.svg" width="64" height="64" alt="Astraxis icon">
  <h3 align="center">Astraxis</h3>
  <img src="https://img.shields.io/github/license/a1ive/Astraxis" alt="License">
  <img src="https://img.shields.io/github/actions/workflow/status/a1ive/Astraxis/build.yml" alt="Build status">
</div>
<br />

Astraxis is an astronomy visualizer that uses orbital data and physical simulations to show planets, spacecraft and stars in motion. Leave it running in the background when you want something to stare at.

## Scenes

| Scene | What you see | Motion source |
|---|---|---|
| Solar System (default) | The Sun, planets, major moons, dwarf planets and largest asteroids | JPL Horizons ephemerides, JPL mean orbital elements, JPL Small-Body Database |
| Jupiter | Jupiter, its faint rings, the Galilean and four inner moons, and the orbits of Galileo (1995 to 2003) and Juno (2016 onward) | JPL Horizons ephemerides, JPL mean orbital elements |
| JWST | JWST's halo orbit around Sun-Earth L2 | JPL Horizons ephemerides |
| Earth-Moon | Artemis II's free-return flight, Artemis I's distant retrograde orbit and CAPSTONE's near-rectilinear halo orbit, in the Earth-Moon rotating frame | JPL Horizons ephemerides |
| Saturn | Saturn's rings (Cassini radio-occultation optical depths), the seven major moons, Cassini's tour and Huygens' descent to Titan | JPL Horizons ephemerides, JPL mean orbital elements (phases fitted to Horizons) |
| Parker Solar Probe | Seven Venus gravity assists lower the perihelion to 9.86 solar radii. The trajectory forms petals in the Sun-Venus rotating frame | JPL Horizons ephemerides |
| Alpha Centauri | A, B and Proxima | Newtonian N-body integration |
| Sgr A* | The Galactic Centre black hole and the star S2 | Kerr geodesics + GPU ray-traced black hole |
| TRAPPIST-1 | Seven Earth-sized planets in a resonant chain around an ultracool dwarf | Newtonian N-body integration from transit-timing fits |
| Kepler-223 | Four sub-Neptunes in a 3:4:6:8 resonant chain | Newtonian N-body integration from transit-timing fits |
| Kepler-47 | Three planets orbiting an eclipsing binary: a Sun-like star and a red dwarf | Newtonian N-body integration from a photodynamical fit |
| Kepler-64 / PH1 | A Neptune-sized planet orbiting an eclipsing binary, with another pair of stars about 1,300 au away | Newtonian N-body integration for the inner binary and planet; assumed Keplerian orbits for the distant pair |
| TIC 168789840 | Six stars in three eclipsing binaries. Two pairs orbit each other, with the third pair about 250 au away | Nested Keplerian orbits from eclipse observations, with an assumed outer orbit |
| PSR B1620-26 | A pulsar, a white dwarf and a planet in the globular cluster M4. The sky shows M4's stars; the pulsar's radio beams are drawn in blue and slowed down | Nested Keplerian orbits from pulsar timing |

Define scenes in `assets/scenes/*.toml`. You can add a scene by writing a TOML file. The scene files cite sources for physical constants and orbital elements. Each `assets/*/SOURCES.md` lists the data licenses.

## Building

CMake downloads pinned versions of SDL3, glm, Dear ImGui, toml++, stb_image and Microsoft's DirectX Shader Compiler. The shader compiler builds HLSL shaders as DXIL and SPIR-V.

### Windows (Direct3D 12)

Use Visual Studio 2026 with its bundled CMake.

```bash
cmake --preset win-msvc
cmake --build --preset win-msvc-debug
```

Run the app with `build/win-msvc/Debug/astraxis.exe` and the unit tests with `build/win-msvc/Debug/astraxis_tests.exe`.

Windows MSVC builds use [VC-LTL5 v5.3.1](https://github.com/Chuyu-Team/VC-LTL5/tree/v5.3.1) by default for non-Debug configurations. VC-LTL uses the Windows 10 system `ucrtbase.dll` to reduce the runtime footprint. CMake downloads the pinned binary package and checks its SHA-256. To build the CI Release configuration locally:

```bash
cmake --preset win-msvc -B build/win-vcltl -DASTRAXIS_USE_VC_LTL=ON -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded
cmake --build build/win-vcltl --config Release
ctest --test-dir build/win-vcltl -C Release --output-on-failure
```

`ASTRAXIS_USE_VC_LTL` defaults to `ON` for Windows MSVC and `OFF` elsewhere. Debug builds use the standard CRT with Debug heap support, even when the option is `ON`. Set `-DASTRAXIS_USE_VC_LTL=OFF` to disable VC-LTL for non-Debug builds too.

CMake caches keep their saved setting. To enable VC-LTL in an existing build directory, pass `-DASTRAXIS_USE_VC_LTL=ON`. VC-LTL requires Windows and MSVC. The app still requires a Windows version that supports its graphics backend.

### Linux x64 (Vulkan)

Use GCC or Clang with C++20 support, CMake 3.28 or later, Ninja and the X11/Wayland development headers required by SDL3. On Ubuntu:

```bash
sudo apt install ninja-build libx11-dev libxext-dev libxrandr-dev libxcursor-dev libxi-dev libxss-dev libxtst-dev libxfixes-dev libxkbcommon-dev libwayland-dev libdecor-0-dev libdrm-dev libgbm-dev libegl-dev
cmake --preset linux
cmake --build --preset linux-release
```

Run the app with `build/linux/Release/astraxis` and the unit tests with `build/linux/Release/astraxis_tests`. The app needs a Vulkan driver.

The Windows build includes SPIR-V shaders too. Set `SDL_GPU_DRIVER=vulkan` to run it on Vulkan.

## Command line

```bash
astraxis --scene sgr_a --event 1
```

`--scene` selects a file from `assets/scenes/`. Use the file name without `.toml`; the default is `solar_system`. `--event` jumps to an entry in the scene's Events list, numbered from 1.

## Controls

| Input | Action |
|---|---|
| Left/right drag | Rotate camera |
| Mouse wheel | Zoom |
| `1` to `9` | Focus on a body |
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

Third-party libraries and data have their own licenses. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
