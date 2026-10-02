# Astraxis

A quiet, physically grounded astronomy visualizer meant to be left running in the background — something to stare at.

## Scenes

| Scene | What you see | Motion source |
|---|---|---|
| Jupiter | Jupiter and the four Galilean moons | JPL mean orbital elements |
| Voyager | The Sun, planets and both Voyager probes from 1977 on | JPL Horizons ephemerides |
| JWST | JWST's halo orbit around Sun–Earth L2 | JPL Horizons ephemerides |
| Alpha Centauri | A, B and Proxima | Newtonian N-body integration |
| Sgr A* | The Galactic Centre black hole and the star S2 | Kerr geodesics + GPU ray-traced black hole |
| TRAPPIST-1 | Seven Earth-sized planets in a resonant chain around an ultracool dwarf | Newtonian N-body integration from transit-timing fits |
| Kepler-223 | Four sub-Neptunes in a 3:4:6:8 resonant chain | Newtonian N-body integration from transit-timing fits |

Scenes are described in `assets/scenes/*.toml`; adding one needs no code changes. Physical constants and orbital elements cite their sources in the scene files, and data licenses are listed in each `assets/*/SOURCES.md`.

## Building

Requirements: Windows, Visual Studio 2026 (with its bundled CMake), and the Windows SDK (for `dxc`). All other dependencies (SDL3, glm, Dear ImGui, toml++, stb_image) are fetched by CMake at pinned versions.

```bash
cmake --preset win-msvc
cmake --build --preset win-msvc-debug
```

Run `build/win-msvc/Debug/astraxis.exe`. Unit tests: `build/win-msvc/Debug/astraxis_tests.exe`.

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
