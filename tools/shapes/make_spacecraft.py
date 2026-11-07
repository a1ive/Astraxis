"""Convert spacecraft models (NASA's GLB and STL, Gaia Sky's OBJ, this project's Orion) to AXMESH2
meshes with a texture atlas.

Run (the decoder and images need packages; the models are not checked in, the argument is
ref/ with a clone of github.com/nasa/NASA-3D-Resources and the models in probes/):
    python -m venv --system-site-packages <venv>
    <venv>/Scripts/python -m pip install DracoPy      (Apache-2.0; only this tool uses it)
    <venv>/Scripts/python -E -P tools/shapes/make_spacecraft.py ref [name ...]

Needs numpy and Pillow. Writes assets/shapes/<name>.mesh and assets/textures/<name>.<ext>.

The body frame: +z is the axis the scene points (`point_at`: the Sun, or the Earth for
an antenna), x and y are chosen below. One texture holds everything: the model's images
and, for materials with a plain color, a small palette cell (a triangle's u and v are then
constant, so the renderer samples mip 0 and neighbors do not bleed in).

The vertex "albedo" channel holds the sunlight reaching the vertex with the Sun on +z, by
a ray along +z against the whole model: 1 where it is free, 0 where another part is in the
way. This is exact while the Sun stays on +z, which is why only spacecraft that keep one
face to the Sun are baked this way; `bake_sun: False` leaves every vertex at 1.

An STL has no materials: `parts` names them, per connected piece, and `colors` gives each one's.
Parts printed separately (ISEE-3) are put together by `assembly`.
"""
import argparse
import io
import json
import struct
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).resolve().parent))  # -P leaves the script's directory out
from axmesh import write_mesh

ROOT = Path(__file__).resolve().parents[2]
PALETTE_CELL = 8   # px
GUTTER = 4         # px around each image, copied from its edge
RAY_START_KM = 3e-6


CRAFT = {
    # The dish (3.81 m in the model) is 3.7 m (the scene's [VOY]). The model's +y is the way
    # the dish faces.
    'voyager': {
        'file': 'NASA-3D-Resources/3D Models/Voyager Probe (B)/Voyager Probe (B).glb',
        'origin_m': (0.0, 0.65, 0.0),
        'scale': 3.7 / 3.81,
        'body_z': (0, 1, 0),
        'body_x': (-1, 0, 0),   # the magnetometer boom toward the body's +y
        'image_scale': {'tex_02_AO.png': 0.5},
        # The NASA logotype on the bus's cream shield: it is not in the public domain
        # (NASA's media guidelines), so it is painted over with the shield's color.
        'erase_px': {'Voyager_tex_01.ps': [(828, 440, 882, 480)]},
        'atlas_width': 1600,
        'format': 'jpg',
    },
    # The model's sunward is -y (the sunshield is under the telescope), and the telescope
    # looks along its +z.
    'jwst': {
        'file': 'NASA-3D-Resources/3D Models/James Webb Space Telescope (B)/James Webb Space Telescope (B).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (0, -1, 0),
        'body_x': (0, 0, 1),
        'merge_m': 0.02,   # 54,000 of its 99,000 triangles are small brackets
        'atlas_width': 128,
        'format': 'png',
    },
    # The model's sunward is -z (the heat shield), and the magnetometer boom goes off along its +z.
    'parker': {
        'file': 'NASA-3D-Resources/3D Models/Parker Solar Probe/Parker Solar Probe.glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (0, 0, -1),
        'body_x': (1, 0, 0),
        'image_scale': {'PSP_Spec': 0.5},
        'atlas_width': 1300,
        'format': 'jpg',
    },
    # The model's +y is the way the dish faces, and it is to scale (a 4 m dish, 6.7 m high).
    # Huygens leaves in 2004 and is its own body in the scene, so this is the file without it.
    'cassini': {
        'file': 'NASA-3D-Resources/3D Models/Cassini-Huygens (A)/Cassini-Huygens (A) (without Hyugens).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (0, 1, 0),
        'body_x': (1, 0, 0),    # the magnetometer boom toward the body's +y
        'image_scale': {'dish_AO.png': 0.5},
        'atlas_width': 1600,
        'format': 'jpg',
    },
    # The model's +y is the way the dish and the solar panels face, and it is to scale (20 m across).
    'juno': {
        'file': 'NASA-3D-Resources/3D Models/Juno (A)/Juno (A).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (0, 1, 0),
        'body_x': (1, 0, 0),
        'image_scale': {'Juno_Texture_Solar': 0.5, 'black_krinkle': 0.5},
        'atlas_width': 2100,
        'format': 'jpg',
    },
    # The model's unit is not stated; 146.7 units from the low-gain antenna to the probe are the
    # 5.3 m of the scene's [GLF]. Its +z is the way the antenna faces, and the boom goes off along
    # +y. The colors in the file are false (a bright purple named "black"): dark blankets and
    # gold. The furled antenna's mesh is half transparent in the file, and opaque here.
    'galileo': {
        'file': 'NASA-3D-Resources/3D Models/Galileo/Galileo.glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 5.3 / 146.7,
        'body_z': (0, 0, 1),
        'body_x': (1, 0, 0),
        'colors': {'bus:black': (92, 92, 98), 'rtga_top:gold_fo': (176, 128, 58), 'ribs:gold_mesh': (196, 150, 70)},
        'subdivide_m': 0.25,    # the sunshield is one fan from its center
        'atlas_width': 64,
        'format': 'png',
    },
    # To scale (a 2.74 m dish). The dish faces the model's -z, and the boom goes off along its -y.
    'pioneer10': {
        'file': 'NASA-3D-Resources/3D Models/Pioneer 10/Pioneer 10.glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (0, 0, -1),
        'body_x': (1, 0, 0),    # the boom toward the body's +y
        'image_scale': {'dis_AO': 0.5, 'other_AO.png': 0.5, 'RTG_AO_new': 0.5, 'Pioneer_tex_01': 0.25,
                        'foil_gold_ramp': 0.25, 'foil_silver_ramp': 0.25, 'Tex_black.png': 0.125},
        'atlas_width': 1600,
        'format': 'jpg',
    },
    # Telescopes. Chandra, Spitzer and Kepler are to scale as they are. +z below is the side that
    # keeps to the Sun: the solar panels, or the shield they make.
    'chandra': {
        'file': 'NASA-3D-Resources/3D Models/Chandra X-ray Observatory/Chandra X-ray Observatory.glb',
        'origin_m': (3.5, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (-0.64, -0.77, 0.0),   # the panels' normal, in the model
        'body_x': (0, 0, -1),            # the telescope toward the body's +y
        'merge_m': 0.015,
        'atlas_width': 1100,
        'format': 'jpg',
    },
    # An early design, with the door open; 411.9 units are the 12.7 m of the scene's [RMS]. The
    # telescope looks along the model's -x, with the solar panels on its +y side.
    'roman': {
        'file': 'NASA-3D-Resources/3D Models/Nancy Grace Roman Space Telescope (B)/Nancy Grace Roman Space Telescope (B).glb',
        'origin_m': (50.0, 0.0, 0.0),
        'scale': 12.7 / 411.9,
        'body_z': (0, 1, 0),
        'body_x': (0, 0, -1),   # the door toward the body's +y
        'merge_m': 0.04,
        'atlas_width': 512,
        'format': 'png',
    },
    'spitzer': {
        'file': 'NASA-3D-Resources/3D Models/Spitzer Space Telescope/Spitzer Space Telescope.glb',
        'origin_m': (0.0, 0.7, 0.0),
        'scale': 1.0,
        'body_z': (0, 0, -1),   # the solar panels and their shield
        'body_x': (-1, 0, 0),   # the telescope toward the body's +y
        'merge_m': 0.015,
        'atlas_width': 512,
        'format': 'png',
    },
    'kepler': {
        'file': 'NASA-3D-Resources/3D Models/Kepler (A)/Kepler (A).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (-1, 0, 0),   # the solar panels wrap around this side
        'body_x': (0, 0, 1),    # the aperture toward the body's +y
        'image_scale': {'r_AO': 0.5, 'dish_AO': 0.5, 'foil_gold_ramp.pn': 0.25, 'foil_silver_ramp': 0.25},
        'atlas_width': 1600,
        'format': 'jpg',
    },
    # 37.1 units across the panels are the 3.9 m of the scene's Horizons notes. The cameras look along
    # the model's +y, at right angles to the Sun.
    'tess': {
        'file': 'NASA-3D-Resources/3D Models/Transiting Exoplanet Survey Satellite (TESS) (A)/Transiting Exoplanet Survey Satellite (TESS) (A).glb',
        'origin_m': (0.0, -1.0, 0.0),
        'scale': 3.9 / 37.1,
        'body_z': (0, 0, 1),
        'body_x': (1, 0, 0),
        'merge_m': 0.01,
        'atlas_width': 256,
        'format': 'png',
    },
    # 278.3 units across the solar array are the 5.1 m of the scene's [WMF]. The array faces the
    # Sun on the model's -y, and shades the telescopes above it.
    'wmap': {
        'file': 'NASA-3D-Resources/3D Models/Wilkinson Microwave Anisotropy Probe (WMAP)/Wilkinson Microwave Anisotropy Probe (WMAP).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 5.1 / 278.3,
        'body_z': (0, -1, 0),
        'body_x': (1, 0, 0),
        'merge_m': 0.02,   # its 112,000 triangles are mostly brackets
        'atlas_width': 512,
        'format': 'png',
    },
    # Mars. To scale; the solar panels face the model's +y (MAVEN, with its dish) or -y (MRO).
    'maven': {
        'file': 'NASA-3D-Resources/3D Models/Mars Atmosphere and Volatile EvolutioN (MAVEN) (B)/Mars Atmosphere and Volatile EvolutioN (MAVEN) (B).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (0, 1, 0),
        'body_x': (1, 0, 0),
        'image_default_scale': 0.5,
        'image_scale': {'tex_01': 1.0, 'foil_gold_ramp': 0.25, 'foil_silver_ramp': 0.25},
        'merge_m': 0.01,
        'atlas_width': 2100,
        'format': 'jpg',
    },
    'mro': {
        'file': 'NASA-3D-Resources/3D Models/Mars Reconnaissance Orbiter (MRO) (C)/Mars Reconnaissance Orbiter (MRO) (C).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (0, -1, 0),
        'body_x': (1, 0, 0),
        'image_scale': {'foil_gold_ramp.pn': 0.25, 'foil_silver_ramp': 0.25},
        'atlas_width': 1600,
        'format': 'jpg',
    },
    # tex_01 has the HEND instrument's emblem and label, and tex_boom the Voyager texture with its
    # NASA logotype: painted over.
    'odyssey': {
        'file': 'NASA-3D-Resources/3D Models/Mars Odyssey/Mars Odyssey.glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (0.0, -0.19, 0.98),   # the solar panels' normal, in the model
        'body_x': (1, 0, 0),
        'image_default_scale': 0.5,
        'image_scale': {'tex_02': 1.0, 'foil_silver_ramp.': 0.25, 'foil_gold': 0.25},
        'erase_px': {'tex_01': [(8, 660, 250, 980)], 'tex_boom': [(828, 440, 882, 480)]},
        'atlas_width': 2100,
        'format': 'jpg',
    },
    # The sunshade faces the model's +x, and the magnetometer boom goes off along its -x.
    'messenger': {
        'file': 'NASA-3D-Resources/3D Models/MErcury Surface, Space ENvironment, GEochemistry, and Ranging (MESSENGER)/MErcury Surface, Space ENvironment, GEochemistry, and Ranging (MESSENGER).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (1, 0, 0),
        'body_x': (0, 0, -1),
        'image_default_scale': 0.5,
        'image_scale': {'Messenger_tex_01.': 1.0, 'foil_gold_ramp.pn': 0.25, 'foil_silver_ramp.': 0.25, 'Tex_black.png': 0.125, 'Tex_dark grey.png': 0.125},
        'atlas_width': 2100,
        'format': 'jpg',
    },
    # The body is 28.3 units across, and 1.85 m ([LDF] in the scene). Solar panels all around.
    'ladee': {
        'file': 'NASA-3D-Resources/3D Models/Lunar Atmosphere and Dust Environment Explorer (LADEE)/Lunar Atmosphere and Dust Environment Explorer (LADEE).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.85 / 28.3,
        'body_z': (1, 0, 0),
        'body_x': (0, 0, -1),
        'merge_m': 0.02,   # small parts: 49,500 triangles for a 2 m spacecraft
        'atlas_width': 512,
        'format': 'png',
    },
    # ARTEMIS is two of the THEMIS spacecraft. To scale; solar panels all around.
    'themis': {
        'file': 'NASA-3D-Resources/3D Models/Time History of Events and Macroscale Interactions During Substorms (THEMIS)/Time History of Events and Macroscale Interactions During Substorms (THEMIS).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (0, 0, -1),
        'body_x': (-1, 0, 0),
        'image_scale': {'foil_gold_ramp.pn': 0.25, 'foil_silver_ramp.': 0.25},
        'atlas_width': 1100,
        'format': 'jpg',
    },
    # To scale (32 m across); the solar panels face the model's +z.
    'rosetta': {
        'file': 'NASA-3D-Resources/3D Models/Rosetta/Rosetta.glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 1.0,
        'body_z': (0, 0, 1),
        'body_x': (1, 0, 0),
        'image_default_scale': 0.5,
        'image_scale': {'tex_01': 1.0, 'foil_gold_ramp.pn': 0.25, 'foil_silver_ramp': 0.25},
        'atlas_width': 1700,
        'format': 'jpg',
    },
    # The body is 4.52 wide in the model, and 4.3 m ([HBN] in the scene), which makes the tube the
    # 13 m it should be. The telescope looks along the model's +y (the door is open). The solar
    # cells, the dark side of the panels, face its -z ([HBP]: the light side, away from the Sun, is
    # the white radiator).
    # The texture has a banner with the NASA and ESA logotypes, neither in the public domain,
    # used only by the NASA_logo decal: the decal is left out and the banner painted over.
    'hubble': {
        'file': 'probes/Hubble-1.glb',
        'origin_m': (0.0, 0.6, 0.0),
        'scale': 4.3 / 4.52,
        'body_z': (0, 0, -1),
        'body_x': (-1, 0, 0),   # the door toward the body's +y
        'skip_materials': ['NASA_logo'],
        'erase_polygon': {'hubble_tex_01.jpg': [(95, 105), (125, 90), (200, 55), (340, 25), (420, 12), (520, 8),
                                                (690, 14), (850, 40), (1000, 95), (1012, 362), (930, 325), (770, 280),
                                                (690, 268), (520, 240), (430, 265), (250, 310), (70, 382)]},
        'image_scale': {'hubble_AO_1.jpg': 0.5, 'hubble_AO_light.jpg': 0.5},
        'atlas_width': 2100,
        'format': 'jpg',
    },
    # The dish (1.888 m in the model) is 2.1 m (the scene's [NHS]). The model's +y is the way the
    # dish faces, and the generator goes off along its +z.
    'new_horizons': {
        'file': 'probes/New_Horizons.glb',
        'origin_m': (0.0, 0.2, 0.0),
        'scale': 2.1 / 1.888,
        'body_z': (0, 1, 0),
        'body_x': (-1, 0, 0),   # the generator toward the body's +y
        # The foil is a metal in the file, near pure yellow where it is lit as a matte surface here.
        'colors': {'foil_gold': (186, 140, 58)},
        'image_scale': {'newhorizons_AO_dish.png': 0.5, 'newhorizons_tex_02.png': 0.5, 'newhorizons_AO_1.png': 0.25},
        'atlas_width': 1600,
        'format': 'jpg',
    },
    # To scale (17.5 m across the panels, the scene's Horizons notes). The solar cells face the
    # model's +z and the dish its -y; Schiaparelli's mounting is on its +y. The metals are near pure
    # white or yellow lit as matte surfaces here, and the blankets (a texture of black with a few
    # gold specks) a silhouette: replaced.
    'tgo': {
        'file': 'probes/TGO.glb',
        'origin_m': (0.0, -0.4, 0.0),
        'scale': 1.0,
        'body_z': (0, 0, 1),
        'body_x': (1, 0, 0),
        'colors': {'black_krinkle_dots': (118, 118, 122), 'black krinkle': (118, 118, 122),
                   'foil gold': (186, 140, 58), 'array gold trim': (186, 140, 58),
                   'foil silver': (190, 190, 195), 'foil silver antenna': (205, 205, 208),
                   'mirror chrome': (125, 130, 140), 'aluminum': (170, 170, 175)},
        'atlas_width': 1100,
        'format': 'jpg',
    },
    # To scale (30.5 m across the arrays, the scene's Horizons notes). The solar cells face the
    # model's +z; the dish faces its -y, the magnetometer boom goes off toward -y and -z, and the
    # radar antennas lie along z.
    'clipper': {
        'file': 'probes/clipper_spacecraft.glb',
        'origin_m': (0.0, 0.0, 1.76),
        'scale': 1.0,
        'body_z': (0, 0, 1),
        'body_x': (1, 0, 0),
        'merge_m': 0.03,
        'merge_uv': 128,   # the textures are mostly flat color, so nearby texels can merge
        # Each 0.56 x 0.79 m rectangle of the arrays repeats one solar cell 7 x 18 times (19 rows,
        # since the rows do not start on a whole tile), and an atlas cannot repeat. Columns, rows,
        # pixels per cell, and a gain: the cell is near black, and would not show against the sky.
        'tile_images': {'clipper_solar_cells_albedo': (7, 19, 24, 2.6)},
        'atlas_width': 2100,
        'format': 'jpg',
    },
    # An STL with no colors, and 40.3 units across are the 19 m of the scene; the bus comes out 1.2 m
    # wide (APL: 1.2 x 1.3 x 1.3 m). The two wings are the same shape turned 180 deg about the
    # telescope's axis, so which face has the cells is a choice: both on the model's +y here. DRACO
    # looks along +z, and the ion engine along -z. The colors are illustrative.
    'dart': {
        'file': 'probes/Double Asteroid Redirection Test (DART).stl',
        'origin_m': (0.0, 1.77, -0.49),
        'scale': 19.0 / 40.34,
        'body_z': (0, 1, 0),
        'body_x': (-1, 0, 0),   # the telescope toward the body's +y
        # Material, then the box a piece has to fit in and the size it has to reach (model units).
        'parts': [('cells', (2.6, 3.1, -2.2), (20.0, 3.19, 1.2), 2.0),
                  ('cells', (-20.0, 0.07, -2.2), (-2.6, 0.15, 1.2), 2.0),
                  ('boom', (1.85, -1.0, -3.0), (20.2, 4.0, 2.0), 0.0),
                  ('boom', (-20.2, -1.0, -3.0), (-1.85, 4.0, 2.0), 0.0),
                  ('radiator', (-1.3, 0.4, -1.8), (-1.2, 3.1, 0.8), 1.0),
                  ('radiator', (1.2, 0.4, -1.8), (1.3, 3.1, 0.8), 1.0),
                  ('blanket', (-1.3, 0.3, -2.0), (1.3, 3.2, 1.01), 2.0),
                  ('dish', (-2.1, 2.8, -1.3), (-1.0, 4.0, 0.0), 0.5),
                  ('camera', (-0.6, 1.2, 1.6), (0.6, 2.35, 2.42), 0.25),
                  ('engine', (-0.73, 1.04, -2.51), (0.73, 2.49, -0.7), 0.3)],
        'fallback_part': 'fitting',
        'back_faces': {'cells': 'array_back'},
        'colors': {'cells': (28, 32, 62), 'array_back': (170, 170, 165), 'boom': (150, 150, 155),
                   'radiator': (205, 205, 200), 'blanket': (128, 128, 134), 'dish': (210, 210, 205),
                   'camera': (190, 190, 195), 'engine': (70, 70, 74), 'fitting': (165, 165, 170)},
        'merge_m': 0.03,
        'atlas_width': 64,
        'format': 'png',
    },
    # The entry assembly, heat shield and back cover. 2.615 m across in the model, and 2.7 m (the
    # scene's). The shield's nose is on the model's +x. The scene points it at Titan, not the Sun, so
    # nothing is baked.
    'huygens': {
        'file': 'probes/Cassini-Huygens (A) (without Cassini).glb',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 2.7 / 2.615,
        'body_z': (1, 0, 0),
        'body_x': (0, 1, 0),
        'image_scale': {'tex_01_': 0.5},
        'bake_sun': False,
        'atlas_width': 530,
        'format': 'jpg',
    },
    # Gaia Sky's models (OBJ with MTL). Gaia: the sunshield is 1 unit across in the model and 10.2 m
    # deployed (ESA); it faces the model's -y, under the payload. Many textures repeat several
    # times across a triangle, which an atlas cannot do: those are plain colors here (each its
    # texture's average), and the solar cells, the one pattern worth keeping, are laid out in a grid
    # first (the panels span up to 7 tiles).
    'gaia': {
        'file': 'probes/mission-gaia/gaia.obj',
        'origin_m': (0.0, 0.0, 0.0),
        'scale': 10.2,
        'body_z': (0, -1, 0),
        'body_x': (1, 0, 0),
        'colors': {'GaiaMetal': (76, 76, 78), 'DarkMetal': (30, 30, 32), 'DarkM_darkmetal': (30, 30, 32),
                   'BelowTTS': (60, 60, 62), 'Bottom': (40, 40, 42), 'Circuit1': (40, 40, 42),
                   'Circuit2': (62, 62, 62), 'DarkM': (48, 48, 48), 'CCD': (40, 44, 90),
                   'Gold': (186, 140, 58), 'AluminumGold': (186, 140, 58), 'GoldFoil': (186, 140, 58)},
        'tile_images': {'solar-panel-small.jpg': (8, 8, 32, 1.0)},
        'atlas_width': 800,
        'format': 'jpg',
    },
    # ESA's Euclid as Gaia Sky has it. 0.869 units from the bottom of the service module to the top
    # of the sunshield, without the dish and thrusters below, are ESA's 4.7 m: it is then 3.5 m
    # across (ESA: 3.7 m). The sunshield faces the model's +z, and the telescope looks along +y.
    'euclid': {
        'file': 'probes/spacecraft-euclid/euclid.obj',
        'origin_m': (-0.015, 0.066, 0.0),
        'scale': 4.7 / 0.869,
        'body_z': (0, 0, 1),
        'body_x': (1, 0, 0),    # the telescope toward the body's +y
        'image_scale': {'Euclid2_Base_Color.png': 0.25},
        'merge_m': 0.02,
        'merge_uv': 256,
        'atlas_width': 1040,
        'format': 'jpg',
    },
    # Orion with the European Service Module and its four arrays, as built for this project in
    # Blender (ref/artemis_orion/build_orion.py, from NASA's and Airbus's published dimensions; no
    # third-party mesh). Metres, +z from the engine toward the crew module's docking ring; the
    # arrays lie flat across +z with their cells on the +z side, so the scenes keep +z to the Sun
    # and the Sun is baked. The origin is halfway up the 7.4 m stack.
    'orion': {
        'file': 'artemis_orion/artemis_orion.obj',
        'origin_m': (0.0, 0.0, 3.7),
        'scale': 1.0,
        'body_z': (0, 0, 1),
        'body_x': (1, 0, 0),
        'merge_m': 0.02,
        'atlas_width': 64,
        'format': 'png',
    },
    # ISEE-3 from NASA's printable parts (inches): the 16-sided drum, the booms and wires each in
    # its own file, put together here after the JHU/APL flight configuration (eoPortal, fig. 3) and
    # the NASA program document (spacecollege.org/?p=245): the drum's y is the spin axis; four rigid
    # booms and four wire antennas leave from the band at mid-height (y 27.9 .. 30.2), the
    # experiment booms (search coil and magnetometer, 302 cm) where the drum has hinge brackets
    # (model +-z), the inertia booms (114 cm) at right angles, the wires (49 m) between them, and
    # the axial antennas (701 cm) along the spin axis. Thruster.stl, a loose part that fits none of
    # the drum's pods, is left out. The solar array covers the side except the
    # band and two rings. The spin axis is kept perpendicular to the ecliptic, pointing north: the
    # body's +y, with +z (a side) to the Sun. It spins every 3 s, so nothing is baked; the colors
    # are illustrative.
    'isee3': {
        'file': 'probes/isee3',
        'assembly': [
            {'file': 'Core 2.stl', 'material': 'structure',
             'bands': [(lo, hi, 33.0, 34.9, 'cells') for lo, hi in ((-0.01, 15.14), (19.85, 42.21), (47.3, 62.07))]
                      + [(27.87, 30.25, 33.0, 34.9, 'structure')],
             'caps': [(0.0, 34.79, 16, 'structure'), (62.06, 34.79, 16, 'structure')]},   # printed open
            {'file': 'Exboom1.stl', 'material': 'boom', 'at': (7.15, 29.06, 34.5), 'toward': (0, 0, 1)},
            {'file': 'Exboom2.stl', 'material': 'boom', 'at': (-7.15, 29.06, -34.5), 'toward': (0, 0, -1)},
            {'file': 'Inertia boom.stl', 'material': 'boom', 'at': (34.6, 29.06, 0.0), 'toward': (1, 0, 0)},
            {'file': 'Inertia boom.stl', 'material': 'boom', 'at': (-34.6, 29.06, 0.0), 'toward': (-1, 0, 0)},
            *({'file': 'Rad antenna.stl', 'material': 'wire', 'at': (24.4 * sx, 29.06, 24.4 * sz),
               'toward': (sx, 0, sz), 'length': 4900.0 / 2.54} for sx, sz in ((1, 1), (1, -1), (-1, 1), (-1, -1))),
            {'file': 'Rad antenna.stl', 'material': 'wire', 'at': (0.0, 62.06, 0.0), 'toward': (0, 1, 0),
             'up': (0, 0, 1), 'length': 701.0 / 2.54},
            {'file': 'Rad antenna.stl', 'material': 'wire', 'at': (0.0, -1.55, 0.0), 'toward': (0, -1, 0),
             'up': (0, 0, 1), 'length': 701.0 / 2.54},
        ],
        'origin_m': (0.0, 31.03, 0.0),
        'scale': 0.0254,
        'body_z': (1, 0, 0),
        'body_x': (0, 0, -1),   # the spin axis (model +y) is the body's +y
        'colors': {'cells': (38, 42, 78), 'structure': (178, 178, 172), 'boom': (205, 205, 200),
                   'wire': (160, 150, 120)},
        'merge_m': 0.01,
        'merge_skip': ['wire'],
        'bake_sun': False,
        'atlas_width': 64,
        'format': 'png',
    },
}


def linear_to_srgb(c):
    c = np.clip(np.asarray(c, dtype=np.float64), 0.0, 1.0)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1.0 / 2.4) - 0.055)


def node_matrix(node):
    if 'matrix' in node:
        return np.array(node['matrix'], dtype=np.float64).reshape(4, 4).T
    x, y, z, w = node.get('rotation', (0, 0, 0, 1))
    r = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                  [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                  [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    m = np.eye(4)
    m[:3, :3] = r * np.asarray(node.get('scale', (1, 1, 1)), dtype=np.float64)
    m[:3, 3] = node.get('translation', (0, 0, 0))
    return m


class Glb:
    plain = ()

    def __init__(self, path):
        with open(path, 'rb') as f:
            if f.read(4) != b'glTF':
                raise SystemExit(f'{path}: not a GLB')
            f.read(8)
            length, _ = struct.unpack('<II', f.read(8))
            self.json = json.loads(f.read(length))
            length, _ = struct.unpack('<II', f.read(8))
            self.bin = f.read(length)

    def view(self, index):
        v = self.json['bufferViews'][index]
        start = v.get('byteOffset', 0)
        return self.bin[start:start + v['byteLength']]

    def accessor(self, index):
        """An uncompressed attribute or index list, as doubles, one row per element."""
        info = self.json['accessors'][index]
        view = self.json['bufferViews'][info['bufferView']]
        kind = np.dtype({5126: '<f4', 5125: '<u4', 5123: '<u2', 5121: '<u1'}[info['componentType']])
        size = {'SCALAR': 1, 'VEC2': 2, 'VEC3': 3}[info['type']]
        if view.get('byteStride', 0) not in (0, kind.itemsize * size):
            raise SystemExit('interleaved buffers are not handled')
        start = view.get('byteOffset', 0) + info.get('byteOffset', 0)
        return np.frombuffer(self.bin, kind, info['count'] * size, start).reshape(-1, size).astype(np.float64)

    def decode(self, prim):
        """(positions, normals, uv or None, triangles) of a primitive, Draco-compressed or not."""
        if 'KHR_draco_mesh_compression' in prim.get('extensions', {}):
            import DracoPy   # only the files from the repository need it
            draco = prim['extensions']['KHR_draco_mesh_compression']
            mesh = DracoPy.decode(self.view(draco['bufferView']))
            uv = None if mesh.tex_coord is None or len(mesh.tex_coord) == 0 else \
                np.asarray(mesh.tex_coord, dtype=np.float64).reshape(-1, 2)
            return (np.asarray(mesh.points, dtype=np.float64).reshape(-1, 3),
                    np.asarray(mesh.normals, dtype=np.float64).reshape(-1, 3), uv,
                    np.asarray(mesh.faces, dtype=np.int64).reshape(-1, 3))
        attributes = prim['attributes']
        uv = self.accessor(attributes['TEXCOORD_0']) if 'TEXCOORD_0' in attributes else None
        return (self.accessor(attributes['POSITION']), self.accessor(attributes['NORMAL']), uv,
                self.accessor(prim['indices']).astype(np.int64).reshape(-1, 3))

    def primitives(self):
        """(material, world positions, world normals, uv or None, triangles) of each primitive."""
        out = []

        def walk(index, parent):
            node = self.json['nodes'][index]
            matrix = parent @ node_matrix(node)
            # '_root' is a unit cube Blender's exporter left as a placeholder.
            if 'mesh' in node and node.get('name') != '_root':
                for prim in self.json['meshes'][node['mesh']]['primitives']:
                    if self.invisible(prim.get('material')):
                        continue
                    points, normals, uv, faces = self.decode(prim)
                    linear = matrix[:3, :3]
                    world = points @ linear.T + matrix[:3, 3]
                    turned = normals @ np.linalg.inv(linear)
                    turned /= np.maximum(np.linalg.norm(turned, axis=1, keepdims=True), 1e-12)
                    out.append((prim.get('material'), world, turned, uv, faces))
            for child in node.get('children', []):
                walk(child, matrix)

        for root in self.json['scenes'][self.json.get('scene', 0)]['nodes']:
            walk(root, np.eye(4))
        return out

    def invisible(self, material):
        """Fully transparent: Juno has copies of faces in such materials, which would flicker here."""
        if material is None:
            return False
        info = self.json['materials'][material]
        factor = info.get('pbrMetallicRoughness', {}).get('baseColorFactor', (1, 1, 1, 1))
        return info.get('alphaMode') in ('MASK', 'BLEND') and factor[3] == 0.0

    def material_image(self, material):
        """Index of the JSON image of a material's base color (the PNG, not the WebP copy), or None:
        also None for a material the spec gives a plain color instead (`plain`)."""
        if self.json['materials'][material].get('name') in self.plain:
            return None
        info = self.json['materials'][material].get('pbrMetallicRoughness', {}).get('baseColorTexture')
        return None if info is None else self.json['textures'][info['index']]['source']

    def material_color(self, material, overrides):
        """sRGB bytes of a plain material; glTF stores the factor in linear light."""
        if material is None:
            return (200, 200, 200)
        if self.json['materials'][material].get('name') in overrides:
            return overrides[self.json['materials'][material]['name']]
        factor = self.json['materials'][material].get('pbrMetallicRoughness', {}).get('baseColorFactor', (1, 1, 1, 1))
        return tuple(int(round(255.0 * c)) for c in linear_to_srgb(factor[:3]))


class Stl:
    """A binary STL, with the interface `build` uses of Glb: one plain material per name in the
    spec's `parts`. The file keeps three vertices per triangle: they are welded to find the connected
    pieces, and the normals are smoothed across edges under 40 deg, as a shaded CAD model looks."""
    CREASE_COS = np.cos(np.radians(40.0))

    RECORD = np.dtype([('normal', '<f4', 3), ('corners', '<f4', (3, 3)), ('attribute', '<u2')])

    def __init__(self, path, spec):
        if 'assembly' in spec:
            corners, material_names = self.assemble(path, spec)
        else:
            corners = np.fromfile(path, dtype=self.RECORD, offset=84)['corners'].astype(np.float64)
        corners = corners.reshape(-1, 3)
        _, first, inverse = np.unique(np.round(corners * 1e4).astype(np.int64), axis=0, return_index=True, return_inverse=True)
        self.positions = corners[first]
        self.triangles = inverse.reshape(-1, 3)
        a, b, c = (self.positions[self.triangles[:, k]] for k in range(3))
        self.faces_normal = np.cross(b - a, c - a)
        self.corner_normals = self.smooth_normals()
        unit = self.faces_normal / np.maximum(np.linalg.norm(self.faces_normal, axis=1, keepdims=True), 1e-30)

        if 'assembly' in spec:
            self.names = list(dict.fromkeys(material_names))
            self.json = {'materials': [{'name': name} for name in self.names]}
            self.material = np.array([self.names.index(name) for name in material_names])
            return
        names = [part[0] for part in spec['parts']] + [spec['fallback_part']] + list(spec.get('back_faces', {}).values())
        self.names = list(dict.fromkeys(names))
        self.json = {'materials': [{'name': name} for name in self.names]}
        pieces = self.connected_pieces()
        low = np.full((pieces.max() + 1, 3), np.inf)
        high = np.full((pieces.max() + 1, 3), -np.inf)
        np.minimum.at(low, pieces, self.positions)
        np.maximum.at(high, pieces, self.positions)
        piece_material = np.full(len(low), self.names.index(spec['fallback_part']))
        for name, box_low, box_high, size in reversed(spec['parts']):
            fits = np.all(low >= box_low, axis=1) & np.all(high <= box_high, axis=1) & ((high - low).max(axis=1) >= size)
            piece_material[fits] = self.names.index(name)
        self.material = piece_material[pieces[self.triangles[:, 0]]]
        sunward = np.array(spec['body_z'], dtype=np.float64)
        sunward /= np.linalg.norm(sunward)
        for name, back in spec.get('back_faces', {}).items():
            away = (self.material == self.names.index(name)) & (unit @ sunward <= 0.5)
            self.material[away] = self.names.index(back)

    @classmethod
    def assemble(cls, folder, spec):
        """Parts printed separately (each in its own file and frame), put together: per entry of
        `assembly`, the part's +z (along a boom) turned to `toward` with its +y as near `up` as it
        goes, scaled along +z to `length` (model units) and its origin moved to `at`; without `at`
        the part is already in place.
        `bands` gives the side faces of a part with all three corners in a y range and a range of
        radius from the y axis another material (a drum's solar array), in order, later ones winning;
        `caps` closes an open drum with flat polygons (y, circumradius, sides, material). One material name per triangle."""
        corners, names = [], []
        for entry in spec['assembly']:
            raw = np.fromfile(folder / entry['file'], dtype=cls.RECORD, offset=84)['corners'].astype(np.float64)
            if 'at' in entry:
                if 'length' in entry:
                    raw[:, :, 2] *= entry['length'] / raw[:, :, 2].max()
                z = np.array(entry['toward'], dtype=np.float64)
                z /= np.linalg.norm(z)
                up = np.array(entry.get('up', (0, 1, 0)), dtype=np.float64)
                y = up - (up @ z) * z
                y /= np.linalg.norm(y)
                x = np.cross(y, z)
                raw = raw @ np.array([x, y, z]) + np.array(entry['at'], dtype=np.float64)
            material = np.full(len(raw), entry['material'], dtype=object)
            radius = np.hypot(raw[:, :, 0], raw[:, :, 2])
            normal = np.cross(raw[:, 1] - raw[:, 0], raw[:, 2] - raw[:, 0])
            side = np.abs(normal[:, 1]) < 0.5 * np.linalg.norm(normal, axis=1)
            for low, high, min_radius, max_radius, name in entry.get('bands', []):
                inside = (raw[:, :, 1] >= low) & (raw[:, :, 1] <= high) & (radius >= min_radius) & (radius <= max_radius)
                material[inside.all(axis=1) & side] = name
            corners.append(raw)
            names.extend(material.tolist())
            for y, radius, sides, name in entry.get('caps', []):
                angle = np.radians(360.0 / sides) * np.arange(sides + 1)
                rim = np.stack([radius * np.cos(angle), np.full(sides + 1, y), radius * np.sin(angle)], axis=1)
                fan = np.stack([np.tile([0.0, y, 0.0], (sides, 1)), rim[:-1], rim[1:]], axis=1)
                corners.append(fan)
                names.extend([name] * sides)
        return np.concatenate(corners), names

    def connected_pieces(self):
        parent = list(range(len(self.positions)))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        for a, b, c in self.triangles.tolist():
            ra = find(a)
            parent[find(b)] = ra
            parent[find(c)] = ra
        return np.unique([find(a) for a in range(len(parent))], return_inverse=True)[1]

    def smooth_normals(self):
        """Per corner of each triangle: the area-weighted normal of the faces at that vertex within
        the crease angle of this one."""
        count = len(self.triangles)
        vertex = self.triangles.reshape(-1)
        order = np.argsort(vertex, kind='stable')
        edges = np.flatnonzero(np.diff(vertex[order])) + 1
        unit = self.faces_normal / np.maximum(np.linalg.norm(self.faces_normal, axis=1, keepdims=True), 1e-30)
        out = np.zeros((3 * count, 3))
        for group in np.split(order, edges):
            faces = group // 3
            similar = unit[faces] @ unit[faces].T > self.CREASE_COS
            out[group] = similar @ self.faces_normal[faces]
        out /= np.maximum(np.linalg.norm(out, axis=1, keepdims=True), 1e-30)
        return out.reshape(-1, 3, 3)

    def primitives(self):
        out = []
        for index in range(len(self.names)):
            chosen = np.flatnonzero(self.material == index)
            if len(chosen) == 0:
                continue
            # Three vertices per triangle: the weld in `build` joins those with the same normal.
            points = self.positions[self.triangles[chosen]].reshape(-1, 3)
            normals = self.corner_normals[chosen].reshape(-1, 3)
            out.append((index, points, normals, None, np.arange(len(points)).reshape(-1, 3)))
        return out

    def material_image(self, material):
        return None

    def material_color(self, material, overrides):
        return overrides[self.names[material]]


class Obj:
    """A Wavefront OBJ with its MTL, with the interface `build` uses of Glb: one primitive per
    material, the images read from the files `map_Kd` names. OBJ texture coordinates start at the
    bottom of the image (glTF's at the top), so v is flipped. Kd is taken as linear, as Blender's
    exporter writes it. Faces with more than three corners are fanned."""
    plain = ()

    def __init__(self, path):
        self.folder = path.parent
        positions, uvs, normals = [], [], []
        faces = {}   # material -> list of corner triples (v, vt, vn)
        material = None
        library = None
        with open(path, encoding='utf-8', errors='replace') as f:
            for line in f:
                head, _, rest = line.strip().partition(' ')
                if head == 'v':
                    positions.append([float(x) for x in rest.split()[:3]])
                elif head == 'vt':
                    uvs.append([float(x) for x in rest.split()[:2]])
                elif head == 'vn':
                    normals.append([float(x) for x in rest.split()[:3]])
                elif head == 'usemtl':
                    material = rest.strip()
                elif head == 'mtllib':
                    library = rest.strip()
                elif head == 'f':
                    corners = []
                    for token in rest.split():
                        parts = (token.split('/') + ['', ''])[:3]
                        corners.append(tuple(int(p) if p else 0 for p in parts))
                    target = faces.setdefault(material, [])
                    for k in range(1, len(corners) - 1):
                        target.append((corners[0], corners[k], corners[k + 1]))
        self.positions = np.array(positions, dtype=np.float64)
        self.uvs = np.array(uvs, dtype=np.float64).reshape(-1, 2)
        self.normals = np.array(normals, dtype=np.float64).reshape(-1, 3)
        self.faces = {name: np.array(corners, dtype=np.int64) for name, corners in faces.items()}

        colors, images = {}, {}
        current = None
        if library is not None:
            with open(self.folder / library, encoding='utf-8', errors='replace') as f:
                for line in f:
                    head, _, rest = line.strip().partition(' ')
                    if head == 'newmtl':
                        current = rest.strip()
                    elif head == 'Kd' and current is not None:
                        colors[current] = [float(x) for x in rest.split()[:3]]
                    elif head == 'map_Kd' and current is not None:
                        images[current] = rest.split()[-1]   # options come before the file name
        self.names = list(self.faces)
        files = sorted(set(images[name] for name in self.names if name in images))
        self.json = {'materials': [{'name': name} for name in self.names],
                     'images': [{'name': file, 'bufferView': k} for k, file in enumerate(files)]}
        self.image_of = {name: files.index(images[name]) for name in self.names if name in images}
        self.colors = colors

    def view(self, index):
        return (self.folder / self.json['images'][index]['name']).read_bytes()

    @staticmethod
    def resolve(index, count):
        """OBJ indices start at 1; negative ones count back from the end."""
        return np.where(index > 0, index - 1, count + index)

    def primitives(self):
        out = []
        for material, name in enumerate(self.names):
            corners = self.faces[name].reshape(-1, 3)   # rows of (v, vt, vn)
            points = self.positions[self.resolve(corners[:, 0], len(self.positions))]
            uv = None
            if name in self.image_of and corners[:, 1].all():
                uv = self.uvs[self.resolve(corners[:, 1], len(self.uvs))] * np.array([1.0, -1.0]) + np.array([0.0, 1.0])
            if corners[:, 2].all():
                normals = self.normals[self.resolve(corners[:, 2], len(self.normals))]
            else:
                a, b, c = points[0::3], points[1::3], points[2::3]
                normals = np.repeat(np.cross(b - a, c - a), 3, axis=0)
            normals = normals / np.maximum(np.linalg.norm(normals, axis=1, keepdims=True), 1e-30)
            out.append((material, points, normals, uv, np.arange(len(points)).reshape(-1, 3)))
        return out

    def material_image(self, material):
        name = self.names[material]
        return None if name in self.plain else self.image_of.get(name)

    def material_color(self, material, overrides):
        name = self.names[material]
        if name in overrides:
            return overrides[name]
        return tuple(int(round(255.0 * c)) for c in linear_to_srgb(self.colors.get(name, (0.8, 0.8, 0.8))))


class Atlas:
    """Images on shelves, then the palette cells."""

    def __init__(self, width):
        self.width = width
        self.images = []  # (key, RGB array)
        self.colors = []
        self.placed = {}
        self.cells = {}
        self.pixels = None

    def add_image(self, key, array):
        self.images.append((key, array))

    def add_color(self, color):
        if color not in self.colors:
            self.colors.append(color)

    def build(self):
        x, y, shelf = 0, 0, 0
        rects = []
        for key, array in sorted(self.images, key=lambda item: -item[1].shape[0]):
            h, w = array.shape[:2]
            if x + w + 2 * GUTTER > self.width:
                x, y = 0, y + shelf
                shelf = 0
            rects.append((key, array, x + GUTTER, y + GUTTER))
            x += w + 2 * GUTTER
            shelf = max(shelf, h + 2 * GUTTER)
        y += shelf
        columns = max(1, self.width // PALETTE_CELL)
        rows = (len(self.colors) + columns - 1) // columns
        height = y + rows * PALETTE_CELL
        pixels = np.zeros((height, self.width, 3), dtype=np.uint8)
        for key, array, px, py in rects:
            h, w = array.shape[:2]
            ys = np.clip(np.arange(-GUTTER, h + GUTTER), 0, h - 1)
            xs = np.clip(np.arange(-GUTTER, w + GUTTER), 0, w - 1)
            pixels[py - GUTTER:py + h + GUTTER, px - GUTTER:px + w + GUTTER] = array[np.ix_(ys, xs)]
            self.placed[key] = (px, py, w, h)
        for k, color in enumerate(self.colors):
            cx, cy = (k % columns) * PALETTE_CELL, y + (k // columns) * PALETTE_CELL
            pixels[cy:cy + PALETTE_CELL, cx:cx + PALETTE_CELL] = color
            self.cells[color] = ((cx + PALETTE_CELL / 2) / self.width, (cy + PALETTE_CELL / 2) / height)
        self.pixels = pixels

    def image_uv(self, key, uv):
        px, py, w, h = self.placed[key]
        H, W = self.pixels.shape[:2]
        inset = 0.5 / np.array([w, h])
        uv = np.clip(uv, inset, 1.0 - inset)
        return (uv * np.array([w, h]) + np.array([px, py])) / np.array([W, H])


def subdivide(positions, uv, triangles, limit_km):
    """Bisect the longest edge of each triangle until no edge is over the limit. An edge's
    midpoint is shared by the triangles on it, so that no cracks open. The sunlight is baked on
    the vertices: a big low-poly fan with one shaded vertex at its center would go dark all over."""
    positions, uv = positions.tolist(), uv.tolist()
    midpoints = {}
    out, stack = [], [tuple(t) for t in triangles.tolist()]
    while stack:
        a, b, c = stack.pop()
        edges = [(a, b, c), (b, c, a), (c, a, b)]   # edge p-q and the opposite corner r
        p, q, r = max(edges, key=lambda e: sum((x - y) ** 2 for x, y in zip(positions[e[0]], positions[e[1]])))
        length = sum((x - y) ** 2 for x, y in zip(positions[p], positions[q])) ** 0.5
        if length <= limit_km:
            out.append((a, b, c))
            continue
        key = (min(p, q), max(p, q))
        if key not in midpoints:
            midpoints[key] = len(positions)
            positions.append([(x + y) / 2.0 for x, y in zip(positions[p], positions[q])])
            uv.append([(x + y) / 2.0 for x, y in zip(uv[p], uv[q])])
        m = midpoints[key]
        stack.append((p, m, r))
        stack.append((m, q, r))
    return np.array(positions), np.array(uv), np.array(out, dtype=np.int64)


def vertex_normals(positions, triangles):
    """As the loader computes them: the area-weighted average of the faces sharing a vertex."""
    a, b, c = (positions[triangles[:, k]] for k in range(3))
    faces = np.cross(b - a, c - a)
    total = np.zeros_like(positions)
    for k in range(3):
        np.add.at(total, triangles[:, k], faces)
    return total / np.maximum(np.linalg.norm(total, axis=1, keepdims=True), 1e-30)


def sun_visibility(positions, normals, triangles):
    """1 or 0 per vertex: a ray toward +z from just above the vertex, on the side of its surface
    that faces the Sun. The mesh is drawn from both sides, so a sheet whose normal points away from
    the Sun still has a side in the Sun. Triangles are binned by xy on a grid."""
    a, b, c = (positions[triangles[:, k]] for k in range(3))
    lo = positions[:, :2].min(0)
    size = (positions[:, :2].max(0) - lo).max() / 512.0
    cell = lambda p: np.floor((p - lo) / size).astype(np.int64)
    low = cell(np.minimum(np.minimum(a, b), c)[:, :2])
    high = cell(np.maximum(np.maximum(a, b), c)[:, :2])
    grid = {}
    for t in range(len(triangles)):
        for i in range(low[t, 0], high[t, 0] + 1):
            for j in range(low[t, 1], high[t, 1] + 1):
                grid.setdefault((i, j), []).append(t)
    grid = {key: np.array(value) for key, value in grid.items()}
    edge1, edge2 = b - a, c - a
    visible = np.ones(len(positions), dtype=np.float32)
    starts = positions + normals * np.sign(normals[:, 2:3]) * RAY_START_KM
    for v in range(len(positions)):
        p = starts[v]
        candidates = grid.get(tuple(cell(p[:2])))
        if candidates is None:
            continue
        # Moller-Trumbore with the ray direction +z.
        e1, e2 = edge1[candidates], edge2[candidates]
        det = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]
        usable = np.abs(det) > 1e-20
        s = p - a[candidates]
        with np.errstate(divide='ignore', invalid='ignore'):
            u = (s[:, 0] * e2[:, 1] - s[:, 1] * e2[:, 0]) / det
            w = (s[:, 1] * e1[:, 0] - s[:, 0] * e1[:, 1]) / det
            q = np.cross(s, e1)
            t = np.einsum('ij,ij->i', e2, q) / -det
            hit = usable & (u >= 0) & (w >= 0) & (u + w <= 1) & (t > RAY_START_KM)
        if hit.any():
            visible[v] = 0.0
    return visible


def build(name, spec, models):
    path = models / spec['file']
    glb = Stl(path, spec) if path.suffix == '.stl' or 'assembly' in spec else Obj(path) if path.suffix == '.obj' else Glb(path)
    glb.plain = tuple(spec.get('colors', {}))
    skip = spec.get('skip_materials', [])
    prims = [p for p in glb.primitives()
             if p[0] is None or glb.json['materials'][p[0]].get('name') not in skip]
    used = {glb.material_image(m) for m, *_ in prims if m is not None} - {None}

    atlas = Atlas(spec['atlas_width'])
    copies = {}    # a file can store the same image several times (Hubble's three of one)
    same_as = {}
    for index in sorted(used):
        info = glb.json['images'][index]
        image = Image.open(io.BytesIO(glb.view(info['bufferView']))).convert('RGB')
        scale = spec.get('image_scale', {}).get(info.get('name'), spec.get('image_default_scale', 1.0))
        pixels = np.array(image)
        for x0, y0, x1, y1 in spec.get('erase_px', {}).get(info.get('name'), []):
            # The median of the frame just outside the rectangle.
            frame = np.concatenate([pixels[y0 - 3:y0, x0:x1].reshape(-1, 3), pixels[y1:y1 + 3, x0:x1].reshape(-1, 3),
                                    pixels[y0:y1, x0 - 3:x0].reshape(-1, 3), pixels[y0:y1, x1:x1 + 3].reshape(-1, 3)])
            pixels[y0:y1, x0:x1] = np.median(frame, axis=0)
        if info.get('name') in spec.get('erase_polygon', {}):
            canvas = Image.fromarray(pixels)
            ImageDraw.Draw(canvas).polygon(spec['erase_polygon'][info['name']], fill=(130, 138, 152))
            pixels = np.array(canvas)
        if scale != 1.0:
            size = (int(pixels.shape[1] * scale), int(pixels.shape[0] * scale))
            pixels = np.array(Image.fromarray(pixels).resize(size, Image.LANCZOS))
        if info.get('name') in spec.get('tile_images', {}):
            columns, rows, cell, gain = spec['tile_images'][info['name']]
            tile = np.array(Image.fromarray(pixels).resize((cell, cell), Image.LANCZOS)) * gain
            pixels = np.tile(np.minimum(tile, 255.0).astype(np.uint8), (rows, columns, 1))
        first = copies.setdefault(pixels.tobytes(), index)
        if first != index:
            same_as[index] = first
            continue
        atlas.add_image(index, pixels)
    for material, *_ in prims:
        if material is None or glb.material_image(material) is None:
            atlas.add_color(glb.material_color(material, spec.get('colors', {})))
    atlas.build()
    for index, first in same_as.items():
        atlas.placed[index] = atlas.placed[first]

    z = np.array(spec['body_z'], dtype=np.float64)
    z /= np.linalg.norm(z)
    x = np.array(spec['body_x'], dtype=np.float64)
    to_body = np.array([x, np.cross(z, x), z])
    origin = np.array(spec['origin_m'], dtype=np.float64)

    all_positions, all_uv, all_triangles = [], [], []
    count = 0
    wrapped = 0
    grid_km = spec.get('merge_m', 0.0) / 1000.0
    for material, world, normals, uv, faces in prims:
        positions = (world - origin) @ to_body.T * spec['scale'] / 1000.0
        normals = normals @ to_body.T
        image = None if material is None else glb.material_image(material)
        if image is None:
            color = glb.material_color(material, spec.get('colors', {}))
            uv = np.tile(np.array(atlas.cells[color]), (len(positions), 1))
        else:
            # Triangles reaching outside 0..1 get their own vertices, moved by whole tiles.
            uv = uv.copy()
            corners = uv[faces]
            out = ((corners < -1e-4) | (corners > 1 + 1e-4)).any(axis=(1, 2))
            if out.any():
                wrapped += int(out.sum())
                ids = faces[out]
                shift = np.floor(corners[out].min(axis=1, keepdims=True))
                start = len(positions)
                positions = np.concatenate([positions, positions[ids].reshape(-1, 3)])
                normals = np.concatenate([normals, normals[ids].reshape(-1, 3)])
                uv = np.concatenate([uv, (corners[out] - shift).reshape(-1, 2)])
                faces = faces.copy()
                faces[out] = np.arange(start, start + 3 * len(ids)).reshape(-1, 3)
            tiles = spec.get('tile_images', {}).get(glb.json['images'][image].get('name'))
            if tiles is not None:
                uv = uv / np.array(tiles[:2], dtype=np.float64)
            uv = atlas.image_uv(image, uv)

        # Weld: the files keep up to three vertices per triangle. A vertex stays separate where
        # the authored normal or the texture differs, so hard edges stay hard.
        key = np.concatenate([np.round(positions * 1e7), np.round(normals * 1e2), np.round(uv * 1e5)], axis=1)
        _, first, inverse = np.unique(key.astype(np.int64), axis=0, return_index=True, return_inverse=True)
        inverse = inverse.reshape(-1)
        positions, normals, uv = positions[first], normals[first], uv[first]
        faces = inverse[faces]

        # Parts smaller than the grid are not seen from any distance the camera goes to: merge
        # vertices in one cell (with about the same normal, so that big flat faces stay flat).
        if grid_km > 0.0 and (material is None or glb.json['materials'][material].get('name') not in spec.get('merge_skip', [])):
            key = np.concatenate([np.floor(positions / grid_km), np.round(normals * 2.0), np.round(uv * spec.get('merge_uv', 1e5))], axis=1)
            _, first, inverse = np.unique(key.astype(np.int64), axis=0, return_index=True, return_inverse=True)
            inverse = inverse.reshape(-1)
            positions, normals, uv = positions[first], normals[first], uv[first]
            faces = inverse[faces]

        a, b, c = (positions[faces[:, k]] for k in range(3))
        cross = np.cross(b - a, c - a)
        keep = np.linalg.norm(cross, axis=1) > 1e-16
        faces, cross = faces[keep], cross[keep]
        if len(faces) == 0:
            continue
        # A few faces are wound against their authored normals: the normals decide.
        wrong = np.einsum('ij,ij->i', cross, normals[faces].mean(axis=1)) < 0.0
        faces[wrong] = faces[wrong][:, [0, 2, 1]]

        all_positions.append(positions)
        all_uv.append(uv)
        all_triangles.append(faces + count)
        count += len(positions)

    positions = np.concatenate(all_positions)
    uv = np.concatenate(all_uv)
    triangles = np.concatenate(all_triangles)
    if 'subdivide_m' in spec:
        positions, uv, triangles = subdivide(positions, uv, triangles, spec['subdivide_m'] / 1000.0)

    # Merging and degenerate faces leave vertices no triangle uses.
    used, inverse = np.unique(triangles, return_inverse=True)
    positions, uv, triangles = positions[used], uv[used], inverse.reshape(-1, 3)
    normals = vertex_normals(positions, triangles)
    visible = sun_visibility(positions, normals, triangles) if spec.get('bake_sun', True) \
        else np.ones(len(positions), dtype=np.float32)

    print(f'{name}: {len(positions)} vertices, {len(triangles)} triangles '
          f'({wrapped} wrapped); '
          f'atlas {atlas.pixels.shape[1]}x{atlas.pixels.shape[0]}, '
          f'{100.0 * float((visible == 0).sum()) / len(visible):.0f}% of the vertices shadowed')
    extent = np.linalg.norm(positions, axis=1).max() * 1000.0
    print(f'  bounds (m) {np.round(positions.min(0) * 1000, 2)} .. {np.round(positions.max(0) * 1000, 2)}, '
          f'farthest {extent:.2f} m from the origin')

    write_mesh(str(ROOT / 'assets' / 'shapes' / f'{name}.mesh'), positions.tolist(), visible.tolist(),
               triangles.tolist(), map_u=uv[:, 0].tolist(), hard_edges=True, map_v=uv[:, 1].tolist())
    target = ROOT / 'assets' / 'textures' / f'{name}.{spec["format"]}'
    image = Image.fromarray(atlas.pixels)
    if spec['format'] == 'jpg':
        image.save(target, quality=90, subsampling=0, optimize=True)
    else:
        image.save(target, optimize=True)
    print(f'-> {target} ({target.stat().st_size // 1024} KB)')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('models', type=Path, help='the directory the `file` entries start from')
    parser.add_argument('names', nargs='*', help=f'of {sorted(CRAFT)}; all by default')
    args = parser.parse_args()
    for name in args.names or sorted(CRAFT):
        build(name, CRAFT[name], args.models)


if __name__ == '__main__':
    main()
