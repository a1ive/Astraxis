"""Build the original Endurance tribute in Blender; see assets/shapes/SOURCES.md.

blender -b --factory-startup -P tools/shapes/make_endurance.py
blender -b ref/Endurance/Endurance.blend -P tools/shapes/make_endurance.py -- --export-only
"""
import math
import sys
from pathlib import Path

import bpy
from mathutils import Vector

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from axmesh import volume, write_mesh

PALETTE = {'Ceramic': .82, 'Seam': .38, 'Graphite': .105, 'Glass': .045,
           'Titanium': .46, 'Thermal blanket': .66, 'Recess': .07}


def circle(radius, count):
    return [(radius * math.cos(math.tau * index / count), radius * math.sin(math.tau * index / count))
            for index in range(count)]


def plate(width, length, cut):
    half, side = width / 2, length / 2
    return [(-half+cut, -side), (half-cut, -side), (half, -side+cut), (half, side-cut),
            (half-cut, side), (-half+cut, side), (-half, side-cut), (-half, -side+cut)]


class Part:
    def __init__(self):
        self.points, self.faces, self.materials = [], [], []

    def loft(self, layers, material, bands=(), closed=False):
        base, count = len(self.points), len(layers[0])
        for layer in layers:
            self.points.extend(layer)
        if not closed:
            self.faces.append(tuple(range(base+count-1, base-1, -1)))
            self.materials.append(material)
        for step in range(len(layers) if closed else len(layers)-1):
            bottom = base + step*count
            top = base + ((step+1) % len(layers))*count
            for index in range(count):
                after = (index+1) % count
                self.faces.append((bottom+index, bottom+after, top+after, top+index))
                self.materials.append(bands[step] if bands else material)
        if not closed:
            top = base + (len(layers)-1)*count
            self.faces.append(tuple(range(top, top+count)))
            self.materials.append(material)

    def stack(self, center, outline, levels, material, bands=()):
        self.loft([[(center[0]+scale*xx, center[1]+scale*yy, height) for xx, yy in outline]
                   for height, scale in levels], material, bands)

    def tile(self, xx, yy, width, length, low, high, material, cut=.15):
        self.stack((xx, yy), plate(width, length, cut), [(low, 1), (high, 1)], material)

    def beam(self, start, end, radius, material, sides=8):
        axis = Vector(end) - Vector(start)
        across = axis.orthogonal().normalized()
        up = axis.normalized().cross(across)
        self.loft([[Vector(point)+xx*across+yy*up for xx, yy in circle(radius, sides)]
                   for point in (start, end)], material)

    def finish(self, collection, name, angle=0):
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(self.points, [], self.faces)
        mesh.update()
        for material in PALETTE:
            mesh.materials.append(bpy.data.materials[material])
        indices = {material: index for index, material in enumerate(PALETTE)}
        mesh.polygons.foreach_set('material_index', [indices[name] for name in self.materials])
        obj = bpy.data.objects.new(name, mesh)
        obj.rotation_euler.z = angle
        collection.objects.link(obj)


def make_module(rotor, station):
    module = Part()
    family = station % 3
    angle = math.tau * station / 12
    outline = plate(10, 11.0 if family == 1 else 10.4, 1.45)
    module.stack((28, 0), outline, [(-3.75, .8), (-3.1, 1), (-.48, 1), (-.42, .965),
                 (.42, .965), (.48, 1), (3.1, 1), (3.75, .8)], 'Ceramic',
                 ['Thermal blanket', 'Ceramic', 'Seam', 'Graphite', 'Seam', 'Ceramic', 'Thermal blanket'])
    module.tile(33.7, 0, 2.4, 5.2, -1.7, 1.7, 'Graphite', .45)
    for offset in (-1.45, 0, 1.45):
        module.tile(34.55, offset, .8, .72, -1.52, 1.52, 'Titanium')
    for sign in (-1, 1):
        def panel(xx, yy, width, length, base, top, material, cut=.15):
            module.tile(xx, yy, width, length, *sorted((sign*base, sign*top)), material, cut)
        if family == 0:
            panel(28, 0, 7.4, 7.8, 3.7, 3.98, 'Graphite', .75)
            for offset in (-2.25, 2.25):
                levels = [(3.94, 1.48), (4.45, 1.12), (5.38, 1.66), (5.38, 1.30), (4.47, .75)]
                layers = [[(28+xx, offset+yy, sign*height) for xx, yy in circle(radius, 16)]
                          for height, radius in levels]
                module.loft([layer if sign > 0 else layer[::-1] for layer in layers], 'Recess',
                            ['Titanium', 'Graphite', 'Titanium', 'Recess'])
        elif family == 1:
            panel(28, 0, 6.9, 8.25, 3.7, 4.08, 'Graphite', .65)
            for offset in (-1.65, 1.65):
                module.stack((28+offset, 0), plate(2.96, 7.75, .56),
                             sorted([(sign*4.02, 1), (sign*4.68, .88)]), 'Ceramic')
                for seam in (-1.85, 1.85):
                    panel(28+offset, seam, 2.65, .35, 4.4, 4.7, 'Seam', .09)
        elif station == 2:
            module.stack((27.3, 0), plate(7.5, 8.4, 1.3), sorted([(sign*3.7, 1), (sign*4.2, 1),
                         (sign*5.55, .78), (sign*5.8, .76)]), 'Ceramic',
                         ['Ceramic', 'Glass', 'Titanium'] if sign > 0 else ['Titanium', 'Glass', 'Ceramic'])
            panel(28, 0, 3.8, 4.8, 5.78, 5.92, 'Thermal blanket', .55)
        else:
            panel(28, 0, 7.2, 8, 3.7, 4.05, 'Graphite', .6)
            for offset in (-2.8, -1.4, 0, 1.4, 2.8):
                panel(27.3, offset, 5.25, .64, 4, 4.47, 'Titanium')
            panel(30.9, 0, .62, 6.9, 4, 4.6, 'Ceramic')
    kind = 'Command' if station == 2 else ('Propulsion', 'Habitat', 'Radiator')[family]
    module.finish(rotor, f'{station+1:02d} {kind}', angle)
    tunnel = Part()
    middle = 28 * math.cos(math.pi / 12)
    for offset in (-1.75, 0, 1.75):
        radius = 1.25 if offset == 0 else .31
        tunnel.beam((middle+offset, -4, 0), (middle+offset, 4, 0), radius, 'Titanium')
    for offset in (-2, -.7, .7, 2):
        tunnel.beam((middle, offset-.18, 0), (middle, offset+.18, 0), 1.42, 'Graphite')
    tunnel.finish(rotor, f'{station+1:02d} Transfer tunnel', angle+math.pi/12)


def make_spokes(rotor):
    bearing = Part()
    bearing.loft([[(xx, yy, height) for xx, yy in circle(radius, 32)] for height, radius in
                  [(-1.6, 5.2), (1.6, 5.2), (1.6, 3.82), (-1.6, 3.82)]], 'Titanium', closed=True)
    bearing.finish(rotor, 'Rotary bearing')
    for index in range(3):
        truss = Part()
        for side in (-1, 1):
            truss.beam((4.65, side*1.1, 0), (23.7, side*2.15, 0), .42, 'Ceramic')
        for bay in range(5):
            side = (-1)**bay
            truss.beam((4.8+bay*3.7, side*(1.1+bay*.21), 0),
                       (8.5+bay*3.7, -side*(1.31+bay*.21), 0), .28, 'Titanium')
        truss.finish(rotor, f'Truss {index+1}', math.tau * index / 3)


def make_hull(hull):
    spine = Part()
    spine.stack((0, 0), circle(1, 16), [(-10.2, 3.6), (-9.7, 3.6), (-9.7, 3), (-6, 3),
                (-5.5, 3.45), (5.5, 3.45), (6, 3), (9.7, 3), (9.7, 3.6), (10.2, 3.6)],
                'Titanium', ['Ceramic', 'Seam', 'Titanium', 'Seam', 'Ceramic', 'Seam',
                            'Titanium', 'Seam', 'Ceramic'])
    for side in (-1, 1):
        spine.loft([[(xx, yy, side*height) for xx, yy in (circle(radius, 24) if side > 0 else
                    circle(radius, 24)[::-1])] for height, radius in [(10.19, 2.75), (10.6, 2.75),
                    (10.6, 2.12), (10.25, 2.12)]], 'Graphite', ['Seam', 'Ceramic', 'Recess'])
        spine.tile(0, 0, 2.3, .5, *sorted((side*10.24, side*10.34)), 'Titanium')
        spine.beam((0, 0, side*6.95), (side*7.5, 0, side*6.95), .98, 'Titanium')
    spine.finish(hull, 'Axial docking spine')
    ranger = Part()
    plan = [(-2.1,-6.6), (2.1,-6.6), (4.85,-4.6), (3.85,.1), (1.2,7.8),
            (-1.2,7.8), (-3.85,.1), (-4.85,-4.6)]
    ranger.loft([[(7.7+xx*scale, yy*scale+shift, height+.09*yy) for xx, yy in plan]
                 for height, scale, shift in [(5.9, .86, 0), (6.55, 1, 0), (7.15, 1, 0), (8.65, .61, -1.65)]],
                'Ceramic', ['Titanium', 'Graphite', 'Ceramic'])
    canopy = [(-1.7,-1.9), (1.7,-1.9), (2.0,-.8), (.9,2.7), (-.9,2.7), (-2.0,-.8)]
    ranger.loft([[(7.7+xx*scale, yy*scale+1, height+.09*(yy*scale+2.65)) for xx, yy in canopy]
                 for scale, height in [(1, 8.58), (.86, 9.08)]], 'Glass')
    for offset in (-1.4, 1.4):
        for ypos in (-4.65, -3.9, -3.15):
            ranger.tile(7.7+offset, ypos, 1.65, .33, 8.25, 8.72+.09*(ypos+1.65), 'Graphite', .08)
    ranger.finish(hull, 'Ranger - lifting body')
    lander = Part()
    lander.stack((-7.9, 0), plate(5.35, 12.6, 1.25), [(-9.25,.76), (-8,1), (-5.8,1), (-5.25,.82)], 'Ceramic',
                 ['Ceramic', 'Titanium', 'Thermal blanket'])
    for offset in (-3.6, 3.6):
        lander.stack((-7.9+offset, -.4), plate(2.1, 9.9, .65), [(-9.55,.77), (-8.75,1), (-6.4,1)], 'Graphite',
                     ['Thermal blanket', 'Graphite'])
    lander.tile(-7.9, 3, 3.5, 2.45, -9.42, -9.15, 'Glass', .45)
    lander.finish(hull, 'Lander - twin side pods')


def export(collection):
    points, shades, triangles = [], [], []
    area_sum = shade_sum = 0.0
    graph = bpy.context.evaluated_depsgraph_get()
    for obj in sorted(collection.all_objects, key=lambda obj: obj.name):
        if obj.type != 'MESH':
            continue
        evaluated = obj.evaluated_get(graph)
        mesh = evaluated.to_mesh()
        mesh.calc_loop_triangles()
        indices = {}
        for triangle in mesh.loop_triangles:
            material = evaluated.material_slots[triangle.material_index].material
            shade = sum(material.diffuse_color[:3]) / 3
            face = []
            for vertex in triangle.vertices:
                key = (triangle.polygon_index, vertex)
                if key not in indices:
                    indices[key] = len(points)
                    points.append((evaluated.matrix_world @ mesh.vertices[vertex].co) * .001)
                    shades.append(shade)
                face.append(indices[key])
            if evaluated.matrix_world.is_negative:
                face.reverse()
            first, second, third = (points[index] for index in face)
            area = (second-first).cross(third-first).length / 2
            assert area > 0, obj.name
            area_sum += area
            shade_sum += area * shade
            triangles.append(face)
        evaluated.to_mesh_clear()
    assert volume(points, triangles) > 0, collection.name
    return points, shades, triangles, shade_sum, area_sum


def main():
    if '--export-only' not in sys.argv:
        bpy.ops.wm.read_factory_settings(use_empty=True)
        scene = bpy.context.scene
        scene.unit_settings.system = 'METRIC'
        for name, shade in PALETTE.items():
            bpy.data.materials.new(name).diffuse_color = (shade, shade, shade, 1)
        hull, rotor = [bpy.data.collections.new(name) for name in ('Hull', 'Rotor')]
        for collection in (hull, rotor):
            scene.collection.children.link(collection)
        make_hull(hull)
        make_spokes(rotor)
        for station in range(12):
            make_module(rotor, station)
        for area in bpy.context.screen.areas:
            if area.type == 'VIEW_3D':
                space = area.spaces.active
                space.overlay.show_overlays = False
                space.shading.color_type = 'MATERIAL'
                space.shading.show_cavity = True
                space.region_3d.view_distance = 97
                space.region_3d.view_rotation = Vector((1, -1, 1.7)).to_track_quat('Z', 'Y')
        path = ROOT / 'ref/Endurance/Endurance.blend'
        path.parent.mkdir(parents=True, exist_ok=True)
        bpy.ops.wm.save_as_mainfile(filepath=str(path))
    bpy.context.view_layer.update()
    meshes = {part: export(bpy.data.collections[name]) for part, name in [('hull', 'Hull'), ('ring', 'Rotor')]}
    mean = sum(mesh[3] for mesh in meshes.values()) / sum(mesh[4] for mesh in meshes.values())
    for part, (points, shades, triangles, _, _) in meshes.items():
        write_mesh(str(ROOT / f'assets/shapes/endurance_{part}.mesh'), points, [shade/mean for shade in shades],
                   triangles, hard_edges=True)
        print(f'{part}: {len(triangles)} triangles')


if __name__ == '__main__':
    main()
