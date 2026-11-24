"""Build the editable ISV Venture Star and export the OBJ for make_spacecraft.py."""
import math
import sys
from pathlib import Path

import bpy
from mathutils import Vector

WORK = Path(__file__).resolve().parents[2] / 'ref/isv'
# Linear RGB. The OBJ's Kd values also supply the runtime palette.
COLORS = {
    'Ceramic': (.77, .81, .80), 'Titanium': (.46, .52, .55),
    'Graphite': (.067, .085, .10), 'Insulation': (.52, .30, .10),
    'Copper': (.63, .19, .089), 'Radiator': (.39, .057, .032),
    'Hot core': (.86, .22, .085), 'Glazing': (.05, .17, .23),
}


def collection(name):
    group = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(group)
    return group


def loft(group, name, sections, material, bands=()):
    """Join equally sized, counterclockwise sections with flat panels."""
    count = len(sections[0])
    ends = (len(sections) - 1) * count
    faces = [tuple(reversed(range(count)))]
    for base in range(0, ends, count):
        for i in range(count):
            j = (i + 1) % count
            faces.append((base+i, base+j, base+j+count, base+i+count))
    faces.append(tuple(range(ends, ends+count)))
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata([p for section in sections for p in section], [], faces)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    group.objects.link(obj)
    for mat in dict.fromkeys((material, *bands)):
        mesh.materials.append(bpy.data.materials[mat])
    for i, mat in enumerate(bands):
        index = mesh.materials.find(mat)
        for face in range(1+i*count, 1+(i+1)*count):
            mesh.polygons[face].material_index = index
    return obj


def spindle(group, name, start, end, profile, material, sides=12, bands=()):
    """Surface of revolution, with stations in metres from start along its axis."""
    axis = (Vector(end) - Vector(start)).normalized()
    side = axis.orthogonal().normalized()
    up = axis.cross(side)
    ring = [side * math.cos(math.tau*i/sides) + up * math.sin(math.tau*i/sides)
            for i in range(sides)]
    return loft(group, name, [[Vector(start) + axis*x + ray*r for ray in ring]
                              for x, r in profile], material, bands)


def strut(group, name, a, b, radius, mat='Titanium'):
    spindle(group, name, a, b, [(0, radius), ((Vector(b)-Vector(a)).length, radius)], mat, 6)


def shell(group, name, x, y, z, sections, mat, bands=()):
    # Chamfered rectangular cross sections in the yz plane.
    outline = [(-.72,-1), (.72,-1), (1,-.72), (1,.72), (.72,1), (-.72,1), (-1,.72), (-1,-.72)]
    return loft(group, name, [[(x+dx, y+u*w, z+v*h) for u,v in outline]
                              for dx,w,h in sections], mat, bands)


def backbone():
    group = collection('Tapered tension truss')
    stations = []
    for i in range(20):
        x = -475 + 57*i
        r = max(8, 38 - 2.3*i)
        stations.append([(x, r*math.cos(a), r*math.sin(a)) for a in (0, math.tau/3, 2*math.tau/3)])
    for i, station in enumerate(stations):
        for j, corner in enumerate(station):
            k = (j+1) % 3
            strut(group, 'Transverse frame', corner, station[k], 1.6)
            if i+1 < len(stations):
                strut(group, 'Longeron', corner, stations[i+1][j], 2.7, 'Ceramic')
                strut(group, 'Diagonal web', corner, stations[i+1][k], 1.65)


def propulsion():
    group = collection('Tractor drives and heat rejection')
    for sign in (-1, 1):
        y = sign*102
        a, b = (-750,y,0), (-400,y,0)
        spindle(group, 'Radiator manifold', a, b, [(0,17),(12,25),(238,32),(263,29)],
                'Graphite', 16, ['Titanium','Hot core','Graphite'])
        for i in range(31):
            x = 6+i*8
            r = 32 + 10*i/30
            spindle(group, 'Radiator lamella', a, b, [(x,r-3),(x+1,r),(x+4,r),(x+5,r-3)],
                    'Radiator', 16, ['Copper','Graphite','Radiator'])
        for z in (-20,20):
            strut(group, 'Coolant return', (-744,y,z), (-499,y,z*1.36), 2.1, 'Copper')
        for angle in range(0,360,45):
            c, s = math.cos(math.radians(angle)), math.sin(math.radians(angle))
            strut(group, 'Longitudinal header', (-746,y+32*c,32*s), (-502,y+42*c,42*s), 2.0, 'Graphite')
        # Nozzle axis cants outboard away from the central truss.
        a, b = (-496,y,0), (-340,y+sign*41,0)
        spindle(group, 'Tractor engine', a, b,
                [(0,30),(17,32),(28,24),(69,23),(84,16),(109,13),(143,22),(143,18),(119,9)],
                'Graphite', 20, ['Titanium','Ceramic','Graphite','Titanium','Copper','Titanium','Ceramic','Graphite'])
        for x in (33,43,53,63):
            spindle(group, 'Field coil', a, b, [(x,24),(x+2,27),(x+6,27),(x+8,24)], 'Copper', 16)
        for z in (-18,18):
            strut(group, 'Engine yoke', (-520,y,z), (-458,0,0), 7)
            strut(group, 'Tank brace', (-390,sign*31,z), (-495,y,z), 4)
        for x,r in [(-468,42),(-388,35)]:
            p = (x-r,sign*43,0)
            profile = [(r*(1-math.cos(math.pi*i/12)), max(0.8,r*math.sin(math.pi*i/12))) for i in range(13)]
            spindle(group, 'Banded propellant tank', p, (x,sign*43,0), profile, 'Titanium', 24,
                    ['Insulation' if i in (2,8) else 'Ceramic' for i in range(12)])
    strut(group, 'Axial service mast', (-735,0,0), (-465,0,0), 3.8)


def payload():
    group = collection('Cargo, cryovaults and docked shuttle')
    for i in range(4):
        x = 238+36*i
        for sign in (-1,1):
            y = sign*15
            shell(group, 'Freight pallet', x,y,0, [(0,12,13),(3,14,15),(29,14,15),(32,12,13)],
                  'Insulation', ['Titanium','Insulation','Titanium'])
            for offset in (7,22):
                shell(group, 'Pallet strap', x+offset,y,0, [(0,14.5,15.5),(2.5,14.5,15.5)], 'Graphite')
    for y in (-19,19):
        shell(group, 'Cryovault jacket', 394,y,0, [(0,8,10),(11,17,18),(72,17,18),(83,8,10)],
              'Ceramic', ['Titanium','Ceramic','Titanium'])
        for x in (417,435,453):
            shell(group, 'Pressure hoop', x,y,0, [(0,17.5,18.5),(3,17.5,18.5)], 'Titanium')
    shell(group, 'Transfer hub', 510,0,0, [(0,12,13),(9,23,22),(36,23,22),(46,14,13)],
          'Titanium', ['Graphite','Ceramic','Titanium'])
    for sign in (-1,1):
        strut(group, 'Habitat arm', (531,0,0), (531,sign*82,sign*40), 5)
        shell(group, 'Habitat gondola', 507,sign*80,sign*39,
              [(0,9,7),(7,15,11),(42,15,11),(49,10,7)], 'Ceramic', ['Graphite','Ceramic','Titanium'])
        for x in (520,540):
            shell(group, 'Gondola belt', x,sign*80,sign*39, [(0,15.6,11.6),(3,15.6,11.6)], 'Graphite')
    strut(group, 'Docking collar', (439,0,10), (439,0,32), 6)
    shell(group, 'Valkyrie fuselage', 403,0,36,
          [(0,11,5),(9,13,7),(43,9,8),(59,6,5),(72,1.4,1.5)], 'Ceramic',
          ['Graphite','Ceramic','Glazing','Titanium'])
    for sign in (-1,1):
        # Planform in the xy plane, with the narrow end swept aft.
        wing = [(409,sign*7), (401,sign*38), (413,sign*37), (449,sign*8)]
        if sign > 0:
            wing.reverse()
        loft(group, 'Shuttle wing', [[(x,y,z) for x,y in wing] for z in (34,37)], 'Titanium')


def shields():
    group = collection('Interstellar debris shields')
    strut(group, 'Shield spine', (552,0,0), (747,0,0), 3.6)
    for i,(x,r) in enumerate([(627,66),(670,78),(711,88),(746,98)]):
        # Separate clipped plates and radial ribs read edge-on as well as face-on.
        shell(group, 'Debris shield '+str(i+1), x,0,0, [(0,r*.97,r*.78),(2,r,r*.81),(4,r,r*.81)],
              'Graphite', ['Titanium','Ceramic'])
        for y,z in ((.70,.74),(-.70,.74),(-.70,-.74),(.70,-.74)):
            strut(group, 'Shield stay', (x-24,0,0), (x,y*r,z*r), 1.45)
            strut(group, 'Face rib', (x+4,0,0), (x+4,y*r,z*r), 1.1, 'Titanium')
        for y in (-.5,0,.5):
            strut(group, 'Panel joint', (x,y*r,-.77*r), (x,y*r,.77*r), .85)
        for z in (-.38,0,.38):
            strut(group, 'Panel joint', (x,-.97*r,z*r), (x,.97*r,z*r), .85)


def main():
    if '--export-only' not in sys.argv:
        bpy.ops.wm.read_factory_settings(use_empty=True)
        bpy.context.scene.unit_settings.system = 'METRIC'
        for name, rgb in COLORS.items():
            mat = bpy.data.materials.new(name)
            mat.diffuse_color = (*rgb, 1)
            mat.node_tree.nodes.get('Principled BSDF').inputs['Base Color'].default_value = (*rgb, 1)
        backbone()
        propulsion()
        payload()
        shields()
        for area in bpy.context.screen.areas:
            if area.type == 'VIEW_3D':
                space = area.spaces.active
                space.clip_end = 10000
                space.overlay.show_overlays = False
                space.shading.color_type = 'MATERIAL'
                space.region_3d.view_distance = 2100
                space.region_3d.view_rotation = Vector((.52,-1,.67)).to_track_quat('Z','Y')
        WORK.mkdir(parents=True, exist_ok=True)
        bpy.ops.wm.save_as_mainfile(filepath=str(WORK / 'Venture Star.blend'))
    bpy.ops.wm.obj_export(filepath=str(WORK / 'venture_star.obj'), forward_axis='Y', up_axis='Z',
                          export_triangulated_mesh=True)


if __name__ == '__main__':
    main()
