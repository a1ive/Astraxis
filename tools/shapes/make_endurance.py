"""Generate the original Endurance tribute hull and 12-module rotor (metres).

Run: python tools/shapes/make_endurance.py
Uses only the standard library; output is AXMESH2 in kilometres, with hard edges
(each flat face has its own vertices; tube walls share theirs and look round).
Design dimensions and materials are illustrative, not film production data.

The rotor turns about +z (the scene's [bodies.rotor]); the hull does not. The
floor of the modules (the inside of their outer wall) is at FLOOR_RADIUS_M, where
the scene's rotor period gives 1 g (tests).
"""
from pathlib import Path
from math import sin, cos, pi, sqrt
from axmesh import write_mesh, volume

ROOT = Path(__file__).resolve().parents[2]
MODULES = 12
RADIUS_M = 28.0        # module centres
FLOOR_RADIUS_M = 32.5  # modules reach out to 33 m; 0.5 m outer wall


class Mesh:
    def __init__(self):
        self.positions, self.albedo, self.triangles = [], [], []

    def vertex(self, p, shade):
        self.positions.append(tuple(c / 1000 for c in p))
        self.albedo.append(shade)
        return len(self.positions) - 1

    def face(self, points, shade):
        """A flat convex polygon, counter-clockwise seen from outside."""
        start = len(self.positions)
        for p in points:
            self.vertex(p, shade)
        self.triangles.extend((start, start + k, start + k + 1) for k in range(1, len(points) - 1))

    def prism(self, outline, z0, z1, shade, angle=0, center=(0, 0)):
        def at(p, z):
            x, y = p
            return (center[0] + x*cos(angle) - y*sin(angle),
                    center[1] + x*sin(angle) + y*cos(angle), z)
        bottom, top = [at(p, z0) for p in outline], [at(p, z1) for p in outline]
        self.face(bottom[::-1], shade)
        self.face(top, shade)
        for k in range(len(outline)):
            j = (k+1) % len(outline)
            self.face([bottom[k], bottom[j], top[j], top[k]], shade)

    def box(self, x, y, z0, z1, shade, angle=0, center=(0, 0), bevel=0):
        a, b, d = x/2, y/2, bevel
        outline = [(-a+d,-b),(a-d,-b),(a,-b+d),(a,b-d),
                   (a-d,b),(-a+d,b),(-a,b-d),(-a,-b+d)] if d else [(-a,-b),(a,-b),(a,b),(-a,b)]
        self.prism(outline, z0, z1, shade, angle, center)

    def wall(self, r, z0, z1, shade, segments, outward):
        """A smooth cylinder wall: the segments share their vertices."""
        def ring(z):
            return [self.vertex((r*cos(2*pi*k/segments), r*sin(2*pi*k/segments), z), shade)
                    for k in range(segments)]
        lo, hi = ring(z0), ring(z1)
        for k in range(segments):
            j = (k+1) % segments
            if outward:
                self.triangles += [(lo[k], lo[j], hi[j]), (lo[k], hi[j], hi[k])]
            else:
                self.triangles += [(lo[j], lo[k], hi[k]), (lo[j], hi[k], hi[j])]

    def tube(self, inner, outer, z0, z1, shade, segments=48):
        self.wall(outer, z0, z1, shade, segments, True)
        self.wall(inner, z0, z1, shade, segments, False)
        for k in range(segments):
            a, b = 2*pi*k/segments, 2*pi*(k+1)/segments
            def p(r, t, z): return (r*cos(t), r*sin(t), z)
            self.face([p(inner,a,z1),p(outer,a,z1),p(outer,b,z1),p(inner,b,z1)],shade)
            self.face([p(inner,b,z0),p(outer,b,z0),p(outer,a,z0),p(inner,a,z0)],shade)

    def area_shade(self):
        """Area-weighted sum of shade and the total area."""
        shade = area = 0.0
        for t in self.triangles:
            a, b, c = (self.positions[i] for i in t)
            u = [b[k]-a[k] for k in range(3)]
            v = [c[k]-a[k] for k in range(3)]
            n = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
            s = sqrt(sum(x*x for x in n)) / 2
            shade += s * sum(self.albedo[i] for i in t) / 3
            area += s
        return shade, area

    def save(self, name, mean_shade):
        assert volume(self.positions, self.triangles) > 0
        # The file holds albedo relative to the mean (the scene's color).
        albedo = [s / mean_shade for s in self.albedo]
        write_mesh(str(ROOT / 'assets/shapes' / name), self.positions, albedo, self.triangles, hard_edges=True)


def main():
    ring = Mesh()
    ring.tube(3.9, 5.2, -1.6, 1.6, .35)
    ring.tube(26.9, 29.1, -1.3, 1.3, .28, 96)
    for k in range(MODULES):
        a = 2*pi*k/MODULES
        c = (RADIUS_M*cos(a), RADIUS_M*sin(a))
        # Radial x / tangential y: separate chamfered pressure modules; module 0
        # is darker, so that the turning shows.
        ring.box(10, 9.2, -3.5, 3.5, .87 if k != 0 else .57, a, c, .9)
        for z0,z1 in [(-3.7,-3.51),(3.51,3.7)]:
            ring.box(6.8, 6.2, z0, z1, .22, a, c, .4)
            for y in [-2.1, 0, 2.1]:
                panel = (c[0]-y*sin(a), c[1]+y*cos(a))
                ring.box(6.5, .18, z0-.03, z1+.03, .64, a, panel)
        # Raised equipment pallet on each module's outside edge.
        outer = (34*cos(a),34*sin(a))
        ring.box(2, 5.8, -1.9, 1.9, .48, a, outer, .3)
    for k in range(3):
        a = 2*pi*k/3
        ring.box(22.0, 1.5, -1, 1, .58, a, (15.5*cos(a),15.5*sin(a)), .2)

    hull = Mesh()
    hull.prism([(3.4*cos(2*pi*k/8),3.4*sin(2*pi*k/8)) for k in range(8)], -9, 9, .82)
    for z0,z1 in [(-10.2,-8.8),(8.8,10.2)]:
        hull.tube(2.0, 3.7, z0, z1, .45, 24)
    hull.box(5.6, 5.6, 9.0, 9.2, .12, bevel=.7)
    # Stationary service frame behind the rotor plane; asymmetric equipment
    # makes the fixed attitude visible as the modules turn.
    hull.box(25, 1.3, -7.4, -6.1, .52, bevel=.2)
    for x in [-10,10]:
        hull.box(4.2, 7, -9, -5, .73, center=(x,0), bevel=.6)
        hull.box(3, 5, -5.0, -4.8, .18, center=(x,0), bevel=.2)
    hull.box(1, 11, -6.9, -5.9, .5, center=(0,7))
    hull.box(4, 3, -7.4, -5.4, .9, center=(0,12), bevel=.5)

    # One reference for both parts (they share the body's color).
    sums = [m.area_shade() for m in (ring, hull)]
    mean_shade = sum(s for s, _ in sums) / sum(a for _, a in sums)
    print(f'mean shade {mean_shade:.3f}')
    ring.save('endurance_ring.mesh', mean_shade)
    hull.save('endurance_hull.mesh', mean_shade)


if __name__ == '__main__':
    main()
