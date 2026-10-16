"""Minimal reader for the SPK files of Philae's trajectory (NAIF DAF/SPK).

Reads segment types 5 (discrete states) and 18 (Hermite or Lagrange packets),
the two used by ESA's lander kernels, and evaluates them by cubic Hermite
interpolation between neighbouring states: their states are 1-2 s apart, far
denser than the motion needs. Times are ET (TDB seconds since J2000), positions
km, velocities km/s, in the segment's frame relative to its center.

Only the Python standard library is used.
"""

import bisect

from daf import Daf


class Segment:
    def __init__(self, target, center, frame, kind, epochs, states):
        self.target = target
        self.center = center
        self.frame = frame
        self.kind = kind
        self.epochs = epochs  # ascending
        self.states = states  # (x, y, z, vx, vy, vz) per epoch

    def start(self):
        return self.epochs[0]

    def end(self):
        return self.epochs[-1]


class Spk:
    def __init__(self, path):
        daf = Daf(path, b'DAF/SPK')
        self.comments = daf.comments
        self.segments = []
        for _, ic, words in daf.segments:
            target, center, frame, kind = ic[:4]
            self.segments.append(self._segment(target, center, frame, kind, words))
        self.segments.sort(key=lambda s: s.start())

    @staticmethod
    def _segment(target, center, frame, kind, w):
        if kind == 5:
            n = int(w[-1])
            states = [w[6 * i:6 * i + 6] for i in range(n)]
            epochs = list(w[6 * n:7 * n])
        elif kind == 18:
            subtype, _, n = (int(x) for x in w[-3:])
            size = {0: 12, 1: 6}[subtype]
            packets = [w[size * i:size * i + size] for i in range(n)]
            # Subtype 0: position, its derivative, velocity, its derivative;
            # subtype 1: position, velocity.
            states = [p[0:3] + p[6:9] if subtype == 0 else p for p in packets]
            epochs = list(w[size * n:size * n + n])
        else:
            raise ValueError('SPK segment type %d not supported' % kind)
        return Segment(target, center, frame, kind, epochs, states)

    def span(self):
        return self.segments[0].start(), self.segments[-1].end()

    def state(self, t):
        """(position, velocity) at t, by cubic Hermite interpolation; None outside the segments."""
        for seg in self.segments:
            if seg.start() <= t <= seg.end():
                break
        else:
            return None
        i = min(max(bisect.bisect_right(seg.epochs, t) - 1, 0), len(seg.epochs) - 2)
        t0, t1 = seg.epochs[i], seg.epochs[i + 1]
        a, b = seg.states[i], seg.states[i + 1]
        h = t1 - t0
        s = (t - t0) / h
        h00, h10, h01, h11 = 2 * s ** 3 - 3 * s ** 2 + 1, s ** 3 - 2 * s ** 2 + s, -2 * s ** 3 + 3 * s ** 2, s ** 3 - s ** 2
        d00, d10, d01, d11 = 6 * s ** 2 - 6 * s, 3 * s ** 2 - 4 * s + 1, -6 * s ** 2 + 6 * s, 3 * s ** 2 - 2 * s
        p = tuple(h00 * a[k] + h10 * h * a[k + 3] + h01 * b[k] + h11 * h * b[k + 3] for k in range(3))
        v = tuple((d00 * a[k] + d01 * b[k]) / h + d10 * a[k + 3] + d11 * b[k + 3] for k in range(3))
        return p, v
