"""Minimal reader for the comet attitude kernels of Rosetta (NAIF DAF/CK).

Reads CK segment type 5 (discrete quaternions, the type of ESA's CATT kernels)
and returns its samples as they are: the kernels hold a quaternion every
10 minutes, far denser than the slow change of the comet's spin needs. The
quaternions are SPICE's (cos(a/2), sin(a/2) axis), rotating vectors from the
segment's reference frame to the body frame. Times are encoded SCLK ticks.

Only the Python standard library is used.
"""

from daf import Daf

# Packet size of each type 5 subtype: quaternion with its derivative, angular
# velocity and its derivative (0, 2); quaternion (1); quaternion and angular
# velocity (3).
PACKET_SIZE = {0: 14, 1: 4, 2: 14, 3: 7}


class Segment:
    def __init__(self, instrument, frame, epochs, quaternions):
        self.instrument = instrument
        self.frame = frame
        self.epochs = epochs            # SCLK ticks, ascending
        self.quaternions = quaternions  # (cos, x sin, y sin, z sin) per epoch


class Ck:
    def __init__(self, path):
        daf = Daf(path, b'DAF/CK')
        self.comments = daf.comments
        self.segments = []
        for _, ic, w in daf.segments:
            instrument, frame, kind = ic[:3]
            if kind != 5:
                raise ValueError('%s: CK segment type %d not supported' % (path, kind))
            # Trailer: seconds per tick, subtype, window size, interpolation
            # intervals, packets.
            subtype, n = int(w[-4]), int(w[-1])
            size = PACKET_SIZE[subtype]
            quaternions = [tuple(w[size * i:size * i + 4]) for i in range(n)]
            epochs = list(w[size * n:size * n + n])
            self.segments.append(Segment(instrument, frame, epochs, quaternions))


def rotation(q):
    """The matrix of a SPICE quaternion: rows are the body axes in the reference frame."""
    s, x, y, z = q
    n = (s * s + x * x + y * y + z * z) ** 0.5
    s, x, y, z = s / n, x / n, y / n, z / n
    return ((1 - 2 * (y * y + z * z), 2 * (x * y - s * z), 2 * (x * z + s * y)),
            (2 * (x * y + s * z), 1 - 2 * (x * x + z * z), 2 * (y * z - s * x)),
            (2 * (x * z - s * y), 2 * (y * z + s * x), 1 - 2 * (x * x + y * y)))
